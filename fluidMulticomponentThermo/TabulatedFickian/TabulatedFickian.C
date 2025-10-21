/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "TabulatedFickian.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "Function2Evaluate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template <class BasicThermophysicalTransportModel>
void TabulatedFickian<BasicThermophysicalTransportModel>::updateDm()
{
    const PtrList<volScalarField> &Y = this->thermo().Y();
    const volScalarField &p = this->thermo().p();
    const volScalarField &T = this->thermo().T();
    const volScalarField Wm(this->thermo().W());
    volScalarField sumXbyD(
        volScalarField::New("sumXbyD",
                            T.mesh(),
                            dimless / dimKinematicViscosity / Wm.dimensions()));
    volScalarField Dij(IOobject("Dij",
                                T.mesh().time().name(),
                                T.mesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE,
                                false),
                       T.mesh(),
                       dimKinematicViscosity);
    forAll(Dm_, i)
    {
        sumXbyD = Zero;

        forAll(Y, j)
        {
            forAll(Dij, cells)
            {
                Dij[cells] = getDij(DpolyCoeffs_[i][j], p[cells], T[cells]);
            }

            forAll(Dij.boundaryField(), patchi)
            {
                forAll(Dij.boundaryField()[patchi], facei)
                {
                    Dij.boundaryFieldRef()[patchi][facei] =
                        getDij(DpolyCoeffs_[i][j],
                               p.boundaryField()[patchi][facei],
                               T.boundaryField()[patchi][facei]);
                }
            }

            if (j != i)
            {
                sumXbyD += Y[j] / (this->thermo().Wi(j) * Dij);
            }
        }

        Dm_.set(
            i,
            (1 / Wm - Y[i] / this->thermo().Wi(i)) /
                max(sumXbyD, dimensionedScalar(sumXbyD.dimensions(), small)));
    }
}

template <class BasicThermophysicalTransportModel>
const PtrList<volScalarField>
    &TabulatedFickian<BasicThermophysicalTransportModel>::Dm() const
{
    if (!Dm_.size())
    {
        updateDm();
    }

    return Dm_;
}

template <class BasicThermophysicalTransportModel>
TabulatedFickian<BasicThermophysicalTransportModel>::TabulatedFickian(
    const word &type,
    const momentumTransportModel &momentumTransport,
    const thermoModel &thermo) :
    BasicThermophysicalTransportModel(type, momentumTransport, thermo),
    TopoChangeableMeshObject(*this),
    DTFuncs_(this->coeffDict().found("DT") ? this->thermo().species().size()
                                           : 0),
    Dm_(this->thermo().species().size())
{
    // Construct DpolyCoeffs_
    DpolyCoeffs_.setSize(this->thermo().species().size());
    forAll(DpolyCoeffs_, i)
    {
        DpolyCoeffs_[i].setSize(this->thermo().species().size());
        forAll(DpolyCoeffs_[i], j)
        {
            DpolyCoeffs_[i][j] = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * *
// //

template <class BasicThermophysicalTransportModel>
bool TabulatedFickian<BasicThermophysicalTransportModel>::read()
{

    if (BasicThermophysicalTransportModel::read())
    {
        const speciesTable &species = this->thermo().species();

        const dictionary &coeffDict = this->coeffDict();

        // coeffDict.lookup("mixtureDiffusionMode")
        //     >> mixtureDiffusionMode_;

        const dictionary &Ddict = coeffDict.subDict("DPolyCoeffs");

        forAll(species, si)
        {
            forAll(species, sj)
            {
                word nameij = species[si] + '-' + species[sj];
                word nameji = species[sj] + '-' + species[si];

                // Read the array of specie binary mass diffusion
                // coefficient polynomials
                if (Ddict.found(nameij))
                {
                    FixedList<scalar, 5> polyCoeffs(
                        Ddict.subDict(nameij).lookup("polyCoeffs"));
                    DpolyCoeffs_[si][sj] = polyCoeffs;
                    DpolyCoeffs_[sj][si] = polyCoeffs;
                }
                else if (Ddict.found(nameji))
                {
                    FixedList<scalar, 5> polyCoeffs(
                        Ddict.subDict(nameji).lookup("polyCoeffs"));
                    DpolyCoeffs_[si][sj] = polyCoeffs;
                    DpolyCoeffs_[sj][si] = polyCoeffs;
                }
                else
                {
                    FatalIOErrorInFunction(Ddict)
                        << "Binary mass diffusion polynomials coefficient "
                           "for "
                           "pair "
                        << nameij << " or " << nameji << " not provided"
                        << exit(FatalIOError);
                }
            }
        }

        Info << "Reading binary diffusion coefficients complete" << endl;

        // Optionally read the List of specie Soret thermal diffusion
        // coefficient functions
        if (coeffDict.found("DT"))
        {
            const dictionary &DTdict = coeffDict.subDict("DT");

            forAll(species, i)
            {
                DTFuncs_.set(i,
                             Function2<scalar>::New(species[i],
                                                    dimPressure,
                                                    dimTemperature,
                                                    dimKinematicViscosity,
                                                    DTdict)
                                 .ptr());
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


template <class BasicThermophysicalTransportModel>
tmp<volScalarField> TabulatedFickian<BasicThermophysicalTransportModel>::DEff(
    const volScalarField &Yi) const
{
    return volScalarField::New(
        "DEff",
        this->momentumTransport().rho() * Dm_[this->thermo().specieIndex(Yi)]);
}


template <class BasicThermophysicalTransportModel>
tmp<scalarField> TabulatedFickian<BasicThermophysicalTransportModel>::DEff(
    const volScalarField &Yi, const label patchi) const
{
    return this->momentumTransport().rho().boundaryField()[patchi] *
           Dm_[this->thermo().specieIndex(Yi)].boundaryField()[patchi];
}


template <class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> TabulatedFickian<BasicThermophysicalTransportModel>::q()
    const
{
    tmp<surfaceScalarField> tmpq(surfaceScalarField::New(
        IOobject::groupName("q",
                            this->momentumTransport().alphaRhoPhi().group()),
        -fvc::interpolate(this->alpha() * this->kappaEff()) *
            fvc::snGrad(this->thermo().T())));

    const PtrList<volScalarField> &Y = this->thermo().Y();
    const volScalarField &p = this->thermo().p();
    const volScalarField &T = this->thermo().T();

    if (Y.size())
    {
        surfaceScalarField sumJ(surfaceScalarField::New(
            "sumJ",
            Y[0].mesh(),
            dimensionedScalar(dimMass / dimArea / dimTime, 0)));

        surfaceScalarField sumJh(surfaceScalarField::New(
            "sumJh",
            Y[0].mesh(),
            dimensionedScalar(sumJ.dimensions() * dimEnergy / dimMass, 0)));

        forAll(Y, i)
        {
            if (i != this->thermo().defaultSpecie())
            {
                const volScalarField hi(this->thermo().hsi(i, p, T));

                const surfaceScalarField ji(this->j(Y[i]));
                sumJ += ji;

                sumJh += ji * fvc::interpolate(hi);
            }
        }

        {
            const label i = this->thermo().defaultSpecie();

            const volScalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ * fvc::interpolate(hi);
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}

template <class BasicThermophysicalTransportModel>
tmp<scalarField> TabulatedFickian<BasicThermophysicalTransportModel>::q(
    const label patchi) const
{
    tmp<scalarField> tmpq(
        -(this->alpha().boundaryField()[patchi] * this->kappaEff(patchi) *
          this->thermo().T().boundaryField()[patchi].snGrad()));

    const PtrList<volScalarField> &Y = this->thermo().Y();
    const scalarField &p = this->thermo().p().boundaryField()[patchi];
    const scalarField &T = this->thermo().T().boundaryField()[patchi];

    if (Y.size())
    {
        scalarField sumJ(tmpq->size(), scalar(0));
        scalarField sumJh(tmpq->size(), scalar(0));

        forAll(Y, i)
        {
            if (i != this->thermo().defaultSpecie())
            {
                const scalarField hi(this->thermo().hsi(i, p, T));

                const scalarField ji(this->j(Y[i], patchi));

                sumJ += ji;

                sumJh += ji * hi;
            }
        }

        {
            const label i = this->thermo().defaultSpecie();

            const scalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ * hi;
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}

template <class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> TabulatedFickian<BasicThermophysicalTransportModel>::divq(
    volScalarField &he) const
{
    tmp<fvScalarMatrix> tmpDivq(fvm::Su(
        -fvc::laplacian(this->alpha() * this->kappaEff(), this->thermo().T()),
        he));

    const PtrList<volScalarField> &Y = this->thermo().Y();
    const volScalarField &p = this->thermo().p();
    const volScalarField &T = this->thermo().T();

    tmpDivq.ref() -=
        fvm::laplacianCorrection(this->alpha() * this->alphaEff(), he);

    surfaceScalarField sumJ(surfaceScalarField::New(
        "sumJ", he.mesh(), dimensionedScalar(dimMass / dimArea / dimTime, 0)));

    surfaceScalarField sumJh(surfaceScalarField::New(
        "sumJh",
        he.mesh(),
        dimensionedScalar(sumJ.dimensions() * he.dimensions(), 0)));

    forAll(Y, i)
    {
        if (i != this->thermo().defaultSpecie())
        {
            const volScalarField hi(this->thermo().hsi(i, p, T));

            const surfaceScalarField ji(this->j(Y[i]));
            sumJ += ji;

            sumJh += ji * fvc::interpolate(hi);
        }
    }

    {
        const label i = this->thermo().defaultSpecie();

        const volScalarField hi(this->thermo().hsi(i, p, T));

        sumJh -= sumJ * fvc::interpolate(hi);
    }

    tmpDivq.ref() += fvc::div(sumJh * he.mesh().magSf());

    return tmpDivq;
}


template <class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> TabulatedFickian<BasicThermophysicalTransportModel>::j(
    const volScalarField &Yi) const
{
    if (DTFuncs_.size())
    {
        const volScalarField &p = this->thermo().p();
        const volScalarField &T = this->thermo().T();

        return BasicThermophysicalTransportModel::j(Yi) -
               fvc::interpolate(
                   evaluate(DTFuncs_[this->thermo().specieIndex(Yi)],
                            dimDynamicViscosity,
                            p,
                            T)) *
                   fvc::snGrad(T) / fvc::interpolate(T);
    }
    else
    {
        return BasicThermophysicalTransportModel::j(Yi);
    }
}


template <class BasicThermophysicalTransportModel>
tmp<scalarField> TabulatedFickian<BasicThermophysicalTransportModel>::j(
    const volScalarField &Yi, const label patchi) const
{
    if (DTFuncs_.size())
    {
        const scalarField &p = this->thermo().p().boundaryField()[patchi];
        const scalarField &T = this->thermo().T().boundaryField()[patchi];

        return BasicThermophysicalTransportModel::j(Yi, patchi) -
               DTFuncs_[this->thermo().specieIndex(Yi)].value(p, T) *
                   this->thermo().T().boundaryField()[patchi].snGrad() / T;
    }
    else
    {
        return BasicThermophysicalTransportModel::j(Yi, patchi);
    }
}


template <class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> TabulatedFickian<BasicThermophysicalTransportModel>::divj(
    volScalarField &Yi) const
{
    if (DTFuncs_.size())
    {
        const volScalarField &p = this->thermo().p();
        const volScalarField &T = this->thermo().T();

        return BasicThermophysicalTransportModel::divj(Yi) -
               fvc::div(fvc::interpolate(
                            evaluate(DTFuncs_[this->thermo().specieIndex(Yi)],
                                     dimDynamicViscosity,
                                     p,
                                     T)) *
                        fvc::snGrad(T) / fvc::interpolate(T) *
                        T.mesh().magSf());
    }
    else
    {
        return BasicThermophysicalTransportModel::divj(Yi);
    }
}


template <class BasicThermophysicalTransportModel>
void TabulatedFickian<BasicThermophysicalTransportModel>::predict()
{
    BasicThermophysicalTransportModel::predict();

    updateDm();
}

template <class BasicThermophysicalTransportModel>
inline scalar TabulatedFickian<BasicThermophysicalTransportModel>::getDij(
    const FixedList<scalar, 5> &coeffs, const scalar &p, const scalar &T) const
{
    scalar logT = log(T);
    scalar sumupValue = 0.0;

    sumupValue =
        coeffs[0] +
        logT * (coeffs[1] +
                logT * (coeffs[2] + logT * (coeffs[3] + logT * coeffs[4])));
    scalar value = sumupValue * pow(T, 1.5) / p;
    return value;
}


template <class BasicThermophysicalTransportModel>
bool TabulatedFickian<BasicThermophysicalTransportModel>::movePoints()
{
    return true;
}


template <class BasicThermophysicalTransportModel>
void TabulatedFickian<BasicThermophysicalTransportModel>::topoChange(
    const polyTopoChangeMap &map)
{
    // Delete the cached Dm, will be re-created in predict
    Dm_.clear();
}


template <class BasicThermophysicalTransportModel>
void TabulatedFickian<BasicThermophysicalTransportModel>::mapMesh(
    const polyMeshMap &map)
{
    // Delete the cached Dm, will be re-created in predict
    Dm_.clear();
}


template <class BasicThermophysicalTransportModel>
void TabulatedFickian<BasicThermophysicalTransportModel>::distribute(
    const polyDistributionMap &map)
{
    // Delete the cached Dm, will be re-created in predict
    Dm_.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// //

}  // End namespace Foam

// ************************************************************************* //
