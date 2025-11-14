#include "schmidtNumberFourier.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

template <class BasicThermophysicalTransportModel>
schmidtNumberFourier<BasicThermophysicalTransportModel>::schmidtNumberFourier(
    const momentumTransportModel &momentumTransport,
    const thermoModel &thermo) :
    laminarThermophysicalTransportModel<BasicThermophysicalTransportModel>(
        typeName, momentumTransport, thermo),
    Sc_("Sc", dimless, 0.0)
{
    if (!this->read())
    {
        FatalErrorInFunction << "Cannot read " << this->type()
                             << " transport properties" << exit(FatalError);
    }
}

template <class BasicThermophysicalTransportModel>
bool schmidtNumberFourier<BasicThermophysicalTransportModel>::read()
{
    if (laminarThermophysicalTransportModel<
            BasicThermophysicalTransportModel>::read())
    {
        this->coeffDict().lookup("Sc") >> Sc_;
        Info << "    Reading " << this->type() << " properties:" << nl
             << "        Sc = " << Sc_ << nl << endl;
        return true;
    }

    return false;
}

template <class BasicThermophysicalTransportModel>
tmp<volScalarField>
schmidtNumberFourier<BasicThermophysicalTransportModel>::DEff(
    const volScalarField &Yi) const
{
    return volScalarField::New("DEff", this->thermo().mu() / Sc_);
}

template <class BasicThermophysicalTransportModel>
tmp<scalarField> schmidtNumberFourier<BasicThermophysicalTransportModel>::DEff(
    const volScalarField &Yi, const label patchi) const
{
    return this->thermo().mu().boundaryField()[patchi] / Sc_.value();
}

template <class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
schmidtNumberFourier<BasicThermophysicalTransportModel>::q() const
{
    // Heat flux calculation remains the same (Fourier's Law)
    return surfaceScalarField::New(
        IOobject::groupName("q",
                            this->momentumTransport().alphaRhoPhi().group()),
        -fvc::interpolate(this->alpha() * this->thermo().kappa() /
                          this->thermo().Cpv()) *
            fvc::snGrad(this->thermo().he()));
}


template <class BasicThermophysicalTransportModel>
tmp<scalarField> schmidtNumberFourier<BasicThermophysicalTransportModel>::q(
    const label patchi) const
{
    return -(this->alpha().boundaryField()[patchi] *
             this->thermo().kappa().boundaryField()[patchi] /
             this->thermo().Cpv()[patchi] *
             this->thermo().he().boundaryField()[patchi].snGrad());
}

template <class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
schmidtNumberFourier<BasicThermophysicalTransportModel>::divq(
    volScalarField &he) const
{
    volScalarField alphahe("alphahe",
                           this->thermo().kappa() / this->thermo().Cpv());

    return -fvm::laplacian(this->alpha() * alphahe, he);
}

template <class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
schmidtNumberFourier<BasicThermophysicalTransportModel>::j(
    const volScalarField &Yi) const
{
    return surfaceScalarField::New(
        IOobject::groupName("j(" + Yi.name() + ')',
                            this->momentumTransport().alphaRhoPhi().group()),
        -fvc::interpolate(this->alpha() * this->DEff(Yi)) * fvc::snGrad(Yi));
}

template <class BasicThermophysicalTransportModel>
tmp<scalarField> schmidtNumberFourier<BasicThermophysicalTransportModel>::j(
    const volScalarField &Yi, const label patchi) const
{
    return -(this->alpha().boundaryField()[patchi] * this->DEff(Yi, patchi) *
             Yi.boundaryField()[patchi].snGrad());
}

template <class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
schmidtNumberFourier<BasicThermophysicalTransportModel>::divj(
    volScalarField &Yi) const
{
    return -fvm::laplacian(this->alpha() * this->DEff(Yi), Yi);
}

}  // End namespace laminarThermophysicalTransportModels
}  // End namespace Foam