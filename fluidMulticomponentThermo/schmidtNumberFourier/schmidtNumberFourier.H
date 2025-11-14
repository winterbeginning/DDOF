#ifndef schmidtNumberFourier_H
#define schmidtNumberFourier_H

#include "laminarThermophysicalTransportModel.H"
#include "dimensionedScalar.H"

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

template <class BasicThermophysicalTransportModel>
class schmidtNumberFourier : public laminarThermophysicalTransportModel<
                                 BasicThermophysicalTransportModel>
{
    //- Schmidt number
    dimensionedScalar Sc_;

public:
    typedef typename BasicThermophysicalTransportModel::alphaField alphaField;

    typedef typename BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename BasicThermophysicalTransportModel::thermoModel thermoModel;


    //- Runtime type information
    TypeName("schmidtNumberFourier");


    // Constructors

    //- Construct from a momentum transport model and a thermo model
    schmidtNumberFourier(const momentumTransportModel &momentumTransport,
                         const thermoModel &thermo);


    //- Destructor
    virtual ~schmidtNumberFourier()
    {
    }


    // Member Functions

    //- Read thermophysicalTransport dictionary
    virtual bool read();

    //- Effective mass diffusion coefficient [kg/m/s]
    //  Calculated as DEff = mu/Sc
    virtual tmp<volScalarField> DEff(const volScalarField &Yi) const;

    //- Effective mass diffusion coefficient for patch [kg/m/s]
    virtual tmp<scalarField> DEff(const volScalarField &Yi,
                                  const label patchi) const;

    //- Return the heat flux [W/m^2]
    virtual tmp<surfaceScalarField> q() const;

    //- Return the heat flux for patch [W/m^2]
    virtual tmp<scalarField> q(const label patchi) const;

    //- Return the source term for the energy equation
    virtual tmp<fvScalarMatrix> divq(volScalarField &he) const;

    //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
    virtual tmp<surfaceScalarField> j(const volScalarField &Yi) const;

    //- Return the specie flux for patch [kg/m^2/s]
    virtual tmp<scalarField> j(const volScalarField &Yi,
                               const label patchi) const;

    //- Return the source term for the given specie mass-fraction equation
    virtual tmp<fvScalarMatrix> divj(volScalarField &Yi) const;
};

}  // End namespace laminarThermophysicalTransportModels
}  // End namespace Foam

#include "schmidtNumberFourier.C"

#endif