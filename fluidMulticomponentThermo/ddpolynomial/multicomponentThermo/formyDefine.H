#ifndef formyDefine_H
#define formyDefine_H

#include "forGases.H"
#include "forLiquids.H"
#include "ddpolynomialTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define formyGasTransports(Macro, Args...) \
    forGasEnergiesAndThermos(ddpolynomialTransport, Macro, Args)

#define formyGases(Macro, Args...) formyGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define formyCoeffGasTransports(Macro, Args...) \
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, Macro, Args)

#define formyCoeffGases(Macro, Args...) formyCoeffGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define formyLiquidTransports(Macro, Args...) \
    forLiquidEnergiesAndThermos(ddpolynomialTransport, Macro, Args);

#define formyLiquids(Macro, Args...) formyLiquidTransports(Macro, Args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define formyCoeffLiquidTransports(Macro, Args...) \
    forCoeffLiquidEnergiesAndThermos(ddpolynomialTransport, Macro, Args);

#define formyCoeffLiquids(Macro, Args...) \
    formyCoeffLiquidTransports(Macro, Args)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//


#endif