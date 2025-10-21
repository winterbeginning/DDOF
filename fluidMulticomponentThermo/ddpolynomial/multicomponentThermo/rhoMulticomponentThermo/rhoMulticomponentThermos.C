/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
#include "rhoFluidMulticomponentThermo.H"

// --- 包含所有需要的组件头文件 ---
#include "coefficientMulticomponentMixture.H"
#include "coefficientWilkeMulticomponentMixture.H"
#include "valueMulticomponentMixture.H"
#include "singleComponentMixture.H"

// 包含核心的宏定义文件
#include "forGases.H"
#include "makeFluidMulticomponentThermo.H"

// --- 包含你的自定义输运模型 ---
#include "ddpolynomialTransport.H"

#include "makeChemistrySolver.H"
#include "chemistryModel.H"
#include "ode.H"

#include "forLiquids.H"

#include "makeChemistryReductionMethod.H"

#include "noChemistryReduction.H"
#include "DAC.H"
#include "DRG.H"
#include "DRGEP.H"
#include "EFA.H"
#include "PFA.H"

#include "makeReaction.H"

#include "ArrheniusReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "FallOffReactionRate.H"
#include "ChemicallyActivatedReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"

#include "MichaelisMentenReactionRate.H"
#include "LangmuirHinshelwoodReactionRate.H"
#include "fluxLimitedLangmuirHinshelwoodReactionRate.H"
#include "surfaceArrheniusReactionRate.H"

#include "formyDefine.H"

namespace Foam
{

formyCoeffGases(makeFluidMulticomponentThermos,
                rhoFluidThermo,
                rhoFluidMulticomponentThermo,
                coefficientMulticomponentMixture);
formyCoeffGases(makeFluidMulticomponentThermos,
                rhoFluidThermo,
                rhoFluidMulticomponentThermo,
                coefficientWilkeMulticomponentMixture);
formyGases(makeFluidMulticomponentThermo,
           rhoFluidMulticomponentThermo,
           singleComponentMixture);

formyCoeffLiquids(makeFluidMulticomponentThermos,
                  rhoFluidThermo,
                  rhoFluidMulticomponentThermo,
                  coefficientMulticomponentMixture);
formyLiquids(makeFluidMulticomponentThermos,
             rhoFluidThermo,
             rhoFluidMulticomponentThermo,
             valueMulticomponentMixture);
formyLiquids(makeFluidMulticomponentThermo,
             rhoFluidMulticomponentThermo,
             singleComponentMixture);

// --- 2. 注册化学求解器 (e.g., ode) ---
formyCoeffGases(defineChemistrySolvers, nullArg);
// formyGases(defineChemistrySolvers, nullArg);

formyCoeffGases(makeChemistrySolvers, ode);
// formyGases(makeChemistrySolvers, ode);

// --- 3. 注册化学简化方法 ---
formyCoeffGases(defineChemistryReductionMethod, nullArg);
formyCoeffGases(makeChemistryReductionMethod, none);
formyCoeffGases(makeChemistryReductionMethod, DAC);
formyCoeffGases(makeChemistryReductionMethod, DRG);
formyCoeffGases(makeChemistryReductionMethod, DRGEP);
formyCoeffGases(makeChemistryReductionMethod, EFA);
formyCoeffGases(makeChemistryReductionMethod, PFA);

// --- 4. 注册所有化学反应类型 ---
formyCoeffGases(defineReaction, nullArg);
formyCoeffLiquids(defineReaction, nullArg);

// Irreversible/reversible/non-equilibrium-reversible reactions
formyCoeffGases(makeIRNReactions, ArrheniusReactionRate);
formyCoeffLiquids(makeIRNReactions, ArrheniusReactionRate);
formyCoeffGases(makeIRNReactions, LandauTellerReactionRate);
formyCoeffLiquids(makeIRNReactions, LandauTellerReactionRate);
formyCoeffGases(makeIRNReactions, thirdBodyArrheniusReactionRate);
formyCoeffLiquids(makeIRNReactions, thirdBodyArrheniusReactionRate);

// Irreversible/reversible reactions
formyCoeffGases(makeIRReactions, JanevReactionRate);
formyCoeffLiquids(makeIRReactions, JanevReactionRate);
formyCoeffGases(makeIRReactions, powerSeriesReactionRate);
formyCoeffLiquids(makeIRReactions, powerSeriesReactionRate);

// Irreversible/reversible fall-off reactions
formyCoeffGases(makeIRTemplate2Reactions,
                FallOffReactionRate,
                ArrheniusReactionRate,
                LindemannFallOffFunction);
formyCoeffLiquids(makeIRTemplate2Reactions,
                  FallOffReactionRate,
                  ArrheniusReactionRate,
                  LindemannFallOffFunction);
formyCoeffGases(makeIRTemplate2Reactions,
                FallOffReactionRate,
                ArrheniusReactionRate,
                TroeFallOffFunction);
formyCoeffLiquids(makeIRTemplate2Reactions,
                  FallOffReactionRate,
                  ArrheniusReactionRate,
                  TroeFallOffFunction);
formyCoeffGases(makeIRTemplate2Reactions,
                FallOffReactionRate,
                ArrheniusReactionRate,
                SRIFallOffFunction);
formyCoeffLiquids(makeIRTemplate2Reactions,
                  FallOffReactionRate,
                  ArrheniusReactionRate,
                  SRIFallOffFunction);

// Irreversible/reversible chemically activated reactions
formyCoeffGases(makeIRTemplate2Reactions,
                ChemicallyActivatedReactionRate,
                ArrheniusReactionRate,
                LindemannFallOffFunction);
formyCoeffLiquids(makeIRTemplate2Reactions,
                  ChemicallyActivatedReactionRate,
                  ArrheniusReactionRate,
                  LindemannFallOffFunction);
formyCoeffGases(makeIRTemplate2Reactions,
                ChemicallyActivatedReactionRate,
                ArrheniusReactionRate,
                TroeFallOffFunction);
formyCoeffLiquids(makeIRTemplate2Reactions,
                  ChemicallyActivatedReactionRate,
                  ArrheniusReactionRate,
                  TroeFallOffFunction);
formyCoeffGases(makeIRTemplate2Reactions,
                ChemicallyActivatedReactionRate,
                ArrheniusReactionRate,
                SRIFallOffFunction);
formyCoeffLiquids(makeIRTemplate2Reactions,
                  ChemicallyActivatedReactionRate,
                  ArrheniusReactionRate,
                  SRIFallOffFunction);

// Michaelis-Menten Reactions
formyCoeffLiquids(makeIReactions, MichaelisMentenReactionRate);

// Langmuir-Hinshelwood Reactions
formyCoeffGases(makeIRReactions, LangmuirHinshelwoodReactionRate);
formyCoeffLiquids(makeIRReactions, LangmuirHinshelwoodReactionRate);

// Flux-limited Langmuir-Hinshelwood Reactions
formyCoeffGases(makeGeneralReaction,
                IrreversibleReaction,
                fluxLimitedLangmuirHinshelwoodReactionRate);

// Surface-Arrhenius Reactions
formyCoeffGases(makeGeneralReaction,
                IrreversibleReaction,
                surfaceArrheniusReactionRate);
formyCoeffLiquids(makeGeneralReaction,
                  IrreversibleReaction,
                  surfaceArrheniusReactionRate);

}  // namespace Foam

// ************************************************************************* //
