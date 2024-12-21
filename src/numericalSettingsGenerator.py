"""
-------------------------------------------------------------------------------
  ***    *     *  ******   *******  ******    *****     ***    *     *  ******   
 *   *   **   **  *     *  *        *     *  *     *   *   *   **    *  *     *  
*     *  * * * *  *     *  *        *     *  *        *     *  * *   *  *     *  
*******  *  *  *  ******   ****     ******    *****   *******  *  *  *  *     *  
*     *  *     *  *        *        *   *          *  *     *  *   * *  *     *  
*     *  *     *  *        *        *    *   *     *  *     *  *    **  *     *  
*     *  *     *  *        *******  *     *   *****   *     *  *     *  ******   
-------------------------------------------------------------------------------
 * AmpersandCFD is a minimalist streamlined OpenFOAM generation tool.
 * Copyright (c) 2024 THAW TAR
 * All rights reserved.
 *
 * This software is licensed under the GNU General Public License version 3 (GPL-3.0).
 * You may obtain a copy of the license at https://www.gnu.org/licenses/gpl-3.0.en.html
 */
"""

from src.models.settings import NumericalSettings, SolverSettings
from src.utils.generation import GenerationUtils


def create_algorithmDict(numericalSettings: NumericalSettings):
    algorithmDict = f"""
PIMPLE
{{
    nOuterCorrectors {numericalSettings.pimpleDict.nOuterCorrectors};
    nCorrectors {numericalSettings.pimpleDict.nCorrectors};
    nNonOrthogonalCorrectors {numericalSettings.pimpleDict.nNonOrthogonalCorrectors};
    pRefCell {numericalSettings.pimpleDict.pRefCell};
    pRefValue {numericalSettings.pimpleDict.pRefValue};
    residualControl
    {{
        "(U|k|omega|epsilon|nut)"
        {{
            tolerance {numericalSettings.pimpleDict.residualControl.U};
            relTol 0;
        }}
        p
        {{
            tolerance {numericalSettings.pimpleDict.residualControl.p};
            relTol 0;
        }}
    }}
}}
SIMPLE
{{
    nNonOrthogonalCorrectors {numericalSettings.simpleDict.nNonOrthogonalCorrectors};
    consistent {numericalSettings.simpleDict.consistent};
    residualControl
    {{
        U {numericalSettings.simpleDict.residualControl.U};
        p {numericalSettings.simpleDict.residualControl.p};
        k {numericalSettings.simpleDict.residualControl.k};
        omega {numericalSettings.simpleDict.residualControl.omega};
        epsilon {numericalSettings.simpleDict.residualControl.epsilon};
        nut {numericalSettings.simpleDict.residualControl.nut};
    }}
}}
potentialFlow
{{
    nNonOrthogonalCorrectors {numericalSettings.potentialFlowDict.nonOrthogonalCorrectors};
}}
relaxationFactors
{{
    equations
    {{
        U {numericalSettings.relaxationFactors.U};

        k {numericalSettings.relaxationFactors.k};
        omega {numericalSettings.relaxationFactors.omega};
        epsilon {numericalSettings.relaxationFactors.epsilon};
        nut {numericalSettings.relaxationFactors.nut};
    }}
    fields
    {{
        p {numericalSettings.relaxationFactors.p};
    }}
}}
"""
    return algorithmDict


def create_solverDict(solverSettings: SolverSettings, solverName="U"):
    solverCfg = getattr(solverSettings, solverName)
    solverDict = f""
    if (solverName == "p" or solverName == "Phi" or solverName == "p_rgh"):
        solverDict += f"""
{solverName}
{{
    solver {solverCfg["type"]};
    smoother {solverCfg["smoother"]};
    agglomerator {solverCfg["agglomerator"]};
    nCellsInCoarsestLevel {solverCfg["nCellsInCoarsestLevel"]};
    mergeLevels {solverCfg["mergeLevels"]};
    cacheAgglomeration {solverCfg["cacheAgglomeration"]};
    tolerance {solverCfg["tolerance"]};
    relTol {solverCfg["relTol"]};
    maxIter {solverCfg["maxIter"]};
    nSweeps {solverCfg["nSweeps"]};
    nPreSweeps {solverCfg["nPreSweeps"]};
}}
"""
    else:
        solverDict += f"""
{solverName}
{{
    solver {solverCfg["type"]};
    smoother {solverCfg["smoother"]};
    tolerance {solverCfg["tolerance"]};
    relTol {solverCfg["relTol"]};
    maxIter 100;
}}
"""
    return solverDict


def create_solverFinalDict(solverSettings: SolverSettings, solverName="U"):
    solverDict = f"""
    {solverName}Final
    {{
        ${solverName}
        tolerance {getattr(solverSettings, solverName)["tolerance"]/100.};
        relTol 0;
    }}
    """
    return solverDict


def create_fvSolutionDict(numericalSettings: NumericalSettings, solverSettings: SolverSettings):
    fvSolutionDict = f"""{GenerationUtils.createFoamHeader("dictionary", "fvSolution")}
solvers
{{
    """
    for solver in solverSettings.model_fields.keys():
        fvSolutionDict += create_solverDict(solverSettings, solver)
        fvSolutionDict += create_solverFinalDict(solverSettings, solver)
    fvSolutionDict += f"""
}}
    """
    fvSolutionDict += create_algorithmDict(numericalSettings)
    return fvSolutionDict


def create_fvSchemesDict(numericalSettings: NumericalSettings):
    fvSchemesDict = f"""{GenerationUtils.createFoamHeader("dictionary", "fvSchemes")}
ddtSchemes
{{
    default {numericalSettings.ddtSchemes.default};
}}
gradSchemes
{{
    default {numericalSettings.gradSchemes.default};
    grad(p) {numericalSettings.gradSchemes.grad_p};
    grad(U) {numericalSettings.gradSchemes.grad_U};
}}
divSchemes
{{
    default {numericalSettings.divSchemes.default};
    div(phi,U) {numericalSettings.divSchemes.div_phi_U};
    div(phi,k) {numericalSettings.divSchemes.div_phi_k};
    div(phi,omega) {numericalSettings.divSchemes.div_phi_omega};
    div(phi,epsilon) {numericalSettings.divSchemes.div_phi_epsilon};
    div(phi,nut) {numericalSettings.divSchemes.div_phi_nut};
    div(nuEff*dev(T(grad(U)))) {numericalSettings.divSchemes.div_nuEff_dev_T_grad_U};
}}
laplacianSchemes
{{
    default {numericalSettings.laplacianSchemes.default};
}}
interpolationSchemes
{{
    default {numericalSettings.interpolationSchemes.default};
}}
snGradSchemes
{{
    default {numericalSettings.snGradSchemes.default};
}}
fluxRequired
{{
    default {numericalSettings.fluxRequired.default};
}}
wallDist
{{
    method {numericalSettings.wallDist};
}}
"""
    return fvSchemesDict
