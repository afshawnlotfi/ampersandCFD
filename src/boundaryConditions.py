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

# This script generates the boundary conditions files for an OpenFOAM pimpleFoam simulation.
# The boundary conditions are specified in the meshSettings.yaml file.
# This is an early version of the script and will be updated in the future.
# Brute force writing is used instead of a more elegant solution.
from typing import Optional
from primitives import ampersandPrimitives
from constants import MeshSettings, BoundaryConditions, InletValues
from src.utils.headers import tuple_to_string
from src.utils.dict_generation import DictGenerationUtils


def write_vector_boundary_condition(patch: str = "inlet1", purpose: str = "inlet", property: Optional[tuple] = None) -> str:
    """
    Write a vector boundary condition 
    """
    bc = f"""{patch}
    {{"""
    # if the purpose is an inlet, then the velocity is specified
    if purpose == "inlet":
        # write the velocity
        bc += f"""
        type            fixedValue;
        value           uniform {tuple_to_string(property)};"""
    # if the purpose is an outlet, give an inletOutlet boundary condition
    elif purpose == "outlet":
        # write the pressure
        bc += f"""
        type            inletOutlet;
        inletValue      uniform (0 0 0);
        value           uniform (0 0 0);"""
    # if the purpose is a wall, give a fixedValue boundary condition
    elif purpose == "wall":
        bc += f"""
        type            fixedValue;
        value           uniform (0 0 0);"""
    # if the purpose is a symmetry, give a symmetry boundary condition
    elif purpose == "symmetry":
        bc += f"""
        type            symmetry;"""
    else:
        raise ValueError("Invalid boundary condition type")
    bc += f"""
    }}"""
    return bc


def write_turbulence_boundary_condition(
    patch: str = "inlet1",
    purpose: str = "inlet",
    property: Optional[tuple[float]] = None,
    wallFunction: str = "kqRWallFunction"
) -> str:
    """
    Write a scalar boundary condition
    """
    bc = f"""{patch}
    {{"""
    # if the purpose is an inlet, then the fixedValue is specified
    if purpose == "inlet":
        # write the velocity
        bc += f"""
        type            fixedValue;
        value           uniform {property};"""
    # if the purpose is an outlet, give an inletOutlet boundary condition
    elif purpose == "outlet":
        # write the pressure
        bc += f"""
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;"""
    # if the purpose is a wall, give a fixedValue boundary condition
    elif purpose == "wall":
        bc += f"""
        type            {wallFunction};
        value           $internalField;"""
    # if the purpose is a symmetry, give a symmetry boundary condition
    elif purpose == "symmetry":
        bc += f"""
        type            symmetry;"""
    else:
        raise ValueError("Invalid boundary condition type")
    bc += f"""
    }}"""
    return bc


def write_pressure_boundary_condition(
    patch: str = "inlet1",
    purpose: str = "inlet",
    property: float = 0.0
):
    """
    Write a scalar boundary condition
    """
    bc = f"""{patch}
    {{"""
    # if the purpose is an inlet, then the fixedValue is specified
    if purpose == "inlet":
        # write the velocity
        bc += f"""
        type            zeroGradient;"""
    # if the purpose is an outlet, give an inletOutlet boundary condition
    elif purpose == "outlet":
        # write the pressure
        bc += f"""
        type            fixedValue;
        value           uniform {property};"""  # to define reference pressure
    # if the purpose is a wall, give a fixedValue boundary condition
    elif purpose == "wall":
        bc += f"""
        type            zeroGradient;"""
    # if the purpose is a symmetry, give a symmetry boundary condition
    elif purpose == "symmetry":
        bc += f"""
        type            symmetry;"""
    else:
        raise ValueError("Invalid boundary condition type")
    bc += f"""
    }}"""
    return bc


def create_scalar_file(
    meshSettings: MeshSettings,
    boundaryConditions: BoundaryConditions,
    objName: str = "k",
    dimensions: tuple = (0, 2, -2)
) -> str:
    header = DictGenerationUtils.createFoamHeader(
        className="volScalarField", objectName=objName)
    dims = DictGenerationUtils.createDimensions(
        M=dimensions[0], L=dimensions[1], T=dimensions[2]
    )
    internalField = DictGenerationUtils.createInternalFieldScalar(
        type="uniform", value=0.0
    )
    s_file = f""+header+dims+internalField+"\n"+"""\nboundaryField 
{"""

    if not meshSettings.internalFlow:
        for patchName, bcPatch in meshSettings.bcPatches.items():
            if objName in ["k", "epsilon", "omega"]:
                s_file += write_turbulence_boundary_condition(
                    patch=patchName, purpose=bcPatch.purpose, property=bcPatch.property)
            elif objName == "p":
                s_file += write_pressure_boundary_condition(
                    patch=patchName, purpose=bcPatch.purpose, property=bcPatch.property)

    # If internal flow, set the boundary conditions for STL patches
    for patch_name, patch in meshSettings.geometry.items():
        if patch.type == 'triSurfaceMesh':
            if objName in ["k", "epsilon", "omega"]:
                s_file += write_turbulence_boundary_condition(
                    patch=patch_name, purpose=patch.purpose, property=patch.property)
            elif objName == "p":
                s_file += write_pressure_boundary_condition(
                    patch=patch_name, purpose=patch.purpose, property=patch.property)
    s_file += """
}"""
    return s_file


def create_u_file(meshSettings: MeshSettings, boundaryConditions: BoundaryConditions) -> str:
    header = DictGenerationUtils.createFoamHeader(
        className="volVectorField", objectName="U")
    dims = DictGenerationUtils.createDimensions(M=0, L=1, T=-1)
    internalField = DictGenerationUtils.createInternalFieldVector(
        type="uniform", value=boundaryConditions.velocityInlet.u_value)
    U_file = f""+header+dims+internalField+"\n"+"""\nboundaryField 
{"""

    if not meshSettings.internalFlow:
        for patch_patch, bcPatch in meshSettings.bcPatches.items():
            U_file += write_vector_boundary_condition(
                patch=patch_patch, purpose=bcPatch.purpose, property=bcPatch.property)

    # If internal flow, set the boundary conditions for STL patches
    for patch_patch, patch in meshSettings.geometry.items():
        if patch.type == 'triSurfaceMesh':
            U_file += write_vector_boundary_condition(
                patch=patch_patch, purpose=patch.purpose, property=patch.property)
    U_file += """
}"""
    return U_file


def create_p_file(meshSettings: MeshSettings, boundaryConditions: BoundaryConditions) -> str:
    p_file = create_scalar_file(
        meshSettings, boundaryConditions, objName="p", dimensions=(0, 2, -2))
    return p_file


def create_k_file(meshSettings: MeshSettings, boundaryConditions: BoundaryConditions) -> str:
    k_file = create_scalar_file(
        meshSettings, boundaryConditions, objName="k", dimensions=(0, 2, -2))
    return k_file


def create_epsilon_file(meshSettings: MeshSettings, boundaryConditions: BoundaryConditions) -> str:
    epsilon_file = create_scalar_file(
        meshSettings, boundaryConditions, objName="epsilon", dimensions=(0, 2, -2))
    return epsilon_file


def create_omega_file(meshSettings: MeshSettings, boundaryConditions: BoundaryConditions) -> str:
    omega_file = create_scalar_file(
        meshSettings, boundaryConditions, objName="omega", dimensions=(0, 2, -2))
    return omega_file


def create_nut_file(meshSettings: MeshSettings, boundaryConditions: BoundaryConditions) -> str:
    header = DictGenerationUtils.createFoamHeader(
        className="volScalarField", objectName="nut")
    dims = DictGenerationUtils.createDimensions(M=0, L=2, T=-1)
    internalField = DictGenerationUtils.createInternalFieldScalar(
        type="calculated", value=0.0)
    nut_file = f""+header+dims+internalField+"\n"+"""\nboundaryField 
{"""

    if not meshSettings.internalFlow:
        for patchName, bcPatch in meshSettings.bcPatches.items():
            nut_file += write_turbulence_boundary_condition(
                patch=patchName, purpose=bcPatch.purpose, property=bcPatch.property, wallFunction="nutkWallFunction")

    # If internal flow, set the boundary conditions for STL patches
    for patch_name, patch in meshSettings.geometry.items():
        if patch.type == 'triSurfaceMesh':
            nut_file += write_turbulence_boundary_condition(
                patch=patch_name, purpose=patch.purpose, property=patch.property, wallFunction="nutkWallFunction")
    nut_file += """
}"""
    return nut_file


def update_boundary_conditions(boundaryConditions: BoundaryConditions, inletValues: InletValues) -> BoundaryConditions:
    """
    Update boundary conditions with inlet values.

    Parameters:
    boundaryConditions (BoundaryConditions): Boundary conditions for U, p, k, and omega.
    inletValues (InletValues): Inlet values for U, p, k, and omega.
    """
    boundaryConditions.velocityInlet.u_value = inletValues.U
    boundaryConditions.velocityInlet.p_value = inletValues.p
    boundaryConditions.velocityInlet.k_value = inletValues.k
    boundaryConditions.velocityInlet.omega_value = inletValues.omega
    boundaryConditions.velocityInlet.epsilon_value = inletValues.epsilon
    boundaryConditions.velocityInlet.nut_value = inletValues.nut
    return boundaryConditions


def create_boundary_conditions(meshSettings: MeshSettings, boundaryConditions: BoundaryConditions, nu: float = 1.e-5) -> None:
    """
    Create boundary condition files for an OpenFOAM pimpleFoam simulation.

    Parameters:
    meshSettings (MeshSettings): Mesh settings.
    boundaryConditions (BoundaryConditions): Boundary conditions for U, p, k, and omega.
    """
    u_file = create_u_file(meshSettings, boundaryConditions)
    p_file = create_p_file(meshSettings, boundaryConditions)
    k_file = create_k_file(meshSettings, boundaryConditions)
    omega_file = create_omega_file(meshSettings, boundaryConditions)
    epsilon_file = create_epsilon_file(meshSettings, boundaryConditions)
    nut_file = create_nut_file(meshSettings, boundaryConditions)
    print("Creating boundary conditions files")
    ampersandPrimitives.write_to_file("U", u_file)
    ampersandPrimitives.write_to_file("p", p_file)
    ampersandPrimitives.write_to_file("k", k_file)
    ampersandPrimitives.write_to_file("omega", omega_file)
    ampersandPrimitives.write_to_file("epsilon", epsilon_file)
    ampersandPrimitives.write_to_file("nut", nut_file)
