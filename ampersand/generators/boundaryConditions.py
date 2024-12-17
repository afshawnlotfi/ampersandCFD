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
from primitives import AmpersandUtils
from ampersand.models.settings import Patch, SnappyHexMeshSettings, BoundaryConditions, InletValues, Patch, TriSurfaceMeshGeometry
from ampersand.utils.generation import GenerationUtils


def write_vector_boundary_condition(patch: Patch, patch_name) -> str:
    """
    Write a vector boundary condition 
    """
    assert isinstance(patch.property, (tuple, list)), f"Property '{
        patch_name}' must be a tuple or list"
    bc = f"""{patch_name}
    {{"""
    # if the purpose is an inlet, then the velocity is specified
    if patch.purpose == "inlet":
        # write the velocity
        bc += f"""
        type            fixedValue;
        value           uniform ({patch.property[0]} {patch.property[1]} {patch.property[2]});"""
    # if the purpose is an outlet, give an inletOutlet boundary condition
    elif patch.purpose == "outlet":
        # write the pressure
        bc += f"""
        type            inletOutlet;
        inletValue      uniform (0 0 0);
        value           uniform (0 0 0);"""
    # if the purpose is a wall, give a fixedValue boundary condition
    elif patch.purpose == "wall":
        bc += f"""
        type            fixedValue;
        value           uniform (0 0 0);"""
    # if the purpose is a symmetry, give a symmetry boundary condition
    elif patch.purpose == "symmetry":
        bc += f"""
        type            symmetry;"""
    else:
        raise ValueError("Invalid boundary condition type")
    bc += f"""
    }}"""
    return bc


def create_turbulence_boundary_condition(patch: Patch, patch_name: str, wallFunction: str = "kqRWallFunction") -> str:
    """
    Write a scalar boundary condition
    """
    bc = f"""{patch_name}
    {{"""
    # if the purpose is an inlet, then the fixedValue is specified
    if patch.purpose == "inlet":
        # write the velocity
        bc += f"""
        type            fixedValue;
        value           uniform {patch.property};"""
    # if the purpose is an outlet, give an inletOutlet boundary condition
    elif patch.purpose == "outlet":
        # write the pressure
        bc += f"""
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;"""
    # if the purpose is a wall, give a fixedValue boundary condition
    elif patch.purpose == "wall":
        bc += f"""
        type            {wallFunction};
        value           $internalField;"""
    # if the purpose is a symmetry, give a symmetry boundary condition
    elif patch.purpose == "symmetry":
        bc += f"""
        type            symmetry;"""
    else:
        raise ValueError("Invalid boundary condition type")
    bc += f"""
    }}"""
    return bc


def create_pressure_boundary_condition(patch: Patch, patch_name: str):
    """
    Write a scalar boundary condition
    """
    bc = f"""{patch_name}
    {{"""
    # if the purpose is an inlet, then the fixedValue is specified
    if patch.purpose == "inlet":
        # write the velocity
        bc += f"""
        type            zeroGradient;"""
    # if the purpose is an outlet, give an inletOutlet boundary condition
    elif patch.purpose == "outlet":
        # write the pressure
        bc += f"""
        type            fixedValue;
        value           uniform {property};"""  # to define reference pressure
    # if the purpose is a wall, give a fixedValue boundary condition
    elif patch.purpose == "wall":
        bc += f"""
        type            zeroGradient;"""
    # if the purpose is a symmetry, give a symmetry boundary condition
    elif patch.purpose == "symmetry":
        bc += f"""
        type            symmetry;"""
    else:
        raise ValueError("Invalid boundary condition type")
    bc += f"""
    }}"""
    return bc


def create_scalar_boundary_condition(patch: Patch, patch_name: str, objName: str):
    if objName in ["k", "epsilon", "omega"]:
        return create_turbulence_boundary_condition(patch, patch_name)
    elif objName == "p":
        return create_pressure_boundary_condition(patch, patch_name)
    raise ValueError("Invalid scalar boundary condition")

def create_scalar_file(
    meshSettings: SnappyHexMeshSettings,
    objName: str = "k",
    dimensions: tuple = (0, 2, -2),
    value: float = 0.0,
) -> str:
    header = GenerationUtils.createFoamHeader("volScalarField", objName)
    dims = GenerationUtils.createDimensions(*dimensions)
    internalField = GenerationUtils.createInternalFieldScalar("uniform", value)
    s_file = f"\n{header}{dims}{internalField}\nboundaryField\n{'{'}"

    if not meshSettings.internalFlow:
        for patch_name, patch in meshSettings.patches.items():
            assert isinstance(patch.property, float), f"Property '{patch_name}' must be a float"
            s_file += create_scalar_boundary_condition(patch, patch_name, objName)

    # If internal flow, set the boundary conditions for STL patches
    for geometry_name, geometry in meshSettings.geometry.items():
        if isinstance(geometry, TriSurfaceMeshGeometry):
            s_file += create_scalar_boundary_condition(geometry, geometry_name, objName)

    s_file += "\n}"
    return s_file


def create_u_file(meshSettings: SnappyHexMeshSettings, boundaryConditions: BoundaryConditions) -> str:
    header = GenerationUtils.createFoamHeader(
        className="volVectorField", objectName="U")
    dims = GenerationUtils.createDimensions(M=0, L=1, T=-1)
    internalField = GenerationUtils.createInternalFieldVector(
        type="uniform", value=boundaryConditions.velocityInlet.u_value)
    U_file = f""+header+dims+internalField+"\n"+"""\nboundaryField 
{"""

    if not meshSettings.internalFlow:
        for patch_name, patch in meshSettings.patches.items():
            assert isinstance(patch.property, (tuple, list)), f"Property '{patch_name}' must be a float"

            U_file += write_vector_boundary_condition(patch, patch_name)

    # If internal flow, set the boundary conditions for STL patches
    for geometry_name, geometry in meshSettings.geometry.items():
        if isinstance(geometry, TriSurfaceMeshGeometry):
            U_file += write_vector_boundary_condition(geometry, geometry_name)
    U_file += """
}"""
    return U_file


def create_p_file(meshSettings: SnappyHexMeshSettings) -> str:
    return create_scalar_file(
        meshSettings, objName="p", dimensions=(0, 2, -2)
    )


def create_k_file(meshSettings: SnappyHexMeshSettings) -> str:
    return create_scalar_file(
        meshSettings, objName="k", dimensions=(0, 2, -2)
    )


def create_epsilon_file(meshSettings: SnappyHexMeshSettings) -> str:
    return create_scalar_file(
        meshSettings, objName="epsilon", dimensions=(0, 2, -2)
    )


def create_omega_file(meshSettings: SnappyHexMeshSettings) -> str:
    return create_scalar_file(
        meshSettings, objName="omega", dimensions=(0, 2, -2)
    )


def create_nut_file(meshSettings: SnappyHexMeshSettings) -> str:
    header = GenerationUtils.createFoamHeader(
        className="volScalarField", objectName="nut")
    dims = GenerationUtils.createDimensions(M=0, L=2, T=-1)
    internalField = GenerationUtils.createInternalFieldScalar(
        type="calculated", value=0.0)
    nut_file = f"\n{header}{dims}{internalField}\nboundaryField\n{'{'}"
    if not meshSettings.internalFlow:
        for patch_name, patch in meshSettings.patches.items():
            assert isinstance(patch.property, float), f"Property '{patch_name}' must be a float"
            nut_file += create_turbulence_boundary_condition(patch, patch_name, wallFunction="nutkWallFunction")

    # If internal flow, set the boundary conditions for STL patches
    for geometry_name, geometry in meshSettings.geometry.items():
        if isinstance(geometry, TriSurfaceMeshGeometry):
            nut_file += create_turbulence_boundary_condition(geometry, geometry_name, wallFunction="nutkWallFunction")
    nut_file += "\n}"
    return nut_file


# def update_boundary_conditions(boundaryConditions: BoundaryConditions, inletValues: InletValues) -> BoundaryConditions:
#     """
#     Update boundary conditions with inlet values.

#     Parameters:
#     boundaryConditions (BoundaryConditions): Boundary conditions for U, p, k, and omega.
#     inletValues (InletValues): Inlet values for U, p, k, and omega.
#     """
#     boundaryConditions.velocityInlet.u_value = inletValues.U
#     boundaryConditions.velocityInlet.p_value = inletValues.p
#     boundaryConditions.velocityInlet.k_value = inletValues.k
#     boundaryConditions.velocityInlet.omega_value = inletValues.omega
#     boundaryConditions.velocityInlet.epsilon_value = inletValues.epsilon
#     boundaryConditions.velocityInlet.nut_value = inletValues.nut
#     return boundaryConditions


def create_boundary_conditions(meshSettings: SnappyHexMeshSettings, boundaryConditions: BoundaryConditions, nu: float = 1.e-5) -> None:
    """
    Create boundary condition files for an OpenFOAM pimpleFoam simulation.

    Parameters:
    meshSettings (MeshSettings): Mesh settings.
    boundaryConditions (BoundaryConditions): Boundary conditions for U, p, k, and omega.
    """
    u_file = create_u_file(meshSettings, boundaryConditions)
    p_file = create_p_file(meshSettings)
    k_file = create_k_file(meshSettings)
    omega_file = create_omega_file(meshSettings)
    epsilon_file = create_epsilon_file(meshSettings)
    nut_file = create_nut_file(meshSettings)
    print("Creating boundary conditions files")
    AmpersandUtils.write_to_file("U", u_file)
    AmpersandUtils.write_to_file("p", p_file)
    AmpersandUtils.write_to_file("k", k_file)
    AmpersandUtils.write_to_file("omega", omega_file)
    AmpersandUtils.write_to_file("epsilon", epsilon_file)
    AmpersandUtils.write_to_file("nut", nut_file)
