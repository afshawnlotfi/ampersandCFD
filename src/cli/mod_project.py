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


from src.models.settings import BoundingBox, Domain
from src.primitives import AmpersandDataInput, AmpersandIO, AmpersandUtils
from src.project import AmpersandProject
from src.utils.stl_analysis import StlAnalysis


# A collection of functions that are used to modify the project
class ModProject:
    @staticmethod
    def modify_project(project: AmpersandProject, modification_type: str):
        if modification_type == "Background Mesh":
            ModProject.change_background_mesh(project)
        elif modification_type == "Mesh Point":
            ModProject.change_mesh_point(project)
        elif modification_type == "Add Geometry":
            ModProject.add_geometry(project)
        elif modification_type == "Refinement Levels":
            ModProject.change_refinement_levels(project)
        elif modification_type == "Boundary Conditions":
            ModProject.change_boundary_conditions(project)
        elif modification_type == "Fluid Properties":
            ModProject.change_fluid_properties(project)
        elif modification_type == "Numerical Settings":
            ModProject.change_numerical_settings(project)
        elif modification_type == "Simulation Control Settings":
            ModProject.change_simulation_settings(project)
        elif modification_type == "Turbulence Model":
            ModProject.change_turbulenc_model(project)
        elif modification_type == "Post Processing Settings":
            ModProject.change_post_process_settings(project)
        
        raise ValueError("Invalid option. Aborting operation")



    @staticmethod
    # this is to change the global refinement level of the mesh
    def change_macro_refinement_level(project: AmpersandProject):
        refLevels = ["coarse", "medium", "fine"]
        AmpersandIO.printMessage(
            "Current refinement level: "+refLevels[project.settings.mesh.fineLevel])
        # ampersandIO.printMessage("Refinement level is the number of cells in the smallest direction")
        refinementLevel = AmpersandIO.get_input_int("Enter new refinement level (0:coarse, 1:medium, 2:fine): ")
        if (refinementLevel < 0 or refinementLevel > 2):
            AmpersandIO.printMessage("Invalid refinement level, please enter the value again")
            ModProject.change_refinement_levels(project)
        project.settings.mesh.fineLevel = refinementLevel
        # return project

    @staticmethod
    def change_domain_size(project: AmpersandProject, bounds: BoundingBox):
        AmpersandIO.printMessage(f"Domain size: {bounds.size} m")
        project.settings.mesh.domain = Domain.update(project.settings.mesh.domain, bounds)
        
    @staticmethod
    def change_mesh_size(project: AmpersandProject, cellSize: float):
        minX = project.settings.mesh.domain.minx
        maxX = project.settings.mesh.domain.maxx
        minY = project.settings.mesh.domain.miny
        maxY = project.settings.mesh.domain.maxy
        minZ = project.settings.mesh.domain.minz
        maxZ = project.settings.mesh.domain.maxz
        domain = (minX, maxX, minY, maxY, minZ, maxZ)
        nx, ny, nz = StlAnalysis.calc_nx_ny_nz(domain, cellSize)
        # check if the values are not too large
        if (nx > 500 or ny > 500 or nz > 500):
            AmpersandIO.printMessage(
                "Warning: Mesh is too fine. Consider increasing the cell size"
            )
        project.settings.mesh.domain.nx = nx
        project.settings.mesh.domain.ny = ny
        project.settings.mesh.domain.nz = nz
        # return project


    @staticmethod
    def change_stl_purpose(geometry, meshSettings):
        stlFile = geometry.file
        AmpersandIO.printMessage(f"Current STL file purpose: {geometry.purpose}")
        purpose = AmpersandIO.get_input("Enter new STL file purpose: ")
        geometry.purpose = purpose
        return geometry

    # this will allow the user to change the details of the stl file if necessary
    # TODO: fix this
    @staticmethod
    def change_stl_details(project, stl_file_number=0):
        project.list_stl_files()
        change_purpose = AmpersandIO.get_input("Change any STL files (y/N)?: ")
        if change_purpose.lower() != 'y':
            AmpersandIO.printMessage("No change in STL files properties")
            return 0
        stl_file_number = AmpersandIO.get_input(
            "Enter the number of the file to change purpose: ")
        try:
            stl_file_number = int(stl_file_number)
        except ValueError:
            AmpersandIO.printMessage("Invalid input. Please try again.")
            ModProject.change_stl_details()
            # return -1
        if stl_file_number < 0 or stl_file_number > len(project.stl_files):
            AmpersandIO.printMessage("Invalid input. Please try again.")
            ModProject.change_stl_details()

        stl_file = project.stl_files[stl_file_number]
        stl_name = stl_file.name
        purpose = AmpersandDataInput.get_purpose()
        # self.add_purpose_(stl_name,purpose)
        return 0

    # add purpose to the stl file. currently not used
    @staticmethod
    def add_purpose_(stl_files, stl_name, purpose='wall'):
        AmpersandIO.printMessage(f"Setting purpose of {stl_name} to")
        for stl in stl_files:
            if stl.name == stl_name:
                AmpersandIO.printMessage(
                    f"Setting purpose of {stl_name} to {purpose}")
                stl.purpose = purpose
                return stl_files
        AmpersandIO.printMessage(
            f"STL file {stl_name} not found in the project")
        return -1

    @staticmethod
    def change_stl_refinement_level(project: AmpersandProject, stl_name: str):
        project.change_stl_refinement_level(stl_name)

    # ---------------------------------------------------------------------#
    # The functions called when modifications are to be made project #

    @staticmethod
    def change_background_mesh(project: AmpersandProject):
        AmpersandIO.printMessage("Current background mesh")
        project.summarize_background_mesh()
        # ask whether to change domain size
        change_domain_size = AmpersandIO.get_input_bool(
            "Change domain size (y/N)?: ")
        # ask new domain size
        if change_domain_size:
            bounds = AmpersandDataInput.get_domain_size()
            ModProject.change_domain_size(project, bounds)
            AmpersandIO.printMessage("Domain size changed")
        # ask new cell size
        change_mesh_size = AmpersandIO.get_input_bool(
            "Change cell size (y/N)?: ")
        if change_mesh_size:
            cellSize = AmpersandDataInput.get_cell_size()
            project.settings.mesh.maxCellSize = cellSize
            # calculate new mesh size
            ModProject.change_mesh_size(project, cellSize)
            AmpersandIO.printMessage("Cell size changed")
        if change_domain_size or change_mesh_size:
            project.summarize_background_mesh()
        else:
            AmpersandIO.printMessage("No change in background mesh")

    # TODO: fix this function
    @staticmethod
    def add_geometry(project: AmpersandProject):
        AmpersandIO.printMessage("Adding geometry")

        yN = AmpersandIO.get_input("Add STL file to the project (y/N)?: ")
        while yN.lower() == 'y':
            stl_path = AmpersandIO.get_file( [("STL Geometry", "*.stl"), ("OBJ Geometry", "*.obj")])
            purpose = AmpersandDataInput.get_purpose()
            project.add_stl_file(stl_path, purpose)
            yN = AmpersandIO.get_input("Add another STL file to the project (y/N)?: ")
        AmpersandUtils.list_stl_files(project.stl_files)

    # TODO: fix this function
    @staticmethod
    def change_refinement_levels(project: AmpersandProject):
        AmpersandIO.printMessage("Changing refinement levels")
        # TODO: Implement this function
        AmpersandUtils.list_stl_files(project.stl_files)

        stl_file_number = AmpersandIO.get_input("Enter the number of the file to change refinement level: ")
        try:
            stl_file_number = int(stl_file_number)
        except ValueError:
            AmpersandIO.printMessage("Invalid input. Please try again.")
        if stl_file_number <= 0 or stl_file_number > len(project.stl_files):
            AmpersandIO.printMessage("Invalid input. Please try again.")
        else:
            ModProject.change_stl_refinement_level(project, stl_file_number-1)
        AmpersandUtils.list_stl_files(project.stl_files)
        return 0

    @staticmethod
    def change_mesh_point(project: AmpersandProject):
        AmpersandIO.printMessage("Changing mesh points")
        currentMeshPoint = project.settings.mesh.castellatedMeshControls.locationInMesh
        AmpersandIO.printMessage(
            f"Current mesh points: ({currentMeshPoint[0]},{currentMeshPoint[1]},{currentMeshPoint[2]})")

        x, y, z = AmpersandIO.get_input_vector("Enter new mesh points: ")
        project.settings.mesh.castellatedMeshControls.locationInMesh = (x, y, z)
        AmpersandIO.printMessage(
            f"New mesh points: ({currentMeshPoint[0]},{currentMeshPoint[1]},{currentMeshPoint[2]})")

    # TODO: fix this function
    @staticmethod
    def change_boundary_conditions(project: AmpersandProject):
        AmpersandIO.printMessage("Changing boundary conditions")
        # TODO: Implement this function
        bcs = project.summarize_boundary_conditions()
        # ampersandIO.printMessage("Current boundary conditions")
        # ampersandIO.printMessage(bcs)

        bc_number = AmpersandIO.get_input( "Enter the number of the boundary to change: ")
        try:
            bc_number = int(bc_number)
        except ValueError:
            AmpersandIO.printMessage("Invalid input. Please try again.")
        if bc_number <= 0 or bc_number > len(bcs):
            AmpersandIO.printMessage("Invalid input. Please try again.")
        else:
            bc = bcs[bc_number-1]
            AmpersandIO.printMessage(
                f"Changing boundary condition for patch: {bc}")
            newBcType = AmpersandDataInput.get_boundary_type()
            project.change_boundary_condition(bc, newBcType)

    @staticmethod
    def change_numerical_settings(project: AmpersandProject):
        AmpersandIO.printMessage("Changing numerical settings")
        # TODO: Implement this function

    @staticmethod
    def change_simulation_settings(project: AmpersandProject):
        AmpersandIO.printMessage("Changing simulation settings")
        # TODO: Implement this function

    @staticmethod
    def change_turbulenc_model(project: AmpersandProject):
        AmpersandIO.printMessage("Changing turbulence model")
        # TODO: Implement this function

    @staticmethod
    def change_post_process_settings(project: AmpersandProject):
        AmpersandIO.printMessage("Changing post process settings")
        # TODO: Implement this function

    @staticmethod
    def change_fluid_properties(project: AmpersandProject):
        AmpersandIO.printMessage("Current fluid properties")
        AmpersandIO.printMessage(f"Density: {project.settings.physicalProperties.rho}")
        AmpersandIO.printMessage(f"Kinematic viscosity: {project.settings.physicalProperties.nu}")
        rho = AmpersandIO.get_input_float("Enter new density (kg/m^3): ")
        nu = AmpersandIO.get_input_float(
            "Enter new kinematic viscosity (m^2/s): ")
        # check if the values are valid
        if (rho <= 0 or nu <= 0):
            AmpersandIO.printMessage(
                "Invalid fluid properties, please enter the values again")
            ModProject.change_fluid_properties(project)
        project.settings.physicalProperties.rho = rho
        project.settings.physicalProperties.nu = nu


