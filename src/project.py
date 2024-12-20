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

# backend module for the ampersandCFD project
# Description: This file contains the code for managing project structure and
# generate OpenFOAM files


import os
from pathlib import Path
import shutil
from typing import Optional, Union
from src.primitives import AmpersandPrimitives, AmpersandIO, AmpersandDataInput
from src.constants import meshSettings, physicalProperties, numericalSettings, inletValues
from src.constants import solverSettings, boundaryConditions, simulationSettings, simulationFlowSettings, parallelSettings, postProcessSettings
from src.utils.stl_analysis import StlAnalysis
from src.blockMeshGenerator import generate_blockMeshDict
from src.decomposeParGenerator import createDecomposeParDict
from src.snappyHexMeshGenerator import generate_snappyHexMeshDict
from src.surfaceExtractor import create_surfaceFeatureExtractDict
from src.transportAndTurbulence import create_transportPropertiesDict, create_turbulencePropertiesDict
from src.boundaryConditionsGenerator import create_boundary_conditions
from src.controlDictGenerator import createControlDict
from src.numericalSettingsGenerator import create_fvSchemesDict, create_fvSolutionDict
from src.scriptGenerator import ScriptGenerator
from src.postProcess import postProcess


# from ../constants/constants import meshSettings


class AmpersandProject:  # ampersandProject class to handle the project creation and manipulation
    # this class will contain the methods to handle the logic and program flow
    def __init__(self, project_path: Union[str, Path]):

        self.project_path = Path(project_path)
        self.current_stl_file = None   # current stl file being processed
        self.project_directory_path = None
        self.user_name = None
        self.caseSettings = None
        self.meshSettings = None
        self.physicalProperties = None
        self.numericalSettings = None
        self.simulationSettings = None
        self.solverSettings = None
        self.inletValues = None
        self.boundaryConditions = None
        self.simulationSettings = None
        self.simulationFlowSettings = None
        self.parallelSettings = None
        self.postProcessSettings = None
        self.stl_files = []  # list to store the settings for stl files
        self.stl_names = []  # list to store the names of the stl files
        self.internalFlow = False  # default is external flow
        self.on_ground = False  # default is off the ground
        self.halfModel = False  # default is full model
        self.transient = False  # default is steady state
        self.refinement = 0  # 0: coarse, 1: medium, 2: fine
        self.useFOs = False  # default is not to use function objects
        self.current_modification = None  # current modification to the project settings
        # flag to check if the current working directory is the project directory
        self.minX, self.maxX, self.minY, self.maxY, self.minZ, self.maxZ = -1e-3, 1e-3, -1e-3, 1e-3, -1e-3, 1e-3
    # --------------------------------------------------------------------
    # Methods to handle the project summary and changes
    # --------------------------------------------------------------------

    def summarize_boundary_conditions(self):
        bcs = AmpersandPrimitives.list_boundary_conditions(self.meshSettings)
        return bcs



    def summarize_project(self):
        trueFalse = {True: 'Yes', False: 'No'}
        AmpersandIO.show_title("Project Summary")

        AmpersandIO.printMessage(
            f"Internal Flow: {trueFalse[self.internalFlow]}")
        if (self.internalFlow == False):
            AmpersandIO.printMessage(
                f"On Ground: {trueFalse[self.on_ground]}")
        AmpersandIO.printMessage(
            f"Transient: {trueFalse[self.transient]}")
        self.summarize_background_mesh()
        AmpersandPrimitives.list_stl_files(self.stl_files)
        

    # this will show the details of the background mesh

    def summarize_background_mesh(self):
        minX = self.meshSettings['domain']["minx"]
        maxX = self.meshSettings['domain']["maxx"]
        minY = self.meshSettings['domain']["miny"]
        maxY = self.meshSettings['domain']["maxy"]
        minZ = self.meshSettings['domain']["minz"]
        maxZ = self.meshSettings['domain']["maxz"]
        nx = self.meshSettings['domain']['nx']
        ny = self.meshSettings['domain']['ny']
        nz = self.meshSettings['domain']['nz']
        AmpersandIO.printMessage(f"Domain size:{'X':>10}{'Y':>10}{'Z':>10}")
        AmpersandIO.printMessage(
            f"Min         {minX:>10.3f}{minY:>10.3f}{minZ:>10.3f}")
        AmpersandIO.printMessage(
            f"Max         {maxX:>10.3f}{maxY:>10.3f}{maxZ:>10.3f}")
        AmpersandIO.printMessage(f"Background mesh size: {nx}x{ny}x{nz} cells")
        AmpersandIO.printMessage(
            f"Background cell size: {self.meshSettings['maxCellSize']} m")

    def change_boundary_condition(self, bcName, newBC):
        if not self.internalFlow:  # if it is external flow
            bcPatches = self.meshSettings['patches']
            for aPatch in self.meshSettings['patches']:
                if bcName == aPatch['name']:
                    aPatch['type'] = newBC
                    AmpersandIO.printMessage(
                        f"Boundary condition {bcName} changed to {newBC}")
                    return 0
            if bcName in bcPatches:
                self.meshSettings['patches'][bcName]['type'] = newBC
                self.meshSettings['bcPatches'][bcName]['purpose'] = newBC
                newProperty = self.set_property(newBC)
                self.meshSettings['bcPatches'][bcName]['property'] = newProperty
                AmpersandIO.printMessage(
                    f"Boundary condition {bcName} changed to {newBC}")
                return 0
            else:
                AmpersandIO.printMessage(
                    "Boundary condition not found in the list")
        geometry = self.meshSettings['geometry']
        for stl in geometry:
            if stl['name'] == bcName:
                stl['purpose'] = newBC
                newProperty = self.set_property(newBC)
                self.meshSettings['bcPatches'][bcName]['property'] = newProperty
                AmpersandIO.printMessage(
                    f"Boundary condition of {bcName} changed to {newBC}")
                return 0
        return -1

    def change_stl_refinement_level(self, stl_file_number=0):
        AmpersandIO.printMessage("Changing refinement level")
        refMin = AmpersandIO.get_input_int("Enter new refMin: ")
        refMax = AmpersandIO.get_input_int("Enter new refMax: ")
        self.stl_files[stl_file_number]['refineMin'] = refMin
        self.stl_files[stl_file_number]['refineMax'] = refMax
        # stl_name = project.stl_files[stl_file_number]['name']

        fileFound = False
        for stl in self.meshSettings['geometry']:
            if stl['name'] == self.stl_files[stl_file_number]['name']:
                fileFound = True
                stl['refineMin'] = refMin
                stl['refineMax'] = refMax
                stl['featureLevel'] = refMax
                break
        if not fileFound:
            AmpersandIO.printMessage("STL file not found in the geometry list")

        # return project


    def choose_modification_categorized(self):
        options = ['Mesh', 'Boundary Conditions', 'Fluid Properties', 'Numerical Settings',
                   'Simulation Control Settings', 'Turbulence Model', 'Post Processing Settings']
        current_modification = AmpersandIO.get_option_choice(prompt="Choose any option for project modification: ",
                                                             options=options, title="\nModify Project Settings")
        mesh_options = ['Background Mesh', 'Mesh Point',
                        'Add Geometry', 'Refinement Levels']

        if current_modification < 0 or current_modification > len(options)-1:
            AmpersandIO.printMessage("Invalid option. Aborting operation")
            return -1
        if current_modification == 0:
            self.current_modification = mesh_options[AmpersandIO.get_option_choice(prompt="Choose any option for mesh modification: ",
                                                                                   options=mesh_options, title="\nModify Mesh Settings")]
        else:
            self.current_modification = options[current_modification]



    # write current settings to the project_settings.yaml file inside the project directory
    def write_settings(self):
        settings = {
            'meshSettings': self.meshSettings,
            'physicalProperties': self.physicalProperties,
            'numericalSettings': self.numericalSettings,
            'inletValues': self.inletValues,
            'boundaryConditions': self.boundaryConditions,
            'solverSettings': self.solverSettings,
            'simulationSettings': self.simulationSettings,
            'parallelSettings': self.parallelSettings,
            'simulationFlowSettings': self.simulationFlowSettings,
            'postProcessSettings': self.postProcessSettings
        }
        # print(self.meshSettings)
        AmpersandIO.printMessage("Writing settings to project_settings.yaml")
        AmpersandPrimitives.dict_to_yaml(settings, f'{self.project_path}/project_settings.yaml')

    # If the project is already existing, load the settings from the project_settings.yaml file
    def load_settings(self):
        AmpersandIO.printMessage("Loading project settings")
        settings = AmpersandPrimitives.yaml_to_dict(f"{self.project_path}/project_settings.yaml")
        self.meshSettings = settings['meshSettings']
        self.physicalProperties = settings['physicalProperties']
        self.numericalSettings = settings['numericalSettings']
        self.inletValues = settings['inletValues']
        self.boundaryConditions = settings['boundaryConditions']
        self.solverSettings = settings['solverSettings']
        self.simulationSettings = settings['simulationSettings']
        self.parallelSettings = settings['parallelSettings']
        self.simulationFlowSettings = settings['simulationFlowSettings']
        self.postProcessSettings = settings['postProcessSettings']
        for geometry in self.meshSettings['geometry']:
            if (geometry['type'] == 'triSurfaceMesh'):
                if (geometry['name'] in self.stl_names):
                    AmpersandIO.printMessage(
                        f"STL file {geometry['name']} already exists in the project, skipping the addition")
                else:
                    self.stl_files.append(geometry)
                    self.stl_names.append(geometry['name'])
        # Change project settings based on loaded settings
        self.internalFlow = self.meshSettings['internalFlow']
        self.on_ground = self.meshSettings['onGround']
        self.transient = self.simulationSettings['transient']
        # self.parallel = self.parallelSettings['parallel']
        # self.snap = self.meshSettings['snap']
        self.refinement = self.meshSettings['fineLevel']
        # self.characteristicLength = self.meshSettings['characteristicLength']
        self.useFOs = self.postProcessSettings['FOs']
        # treat bounds as tuple
        self.meshSettings["geometry"] = AmpersandPrimitives.treat_bounds(
            self.meshSettings["geometry"]
        )
        

    def show_settings(self):
        AmpersandIO.printMessage("Project settings")
        AmpersandIO.printMessage("Mesh Settings")
        AmpersandIO.print_dict(self.meshSettings)
        AmpersandIO.printMessage("Physical Properties")
        AmpersandIO.print_dict(self.physicalProperties)
        AmpersandIO.printMessage("Numerical Settings")
        AmpersandIO.print_dict(self.numericalSettings)
        AmpersandIO.printMessage("Inlet Values")
        AmpersandIO.print_dict(self.inletValues)
        AmpersandIO.printMessage("Boundary Conditions")
        AmpersandIO.print_dict(self.boundaryConditions)
        AmpersandIO.printMessage("Solver Settings")
        AmpersandIO.print_dict(self.solverSettings)
        AmpersandIO.printMessage("Simulation Settings")
        AmpersandIO.print_dict(self.simulationSettings)
        AmpersandIO.printMessage("Parallel Settings")
        AmpersandIO.print_dict(self.parallelSettings)
        AmpersandIO.printMessage("Simulation Flow Settings")
        AmpersandIO.print_dict(self.simulationFlowSettings)
        AmpersandIO.printMessage("Post Process Settings")
        AmpersandIO.print_dict(self.postProcessSettings)

    # If the project is not existing, load the default settings
    def load_default_settings(self):
        self.meshSettings = meshSettings
        self.physicalProperties = physicalProperties
        self.numericalSettings = numericalSettings
        self.inletValues = inletValues
        self.boundaryConditions = boundaryConditions
        self.simulationSettings = simulationSettings
        self.solverSettings = solverSettings
        self.parallelSettings = parallelSettings
        self.simulationFlowSettings = simulationFlowSettings
        self.postProcessSettings = postProcessSettings


    # Add a stl file to the project settings (self.meshSettings)
    def add_stl_to_mesh_settings(self, stl_name, refMin=0, refMax=0, featureEdges=True,
                                 featureLevel=1, purpose='wall', property=None, bounds=None, nLayers=3):
        already_exists = False  # flag to check if the stl file already exists in the project
        idx = 0  # index of the stl file in the list
        # Purpose is wall by default
        # Other purposes are patch, refinementRegion, refinementSurface, cellZone, baffles

        if self.refinement == 0:
            nLayers = 3
        elif self.refinement == 1:
            nLayers = 5
        else:
            nLayers = 7

        stl_geometry = {'name': stl_name, 'type': 'triSurfaceMesh', 'purpose': purpose, 'refineMin': refMin, 'refineMax': refMax,
                'featureEdges': featureEdges, 'featureLevel': featureLevel, 'nLayers': nLayers, 'property': property, 'bounds': bounds}

        # this snippet is to prevent the same stl file from being added multiple times
        # Instead of using set, we are using a list to store the stl files
        # To check if the stl file already exists in the project, stl names are compared.
        # list all the stl files in the project
        for stl_name in self.stl_names:
            if stl_name == stl_geometry['name']:
                AmpersandIO.printMessage(
                    f"STL file {stl_name} already exists in the project")
                already_exists = True
                break
            idx += 1
        if not already_exists:
            self.stl_names.append(stl_name)
            self.stl_files.append(stl_geometry)
            self.meshSettings['geometry'].append(stl_geometry)
        else:
            self.stl_files[idx] = stl_geometry
            return 1  # flag to indicate that the stl file is replaced
        return 0


    def set_property(self, purpose='wall'):
        if purpose == 'inlet':
            U = AmpersandDataInput.get_inlet_values()
            property = tuple(U)
            AmpersandIO.printMessage(
                f"Setting property of {purpose} to {property}")
        elif purpose == 'refinementRegion':
            refLevel = AmpersandIO.get_input_int("Enter refinement level: ")
            property = refLevel
        elif purpose == 'cellZone':
            refLevel = AmpersandIO.get_input_int("Enter refinement level: ")
            createPatches = AmpersandIO.get_input_bool(
                "Create patches for this cellZone? (y/N): ")
            # 0 is just a placeholder for listing the patches
            property = (refLevel, createPatches, 0)
        elif purpose == 'refinementSurface':
            refLevel = AmpersandIO.get_input_int("Enter refinement level: ")
            property = refLevel
        else:
            property = None
        return property

    def add_stl_file(self, stl_path: Union[str, Path], purpose='wall'):
        # Convert paths to Path objects
        stl_path = Path(stl_path)
        
        # Validate input file
        if not stl_path.exists():
            raise FileNotFoundError(f"STL file {stl_path} does not exist")
            
        stl_name = stl_path.name
        if stl_name in self.stl_names:
            raise ValueError(f"STL file {stl_name} already exists in project")

        # Get purpose and properties
        property = None if AmpersandIO.GUIMode else self.set_property(purpose)

        # Calculate bounds from STL
        bounds = tuple(StlAnalysis.compute_bounding_box(stl_path))
        AmpersandIO.printMessage(f"Geometry bounds: {bounds}")

        # Skip feature edges for refinement regions
        feature_edges = purpose not in ('refinementRegion', 'refinementSurface')
        
        # Add to mesh settings
        self.add_stl_to_mesh_settings(
            stl_name=stl_name, 
            purpose=purpose,
            property=property,
            featureEdges=feature_edges,
            bounds=bounds
        )

        self.analyze_stl_file(stl_path)


        # Copy STL file to project
        dest_path = self.project_path / "constant" / "triSurface" / stl_name
        try:
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(stl_path, dest_path)
            AmpersandIO.printMessage(f"Copied {stl_name} to {dest_path}")
        except OSError as e:
            raise RuntimeError(f"Failed to copy STL file: {e}")

        # Set solid name in STL file
        try:
            StlAnalysis.set_stl_solid_name(dest_path)
        except Exception as e:
            raise RuntimeError(f"Failed to set STL solid name: {e}")

        self.current_stl_file = dest_path


    def list_stl_paths(self, project_path: Union[str, Path]):
        stl_dir = Path(project_path) / "constant" / "triSurface"
        if not stl_dir.exists():
            return []
        return list(stl_dir.glob("*.stl"))

    # TODO: this function is not used here
    def remove_stl_file(self, stl_file_number=0):
        # self.list_stl_files()
        stl_file_number = AmpersandIO.get_input(
            "Enter the number of the file to remove: ")
        try:
            stl_file_number = int(stl_file_number)
        except ValueError:
            AmpersandIO.printMessage("Invalid input. Aborting operation")
            return -1
        if stl_file_number < 0 or stl_file_number > len(self.stl_files):
            AmpersandIO.printMessage("Invalid file number. Aborting operation")
            return -1
        stl_file = self.stl_files[stl_file_number]
        stl_name = stl_file['name']
        self.stl_files.remove(stl_file)
        self.stl_names.remove(stl_name)
        stl_path = os.path.join(self.project_path, "constant", "triSurface", stl_name)
        try:
            os.remove(stl_path)
        except OSError as error:
            AmpersandIO.printError(error)
            return -1
        return 0

    def set_flow_type(self, is_internal_flow=False):
        self.internalFlow = is_internal_flow
        self.meshSettings['internalFlow'] = is_internal_flow

    def set_half_model(self, is_half_model: bool):
        self.halfModel = is_half_model
        self.meshSettings['halfModel'] = is_half_model

        if is_half_model:
            AmpersandPrimitives.change_patch_type(
                self.meshSettings['patches'], patch_name='back', new_type='symmetry'
            )

    def set_max_domain_size(self, domain_size, nx, ny, nz):
        self.minX = min(domain_size[0], self.minX)
        self.maxX = max(domain_size[1], self.maxX)
        self.minY = min(domain_size[2], self.minY)
        self.maxY = max(domain_size[3], self.maxY)
        self.minZ = min(domain_size[4], self.minZ)
        self.maxZ = max(domain_size[5], self.maxZ)

        self.meshSettings['domain']['nx'] = nx
        self.meshSettings['domain']['ny'] = ny
        self.meshSettings['domain']['nz'] = nz

    def analyze_stl_file(self, stl_path: Union[str, Path]):
        stl_path = Path(stl_path)
        rho = self.physicalProperties['rho']
        nu = self.physicalProperties['nu']
        U = max(self.inletValues['U'])
        ER = self.meshSettings['addLayersControls']['expansionRatio']

        AmpersandIO.printMessage(f"Analyzing {stl_path.name}")
        # stl_path = os.path.join(self.project_path, "constant", "triSurface", stl_name)
        stlBoundingBox = StlAnalysis.compute_bounding_box(stl_path)
        domain_size, nx, ny, nz, refLevel, target_y, nLayers = StlAnalysis.calc_mesh_settings(stlBoundingBox, nu, rho, U=U, maxCellSize=2.0, expansion_ratio=ER,
                                                                                              onGround=self.on_ground, internalFlow=self.internalFlow,
                                                                                              refinement=self.refinement, halfModel=self.halfModel)
        featureLevel = max(refLevel, 1)
        self.meshSettings = StlAnalysis.set_mesh_settings(self.meshSettings, domain_size, nx, ny, nz, refLevel, featureLevel, nLayers)
        self.set_max_domain_size(domain_size, nx, ny, nz)
        self.meshSettings = StlAnalysis.set_mesh_location(self.meshSettings, stl_path, self.internalFlow)
        refinementBoxLevel = max(2, refLevel-3)
        self.meshSettings = StlAnalysis.addRefinementBoxToMesh(self.meshSettings, stl_path, refLevel=refinementBoxLevel, internalFlow=self.internalFlow)
        if (self.internalFlow == False and self.on_ground == True):
            # if the flow is external and the geometry is on the ground, add a ground refinement box
            self.meshSettings = StlAnalysis.addGroundRefinementBoxToMesh(self.meshSettings, stl_path, refinementBoxLevel)
        # set the layer thickness to 0.5 times the cell size
        self.meshSettings = StlAnalysis.set_layer_thickness(self.meshSettings, 0.5)
        # store the background mesh size for future reference
        maxCellSize = abs((domain_size[1]-domain_size[0])/nx)
        self.meshSettings['maxCellSize'] = maxCellSize
        # self.meshSettings = stlAnalysis.set_min_vol(self.meshSettings, minVol)

    # TODO: not used here
    def adjust_domain_size(self):
        # adjust the domain size based on the bounding box of the stl files
        AmpersandIO.printMessage(
            "Adjusting domain size based on the bounding box of the stl files")
        for stl_file in self.stl_files:
            stl_name = stl_file['name']
            stl_path = os.path.join(self.project_path, "constant", "triSurface", stl_name)
            stlBoundingBox = StlAnalysis.compute_bounding_box(stl_path)
            xmin, xmax, ymin, ymax, zmin, zmax = stlBoundingBox
            self.minX = min(xmin, self.minX)
            self.maxX = max(xmax, self.maxX)
            self.minY = min(ymin, self.minY)
            self.maxY = max(ymax, self.maxY)
            self.minZ = min(zmin, self.minZ)
            self.maxZ = max(zmax, self.maxZ)
            # self.meshSettings = stlAnalysis.set_mesh_location(self.meshSettings, stl_path,self.internalFlow)
        # if the flow is internal, the domain size should be adjusted to include the entire geometry

        self.meshSettings['domain']['minX'] = self.minX
        self.meshSettings['domain']['maxX'] = self.maxX
        self.meshSettings['domain']['minY'] = self.minY
        self.meshSettings['domain']['maxY'] = self.maxY
        self.meshSettings['domain']['minZ'] = self.minZ
        self.meshSettings['domain']['maxZ'] = self.maxZ

    def set_inlet_values(self, U: Optional[tuple[float, float, float]] = None):
        if (not self.internalFlow):  # external flow
            assert U is not None, "Inlet velocity is not set, required for external flow"
            self.inletValues['U'] = U
            self.boundaryConditions['velocityInlet']['u_value'] = U
        else:  # internal flow
            # Use inlet values from the stl file
            AmpersandIO.printMessage(
                "Setting inlet values for various inlet boundaries")
            for stl_file in self.stl_files:
                if stl_file['purpose'] == 'inlet':
                    stl_U = list(stl_file['property']) or U
                    self.boundaryConditions['velocityInlet']['u_value'] = stl_U or U
                    self.inletValues['U'] = stl_U or U

    def set_fluid_properties(self, fluid: dict):
        self.physicalProperties['rho'] = fluid['rho']
        self.physicalProperties['nu'] = fluid['nu']

    def set_parallel(self, n_core: int):
        self.parallelSettings['numberOfSubdomains'] = n_core

    # setting the purpose of a patch. Used for setting the boundary conditions
    def set_purpose(self, patch, purpose='wall'):
        purposes = ['wall', 'inlet', 'outlet', 'refinementRegion', 'refinementSurface',
                    'cellZone', 'baffles', 'symmetry', 'cyclic', 'empty',]
        # purposes = ['wall', 'inlet','outlet', 'refinementRegion', 'refinementSurface', 'cellZone', 'baffles']
        if purpose not in purposes:
            AmpersandIO.printMessage(
                "Invalid purpose. Setting purpose to wall")
            purpose = 'wall'
        patch['purpose'] = purpose

    # set the turbulence model for the simulation
    def set_turbulence_model(self, turbulence_model='kOmegaSST'):
        turbulence_model = AmpersandDataInput.choose_turbulence_model()
        self.solverSettings['turbulenceModel'] = turbulence_model

    def set_transient_settings(self, is_transient: bool):
        self.transient = is_transient

        if self.transient:
            AmpersandIO.printMessage("Transient simulation settings")
            self.simulationSettings['transient'] = True
            self.simulationSettings['application'] = 'pimpleFoam'
            self.simulationFlowSettings['solver'] = 'pimpleFoam'
            self.simulationSettings['endTime'] = AmpersandIO.get_input_float(
                "End time: ")
            self.simulationSettings['writeInterval'] = AmpersandIO.get_input_float(
                "Write interval: ")
            self.simulationSettings['deltaT'] = AmpersandIO.get_input_float(
                "Time step: ")
            self.simulationSettings['adjustTimeStep'] = 'no'
            self.simulationSettings['maxCo'] = 0.9
            self.numericalSettings['ddtSchemes']['default'] = 'Euler'
            # if steady state, SIMPLEC is used. If transient, PIMPLE is used
            # for PIMPLE, the relaxation factors are set to 0.7 and p = 0.3
            self.numericalSettings['relaxationFactors']['p'] = 0.3

    def set_on_ground(self, on_ground: bool):
        self.on_ground = on_ground
        self.meshSettings['onGround'] = on_ground
        if on_ground:
            AmpersandPrimitives.change_patch_type(self.meshSettings['patches'], patch_name='ground',
                                                  new_type='wall')

    def set_refinement_level(self, fine_level: int):
        self.refinement = fine_level
        self.meshSettings['fineLevel'] = fine_level


    def set_post_process_settings(self, useFOs: bool):
        self.useFOs = useFOs
        self.postProcessSettings['FOs'] = useFOs

        meshPoint = list(
            self.meshSettings['castellatedMeshControls']['locationInMesh'])
        self.postProcessSettings['massFlow'] = True
        self.postProcessSettings['minMax'] = True
        self.postProcessSettings['yPlus'] = True
        self.postProcessSettings['forces'] = True
        # the default probe location for monitoring of flow variables
        self.postProcessSettings['probeLocations'].append(meshPoint)

    def write_project_files(self):
        if (not self.project_path.exists()):
            raise FileNotFoundError(f"Project not found at: {self.project_path}")


        create_boundary_conditions(self.meshSettings, self.boundaryConditions, f"{self.project_path}/0")

        # go inside the constant directory
        AmpersandIO.printMessage("Creating physical properties and turbulence properties")
        # create transportProperties file
        tranP = create_transportPropertiesDict(self.physicalProperties)
        # create turbulenceProperties file
        turbP = create_turbulencePropertiesDict(self.physicalProperties)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/constant/transportProperties", tranP)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/constant/turbulenceProperties", turbP)

        # create the controlDict file
        AmpersandIO.printMessage("Creating the system files")
        controlDict = createControlDict(self.simulationSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/system/controlDict", controlDict)
        
        blockMeshDict = generate_blockMeshDict(self.meshSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/system/blockMeshDict", blockMeshDict)
        
        snappyHexMeshDict = generate_snappyHexMeshDict(self.meshSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/system/snappyHexMeshDict", snappyHexMeshDict)
        
        surfaceFeatureExtractDict = create_surfaceFeatureExtractDict(self.meshSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/system/surfaceFeatureExtractDict", surfaceFeatureExtractDict)
        
        fvSchemesDict = create_fvSchemesDict(self.numericalSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/system/fvSchemes", fvSchemesDict)
        
        fvSolutionDict = create_fvSolutionDict(self.numericalSettings, self.solverSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/system/fvSolution", fvSolutionDict)
        
        decomposeParDict = createDecomposeParDict(self.parallelSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/system/decomposeParDict", decomposeParDict)
        
        FODict = postProcess.create_FOs(self.meshSettings, self.postProcessSettings, useFOs=self.useFOs)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/system/FOs", FODict)

        # create mesh script
        AmpersandIO.printMessage("Creating scripts for meshing and running the simulation")
        meshScript = ScriptGenerator.generate_mesh_script(self.simulationFlowSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/mesh", meshScript)
        
        # create simulation script
        simulationScript = ScriptGenerator.generate_simulation_script(self.simulationFlowSettings)
        AmpersandPrimitives.write_dict_to_file(f"{self.project_path}/run", simulationScript)
        
        AmpersandPrimitives.crlf_to_LF(f"{self.project_path}/mesh")
        AmpersandPrimitives.crlf_to_LF(f"{self.project_path}/run")
        
        if os.name != 'nt':
            os.chmod(f"{self.project_path}/mesh", 0o755)
            os.chmod(f"{self.project_path}/run", 0o755)
        
        AmpersandIO.printMessage("\n-----------------------------------")
        AmpersandIO.printMessage("Project files created successfully!")
        AmpersandIO.printMessage("-----------------------------------\n")

