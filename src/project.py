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
import shutil
from typing import Optional
from src.primitives import ampersandPrimitives, ampersandIO, AmpersandDataInput
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
    def __init__(self, GUIMode=False, window=None):
        # project path = project_directory_path/user_name/project_name
        self.GUIMode = GUIMode
        if GUIMode and window != None:
            self.window = window
        else:
            self.window = None
        self.current_stl_file = None   # current stl file being processed
        self.project_directory_path = None
        self.project_name = None
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
        self.project_path = None
        self.existing_project = False  # flag to check if the project is already existing
        self.stl_files = []  # list to store the settings for stl files
        self.stl_names = []  # list to store the names of the stl files
        self.internalFlow = False  # default is external flow
        self.on_ground = False  # default is off the ground
        self.halfModel = False  # default is full model
        self.parallel = True  # default is parallel simulation
        self.snap = True  # default is to use snappyHexMesh
        self.transient = False  # default is steady state
        self.refinement = 0  # 0: coarse, 1: medium, 2: fine
        self.characteristicLength = None  # default characteristic length
        self.useFOs = False  # default is not to use function objects
        self.current_modification = None  # current modification to the project settings
        # flag to check if the current working directory is the project directory
        self.inside_project_directory = False
        self.mod_options = ["Background Mesh", "Add Geometry", "Refinement Levels", "Mesh Point", "Boundary Conditions", "Fluid Properties", "Numerical Settings",
                            "Simulation Control Settings", "Turbulence Model", "Post Processing Settings"]
        self.minX, self.maxX, self.minY, self.maxY, self.minZ, self.maxZ = - \
            1e-3, 1e-3, -1e-3, 1e-3, -1e-3, 1e-3
    # --------------------------------------------------------------------
    # Methods to handle the project summary and changes
    # --------------------------------------------------------------------

    def summarize_boundary_conditions(self):
        bcs = ampersandPrimitives.list_boundary_conditions(self.meshSettings)
        return bcs



    def summarize_project(self):
        trueFalse = {True: 'Yes', False: 'No'}
        ampersandIO.show_title("Project Summary")
        # ampersandIO.printMessage(f"Project directory: {self.project_directory_path}")
        ampersandIO.printFormat(
            "Project name", self.project_name, GUIMode=self.GUIMode, window=self.window)
        ampersandIO.printFormat(
            "Project path", self.project_path, GUIMode=self.GUIMode, window=self.window)

        ampersandIO.printMessage(
            f"Internal Flow: {trueFalse[self.internalFlow]}", GUIMode=self.GUIMode, window=self.window)
        if (self.internalFlow == False):
            ampersandIO.printMessage(
                f"On Ground: {trueFalse[self.on_ground]}", GUIMode=self.GUIMode, window=self.window)
        ampersandIO.printMessage(
            f"Transient: {trueFalse[self.transient]}", GUIMode=self.GUIMode, window=self.window)
        self.summarize_background_mesh()
        self.list_stl_files()

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
        ampersandIO.printMessage(f"Domain size:{'X':>10}{'Y':>10}{'Z':>10}")
        ampersandIO.printMessage(
            f"Min         {minX:>10.3f}{minY:>10.3f}{minZ:>10.3f}")
        ampersandIO.printMessage(
            f"Max         {maxX:>10.3f}{maxY:>10.3f}{maxZ:>10.3f}")
        ampersandIO.printMessage(f"Background mesh size: {nx}x{ny}x{nz} cells")
        ampersandIO.printMessage(
            f"Background cell size: {self.meshSettings['maxCellSize']} m")

    def change_boundary_condition(self, bcName, newBC):
        if not self.internalFlow:  # if it is external flow
            bcPatches = self.meshSettings['patches']
            for aPatch in self.meshSettings['patches']:
                if bcName == aPatch['name']:
                    aPatch['type'] = newBC
                    ampersandIO.printMessage(
                        f"Boundary condition {bcName} changed to {newBC}", GUIMode=self.GUIMode, window=self.window)
                    return 0
            if bcName in bcPatches:
                self.meshSettings['patches'][bcName]['type'] = newBC
                self.meshSettings['bcPatches'][bcName]['purpose'] = newBC
                newProperty = self.set_property(newBC)
                self.meshSettings['bcPatches'][bcName]['property'] = newProperty
                ampersandIO.printMessage(
                    f"Boundary condition {bcName} changed to {newBC}", GUIMode=self.GUIMode, window=self.window)
                return 0
            else:
                ampersandIO.printMessage(
                    "Boundary condition not found in the list")
        geometry = self.meshSettings['geometry']
        for stl in geometry:
            if stl['name'] == bcName:
                stl['purpose'] = newBC
                newProperty = self.set_property(newBC)
                self.meshSettings['bcPatches'][bcName]['property'] = newProperty
                ampersandIO.printMessage(
                    f"Boundary condition of {bcName} changed to {newBC}", GUIMode=self.GUIMode, window=self.window)
                return 0
        return -1

    def change_stl_refinement_level(self, stl_file_number=0):
        ampersandIO.printMessage("Changing refinement level")
        refMin = ampersandIO.get_input_int("Enter new refMin: ")
        refMax = ampersandIO.get_input_int("Enter new refMax: ")
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
            ampersandIO.printMessage("STL file not found in the geometry list")

        # return project

    def choose_modification(self):
        current_modification = ampersandIO.get_option_choice(prompt="Choose any option for project modification: ",
                                                             options=self.mod_options, title="\nModify Project Settings")
        self.current_modification = self.mod_options[current_modification]
        ampersandIO.printMessage(
            f"Current modification: {self.current_modification}", GUIMode=self.GUIMode, window=self.window)

    def choose_modification_categorized(self):
        options = ['Mesh', 'Boundary Conditions', 'Fluid Properties', 'Numerical Settings',
                   'Simulation Control Settings', 'Turbulence Model', 'Post Processing Settings']
        current_modification = ampersandIO.get_option_choice(prompt="Choose any option for project modification: ",
                                                             options=options, title="\nModify Project Settings")
        mesh_options = ['Background Mesh', 'Mesh Point',
                        'Add Geometry', 'Refinement Levels']

        if current_modification < 0 or current_modification > len(options)-1:
            ampersandIO.printMessage("Invalid option. Aborting operation")
            return -1
        if current_modification == 0:
            self.current_modification = mesh_options[ampersandIO.get_option_choice(prompt="Choose any option for mesh modification: ",
                                                                                   options=mesh_options, title="\nModify Mesh Settings")]
        else:
            self.current_modification = options[current_modification]

    def remove_duplicate_stl_files(self):
        # detect duplicate dictionaries in the list
        seen = set()
        new_list = []
        for d in self.stl_files:
            t = tuple(d.items())
            if t not in seen:
                seen.add(t)
                new_list.append(d)
        self.stl_files = new_list
        self.meshSettings['geometry'] = ampersandPrimitives.remove_duplicate_dicts(
            self.meshSettings['geometry'])
        # print("stl_files",self.stl_files)
        # print("Mesh settings",self.meshSettings["geometry"])

    def set_project_directory(self, project_directory_path):
        if self.GUIMode:
            stopWhenError = False
        else:
            stopWhenError = True
        if project_directory_path is None:
            if stopWhenError:
                ampersandIO.printMessage(
                    "No directory selected. Aborting project creation.")
                exit()
            else:
                return -1
        # assert os.path.exists(project_directory_path), "The chosen directory does not exist"
        if not os.path.exists(project_directory_path):
            if stopWhenError:
                ampersandIO.printMessage(
                    "The chosen directory does not exist. Aborting project creation.")
                exit()
            else:
                self.project_directory_path = None
                return -1
        self.project_directory_path = project_directory_path

    def set_project_name(self, project_name):
        self.project_name = project_name

    def set_user_name(self, user_name):
        self.user_name = user_name

    # create the project path for the user and project name
    def create_project_path_user(self):
        if not self.project_directory_path:
            ampersandIO.printWarning(
                "No directory selected. Aborting project creation.", GUIMode=self.GUIMode)
            return -1
        self.project_path = os.path.join(
            self.project_directory_path, self.user_name, self.project_name)

    # To create the project path for a new project with the project name
    def create_project_path(self):
        if not self.project_directory_path:
            ampersandIO.printWarning(
                "No directory selected. Aborting project creation.")
            return -1
        self.project_path = os.path.join(
            self.project_directory_path, self.project_name)

    # this is to set the project path if the project is already existing
    # useful for opening existing projects and modifying the settings
    def set_project_path(self, project_path):
        if project_path is None:
            if self.GUIMode == False:
                ampersandIO.printWarning(
                    "No project path selected. Aborting project creation.", GUIMode=self.GUIMode)
            # ampersandIO.printWarning("No project path selected. Aborting project creation/modification.",GUIMode=self.GUIMode)
            return -1
            # exit()
        if os.path.exists(project_path):
            settings_file = os.path.join(project_path, "project_settings.yaml")
            if os.path.exists(settings_file):
                ampersandIO.printMessage(
                    "Project found, loading project settings", GUIMode=self.GUIMode, window=self.window)
                self.existing_project = True
                self.project_path = project_path
                return 0
            else:
                if self.GUIMode == False:
                    ampersandIO.printWarning(
                        "Settings file not found. Please open an Ampersand case directory.", GUIMode=self.GUIMode)
                # ampersandIO.printWarning("Settings file not found. Please open an Ampersand case directory.",GUIMode=self.GUIMode)
                # TO DO: Add the code socket to create a new project here
                return -1
        else:
            if self.GUIMode == False:
                ampersandIO.printWarning(
                    "Project path does not exist. Aborting project creation/opening.", GUIMode=self.GUIMode)
            # ampersandIO.printWarning("Project path does not exist. Aborting project creation/opening.",GUIMode=self.GUIMode)
            return -1

    # check if the project path exists and if the project is already existing
    def check_project_path(self):
        if os.path.exists(self.project_path):
            settings_file = os.path.join(
                self.project_path, "project_settings.yaml")
            if os.path.exists(settings_file):
                ampersandIO.printWarning(
                    "Project already exists, loading project settings", GUIMode=self.GUIMode)
                self.existing_project = True
                return 0
            else:
                self.existing_project = False
                return -1
        else:
            self.existing_project = False
            return -1

    # Wrapper of cwd with error handling
    def go_inside_directory(self):
        try:
            os.chdir(self.project_path)
        except OSError as error:
            ampersandIO.printError(error)
        cwd = os.getcwd()
        ampersandIO.printMessage(
            f"Working directory: {cwd}", GUIMode=self.GUIMode, window=self.window)
        self.inside_project_directory = True

    # Check if the 0 directory exists in the project directory
    def check_0_directory(self):
        if not os.path.exists("0"):
            ampersandIO.printWarning(
                "0 directory does not exist.", GUIMode=self.GUIMode)
            ampersandIO.printMessage(
                "Checking for 0.orig directory", GUIMode=self.GUIMode, window=self.window)
            if os.path.exists("0.orig"):
                ampersandIO.printMessage(
                    "0.orig directory found. Copying to 0 directory", GUIMode=self.GUIMode, window=self.window)
                shutil.copytree("0.orig", "0")
            else:
                ampersandIO.printWarning(
                    "0.orig directory not found. Aborting project creation.", GUIMode=self.GUIMode)
                return -1
        return 0

    # Check if the constant directory exists in the project directory
    def check_constant_directory(self):
        if not os.path.exists("constant"):
            ampersandIO.printWarning(
                "constant directory does not exist.", GUIMode=self.GUIMode)
            # ampersandIO.printError("constant directory is necessary for the project")
            return -1
        return 0

    # Check if the system directory exists in the project directory
    def check_system_directory(self):
        if not os.path.exists("system"):
            ampersandIO.printWarning(
                "system directory does not exist.", GUIMode=self.GUIMode)
            # ampersandIO.printError("system directory is necessary for the project")
            return -1
        return 0

    # Check if the constant/triSurface directory exists in the project directory
    def check_triSurface_directory(self):
        if not os.path.exists("constant/triSurface"):
            ampersandIO.printWarning(
                "triSurface directory does not exist.", GUIMode=self.GUIMode)
            # ampersandIO.printError("triSurface directory is necessary for the project")
            return -1
        # if exists, check if the stl files are present
        stl_files = os.listdir("constant/triSurface")

    # to check whether log files are present in the project directory
    def check_log_files(self):
        log_files = os.listdir()
        if 'log.simpleFoam' in log_files:
            ampersandIO.printMessage(
                "Simulation log file found", GUIMode=self.GUIMode, window=self.window)
            return 1
        if 'log.pimpleFoam' in log_files:
            ampersandIO.printMessage(
                "Simulation log file found", GUIMode=self.GUIMode, window=self.window)
            return 1
        return 0

    # to check whether the U and p files are present in the postProcess directory
    def check_post_process_files(self):
        if (not os.path.exists("postProcessing/probe/0")):
            ampersandIO.printWarning(
                "postProcess directory does not exist", GUIMode=self.GUIMode)
            return 0
        postProcess_files = os.listdir("postProcessing/probe/0")
        if 'U' in postProcess_files and 'p' in postProcess_files:
            ampersandIO.printMessage(
                "U and p files found in postProcess directory", GUIMode=self.GUIMode, window=self.window)
            return 1
        return 0

    def check_forces_files(self):
        if (not os.path.exists("postProcessing/forces/0")):
            ampersandIO.printWarning(
                "forces directory does not exist", GUIMode=self.GUIMode)
            return 0
        forces_files = os.listdir("postProcessing/forces/0")
        if 'force.dat' in forces_files:
            ampersandIO.printMessage("force.dat found in forces directory")
            return 1
        return 0

    # Create the project directory in the specified location.
    # 0, constant, system, constant/triSurface directories are created.

    def create_project(self):
        # check if the project path exists
        if self.project_path is None:
            ampersandIO.printError(
                "No project path selected. Aborting project creation.", GUIMode=self.GUIMode)
            return -1
        if os.path.exists(self.project_path):
            ampersandIO.printWarning(
                "Project already exists. Skipping the creation of directories", GUIMode=self.GUIMode)
            self.existing_project = True
        else:
            ampersandIO.printMessage("Creating project directory")
            try:
                os.makedirs(self.project_path)

            except OSError as error:
                ampersandIO.printError(error)
        try:
            os.chdir(self.project_path)
        except OSError as error:
            ampersandIO.printError(error)
        cwd = os.getcwd()
        ampersandIO.printMessage(
            f"Working directory: {cwd}", GUIMode=self.GUIMode, window=self.window)

        # create 0, constant and system directory
        try:
            os.mkdir("0")
            os.mkdir("constant")
            os.mkdir("system")
            os.mkdir("constant/triSurface")
        except OSError as error:
            ampersandIO.printError(
                "File system already exists. Skipping the creation of directories", GUIMode=self.GUIMode)
            return -1
        return 0  # return 0 if the project is created successfully

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
        ampersandIO.printMessage(
            "Writing settings to project_settings.yaml", GUIMode=self.GUIMode, window=self.window)
        ampersandPrimitives.dict_to_yaml(settings, 'project_settings.yaml')

    # If the project is already existing, load the settings from the project_settings.yaml file
    def load_settings(self):
        ampersandIO.printMessage(
            "Loading project settings", GUIMode=self.GUIMode, window=self.window)
        settings = ampersandPrimitives.yaml_to_dict('project_settings.yaml')
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
                    ampersandIO.printMessage(
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
        project_name = self.project_path.split("/")[-1]
        self.project_name = project_name
        # treat bounds as tuple
        self.meshSettings["geometry"] = ampersandPrimitives.treat_bounds(
            self.meshSettings["geometry"])

    def show_settings(self):
        ampersandIO.printMessage("Project settings")
        ampersandIO.printMessage("Mesh Settings")
        ampersandIO.print_dict(self.meshSettings)
        ampersandIO.printMessage("Physical Properties")
        ampersandIO.print_dict(self.physicalProperties)
        ampersandIO.printMessage("Numerical Settings")
        ampersandIO.print_dict(self.numericalSettings)
        ampersandIO.printMessage("Inlet Values")
        ampersandIO.print_dict(self.inletValues)
        ampersandIO.printMessage("Boundary Conditions")
        ampersandIO.print_dict(self.boundaryConditions)
        ampersandIO.printMessage("Solver Settings")
        ampersandIO.print_dict(self.solverSettings)
        ampersandIO.printMessage("Simulation Settings")
        ampersandIO.print_dict(self.simulationSettings)
        ampersandIO.printMessage("Parallel Settings")
        ampersandIO.print_dict(self.parallelSettings)
        ampersandIO.printMessage("Simulation Flow Settings")
        ampersandIO.print_dict(self.simulationFlowSettings)
        ampersandIO.printMessage("Post Process Settings")
        ampersandIO.print_dict(self.postProcessSettings)

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

    # Create the settings for the project or load the existing settings
    def create_settings(self):
        if self.existing_project:
            ampersandIO.printMessage(
                "Project already exists. Loading project settings")
            try:
                self.load_settings()
            except FileNotFoundError:
                ampersandIO.printMessage(
                    "Settings file not found. Loading default settings")
                self.load_default_settings()
                self.write_settings()
        else:
            self.load_default_settings()
            self.write_settings()

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

        stl_ = {'name': stl_name, 'type': 'triSurfaceMesh', 'purpose': purpose, 'refineMin': refMin, 'refineMax': refMax,
                'featureEdges': featureEdges, 'featureLevel': featureLevel, 'nLayers': nLayers, 'property': property, 'bounds': bounds}

        # this snippet is to prevent the same stl file from being added multiple times
        # Instead of using set, we are using a list to store the stl files
        # To check if the stl file already exists in the project, stl names are compared.
        # list all the stl files in the project
        for stl_name in self.stl_names:
            if stl_name == stl_['name']:
                ampersandIO.printMessage(
                    f"STL file {stl_name} already exists in the project")
                already_exists = True
                break
            idx += 1
        if not already_exists:
            self.stl_names.append(stl_name)
            self.stl_files.append(stl_)
        else:
            self.stl_files[idx] = stl_
            return 1  # flag to indicate that the stl file is replaced
        return 0


    def set_property(self, purpose='wall'):
        if purpose == 'inlet':
            U = AmpersandDataInput.get_inlet_values()
            property = tuple(U)
            ampersandIO.printMessage(
                f"Setting property of {purpose} to {property}")
        elif purpose == 'refinementRegion':
            refLevel = ampersandIO.get_input_int("Enter refinement level: ")
            property = refLevel
        elif purpose == 'cellZone':
            refLevel = ampersandIO.get_input_int("Enter refinement level: ")
            createPatches = ampersandIO.get_input_bool(
                "Create patches for this cellZone? (y/N): ")
            # 0 is just a placeholder for listing the patches
            property = (refLevel, createPatches, 0)
        elif purpose == 'refinementSurface':
            refLevel = ampersandIO.get_input_int("Enter refinement level: ")
            property = refLevel
        else:
            property = None
        return property

    def ask_stl_settings(self, stl_file):
        ampersandIO.printMessage(f"Settings of the {stl_file['name']} file")
        stl_file['refineMin'] = ampersandIO.get_input("Min Refinement: ")
        stl_file['refineMax'] = ampersandIO.get_input("Max Refinement: ")
        featureEdges = ampersandIO.get_input("Refine Feature Edges?: (y/N) ")
        if (featureEdges == 'y'):
            stl_file['featureEdges'] = True
        else:
            stl_file['featureEdges'] = False
        stl_file['featureLevel'] = ampersandIO.get_input("Feature Level: ")
        stl_file['nLayers'] = ampersandIO.get_input("Number of Layers: ")

    def add_stl_to_project(self):
        for stl_file in self.stl_files:
            self.meshSettings['geometry'].append(stl_file)
        self.remove_duplicate_stl_files()

    # to only copy the STL file to the project directory and add it to the STL list
    def add_stl_file(self):
        stl_file = ampersandPrimitives.ask_for_file(
            [("STL Geometry", "*.stl"), ("OBJ Geometry", "*.obj")], self.GUIMode)
        if stl_file is None:
            ampersandIO.printWarning(
                "No file selected. Please select STL file if necessary.", GUIMode=self.GUIMode)
            return -1
        if os.path.exists(stl_file):
            # add the stl file to the project
            # This is a bit confusing.
            # stl_name is the name of the file, stl_file is the path to the file
            file_path_to_token = stl_file.split("/")
            stl_name = file_path_to_token[-1]
            if stl_name in self.stl_names:
                ampersandIO.printWarning(
                    f"STL file {stl_name} already exists in the project", GUIMode=self.GUIMode)
                return -1
            else:  # this is to prevent the bug of having the same file added multiple times
                if self.GUIMode:
                    purpose = "wall"
                    property = None
                else:
                    purpose = AmpersandDataInput.get_purpose()
                    property = self.set_property(purpose)
                bounds = StlAnalysis.compute_bounding_box(stl_file)
                bounds = tuple(bounds)
                ampersandIO.printMessage(
                    f"Bounds of the geometry: {bounds}", GUIMode=self.GUIMode, window=self.window)
                if purpose == 'refinementRegion' or purpose == 'refinementSurface':
                    featureEdges = False
                else:
                    featureEdges = True
                self.add_stl_to_mesh_settings(
                    stl_name, purpose=purpose, property=property, featureEdges=featureEdges, bounds=bounds)
            # this is the path to the constant/triSurface inside project directory where STL will be copied
            stl_path = os.path.join(
                self.project_path, "constant", "triSurface", stl_name)
            try:
                ampersandIO.printMessage(
                    f"Copying {stl_name} to the project directory", GUIMode=self.GUIMode, window=self.window)
                shutil.copy(stl_file, stl_path)
            except OSError as error:
                ampersandIO.printError(error, GUIMode=self.GUIMode)
                return -1
            try:
                StlAnalysis.set_stl_solid_name(stl_path)
            except Exception as error:
                ampersandIO.printError(error, GUIMode=self.GUIMode)
                return -1
        else:
            ampersandIO.printError(
                "File does not exist. Aborting project creation.", GUIMode=self.GUIMode)
            return -1
        self.current_stl_file = stl_path
        return 0

    # this is a wrapper of the primitives
    def list_stl_files(self):
        ampersandPrimitives.list_stl_files(
            self.stl_files, self.GUIMode, self.window)

    def list_stl_paths(self):
        stl_paths = []
        for stl_file in self.stl_files:
            stl_path = os.path.join(
                self.project_path, "constant", "triSurface", stl_file['name'])
            # ampersandIO.printMessage(stl_path)
            stl_paths.append(stl_path)
        return stl_paths

    def remove_stl_file(self, stl_file_number=0):
        # self.list_stl_files()
        stl_file_number = ampersandIO.get_input(
            "Enter the number of the file to remove: ")
        try:
            stl_file_number = int(stl_file_number)
        except ValueError:
            ampersandIO.printMessage("Invalid input. Aborting operation")
            return -1
        if stl_file_number < 0 or stl_file_number > len(self.stl_files):
            ampersandIO.printMessage("Invalid file number. Aborting operation")
            return -1
        stl_file = self.stl_files[stl_file_number]
        stl_name = stl_file['name']
        self.stl_files.remove(stl_file)
        self.stl_names.remove(stl_name)
        stl_path = os.path.join(
            self.project_path, "constant", "triSurface", stl_name)
        try:
            os.remove(stl_path)
        except OSError as error:
            ampersandIO.printError(error)
            return -1
        return 0

    def set_flow_type(self, internalFlow=False):
        self.internalFlow = internalFlow
        self.meshSettings['internalFlow'] = self.internalFlow

    def set_half_model(self, half_model: bool):
        self.halfModel = half_model
        self.meshSettings['halfModel'] = half_model

        if half_model:
            ampersandPrimitives.change_patch_type(
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

    def analyze_stl_file(self, stl_file_number=0):
        rho = self.physicalProperties['rho']
        nu = self.physicalProperties['nu']
        U = max(self.inletValues['U'])
        ER = self.meshSettings['addLayersControls']['expansionRatio']
        try:
            stl_file_number = int(stl_file_number)
        except ValueError:
            ampersandIO.printError(
                "Invalid input. Aborting operation", GUIMode=self.GUIMode)
            return -1
        if stl_file_number < 0 or stl_file_number > len(self.stl_files):
            ampersandIO.printError(
                "Invalid file number. Aborting operation", GUIMode=self.GUIMode)
            return -1
        stl_file = self.stl_files[stl_file_number]
        stl_name = stl_file['name']
        ampersandIO.printMessage(
            f"Analyzing {stl_name}", GUIMode=self.GUIMode, window=self.window)
        stl_path = os.path.join(
            self.project_path, "constant", "triSurface", stl_name)
        stlBoundingBox = StlAnalysis.compute_bounding_box(stl_path)
        domain_size, nx, ny, nz, refLevel, target_y, nLayers = StlAnalysis.calc_mesh_settings(stlBoundingBox, nu, rho, U=U, maxCellSize=2.0, expansion_ratio=ER,
                                                                                              onGround=self.on_ground, internalFlow=self.internalFlow,
                                                                                              refinement=self.refinement, halfModel=self.halfModel,
                                                                                              GUI=self.GUIMode, window=self.window)
        featureLevel = max(refLevel, 1)
        self.meshSettings = StlAnalysis.set_mesh_settings(
            self.meshSettings, domain_size, nx, ny, nz, refLevel, featureLevel, nLayers=nLayers)
        self.set_max_domain_size(domain_size, nx, ny, nz)
        self.meshSettings = StlAnalysis.set_mesh_location(
            self.meshSettings, stl_path, self.internalFlow)
        refinementBoxLevel = max(2, refLevel-3)
        self.meshSettings = StlAnalysis.addRefinementBoxToMesh(
            meshSettings=self.meshSettings, stl_path=stl_path, refLevel=refinementBoxLevel, internalFlow=self.internalFlow)
        if (self.internalFlow == False and self.on_ground == True):
            # if the flow is external and the geometry is on the ground, add a ground refinement box
            self.meshSettings = StlAnalysis.addGroundRefinementBoxToMesh(
                meshSettings=self.meshSettings, stl_path=stl_path, refLevel=refinementBoxLevel)
        # set the layer thickness to 0.5 times the cell size
        self.meshSettings = StlAnalysis.set_layer_thickness(
            self.meshSettings, 0.5)
        # store the background mesh size for future reference
        maxCellSize = abs((domain_size[1]-domain_size[0])/nx)
        self.meshSettings['maxCellSize'] = maxCellSize
        # self.meshSettings = stlAnalysis.set_min_vol(self.meshSettings, minVol)
        return 0

    def adjust_domain_size(self):
        # adjust the domain size based on the bounding box of the stl files
        ampersandIO.printMessage(
            "Adjusting domain size based on the bounding box of the stl files", GUIMode=self.GUIMode, window=self.window)
        for stl_file in self.stl_files:
            stl_name = stl_file['name']
            stl_path = os.path.join(
                self.project_path, "constant", "triSurface", stl_name)
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

    def set_inlet_values(self, U: Optional[tuple[float]] = None):
        if (not self.internalFlow):  # external flow
            assert U is not None, "Inlet velocity is not set, required for external flow"
            self.inletValues['U'] = U
            self.boundaryConditions['velocityInlet']['u_value'] = U
        else:  # internal flow
            # Use inlet values from the stl file
            ampersandIO.printMessage(
                "Setting inlet values for various inlet boundaries")
            for stl_file in self.stl_files:
                if stl_file['purpose'] == 'inlet':
                    stl_U = list(stl_file['property']) or U
                    self.boundaryConditions['velocityInlet']['u_value'] = stl_U or U
                    self.inletValues['U'] = stl_U or U

    def set_fluid_properties(self, fluid: dict):
        self.physicalProperties['rho'] = fluid['rho']
        self.physicalProperties['nu'] = fluid['nu']

    # def set_transient(self):
    #    self.transient = ampersandIO.get_input_bool("Transient simulation (y/N)?: ")

    def set_parallel(self, n_core: int):
        self.parallelSettings['numberOfSubdomains'] = n_core

    # setting the purpose of a patch. Used for setting the boundary conditions
    def set_purpose(self, patch, purpose='wall'):
        purposes = ['wall', 'inlet', 'outlet', 'refinementRegion', 'refinementSurface',
                    'cellZone', 'baffles', 'symmetry', 'cyclic', 'empty',]
        # purposes = ['wall', 'inlet','outlet', 'refinementRegion', 'refinementSurface', 'cellZone', 'baffles']
        if purpose not in purposes:
            ampersandIO.printMessage(
                "Invalid purpose. Setting purpose to wall")
            purpose = 'wall'
        patch['purpose'] = purpose

    # choose turbulence model for the simulation
    def choose_turbulence_model(self):
        turbulence_models = ['kOmegaSST', 'kEpsilon', 'SpalartAllmaras']
        turbulence_model = AmpersandDataInput.get_option_choice(
            "Choose turbulence model: ", turbulence_models)
        self.solverSettings['turbulenceModel'] = turbulence_model

    # set the turbulence model for the simulation
    def set_turbulence_model(self, turbulence_model='kOmegaSST'):
        turbulence_model = AmpersandDataInput.choose_turbulence_model()
        self.solverSettings['turbulenceModel'] = turbulence_model

    def set_transient_settings(self, transient: bool):
        self.transient = transient

        if self.transient:
            ampersandIO.printMessage("Transient simulation settings")
            self.simulationSettings['transient'] = True
            self.simulationSettings['application'] = 'pimpleFoam'
            self.simulationFlowSettings['solver'] = 'pimpleFoam'
            self.simulationSettings['endTime'] = ampersandIO.get_input_float(
                "End time: ")
            self.simulationSettings['writeInterval'] = ampersandIO.get_input_float(
                "Write interval: ")
            self.simulationSettings['deltaT'] = ampersandIO.get_input_float(
                "Time step: ")
            self.simulationSettings['adjustTimeStep'] = 'no'
            self.simulationSettings['maxCo'] = 0.9
            self.numericalSettings['ddtSchemes']['default'] = 'Euler'
            # if steady state, SIMPLEC is used. If transient, PIMPLE is used
            # for PIMPLE, the relaxation factors are set to 0.7 and p = 0.3
            self.numericalSettings['relaxationFactors']['p'] = 0.3

    def set_ground_type(self, on_ground: bool):
        self.on_ground = on_ground
        self.meshSettings['onGround'] = on_ground
        if on_ground:
            ampersandPrimitives.change_patch_type(self.meshSettings['patches'], patch_name='ground',
                                                  new_type='wall')

    def set_refinement_level(self, fine_level: int):
        self.refinement = fine_level
        self.meshSettings['fineLevel'] = fine_level

    def set_global_refinement_level(self, refinement=0):
        self.refinement = refinement
        self.meshSettings['fineLevel'] = refinement

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

    def get_probe_location(self):
        point = postProcess.get_probe_location()
        # for internal flows, the point should be inside stl
        # for external flows, the point should be outside stl
        # TO DO
        self.postProcessSettings['probeLocations'].append(point)

    def create_project_files(self):
        # (meshSettings, physicalProperties, numericalSettings, inletValues, boundaryConditions)=caseSettings
        # check if the current working directory is the project directory
        if not os.path.exists(self.project_path):
            ampersandIO.printMessage(
                "Project directory does not exist. Aborting project creation.")
            return -1
        if os.getcwd() != self.project_path:
            os.chdir(self.project_path)

        # Remove the existing 0.orig directory if it exists.
        # This is to prevent the error of copying the old 0.orig directory to 0 directory
        if os.path.exists("0.orig"):
            shutil.rmtree("0.orig")
        # create the initial conditions file
        ampersandIO.printMessage("Creating boundary conditions")
        # check if the 0 directory exists
        if not os.path.exists("0"):
            # create the 0 directory
            os.mkdir("0")

        create_boundary_conditions(
            self.meshSettings, self.boundaryConditions, "0")

        # go inside the constant directory
        ampersandIO.printMessage(
            "Creating physical properties and turbulence properties", GUIMode=self.GUIMode, window=self.window)
        # create transportProperties file
        tranP = create_transportPropertiesDict(self.physicalProperties)
        # create turbulenceProperties file
        turbP = create_turbulencePropertiesDict(self.physicalProperties)
        ampersandPrimitives.write_dict_to_file(
            "constant/transportProperties", tranP)
        ampersandPrimitives.write_dict_to_file(
            "constant/turbulenceProperties", turbP)

        # create the controlDict file
        ampersandIO.printMessage(
            "Creating the system files", GUIMode=self.GUIMode, window=self.window)
        controlDict = createControlDict(self.simulationSettings)
        ampersandPrimitives.write_dict_to_file(
            "system/controlDict", controlDict)
        blockMeshDict = generate_blockMeshDict(self.meshSettings)
        ampersandPrimitives.write_dict_to_file(
            "system/blockMeshDict", blockMeshDict)
        snappyHexMeshDict = generate_snappyHexMeshDict(self.meshSettings)
        ampersandPrimitives.write_dict_to_file(
            "system/snappyHexMeshDict", snappyHexMeshDict)
        surfaceFeatureExtractDict = create_surfaceFeatureExtractDict(
            self.meshSettings)
        ampersandPrimitives.write_dict_to_file(
            "system/surfaceFeatureExtractDict", surfaceFeatureExtractDict)
        fvSchemesDict = create_fvSchemesDict(self.numericalSettings)
        ampersandPrimitives.write_dict_to_file(
            "system/fvSchemes", fvSchemesDict)
        fvSolutionDict = create_fvSolutionDict(
            self.numericalSettings, self.solverSettings)
        ampersandPrimitives.write_dict_to_file(
            "system/fvSolution", fvSolutionDict)
        decomposeParDict = createDecomposeParDict(self.parallelSettings)
        ampersandPrimitives.write_dict_to_file(
            "system/decomposeParDict", decomposeParDict)
        FODict = postProcess.create_FOs(
            self.meshSettings, self.postProcessSettings, useFOs=self.useFOs)
        ampersandPrimitives.write_dict_to_file("system/FOs", FODict)

        # create mesh script
        ampersandIO.printMessage(
            "Creating scripts for meshing and running the simulation", GUIMode=self.GUIMode, window=self.window)
        meshScript = ScriptGenerator.generate_mesh_script(
            self.simulationFlowSettings)
        ampersandPrimitives.write_dict_to_file("mesh", meshScript)
        # create simulation script
        simulationScript = ScriptGenerator.generate_simulation_script(
            self.simulationFlowSettings)
        ampersandPrimitives.write_dict_to_file("run", simulationScript)
        ampersandPrimitives.crlf_to_LF("mesh")
        ampersandPrimitives.crlf_to_LF("run")
        if os.name != 'nt':
            os.chmod("mesh", 0o755)
            os.chmod("run", 0o755)
        # go back to the main directory
        os.chdir("..")
        ampersandIO.printMessage(
            "\n-----------------------------------", GUIMode=self.GUIMode, window=self.window)
        ampersandIO.printMessage(
            "Project files created successfully!", GUIMode=self.GUIMode, window=self.window)
        ampersandIO.printMessage(
            "-----------------------------------\n", GUIMode=self.GUIMode, window=self.window)
        return 0


if __name__ == '__main__':
    try:
        from src.cli.create_project import create_project
        create_project()
    except KeyboardInterrupt:
        ampersandIO.printMessage(
            "\nKeyboardInterrupt detected! Aborting project creation")
        exit()
    except Exception as error:
        ampersandIO.printError(error)
        exit()
