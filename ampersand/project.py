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
from primitives import AmpersandUtils, AmpersandIO, FluidPhysicalProperties
from ampersand.models.settings import Domain, Patch, ProjectSettings, BCPatchType, PatchPurpose
from ampersand.utils.stlAnalysis import BoundingBox, STLAnalysis
from ampersand.generators.blockMeshGenerator import create_blockMeshDict
from ampersand.generators.decomposeParGenerator import createDecomposeParDict
from ampersand.snappyHexMeshGenerator import create_snappyHexMeshDict
from ampersand.generators.surfaceExtractor import create_surfaceDict
from ampersand.generators.transportAndTurbulence import create_transportPropertiesDict, create_turbulencePropertiesDict
from boundaryConditionsGenerator import create_boundary_conditions
from ampersand.generators.controlDictGenerator import createControlDict
from ampersand.generators.numericalSettingsGenerator import create_fvSchemesDict, create_fvSolutionDict
from ampersand.generators.scriptGenerator import ScriptGenerator
from ampersand.utils.postProcess import PostProcess
from mod_project import ModProject


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
        self.settings = ProjectSettings()
        self.project_path = ""
        self.stl_files = []  # list to store the settings for stl files
        self.stl_names = []  # list to store the names of the stl files

        self.internalFlow = False  # default is external flow
        self.onGround = False  # default is off the ground
        self.halfModel = False  # default is full model
        self.parallel = True  # default is parallel simulation
        self.snap = True  # default is to use snappyHexMesh
        self.transient = False  # default is steady state
        self.useFOs = False  # default is not to use function objects
        self.current_modification = None  # current modification to the project settings
        # flag to check if the current working directory is the project directory
        self.inside_project_directory = False
        self.mod_options = [
            "Background Mesh", 
            "Add Geometry", 
            "Refinement Levels", 
            "Mesh Point", 
            "Boundary Conditions", 
            "Fluid Properties", 
            "Numerical Settings",
            "Simulation Control Settings", 
            "Turbulence Model", 
            "Post Processing Settings"
        ]
        self.minX, self.maxX, self.minY, self.maxY, self.minZ, self.maxZ = - \
            1e-3, 1e-3, -1e-3, 1e-3, -1e-3, 1e-3
        
        self.domain = Domain()
    # --------------------------------------------------------------------
    # Methods to handle the project summary and changes
    # --------------------------------------------------------------------
    # @property
    # def refinement(self):
    #     return self.meshSettings.fineLevel



# Modification methods
# --------------------------------------------------------------------
    def set_inlet_values(self, U: tuple[float, float, float]):
        if(self.internalFlow): # external flow
            AmpersandIO.printMessage("Setting inlet values for various inlet boundaries")
            for stl_file in self.stl_files:
                if stl_file.purpose == 'inlet':
                    U = stl_file.property
                    self.settings.boundaryConditions.velocityInlet.u_value = U
                    self.settings.inletValues.U = U

        else: # internal flow
            self.settings.inletValues.U = U
            self.settings.boundaryConditions.velocityInlet.u_value = U

        
    def set_fluid_properties(self, fluid: FluidPhysicalProperties):
        self.settings.physicalProperties.rho = fluid.rho
        self.settings.physicalProperties.nu = fluid.nu

    # setting the purpose of a patch. Used for setting the boundary conditions
    def set_purpose(self,patch: Patch, purpose: PatchPurpose='wall'):
        purposes = ['wall', 'inlet','outlet', 'refinementRegion', 'refinementSurface', 
                    'cellZone', 'baffles','symmetry','cyclic','empty']
        if purpose not in purposes:
            AmpersandIO.printMessage("Invalid purpose. Setting purpose to wall")
            purpose = 'wall'
        patch.purpose = purpose


    def set_parallel(self, num_subdomains: int):
        self.settings.parallel.numberOfSubdomains = num_subdomains

    # set the turbulence model for the simulation
    def set_turbulence_model(self,turbulence_model='kOmegaSST'):
        self.settings.physicalProperties.turbulenceModel = turbulence_model

    def set_patch_type(self, patch_name: str, new_type: BCPatchType="patch"):
        if (patch_name not in self.settings.mesh.patches):
            raise ValueError("Patch not found")

        self.settings.mesh.patches[patch_name].type = new_type
     
    def set_ground_type(self, ground_type: bool):
        self.settings.mesh.onGround = ground_type
        if ground_type:
            self.set_patch_type(patch_name='ground',  new_type='wall')

    def set_refinement_level(self,refinement: int):
        self.meshSettings.fineLevel = refinement

    def set_flow_type(self,internalFlow=False):
        self.meshSettings.internalFlow = internalFlow

    def set_half_model(self,halfModel=False):
        self.meshSettings.halfModel = halfModel
        if halfModel:
            self.set_patch_type(patch_name='back',new_type='symmetry')


    # TODO: this was always True to begin with
    def set_post_process_settings(self):
        if self.useFOs:
            self.settings.postProcess.FOs = True
        meshPoint = list(self.meshSettings.castellatedMeshControls.locationInMesh)
        self.settings.postProcess.massFlow = True
        self.settings.postProcess.minMax = True
        self.settings.postProcess.yPlus = True
        self.settings.postProcess.forces = True
        # the default probe location for monitoring of flow variables
        self.settings.postProcess.probeLocations.append(meshPoint)
      
# --------------------------------------------------------------------


    def update_domain_size(self, domain_size: BoundingBox, nx: Optional[int] = None, ny: Optional[int] = None, nz: Optional[int] = None):
        self.domain = Domain.update(self.meshSettings.domain, domain_size, nx, ny, nz)

    def update_domain_size_from_stls(self):
            # adjust the domain size based on the bounding box of the stl files
            AmpersandIO.printMessage("Update domain size based on the bounding box of the stl files",GUIMode=self.GUIMode,window=self.window)
            for stl_file in self.stl_files:
                stl_name = stl_file.name
                stl_path = os.path.join(self.project_path, "constant", "triSurface", stl_name)
                stlBoundingBox = STLAnalysis.compute_bounding_box(stl_path)
                self.update_domain_size(stlBoundingBox)
                #self.meshSettings = stlAnalysis.set_mesh_location(self.meshSettings, stl_path,self.internalFlow)
            


    def update_boundary_condition(self, patch_name: str, new_bc_type: BCPatchType, new_property: Optional[tuple] = None):
        if patch_name in self.settings.mesh.patches:
            self.settings.mesh.patches[patch_name].type = new_bc_type
            self.settings.mesh.patches[patch_name].purpose = new_bc_type
            self.settings.mesh.patches[patch_name].property = new_property
            AmpersandIO.printMessage(f"Boundary condition {patch_name} changed to {new_bc_type}", GUIMode=self.GUIMode, window=self.window)
        elif patch_name in self.meshSettings.geometry:
            self.meshSettings.geometry[patch_name].purpose = new_bc_type
            self.meshSettings.geometry[patch_name].property = new_property
            AmpersandIO.printMessage(f"Boundary condition of {patch_name} changed to {new_bc_type}", GUIMode=self.GUIMode, window=self.window)
        else:
            raise ValueError("Boundary condition not found")

    def change_stl_refinement_level(self,stl_file_number=0):
        AmpersandIO.printMessage("Changing refinement level")
        refMin = AmpersandIO.get_input_int("Enter new refMin: ")
        refMax = AmpersandIO.get_input_int("Enter new refMax: ")
        self.stl_files[stl_file_number].refineMin = refMin
        self.stl_files[stl_file_number].refineMax = refMax
        #stl_name = project.stl_files[stl_file_number].name
        
        fileFound = False
        for stl in self.meshSettings.geometry:
            if stl.name == self.stl_files[stl_file_number].name:
                fileFound = True
                stl.refineMin = refMin
                stl.refineMax = refMax
                stl.featureLevel = refMax
                break
        if not fileFound:
            AmpersandIO.printMessage("STL file not found in the geometry list")
        
        #return project

    

    def modify_project(self):
        if self.current_modification=="Background Mesh":
            ModProject.change_background_mesh(self)
        elif self.current_modification=="Mesh Point":
            ModProject.change_mesh_point(self)
        elif self.current_modification=="Add Geometry":
            ModProject.add_geometry(self)
        elif self.current_modification=="Refinement Levels":
            ModProject.change_refinement_levels(self)
        elif self.current_modification=="Boundary Conditions":
            ModProject.change_boundary_conditions(self)
        elif self.current_modification=="Fluid Properties":
            ModProject.change_fluid_properties(self)
        elif self.current_modification=="Numerical Settings":
            ModProject.change_numerical_settings(self)
        elif self.current_modification=="Simulation Control Settings":
            ModProject.change_simulation_settings(self)
        elif self.current_modification=="Turbulence Model":
            ModProject.change_turbulenc_model(self)
        elif self.current_modification=="Post Processing Settings":
            ModProject.change_post_process_settings(self)
        else:
            AmpersandIO.printMessage("Invalid option. Aborting operation")
            return -1
        return False
    
    
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
        self.meshSettings.geometry = AmpersandUtils.remove_duplicate_dicts(self.meshSettings.geometry)
        #print("stl_files",self.stl_files)
        #print("Mesh settings",self.meshSettings.geometry)

    def set_project_directory(self, project_directory_path: str):
        if self.GUIMode:
            stopWhenError = False
        else:
            stopWhenError = True
        if project_directory_path is None:
            if stopWhenError:
                AmpersandIO.printMessage("No directory selected. Aborting project creation.")
                exit()
            else:
                return -1
        #assert os.path.exists(project_directory_path), "The chosen directory does not exist"
        if not os.path.exists(project_directory_path):
            if stopWhenError:
                AmpersandIO.printMessage("The chosen directory does not exist. Aborting project creation.")
                exit()
            else:
                self.project_directory_path = None
                return -1
        self.project_directory_path = project_directory_path
        
    # To create the project path for a new project with the project name
    def create_project_path(self):
        if not self.project_directory_path:
            AmpersandIO.printWarning("No directory selected. Aborting project creation.")
            return -1
        self.project_path = os.path.join(self.project_directory_path, self.project_name)
    
    # this is to set the project path if the project is already existing
    # useful for opening existing projects and modifying the settings
    def set_project_path(self,project_path):
        if project_path is None:
            if self.GUIMode==False:
                AmpersandIO.printWarning("No project path selected. Aborting project creation.",GUIMode=self.GUIMode)   
            #ampersandIO.printWarning("No project path selected. Aborting project creation/modification.",GUIMode=self.GUIMode)
            return -1
            #exit()
        if os.path.exists(project_path):
            settings_file = os.path.join(project_path, "project_settings.yaml")
            if os.path.exists(settings_file):
                AmpersandIO.printMessage("Project found, loading project settings",GUIMode=self.GUIMode,window=self.window)
                self.existing_project = True
                self.project_path = project_path
                return False
            else:
                if self.GUIMode==False:
                    AmpersandIO.printWarning("Settings file not found. Please open an Ampersand case directory.",GUIMode=self.GUIMode)
                #ampersandIO.printWarning("Settings file not found. Please open an Ampersand case directory.",GUIMode=self.GUIMode)
                # TO DO: Add the code socket to create a new project here
                return -1
        else:
            if self.GUIMode==False:
                AmpersandIO.printWarning("Project path does not exist. Aborting project creation/opening.",GUIMode=self.GUIMode)
            #ampersandIO.printWarning("Project path does not exist. Aborting project creation/opening.",GUIMode=self.GUIMode)
            return -1



    # Create the project directory in the specified location.
    # 0, constant, system, constant/triSurface directories are created.
    def create_project(self):
        # check if the project path exists
        if self.project_path is None:
            AmpersandIO.printError("No project path selected. Aborting project creation.",GUIMode=self.GUIMode)
            return -1
        if os.path.exists(self.project_path):
            AmpersandIO.printWarning("Project already exists. Skipping the creation of directories",GUIMode=self.GUIMode)
            self.existing_project = True
        else:
            AmpersandIO.printMessage("Creating project directory")
            try:
                os.makedirs(self.project_path)
                
            except OSError as error:
                AmpersandIO.printError(error)
        try:
            os.chdir(self.project_path)
        except OSError as error:
                AmpersandIO.printError(error)
        cwd = os.getcwd()
        AmpersandIO.printMessage(f"Working directory: {cwd}",GUIMode=self.GUIMode,window=self.window)

        # create 0, constant and system directory
        try:
            os.mkdir("0")
            os.mkdir("constant")
            os.mkdir("system")
            os.mkdir("constant/triSurface")
        except OSError as error:
            AmpersandIO.printError("File system already exists. Skipping the creation of directories",GUIMode=self.GUIMode)   
            return -1
        return False # return False if the project is created successfully 

    # write current settings to the project_settings.yaml file inside the project directory
    def write_settings(self):
        settings = self.settings.model_dump()
        AmpersandIO.printMessage("Writing settings to project_settings.yaml",GUIMode=self.GUIMode,window=self.window)
        AmpersandUtils.dict_to_yaml(settings, 'project_settings.yaml')

    # If the project is already existing, load the settings from the project_settings.yaml file
    def load_settings(self):
        AmpersandIO.printMessage("Loading project settings",GUIMode=self.GUIMode,window=self.window)
        self.settings = ProjectSettings.model_validate(AmpersandUtils.yaml_to_dict('project_settings.yaml'))
        for patch_name, patch in self.meshSettings.geometry.items():
            if(patch.type=='triSurfaceMesh'):
                if(patch_name in self.stl_names):
                    AmpersandIO.printMessage(f"STL file {patch_name} already exists in the project, skipping the addition")
                else:
                    self.stl_files.append(patch)
                    self.stl_names.append(patch_name)

        self.project_name = self.project_path.split("/")[-1]
        # treat bounds as tuple
        self.meshSettings.geometry = AmpersandUtils.treat_bounds(self.meshSettings.geometry)


    def setup_default_settings(self):
        self.settings = ProjectSettings()
        self.write_settings()

    # Create the settings for the project or load the existing settings
    def create_settings(self):
        if self.existing_project:
            AmpersandIO.printMessage("Project already exists. Loading project settings")
            try:
                self.load_settings()
            except FileNotFoundError:
                AmpersandIO.printMessage("Settings file not found. Setting up default settings")
                self.setup_default_settings()
        else:
            self.setup_default_settings()

    # Add a stl file to the project settings (self.meshSettings)
    def add_stl_to_mesh_settings(self, stl_name,refMin=0, refMax=0, featureEdges=True, 
                                 featureLevel=1,purpose='wall',property=None,bounds=None, nLayers=3):
        already_exists = False # flag to check if the stl file already exists in the project
        idx = 0 # index of the stl file in the list
        # Purpose is wall by default
        # Other purposes are patch, refinementRegion, refinementSurface, cellZone, baffles
        
        if self.refinement == 0:
            nLayers = 3
        elif self.refinement == 1:
            nLayers = 5
        else:
            nLayers = 7
        
        stl_ = {'name': stl_name, 'type':'triSurfaceMesh','purpose':purpose, 'refineMin': refMin, 'refineMax': refMax, 
                'featureEdges':featureEdges, 'featureLevel':featureLevel, 'nLayers':nLayers, 'property':property, 'bounds':bounds}
        
        # this snippet is to prevent the same stl file from being added multiple times
        # Instead of using set, we are using a list to store the stl files
        # To check if the stl file already exists in the project, stl names are compared.
        # list all the stl files in the project
        for stl_name in self.stl_names:
            if stl_name == stl_.name:
                AmpersandIO.printMessage(f"STL file {stl_name} already exists in the project")
                already_exists = True
                break
            idx += 1
        if not already_exists:
            self.stl_names.append(stl_name)
            self.stl_files.append(stl_)
        else:    
            self.stl_files[idx] = stl_
            return True # flag to indicate that the stl file is replaced
        return False

    

    def add_stl_to_project(self):
        for stl_file in self.stl_files:
            self.meshSettings.geometry.append(stl_file)
        self.remove_duplicate_stl_files()

    def add_stl_file(self): # to only copy the STL file to the project directory and add it to the STL list
        stl_file = AmpersandUtils.ask_for_file([("STL Geometry", "*.stl"), ("OBJ Geometry", "*.obj")],self.GUIMode)
        if stl_file is None:
            AmpersandIO.printWarning("No file selected. Please select STL file if necessary.",GUIMode=self.GUIMode)
            return -1
        if os.path.exists(stl_file):
            # add the stl file to the project
            # This is a bit confusing. 
            # stl_name is the name of the file, stl_file is the path to the file
            file_path_to_token = stl_file.split("/")
            stl_name = file_path_to_token[-1]
            if stl_name in self.stl_names:
                AmpersandIO.printWarning(f"STL file {stl_name} already exists in the project",GUIMode=self.GUIMode)
                return -1
            else: # this is to prevent the bug of having the same file added multiple times
                if self.GUIMode:
                    purpose = "wall"
                    property = None
                else:
                    purpose = self.ask_purpose()
                    property = self.ask_property(purpose)
                bounds = STLAnalysis.compute_bounding_box(stl_file)
                AmpersandIO.printMessage(f"Bounds of the geometry: {bounds.to_tuple()}",GUIMode=self.GUIMode,window=self.window)
                if purpose == 'refinementRegion' or purpose == 'refinementSurface':
                    featureEdges = False
                else:  
                    featureEdges = True
                self.add_stl_to_mesh_settings(stl_name,purpose=purpose,property=property,featureEdges=featureEdges,bounds=bounds)
            # this is the path to the constant/triSurface inside project directory where STL will be copied
            stl_path = os.path.join(self.project_path, "constant", "triSurface", stl_name)
            try:
                AmpersandIO.printMessage(f"Copying {stl_name} to the project directory",GUIMode=self.GUIMode,window=self.window)
                shutil.copy(stl_file, stl_path)
            except OSError as error:
                AmpersandIO.printError(error,GUIMode=self.GUIMode)
                return -1
            try:
                STLAnalysis.set_stl_solid_name(stl_path)
            except Exception as error:
                AmpersandIO.printError(error,GUIMode=self.GUIMode)
                return -1
        else:
            AmpersandIO.printError("File does not exist. Aborting project creation.",GUIMode=self.GUIMode)
            return -1
        self.current_stl_file = stl_path
        return False
            
    # this is a wrapper of the primitives 
    def list_stl_files(self):
        AmpersandUtils.list_stl_files(self.stl_files,self.GUIMode,self.window)

    def list_stl_paths(self):
        stl_paths = []
        for stl_file in self.stl_files:
            stl_path = os.path.join(self.project_path, "constant", "triSurface", stl_file.name)
            #ampersandIO.printMessage(stl_path)
            stl_paths.append(stl_path)
        return stl_paths


    def remove_stl_file(self):
        #self.list_stl_files()
        stl_file_number = AmpersandIO.get_input("Enter the number of the file to remove: ")
        try:
            stl_file_number = int(stl_file_number)
        except ValueError:
            AmpersandIO.printMessage("Invalid input. Aborting operation")
            return -1
        if stl_file_number < 0 or stl_file_number > len(self.stl_files):
            AmpersandIO.printMessage("Invalid file number. Aborting operation")
            return -1
        stl_file = self.stl_files[stl_file_number]
        stl_name = stl_file.name
        self.stl_files.remove(stl_file)
        self.stl_names.remove(stl_name)
        stl_path = os.path.join(self.project_path, "constant", "triSurface", stl_name)
        try:
            os.remove(stl_path)
        except OSError as error:
            AmpersandIO.printError(error)
            return -1
        return False



    def analyze_stl_file(self,stl_file_number=0):
        rho = self.settings.physicalProperties.rho
        nu = self.settings.physicalProperties.nu
        U = max(self.settings.inletValues.U)
        ER = self.meshSettings.addLayersControls.expansionRatio
        try:
            stl_file_number = int(stl_file_number)
        except ValueError:
            AmpersandIO.printError("Invalid input. Aborting operation",GUIMode=self.GUIMode)
            return -1
        if stl_file_number < 0 or stl_file_number > len(self.stl_files):
            AmpersandIO.printError("Invalid file number. Aborting operation",GUIMode=self.GUIMode)
            return -1
        stl_file = self.stl_files[stl_file_number]
        stl_name = stl_file.name
        AmpersandIO.printMessage(f"Analyzing {stl_name}",GUIMode=self.GUIMode,window=self.window)
        stl_path = os.path.join(self.project_path, "constant", "triSurface", stl_name)
        stlBoundingBox = STLAnalysis.compute_bounding_box(stl_path)
        domain_size, nx, ny, nz, refLevel,target_y,nLayers = STLAnalysis.calc_mesh_settings(stlBoundingBox, nu, rho,U=U,maxCellSize=2.0,expansion_ratio=ER,
                                                                           onGround=self.onGround,internalFlow=self.internalFlow,
                                                                           refinement=self.refinement,halfModel=self.halfModel,
                                                                           GUI=self.GUIMode,window=self.window)
        featureLevel = max(refLevel,1)
        self.meshSettings = STLAnalysis.set_mesh_settings(self.meshSettings, domain_size, nx, ny, nz, refLevel, featureLevel,nLayers=nLayers) 
        self.update_domain_size(domain_size,nx,ny,nz)
        self.meshSettings = STLAnalysis.set_mesh_location(self.meshSettings, stl_path,self.internalFlow)
        refinementBoxLevel = max(2,refLevel-3)
        self.meshSettings = STLAnalysis.addRefinementBoxToMesh(meshSettings=self.meshSettings, stl_path=stl_path,refLevel=refinementBoxLevel,internalFlow=self.internalFlow)
        if(self.internalFlow==False and self.onGround==True):
            # if the flow is external and the geometry is on the ground, add a ground refinement box
            self.meshSettings = STLAnalysis.addGroundRefinementBoxToMesh(meshSettings=self.meshSettings, stl_path=stl_path,refLevel=refinementBoxLevel)
        self.meshSettings = STLAnalysis.set_layer_thickness(self.meshSettings, 0.5) # set the layer thickness to 0.5 times the cell size
        # store the background mesh size for future reference
        maxCellSize = abs((domain_size.maxx-domain_size.minx)/nx)
        self.settings.mesh.maxCellSize = maxCellSize
        #self.meshSettings = stlAnalysis.set_min_vol(self.meshSettings, minVol)
        return False
    

    # def get_probe_location(self):
    #     point = PostProcess.get_probe_location()
    #     # for internal flows, the point should be inside stl
    #     # for external flows, the point should be outside stl
    #     # TO DO
    #     self.postProcessSettings.probeLocations.append(point)


    # def summarize_boundary_conditions(self):
    #     return AmpersandUtils.list_boundary_conditions(self.meshSettings)


    # def summarize_project(self):
    #     trueFalse = {True: 'Yes', False: 'No'}
    #     AmpersandIO.show_title("Project Summary")
    #     #ampersandIO.printMessage(f"Project directory: {self.project_directory_path}")
    #     AmpersandIO.printFormat("Project name", self.project_name, GUIMode=self.GUIMode,window=self.window)
    #     AmpersandIO.printFormat("Project path", self.project_path, GUIMode=self.GUIMode,window=self.window)
        
    #     AmpersandIO.printMessage(f"Internal Flow: {trueFalse[self.internalFlow]}",GUIMode=self.GUIMode,window=self.window)
    #     if(self.internalFlow==False):
    #         AmpersandIO.printMessage(f"On Ground: {trueFalse[self.onGround]}",GUIMode=self.GUIMode,window=self.window)
    #     AmpersandIO.printMessage(f"Transient: {trueFalse[self.transient]}",GUIMode=self.GUIMode,window=self.window)
    #     self.summarize_background_mesh()
    #     self.list_stl_files()
        

    # # this will show the details of the background mesh
    # def summarize_background_mesh(self):
    #     AmpersandIO.printMessage(
    #       f"{self.settings.mesh.domain.__repr__()}\n"
    #       f"Background cell size: {self.settings.mesh.maxCellSize} m"
    #     )



    # def check_project_path(self): # check if the project path exists and if the project is already existing
    #     if os.path.exists(self.project_path):
    #         settings_file = os.path.join(self.project_path, "project_settings.yaml")
    #         if os.path.exists(settings_file):
    #             AmpersandIO.printWarning("Project already exists, loading project settings",GUIMode=self.GUIMode)
    #             self.existing_project = True
    #             return False
    #         else:
    #             self.existing_project = False
    #             return -1
    #     else:
    #         self.existing_project = False
    #         return -1
        


    # # Check if the 0 directory exists in the project directory
    # def check_0_directory(self):
    #     if not os.path.exists("0"):
    #         AmpersandIO.printWarning("0 directory does not exist.",GUIMode=self.GUIMode)
    #         AmpersandIO.printMessage("Checking for 0.orig directory",GUIMode=self.GUIMode,window=self.window)
    #         if os.path.exists("0.orig"):
    #             AmpersandIO.printMessage("0.orig directory found. Copying to 0 directory",GUIMode=self.GUIMode,window=self.window)
    #             shutil.copytree("0.orig", "0")
    #         else:
    #             AmpersandIO.printWarning("0.orig directory not found. Aborting project creation.",GUIMode=self.GUIMode)
    #             return -1
    #     return False
    
    # # Check if the constant directory exists in the project directory
    # def check_constant_directory(self):
    #     if not os.path.exists("constant"):
    #         AmpersandIO.printWarning("constant directory does not exist.",GUIMode=self.GUIMode)
    #         return -1
    #     return False
    
    # # Check if the system directory exists in the project directory
    # def check_system_directory(self):
    #     if not os.path.exists("system"):
    #         AmpersandIO.printWarning("system directory does not exist.",GUIMode=self.GUIMode)
    #         #ampersandIO.printError("system directory is necessary for the project")
    #         return -1
    #     return False
    
    # # Check if the constant/triSurface directory exists in the project directory
    # def check_triSurface_directory(self):
    #     if not os.path.exists("constant/triSurface"):
    #         AmpersandIO.printWarning("triSurface directory does not exist.",GUIMode=self.GUIMode)
    #         #ampersandIO.printError("triSurface directory is necessary for the project")
    #         return -1
    #     # if exists, check if the stl files are present
    #     stl_files = os.listdir("constant/triSurface")

    # # to check whether log files are present in the project directory
    # def check_log_files(self):
    #     log_files = os.listdir()
    #     if 'log.simpleFoam' in log_files:
    #         AmpersandIO.printMessage("Simulation log file found",GUIMode=self.GUIMode,window=self.window)
    #         return True
    #     if 'log.pimpleFoam' in log_files:
    #         AmpersandIO.printMessage("Simulation log file found",GUIMode=self.GUIMode,window=self.window)
    #         return True
    #     return False
    
    # # to check whether the U and p files are present in the postProcess directory
    # def check_post_process_files(self):
    #     if(not os.path.exists("postProcessing/probe/0")):
    #         AmpersandIO.printWarning("postProcess directory does not exist",GUIMode=self.GUIMode)
    #         return False
    #     postProcess_files = os.listdir("postProcessing/probe/0")
    #     if 'U' in postProcess_files and 'p' in postProcess_files:
    #         AmpersandIO.printMessage("U and p files found in postProcess directory",GUIMode=self.GUIMode,window=self.window)
    #         return True
    #     return False
    
    # def check_forces_files(self):
    #     if(not os.path.exists("postProcessing/forces/0")):
    #         AmpersandIO.printWarning("forces directory does not exist", GUIMode=self.GUIMode)
    #         return False
    #     forces_files = os.listdir("postProcessing/forces/0")
    #     if 'force.dat' in forces_files:
    #         AmpersandIO.printMessage("force.dat found in forces directory")
    #         return True
    #     return False






    # Wrapper of cwd with error handling
    # def go_inside_directory(self):
    #     try:
    #         os.chdir(self.project_path)
    #     except OSError as error:
    #             AmpersandIO.printError(error)
    #     cwd = os.getcwd()
    #     AmpersandIO.printMessage(f"Working directory: {cwd}",GUIMode=self.GUIMode,window=self.window)
    #     self.inside_project_directory = True


    #def set_transient(self):
    #    self.transient = ampersandIO.get_input_bool("Transient simulation (y/N)?: ")


        
    # def show_settings(self):
    #     AmpersandIO.printMessage("Project settings")
    #     AmpersandIO.printMessage("Mesh Settings")
    #     AmpersandIO.print_dict(self.meshSettings)
    #     AmpersandIO.printMessage("Physical Properties")
    #     AmpersandIO.print_dict(self.physicalProperties)
    #     AmpersandIO.printMessage("Numerical Settings")
    #     AmpersandIO.print_dict(self.numericalSettings)
    #     AmpersandIO.printMessage("Inlet Values")
    #     AmpersandIO.print_dict(self.inletValues)
    #     AmpersandIO.printMessage("Boundary Conditions")
    #     AmpersandIO.print_dict(self.boundaryConditions)
    #     AmpersandIO.printMessage("Solver Settings")
    #     AmpersandIO.print_dict(self.solverSettings)
    #     AmpersandIO.printMessage("Simulation Settings")
    #     AmpersandIO.print_dict(self.simulationSettings)
    #     AmpersandIO.printMessage("Parallel Settings")
    #     AmpersandIO.print_dict(self.parallelSettings)
    #     AmpersandIO.printMessage("Simulation Flow Settings")
    #     AmpersandIO.print_dict(self.simulationFlowSettings)
    #     AmpersandIO.printMessage("Post Process Settings")
    #     AmpersandIO.print_dict(self.postProcessSettings)


    # def create_project_files(self):
    #     #(meshSettings, physicalProperties, numericalSettings, inletValues, boundaryConditions)=caseSettings
    #     # check if the current working directory is the project directory
    #     if not os.path.exists(self.project_path):
    #         AmpersandIO.printMessage("Project directory does not exist. Aborting project creation.")
    #         return -1
    #     if os.getcwd() != self.project_path:
    #         os.chdir(self.project_path)

    #     # Remove the existing 0.orig directory if it exists.
    #     # This is to prevent the error of copying the old 0.orig directory to 0 directory
    #     if os.path.exists("0.orig"):
    #         shutil.rmtree("0.orig")
    #     # create the initial conditions file
    #     AmpersandIO.printMessage("Creating boundary conditions")
    #     # check if the 0 directory exists
    #     if not os.path.exists("0"):
    #         # create the 0 directory
    #         os.mkdir("0")
    #     # go inside the 0 directory
    #     os.chdir("0")
    #     create_boundary_conditions(self.meshSettings, self.boundaryConditions)    
    #     # go back to the main directory 
    #     os.chdir("..")
    #     # go inside the constant directory
    #     os.chdir("constant")
    #     AmpersandIO.printMessage("Creating physical properties and turbulence properties",GUIMode=self.GUIMode,window=self.window)
    #     # create transportProperties file
    #     tranP = create_transportPropertiesDict(self.physicalProperties)
    #     # create turbulenceProperties file
    #     turbP = create_turbulencePropertiesDict(self.physicalProperties)
    #     AmpersandUtils.write_dict_to_file("transportProperties", tranP)
    #     AmpersandUtils.write_dict_to_file("turbulenceProperties", turbP)
    #     # go back to the main directory
    #     os.chdir("..")
        
    #     # go inside the system directory
    #     os.chdir("system")
    #     # create the controlDict file
    #     AmpersandIO.printMessage("Creating the system files",GUIMode=self.GUIMode,window=self.window)
    #     controlDict = createControlDict(self.simulationSettings)
    #     AmpersandUtils.write_dict_to_file("controlDict", controlDict)
    #     blockMeshDict = create_blockMeshDict(self.meshSettings)
    #     AmpersandUtils.write_dict_to_file("blockMeshDict", blockMeshDict)
    #     snappyHexMeshDict = create_snappyHexMeshDict(self.meshSettings)
    #     AmpersandUtils.write_dict_to_file("snappyHexMeshDict", snappyHexMeshDict)
        
    #     # TODO: fix this for surfaceExtracts
    #     surfaceFeatureExtractDict = create_surfaceDict(self.meshSettings, "surfaceFeatureExtractDict")
    #     AmpersandUtils.write_dict_to_file("surfaceFeatureExtractDict", surfaceFeatureExtractDict)

            

    #     fvSchemesDict = create_fvSchemesDict(self.numericalSettings)
    #     AmpersandUtils.write_dict_to_file("fvSchemes", fvSchemesDict)
    #     fvSolutionDict = create_fvSolutionDict(self.numericalSettings, self.solverSettings)
    #     AmpersandUtils.write_dict_to_file("fvSolution", fvSolutionDict)
    #     decomposeParDict = createDecomposeParDict(self.parallelSettings)
    #     AmpersandUtils.write_dict_to_file("decomposeParDict", decomposeParDict)
    #     FODict = PostProcess.create_FOs(self.meshSettings,self.postProcessSettings,useFOs=self.useFOs)
    #     AmpersandUtils.write_dict_to_file("FOs", FODict)
    #     # go back to the main directory
    #     os.chdir("..")
    #     # create mesh script
    #     AmpersandIO.printMessage("Creating scripts for meshing and running the simulation",GUIMode=self.GUIMode,window=self.window)
    #     meshScript = ScriptGenerator.create_mesh_script(self.simulationFlowSettings)
    #     AmpersandUtils.write_dict_to_file("mesh", meshScript)
    #     # create simulation script
    #     simulationScript = ScriptGenerator.create_simulation_script(self.simulationFlowSettings)
    #     AmpersandUtils.write_dict_to_file("run", simulationScript)
    #     AmpersandUtils.crlf_to_LF("mesh")
    #     AmpersandUtils.crlf_to_LF("run")
    #     if os.name != 'nt':
    #         os.chmod("mesh", 0o755)
    #         os.chmod("run", 0o755)
    #     # go back to the main directory
    #     os.chdir("..")
    #     AmpersandIO.printMessage("\n-----------------------------------",GUIMode=self.GUIMode,window=self.window)
    #     AmpersandIO.printMessage("Project files created successfully!",GUIMode=self.GUIMode,window=self.window)
    #     AmpersandIO.printMessage("-----------------------------------\n",GUIMode=self.GUIMode,window=self.window)
    #     return False



# def main():
#     project = AmpersandProject()
#     # Clear the screen
#     os.system('cls' if os.name == 'nt' else 'clear')
#     AmpersandIO.printMessage(get_ampersand_header())
#     project.set_project_directory(AmpersandUtils.ask_for_directory())
#     project_name = AmpersandIO.get_input("Enter the project name: ")
#     project.set_project_name(project_name)
#     #user_name = input("Enter the user name: ")
#     #project.set_user_name(user_name)
#     project.create_project_path()
#     AmpersandIO.printMessage("Creating the project")
#     AmpersandIO.printMessage(f"Project path: {project.project_path}")
#     #project.project_path = r"C:\Users\Ridwa\Desktop\CFD\ampersandTests\drivAer2"
#     project.create_project()
#     project.create_settings()
#     yN = AmpersandIO.get_input("Add STL file to the project (y/N)?: ")
#     while yN.lower() == 'y':
#         project.add_stl_file()
#         yN = AmpersandIO.get_input("Add another STL file to the project (y/N)?: ")
#     project.add_stl_to_project()
#     # Before creating the project files, the settings are flushed to the project_settings.yaml file
#     project.list_stl_files()
#     project.ask_flow_type()
#     if(project.internalFlow!=True):
#         project.ask_ground_type()
#     if(len(project.stl_files)>0):
#         project.analyze_stl_file()

#     #project.analyze_stl_file()
#     project.write_settings()
#     project.create_project_files()

# if __name__ == '__main__':
#     # Specify the output YAML file
#     try:
#         main()
#     except KeyboardInterrupt:
#         AmpersandIO.printMessage("\nKeyboardInterrupt detected! Aborting project creation")
#         exit()
#     except Exception as error:
#         AmpersandIO.printError(error)
#         exit()
