from ampersand.primitives import AmpersandIO, AmpersandUtils
from ampersand.project import AmpersandProject


class ProjectIOUtils:
    @staticmethod
    def ask_transient_settings(project: AmpersandProject):
        self.ask_transient()
        if self.transient:
            AmpersandIO.printMessage("Transient simulation settings")
            self.simulationSettings.transient = True
            self.simulationSettings.application = 'pimpleFoam'
            self.simulationFlowSettings.solver = 'pimpleFoam'
            self.simulationSettings.endTime = AmpersandIO.get_input_float("End time: ")
            self.simulationSettings.writeInterval = AmpersandIO.get_input_float("Write interval: ")
            self.simulationSettings.deltaT = AmpersandIO.get_input_float("Time step: ")
            self.simulationSettings.adjustTimeStep = 'no'
            self.simulationSettings.maxCo = 0.9
            self.numericalSettings.ddtSchemes.default = 'Euler'
            # if steady state, SIMPLEC is used. If transient, PIMPLE is used
            # for PIMPLE, the relaxation factors are set to 0.7 and p = 0.3
            self.numericalSettings.relaxationFactors.p = 0.3
            
    @staticmethod
    def ask_purpose(project: AmpersandProject):
        purposes = ['wall', 'inlet', 'outlet', 'refinementRegion', 'refinementSurface', 'cellZone', 'baffles','symmetry','cyclic','empty',]
        AmpersandIO.printMessage(f"Enter purpose for this STL geometry")
        AmpersandIO.print_numbered_list(purposes)
        purpose_no = AmpersandIO.get_input_int("Enter purpose number: ")-1
        if(purpose_no < 0 or purpose_no > len(purposes)-1):
                AmpersandIO.printMessage("Invalid purpose number. Setting purpose to wall")
                purpose = 'wall'
        else:
            purpose = purposes[purpose_no]
        return purpose
    
    @staticmethod
    def ask_property(project: AmpersandProject, purpose: PatchPurpose='wall'):
        if purpose == 'inlet':
            U = ampersandDataInput.get_inlet_values()
            property = tuple(U)
            AmpersandIO.printMessage(f"Setting property of {purpose} to {property}")
        elif purpose == 'refinementRegion' :
            refLevel = AmpersandIO.get_input_int("Enter refinement level: ")
            property = refLevel
        elif purpose == 'cellZone':
            refLevel = AmpersandIO.get_input_int("Enter refinement level: ")
            createPatches = AmpersandIO.get_input_bool("Create patches for this cellZone? (y/N): ")
            property = (refLevel, createPatches,0) # 0 is just a placeholder for listing the patches
        elif purpose == 'refinementSurface':
            refLevel = AmpersandIO.get_input_int("Enter refinement level: ")
            property = refLevel
        else:
            property = None
        return property
    @staticmethod 
    def ask_boundary_type(project: AmpersandProject):
        bcTypes = ["inlet","outlet","wall","symmetry","cyclic","empty","movingWall",]
        AmpersandIO.printMessage("List of boundary types")
        AmpersandIO.print_numbered_list(bcTypes)
        bcType = AmpersandIO.get_input_int("Enter the number of the boundary type: ")
        if bcType <= 0 or bcType > len(bcTypes):
            AmpersandIO.printMessage("Invalid boundary type. Setting to wall")
            return "wall"
        return bcTypes[bcType-1]



    @staticmethod
    def ask_stl_settings(project: AmpersandProject,stl_file):
        AmpersandIO.printMessage(f"Settings of the {stl_file.name} file")
        stl_file.refineMin = AmpersandIO.get_input("Min Refinement: ")
        stl_file.refineMax = AmpersandIO.get_input("Max Refinement: ")
        featureEdges = AmpersandIO.get_input("Refine Feature Edges?: (y/N) ")
        if(featureEdges == 'y'):
            stl_file.featureEdges = True
        else:    
            stl_file.featureEdges = False
        stl_file.featureLevel = AmpersandIO.get_input("Feature Level: ")
        stl_file.nLayers = AmpersandIO.get_input("Number of Layers: ")

    @staticmethod
    def assign_ground_type(project: AmpersandProject):
        ground_type = AmpersandIO.get_input_bool("Is the ground touching the body (y/N): ")
        project.set_ground_type(ground_type)

    @staticmethod
    def assign_refinement_level(project: AmpersandProject):
        refLevel = AmpersandIO.get_input_int("Enter the mesh refinement (0: coarse, 1: medium, 2: fine): ")
        if refLevel not in [0, 1, 2]:
            AmpersandIO.printMessage("Invalid mesh refinement level. Defaulting to medium.")
            refLevel = 1
        project.set_refinement_level(refLevel)



    def ask_flow_type(self):
        flow_type = AmpersandIO.get_input("Internal or External Flow (I/E)?: ")
        if flow_type.lower() == 'i':
            self.internalFlow = True
        else:
            self.internalFlow = False
        self.meshSettings.internalFlow = self.internalFlow

    def ask_transient(self):
        transient = AmpersandIO.get_input("Transient or Steady State (T/S)?: ")
        if transient.lower() == 't':
            self.transient = True
        else:
            self.transient = False


    # choose turbulence model for the simulation
    @staticmethod
    def choose_turbulence_model(project: AmpersandProject):
        turbulence_models = ['kOmegaSST', 'kEpsilon', 'SpalartAllmaras']
        turbulence_model = AmpersandDataInput.get_option_choice("Choose turbulence model: ", turbulence_models)
        self.solverSettings.turbulenceModel = turbulence_model

    @staticmethod
    def assign_parallel(project: AmpersandProject):
        n_core = AmpersandIO.get_input_int("Number of cores for parallel simulation: ")


    @staticmethod
    def assign_fluid(project: AmpersandProject):
        fluid = AmpersandDataInput.choose_fluid_properties()
        if fluid == -1:
            fluid = AmpersandDataInput.get_physical_properties()


    @staticmethod
    def ask_half_model(project: AmpersandProject):
        half_model = AmpersandIO.get_input_bool("Half Model (y/N)?: ")
        project.set_half_model(half_model)

    
    def choose_modification(self):
        current_modification = AmpersandIO.get_option_choice(prompt="Choose any option for project modification: ",
                                      options=self.mod_options,title="\nModify Project Settings")
        self.current_modification = self.mod_options[current_modification]
        AmpersandIO.printMessage(f"Current modification: {self.current_modification}",GUIMode=self.GUIMode,window=self.window)
        
    def choose_modification_categorized(self):
        options = ['Mesh','Boundary Conditions','Fluid Properties','Numerical Settings','Simulation Control Settings','Turbulence Model','Post Processing Settings']
        current_modification = AmpersandIO.get_option_choice(prompt="Choose any option for project modification: ",
                                      options=options,title="\nModify Project Settings")
        mesh_options = ['Background Mesh','Mesh Point','Add Geometry','Refinement Levels']
        
        if current_modification < 0 or current_modification > len(options)-1:
            AmpersandIO.printMessage("Invalid option. Aborting operation")
            return -1
        if current_modification == 0:
            self.current_modification = mesh_options[AmpersandIO.get_option_choice(prompt="Choose any option for mesh modification: ",
                                      options=mesh_options,title="\nModify Mesh Settings")]
        else:
            self.current_modification = options[current_modification]
        

