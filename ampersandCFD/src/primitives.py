import os
import yaml
import sys
from tkinter import filedialog, Tk

class ampersandPrimitives:
    def __init__(self):
        pass

    @staticmethod
    def crlf_to_LF(file_path):
        WINDOWS_LINE_ENDING = b'\r\n'
        UNIX_LINE_ENDING = b'\n'
        with open(file_path, 'rb') as f:
            content = f.read()
        content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)
        with open(file_path, 'wb') as f:
            f.write(content)


    @staticmethod
    def ask_for_directory():
        root = Tk()
        root.withdraw()  # Hide the main window
        directory = filedialog.askdirectory(title="Select Project Directory")
        return directory if directory else None
    
    @staticmethod
    def ask_for_file(filetypes=[("STL Geometry", "*.stl")]):
        root = Tk()
        root.withdraw()
        file = filedialog.askopenfilename(title="Select File", filetypes=filetypes)
        return file if file else None

    @staticmethod
    def dict_to_yaml(data, output_file):
        """
        Convert a dictionary to a YAML file.

        Parameters:
        - data (dict): The dictionary to be converted.
        - output_file (str): The name of the output YAML file.
        """
        with open(output_file, 'w') as file:
            yaml.dump(data, file, default_flow_style=False, sort_keys=False)
        #print(f"YAML file '{output_file}' has been created.")


    @staticmethod
    def yaml_to_dict(input_file):
        """
        Read a YAML file and convert it to a dictionary.

        Parameters:
        - input_file (str): The name of the input YAML file.

        Returns:
        - dict: The dictionary representation of the YAML file.
        """
        with open(input_file, 'r') as file:
            data = yaml.safe_load(file)
        return data

    @staticmethod
    # This file contains the basic primitives used in the generation of OpenFOAM casefiles
    def createFoamHeader(className="dictionary",objectName="blockMeshDict"):
        header = f"""/*----------------------------*- AmpersandCFD -*------------------------------*\ 
/*--------------------------------*- C++ -*----------------------------------*\         
| ========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

/*This file is part of OpenFOAM casefiles automatically generated by AmpersandCFD*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {className};
    object      {objectName};
}}"""
        return header

    @staticmethod
    def createDimensions(M=1,L=1,T=1):
        return f"\ndimensions      [{M} {L} {T} 0 0 0 0];"
    
    @staticmethod
    def createInternalFieldScalar(type="uniform",value=0):
        return f"""\ninternalField   {type} {value};"""
    
    @staticmethod
    def createInternalFieldVector(type="uniform",value=[0,0,0]):
        return f"""\ninternalField   {type} ({value[0]} {value[1]} {value[2]});"""

    @staticmethod
    def write_to_file(filename, content):
        with open(filename, 'w') as f:
            f.write(content)
    @staticmethod
    def createScalarFixedValue(patch_name="inlet",value=0):
        return f"""\n{patch_name}
        {{
            type            fixedValue;
            value           uniform {value};
        }};"""

    @staticmethod
    def createScalarZeroGradient(patch_name="inlet"):
        return f"""\n{patch_name}
        {{
            type            zeroGradient;
        }};"""

    @staticmethod
    def createVectorFixedValue(patch_name="inlet",value=[0,0,0]):
        return f"""\n{patch_name}
        {{
            type            fixedValue;
            value           uniform ("{value[0]} {value[1]} {value[2]})";
        }};""" 

    @staticmethod
    def createVectorZeroGradient(patch_name="inlet"):
        return f"""\n{patch_name}
        {{
            type            zeroGradient;
        }};"""
    
    @staticmethod
    def write_dict_to_file(filename, content):
        with open(filename, 'w') as f:
            f.write(content)

    @staticmethod
    # to remove duplicates from a list
    def remove_duplicates(lst):
        return list(set(lst))

class ampersandIO:
    def __init__(self):
        pass

    @staticmethod
    def printMessage(*args):
        print(*args)
    
    @staticmethod
    def printError(*args):
        print(*args, file=sys.stderr)
    
    @staticmethod
    def get_input(prompt):
        return input(prompt)
    
    @staticmethod
    def get_input_int(prompt):
        return int(input(prompt))
    
    @staticmethod  
    def get_input_float(prompt):
        return float(input(prompt))
    
    @staticmethod
    def print_numbered_list(lst):
        for i in range(len(lst)):
            print(f"{i+1}. {lst[i]}")
    
    @staticmethod
    def get_input_vector(prompt):
        return list(map(float, input(prompt).split()))
    
    @staticmethod
    def get_input_bool(prompt):
        return input(prompt).lower() in ['y', 'yes', 'true', '1']
    

class ampersandDataInput:
    def __init__(self):
        pass

    @staticmethod
    def get_inlet_values():
        U = ampersandIO.get_input_vector("Enter the velocity vector at the inlet (m/s): ")
        return U
    
    @staticmethod
    def get_physical_properties():
        rho = ampersandIO.get_input_float("Enter the density of the fluid (kg/m^3): ")
        nu = ampersandIO.get_input_float("Enter the kinematic viscosity of the fluid (m^2/s): ")
        return rho, nu
    
    @staticmethod
    def choose_fluid_properties():
        fluids = {"Air":{'rho':1.225, 'nu':1.5e-5}, "Water":{'rho':1000, 'nu':1e-6}, }
        fluid_names = list(fluids.keys())
        fluid_name = ampersandIO.get_input_int("Choose the fluid properties:\n" + "\n".join([f"{i+1}. {fluid_names[i]}" for i in range(len(fluid_names))]) + "\n")
        if(fluid_name>len(fluids) or fluid_name<=0):
            ampersandIO.printMessage("Invalid fluid choice. Please input fluid properties manually.")
            rho, nu = ampersandDataInput.get_physical_properties()
            return {'rho':rho, 'nu':nu}
        fluid = fluids[fluid_names[fluid_name-1]]
        return fluid
    
    @staticmethod
    def get_mesh_refinement_level():
        refLevel = ampersandIO.get_input_int("Enter the mesh refinement (0: coarse, 1: medium, 2: fine): ")
        if refLevel not in [0,1,2]:
            ampersandIO.printMessage("Invalid mesh refinement level. Defaulting to medium.")
            refLevel = 1
        return refLevel
    


if __name__ == "__main__":
    print(ampersandPrimitives.createFoamHeader(className="dictionary",objectName="snappyHexMeshDict"))
    print(ampersandPrimitives.createDimensions(M=1,L=1,T=1))
    print(ampersandPrimitives.createScalarFixedValue(patch_name="inlet",value=0))
    print(ampersandPrimitives.createScalarZeroGradient(patch_name="inlet"))
    print(ampersandPrimitives.createVectorFixedValue(patch_name="inlet",value=[0,0,0]))



