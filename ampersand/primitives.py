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

import os
from typing import Iterable
import numpy as np
import yaml
import sys
from tkinter import filedialog, Tk
from ampersand.models.settings import SnappyHexMeshSettings
from ampersand.utils.generation import GenerationUtils
from pydantic import BaseModel

try:
    from PySide6.QtWidgets import QMessageBox # type: ignore
    from ampersand.gui.dialogBoxes import inputDialogDriver, vectorInputDialogDriver
except:
    pass


class FluidPhysicalProperties(BaseModel):
    rho: float
    nu: float

FLUID_PYSICAL_PROPERTIES = {
    "Air": FluidPhysicalProperties(rho= 1.225, nu= 1.5e-5),
    "Water": FluidPhysicalProperties(rho= 1000, nu= 1e-6),
}

class AmpersandUtils:
    @staticmethod
    def list_stl_files(stl_files, GUIMode=False, window=None):
        if GUIMode:
            stl_names = [stl_file.name for stl_file in stl_files]
            # window.listWidgetObjList.clear()
            for i in range(len(stl_names)):
                AmpersandIO.printMessage(
                    f"{i+1}. {stl_names[i]}", GUIMode=GUIMode, window=window)
            return 0
        i = 1
        AmpersandIO.show_title("STL Files")

        AmpersandIO.printMessage(f"{'No.':<5}{'Name':<20}{'Purpose':<20}{
                                 'RefineMent':<15}{'Property':<15}")
        for stl_file in stl_files:
            if (stl_file.property == None):
                stl_property = "None"
                if stl_file.purpose == 'wall':
                    stl_property = f"nLayers: {stl_file.nLayers}"
                else:
                    stl_property = "None"
                # stl_file.property = "None"
            elif isinstance(stl_file.property, list):
                stl_property = f"[{stl_file.property[0]} {
                    stl_file.property[1]} {stl_file.property[2]}]"
            elif isinstance(stl_file.property, tuple):
                if stl_file.purpose == 'inlet':
                    stl_property = f"U: [{stl_file.property[0]} {
                        stl_file.property[1]} {stl_file.property[2]}]"
                elif stl_file.purpose == 'cellZone':
                    stl_property = f"Refinement: {stl_file.property[0]}"
                # stl_property = f"[{stl_file.property[0]} {stl_file.property[1]} {stl_file.property[2]}]"
            else:
                stl_property = stl_file.property
            AmpersandIO.printMessage(f"{i:<5}{stl_file.name:<20}{stl_file.purpose:<20}({
                                     stl_file.refineMin} {stl_file.refineMax}{')':<11}{stl_property:<15}")
            i += 1
        AmpersandIO.show_line()
        return 0

    @staticmethod
    def list_boundary_conditions(meshSettings: SnappyHexMeshSettings):
        i = 1
        boundaries = []
        AmpersandIO.show_title("Boundary Conditions")
        AmpersandIO.printMessage(f"{'No.':<5}{'Name':<20}{
                                 'Purpose':<20}{'Value':<15}")
        # for external flows, show the boundary conditions for domain first
        if meshSettings.internalFlow == False:
            for patch in meshSettings.patches:
                if patch.property == None:
                    property = "None"
                elif isinstance(patch.property, list):
                    property = f"[{patch.property[0]} {
                        patch.property[1]} {patch.property[2]}]"
                elif isinstance(patch.property, tuple):
                    property = f"[{patch.property[0]} {
                        patch.property[1]} {patch.property[2]}]"
                else:
                    property = patch.property
                # ampersandIO.printMessage(f"{patch.name}: {patch.purpose}\t{patch.property}")
                AmpersandIO.printMessage(f"{i:<5}{patch.name:<20}{
                                         patch.purpose:<20}{property:<15}")
                i += 1
                boundaries.append(patch.name)
        for patch in meshSettings.geometry:
            if patch.purpose != 'refinementRegion' and patch.purpose != 'refinementSurface':
                # ampersandIO.printMessage(patch)
                if patch.property == None:
                    property = "None"
                elif isinstance(patch.property, list):
                    property = f"[{patch.property[0]} {
                        patch.property[1]} {patch.property[2]}]"
                elif isinstance(patch.property, tuple):
                    property = f"[{patch.property[0]} {
                        patch.property[1]} {patch.property[2]}]"
                else:
                    property = "None"
                AmpersandIO.printMessage(f"{i:<5}{patch.name:<20}{
                                         patch.purpose:<20}{property:<15}")
                i += 1
                boundaries.append(patch.name)
        return boundaries  # return the number of boundarys
        # ampersandIO.printMessage(f"{patch.name}: {patch.purpose}\t{patch.property}")


    @staticmethod
    # Function to recursively convert tuples to lists (or any other conversion)
    def sanitize_yaml(data):
        if isinstance(data, tuple):
            return list(data)
        elif isinstance(data, dict):
            return {k: AmpersandUtils.sanitize_yaml(v) for k, v in data.items()}
        elif isinstance(data, list):
            return [AmpersandUtils.sanitize_yaml(item) for item in data]
        else:
            return data

    # Function to remove duplicates in a YAML file
    @staticmethod
    def remove_duplicate_dicts(dict_list):
        seen = set()
        unique_dicts = []
        for d in dict_list:
            # Convert dictionary to a frozenset of its items to make it hashable
            # print(d)
            dict_tuple = frozenset(d.items())
            if dict_tuple not in seen:
                seen.add(dict_tuple)
                unique_dicts.append(d)
        return unique_dicts

    @staticmethod
    def treat_bounds(geometry):
        for anObject in geometry:
            AmpersandUtils.list_to_tuple_dict(anObject)
            # print(anObject)
        return geometry

    # Function to convert a list to a tuple inside a dictionary
    @staticmethod
    def list_to_tuple_dict(data):
        for key, value in data.items():
            if isinstance(value, dict):
                AmpersandUtils.list_to_tuple_dict(value)
            elif isinstance(value, list):
                data[key] = tuple(value)
        return data

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
    def ask_for_directory(qt=False):
        try:
            if qt:
                from PySide6.QtWidgets import QFileDialog # type: ignore
                directory = QFileDialog.getExistingDirectory(
                    None, "Select Project Directory")
                return directory if directory else None
            else:
                root = Tk()
                root.withdraw()  # Hide the main window
                directory = filedialog.askdirectory(
                    title="Select Project Directory")
                return directory if directory else None
        except:
            return AmpersandIO.get_input("Select Project Directory: ")

    @staticmethod
    def ask_for_file(filetypes=[("STL Geometry", "*.stl")], qt=False):
        try:
            if qt:
                from PySide6.QtWidgets import QFileDialog # type: ignore
                file = QFileDialog.getOpenFileName(
                    None, "Select File", filter="STL Geometry (*.stl)")
                return file[0] if file[0] else None
            else:
                root = Tk()
                root.withdraw()
                file = filedialog.askopenfilename(
                    title="Select File", filetypes=filetypes)
                return file if file else None
        except:
            return AmpersandIO.get_input("Select file: ")

    @staticmethod
    def check_dict(dict_):
        # check every elements of the dictionary and whether there are tuples
        """
        dict_: The dictionary to be checked.
        This function checks every element of the dictionary
        and converts tuples to lists."""
        for key, value in dict_.items():
            if isinstance(value, dict):
                AmpersandUtils.check_dict(value)
            elif isinstance(value, tuple):
                dict_[key] = list(value)
        return dict_

    @staticmethod
    def dict_to_yaml(data, output_file):
        """
        Convert a dictionary to a YAML file.

        Parameters:
        - data (dict): The dictionary to be converted.
        - output_file (str): The name of the output YAML file.
        """
        # data = ampersandPrimitives.check_dict(data)
        data = AmpersandUtils.sanitize_yaml(data)
        with open(output_file, 'w') as file:
            yaml.dump(data, file, default_flow_style=False, sort_keys=False)
        # print(f"YAML file '{output_file}' has been created.")

    @staticmethod
    def yaml_to_dict(input_file):
        """
        Read a YAML file and convert it to a dictionary.

        Parameters:
        - input_file (str): The name of the input YAML file.

        Returns:
        - dict: The dictionary representation of the YAML file.
        """
        try:
            with open(input_file, 'r') as file:
                data = yaml.safe_load(file)
            return data
        except Exception as e:
            print(f"Error reading YAML file: {e}")
            yN = AmpersandIO.get_input_bool("Continue y/N?")
            if yN:
                return None
            else:
                exit()
            return None

    @staticmethod
    def write_to_file(filename, content):
        with open(filename, 'w') as f:
            f.write(content)

    @staticmethod
    def write_dict_to_file(filename, content):
        try:
            with open(filename, 'w') as f:
                f.write(content)
        except Exception as e:
            print(f"Error writing to file: {e}")

    @staticmethod
    def calc_Umag(U: Iterable[float]):
        return sum([u**2 for u in U])**0.5


class AmpersandIO:
    def __init__(self):
        pass

    @staticmethod
    def printMessage(*args, GUIMode=False, window=None):
        if GUIMode and window != None:
            window.updateTerminal(*args)
        else:
            print(*args)

    @staticmethod
    def printWarning(*args, GUIMode=False):
        if GUIMode:
            # ampersandIO.printMessage(*args)
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Warning")
            msg.setInformativeText(*args)
            msg.setWindowTitle("Warning")
            msg.exec_()
        else:
            print(*args)

    @staticmethod
    def printError(*args, GUIMode=False):
        if GUIMode:
            # ampersandIO.printMessage(*args)
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error")
            msg.setInformativeText(*args)
            msg.setWindowTitle("Error")
            msg.exec_()
        else:
            print(*args, file=sys.stderr)

    @staticmethod
    def get_input(prompt, GUIMode=False):
        if GUIMode:
            return inputDialogDriver(prompt)
        else:
            return input(prompt)

    @staticmethod
    def print_dict(data):
        for key, value in data.items():
            print(f"{key}: {value}")

    @staticmethod
    def get_input_int(prompt, GUIMode=False):
        if GUIMode:
            return inputDialogDriver(prompt, input_type=int)
        else:
            try:
                return int(input(prompt))
            except:
                AmpersandIO.printError(
                    "Invalid input. Please enter an integer.")
                return AmpersandIO.get_input_int(prompt)

    @staticmethod
    def get_input_float(prompt, GUIMode=False):
        if GUIMode:
            return inputDialogDriver(prompt, input_type=float)
        else:
            try:
                return float(input(prompt))
            except:
                AmpersandIO.printError("Invalid input. Please enter a number.")
                return AmpersandIO.get_input_float(prompt)

    @staticmethod
    def show_list(lst):
        i = 1
        for item in lst:
            AmpersandIO.printMessage(f"{i}. {item}")

    @staticmethod
    def print_numbered_list(lst):
        for i in range(len(lst)):
            print(f"{i+1}. {lst[i]}")

    @staticmethod
    def get_input_vector(prompt, GUIMode=False):
        if GUIMode:
            return vectorInputDialogDriver(prompt)
        else:
            inp = input(prompt).split()
            # output = [0.,0.,0.]
            # Check if the input is a list of floats
            try:
                vec = tuple(map(float, inp))
                if len(vec) != 3:
                    AmpersandIO.printError(
                        "Invalid input. Please enter 3 numbers.")
                    # Recursively call the function until a valid input is given
                    return AmpersandIO.get_input_vector(prompt)
                return vec
            except:
                AmpersandIO.printError(
                    "Invalid input. Please enter a list of numbers.")
                # Recursively call the function until a valid input is given
                return AmpersandIO.get_input_vector(prompt)
        
        # return list(map(float, input(prompt).split()))

    @staticmethod
    def get_input_bool(prompt):
        try:
            return input(prompt).lower() in ['y', 'yes', 'true', '1']
        except:
            AmpersandIO.printError(
                "Invalid input. Please enter a boolean value.")
            return AmpersandIO.get_input_bool(prompt)



    @staticmethod
    def show_title(title):
        total_len = 60
        half_len = (total_len - len(title))//2
        title = "-"*half_len + title + "-"*half_len
        AmpersandIO.printMessage("\n" + title)

    @staticmethod
    def show_line():
        AmpersandIO.printMessage("-"*60)

    @staticmethod
    def printFormat(item_name, item_value, GUIMode=False, window=None):
        if GUIMode:
            window.updateStatusBar(f"{item_name}: {item_value}")
        else:
            print(f"{item_name:12}\t{item_value}")


class AmpersandDataInput:

    @staticmethod
    def get_option_choice(prompt: str, options: list, title=None):
        if title:
            AmpersandIO.printMessage(title)
        AmpersandIO.print_numbered_list(options)
        choice = AmpersandIO.get_input_int(prompt)
        if choice > len(options) or choice <= 0:
            AmpersandIO.printError(
                "Invalid choice. Please choose from the given options.")
            return AmpersandDataInput.get_option_choice(prompt, options)
        
        return choice-1

    @staticmethod
    def get_inlet_values():
        U = AmpersandIO.get_input_vector(
            "Enter the velocity vector at the inlet (m/s): ")
        return U

    @staticmethod
    def get_physical_properties():
        rho = AmpersandIO.get_input_float(
            "Enter the density of the fluid (kg/m^3): ")
        nu = AmpersandIO.get_input_float(
            "Enter the kinematic viscosity of the fluid (m^2/s): ")
        return FluidPhysicalProperties(rho=rho, nu=nu)

    @staticmethod
    def get_turbulence_model():
        turbulence_models = ['kOmegaSST', 'kEpsilon']
        AmpersandIO.show_title("Turbulence models")
        for i in range(len(turbulence_models)):
            AmpersandIO.printMessage(f"{i+1}. {turbulence_models[i]}")
        turbulence_model = AmpersandIO.get_input_int(
            "Choose the turbulence model: ")
        if turbulence_model > len(turbulence_models) or turbulence_model <= 0:
            AmpersandIO.printError(
                "Invalid turbulence model. Defaulting to kOmegaSST.")
            turbulence_model = 1
        return turbulence_models[turbulence_model-1]

    @staticmethod
    def choose_fluid_properties():

        fluid_names = list(FLUID_PYSICAL_PROPERTIES.keys())
        AmpersandIO.printMessage("Fluid properties")
        AmpersandIO.printMessage("0. Enter fluid properties manually")
        for i in range(len(fluid_names)):
            AmpersandIO.printMessage(f"{i+1}. {fluid_names[i]}")
        fluid_name = AmpersandIO.get_input_int("Choose the fluid properties:")

        if (fluid_name > len(FLUID_PYSICAL_PROPERTIES) or fluid_name <= 0):
            AmpersandIO.printMessage("Please input fluid properties manually.")
            return AmpersandDataInput.get_physical_properties()
        return FLUID_PYSICAL_PROPERTIES[fluid_names[fluid_name-1]]


if __name__ == "__main__":
    print(GenerationUtils.createFoamHeader(
        className="dictionary", objectName="snappyHexMeshDict")
    )
    print(GenerationUtils.createDimensions(M=1, L=1, T=1))
    print(GenerationUtils.createScalarFixedValue(patch_name="inlet", value=0))
    print(GenerationUtils.createScalarZeroGradient(patch_name="inlet"))
    print(GenerationUtils.createVectorFixedValue(
        patch_name="inlet", value=[0, 0, 0]))
