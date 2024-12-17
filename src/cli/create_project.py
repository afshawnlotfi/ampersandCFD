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

from src.project import AmpersandProject
from src.primitives import AmpersandDataInput, ampersandPrimitives, ampersandIO
import os


def create_project():
    project = AmpersandProject()
    # Clear the screen
    os.system('cls' if os.name == 'nt' else 'clear')
    project_directory = ampersandPrimitives.ask_for_directory()
    project.set_project_directory(project_directory)
    
    if project.project_directory_path == None:
        ampersandIO.printMessage("No project directory selected. Exiting...")
        exit()

    project_name = ampersandIO.get_input("Enter the project name: ")
    project.set_project_name(project_name)

    project.create_project_path()
    ampersandIO.printMessage("Creating the project")
    ampersandIO.printMessage(f"Project path: {project.project_path}")
    # project.project_path = r"C:\Users\Ridwa\Desktop\CFD\ampersandTests\drivAer2"
    project.create_project()
    project.create_settings()
    ampersandIO.printMessage("Preparing for mesh generation")


    refinement_level = AmpersandDataInput.get_mesh_refinement_level()
    project.set_refinement_level(refinement_level)

    yN = ampersandIO.get_input("Add STL file to the project (y/N)?: ")
    while yN.lower() == 'y':
        project.add_stl_file()
        yN = ampersandIO.get_input("Add another STL file to the project (y/N)?: ")
    project.add_stl_to_project()

    # Before creating the project files, the settings are flushed to the project_settings.yaml file

    flow_type = ampersandIO.get_input(
        "Internal or External Flow (I/E)?: ").lower() == 'i'
    project.set_flow_type(flow_type)

    if (not project.internalFlow):
        ground_type = ampersandIO.get_input_bool(
            "Is the ground touching the body (y/N): ")
        project.set_ground_type(ground_type)

    ampersandIO.printMessage(
        "Fluid properties and inlet values are necessary for mesh size calculations")

    fluid = AmpersandDataInput.choose_fluid_properties()
    project.set_fluid_properties(fluid)

    U = AmpersandDataInput.get_inlet_values()
    project.set_inlet_values(U)

    transient = ampersandIO.get_input(
        "Transient or Steady State (T/S)?: ").lower() == 't'
    project.set_transient_settings(transient)

    n_core = ampersandIO.get_input_int(
        "Number of cores for parallel simulation: ")
    project.set_parallel(n_core)

    half_model = ampersandIO.get_input_bool("Half Model (y/N)?: ")
    project.set_half_model(half_model)

    if (len(project.stl_files) > 0):
        project.analyze_stl_file()

    useFOs = ampersandIO.get_input_bool(
        "Use function objects for post-processing (y/N)?: ")
    project.set_post_process_settings(useFOs)

    project.summarize_project()

    project.write_settings()
    project.create_project_files()


if __name__ == '__main__':
    # Specify the output YAML file
    try:
        create_project()
    except KeyboardInterrupt:
        ampersandIO.printMessage(
            "\nKeyboardInterrupt detected! Aborting project creation")
        exit()
    except Exception as error:
        ampersandIO.printError(error)
        exit()
