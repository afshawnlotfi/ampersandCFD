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

from project import AmpersandProject
from primitives import AmpersandUtils, AmpersandIO
from ampersand.utils.headers import get_ampersand_header
import os

def create_project():
    project = AmpersandProject()
    # Clear the screen
    os.system('cls' if os.name == 'nt' else 'clear')
    AmpersandIO.printMessage(get_ampersand_header())
    project.set_project_directory(AmpersandUtils.ask_for_directory())
    if project.project_directory_path == None:
        AmpersandIO.printMessage("No project directory selected. Exiting...")
        exit()
    project_name = AmpersandIO.get_input("Enter the project name: ")
    project.set_project_name(project_name)
    #user_name = input("Enter the user name: ")
    #project.set_user_name(user_name)
    project.create_project_path()
    AmpersandIO.printMessage("Creating the project")
    AmpersandIO.printMessage(f"Project path: {project.project_path}")
    #project.project_path = r"C:\Users\Ridwa\Desktop\CFD\ampersandTests\drivAer2"
    project.create_project()
    project.create_settings()
    AmpersandIO.printMessage("Preparing for mesh generation")
    project.ask_refinement_level()
    yN = AmpersandIO.get_input("Add STL file to the project (y/N)?: ")
    while yN.lower() == 'y':
        project.add_stl_file()
        yN = AmpersandIO.get_input("Add another STL file to the project (y/N)?: ")
    project.add_stl_to_project()
    
    # Before creating the project files, the settings are flushed to the project_settings.yaml file
    
    project.ask_flow_type()
    if(project.internalFlow!=True):
        project.ask_ground_type()
    AmpersandIO.printMessage("Fluid properties and inlet values are necessary for mesh size calculations")
    project.set_fluid_properties()
    project.set_inlet_values()
    project.ask_transient_settings()
    project.set_parallel()
    project.ask_half_model()
    
    if(len(project.stl_files)>0):
        project.analyze_stl_file()
    
    project.useFOs = AmpersandIO.get_input_bool("Use function objects for post-processing (y/N)?: ")
    project.set_post_process_settings()
    #project.list_stl_files()
    project.summarize_project()
    #project.analyze_stl_file()
    
    project.write_settings()
    project.create_project_files()

if __name__ == '__main__':
    # Specify the output YAML file
    try:
        create_project()
    except KeyboardInterrupt:
        AmpersandIO.printMessage("\nKeyboardInterrupt detected! Aborting project creation")
        exit()
    except Exception as error:
        AmpersandIO.printError(error)
        exit()