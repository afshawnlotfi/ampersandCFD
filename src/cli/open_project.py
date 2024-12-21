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

from pathlib import Path
from src.cli.mod_project import ModProject
from src.primitives import AmpersandDataInput, AmpersandUtils, AmpersandIO

from src.services.project_service import ProjectService

def open_project():
    AmpersandIO.printMessage("Please select the project directory to open")
    
    parent_directory = AmpersandUtils.ask_for_directory()
    project_name = AmpersandIO.get_input("Enter the project name: ")
    project_path = Path(f"{parent_directory}/{project_name}")


    project = ProjectService.load_project(project_path)

    AmpersandIO.printMessage("Project loaded successfully")
    project.summarize_project()
    modify_project = AmpersandIO.get_input_bool(
        "Do you want to modify the project settings (y/N)?: ")
    project_modified = False  # flag to check if the project has been modified
   
    while modify_project:
        modification_type = AmpersandDataInput.choose_modification_categorized()
        ModProject.modify_project(project, modification_type)
        project.write_settings()
        project_modified = True
        modify_project = AmpersandIO.get_input_bool("Do you want to modify another settings (y/N)?: ")
    # project.choose_modification()
    if project_modified:  # if the project is modified at least once
        AmpersandIO.printMessage(
            "Generating the project files based on the new settings")
        # if everything is successful, write the settings to the project_settings.yaml file
        project.write_settings()
        # Then create the project files with the new settings
        project.write_project_files()
    else:
        AmpersandIO.printMessage(
            "No modifications were made to the project settings")
    return 0


if __name__ == '__main__':
    # Specify the output YAML file
    try:
        open_project()
    except KeyboardInterrupt:
        AmpersandIO.printMessage(
            "\nKeyboardInterrupt detected! Aborting project creation")
        exit()
    # except Exception as error:
    #    ampersandIO.printError(error)
    #    exit()
