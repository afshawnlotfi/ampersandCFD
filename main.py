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
from src.primitives import AmpersandIO
from src.cli.create_project import create_project
from src.cli.open_project import open_project
from src.headers import AMPERSAND_HEADER
from src.utils.watch_sim import watch_sim
import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Ampersand CFD Automation Tool')
    parser.add_argument('--create', action='store_true',
                        help='Create a new project')
    parser.add_argument('--open', action='store_true',
                        help='Open an existing project')
    parser.add_argument('--post', action='store_true',
                        help='Post-process the simulation')
    args = parser.parse_args()

    # Clear the screen
    os.system('cls' if os.name == 'nt' else 'clear')

    AmpersandIO.printMessage(AMPERSAND_HEADER)

    if args.create:
        try:
            create_project()
        except KeyboardInterrupt:
            AmpersandIO.printMessage(
                "\nKeyboardInterrupt detected! Aborting project creation")
            exit()
        except Exception as error:
            AmpersandIO.printError(error)
    elif args.open:
        try:
            open_project()
        except KeyboardInterrupt:
            AmpersandIO.printMessage(
                "\nKeyboardInterrupt detected! Aborting project creation")
            exit()
        except Exception as error:
            AmpersandIO.printError(error)
    elif args.post:
        try:
            watch_sim()
        except KeyboardInterrupt:
            AmpersandIO.printMessage(
                "\nKeyboardInterrupt detected! Aborting project creation")
            exit()
        except Exception as error:
            AmpersandIO.printError(error)
    else:
        AmpersandIO.printMessage(
            "Please specify an action to perform. Use --help for more information.")
        parser.print_help()


if __name__ == '__main__':
    # Specify the output YAML file
    main()
    # open_project()
    # create_project()
