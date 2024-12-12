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

from typing import Literal
from primitives import ampersandPrimitives
from constants import MeshSettings
from src.utils.headers import createFoamHeader


def create_surfaceDict(meshSettings: MeshSettings, objectName: Literal["surfaceFeatureExtractDict", "surfaceFeaturesDict"]) -> str:
    header = createFoamHeader(className="dictionary", objectName=objectName)
    surfaceDict = header
    for anEntry in meshSettings.geometry:
        if anEntry.type == 'triSurfaceMesh':
            surfaceFeature = f"""\n{anEntry.name}
{{
extractionMethod    extractFromSurface;
includedAngle   170;
subsetFeatures
{{
    nonManifoldEdges       no;
    openEdges       yes;
}}
writeObj            yes;
writeSets           no;
}}"""
            surfaceDict += surfaceFeature
    return surfaceDict


if __name__ == "__main__":
    meshSettings = MeshSettings.model_validate(
        ampersandPrimitives.yaml_to_dict('meshSettings.yaml'))
    surfaceFeatureExtractDict = create_surfaceDict(
        meshSettings, "surfaceFeatureExtractDict")
    with open('surfaceFeatureExtractDict', 'w') as file:
        file.write(surfaceFeatureExtractDict)
