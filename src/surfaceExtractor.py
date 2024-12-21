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
from src.models.settings import MeshSettings, TriSurfaceMeshGeometry
from src.utils.generation import GenerationUtils


def create_surfaceDict(meshSettings: MeshSettings, objectName: Literal["surfaceFeatureExtractDict", "surfaceFeaturesDict"]) -> str:
    surfaceDict = GenerationUtils.createFoamHeader("dictionary", objectName)
    for geometry_name, geometry in meshSettings.geometry.items():
        if isinstance(geometry, TriSurfaceMeshGeometry):
            surfaceDict += f"""\n{geometry_name}
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
    return surfaceDict

if __name__ == "__main__":
    from src.primitives import AmpersandUtils

    meshSettings = MeshSettings.model_validate(
        AmpersandUtils.yaml_to_dict('meshSettings.yaml'))
    surfaceFeatureExtractDict = create_surfaceDict(
        meshSettings, "surfaceFeatureExtractDict")
    with open('surfaceFeatureExtractDict', 'w') as file:
        file.write(surfaceFeatureExtractDict)
