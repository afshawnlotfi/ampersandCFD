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

from constants import MeshSettings, meshSettings
from src.primitives import ampersandPrimitives
from src.utils.dict_generation import DictGenerationUtils


def create_blockMeshDict(meshSettings: MeshSettings) -> str:
    header = DictGenerationUtils.createFoamHeader(
        className="dictionary", objectName="blockMeshDict"
    )
    blockMeshDict = header + f"""

// ********* Domain *********
scale {meshSettings.scale};

vertices
(
    ({meshSettings.domain.minx} {meshSettings.domain.miny} {meshSettings.domain.minz})
    ({meshSettings.domain.maxx} {meshSettings.domain.miny} {meshSettings.domain.minz})
    ({meshSettings.domain.maxx} {meshSettings.domain.maxy} {meshSettings.domain.minz})
    ({meshSettings.domain.minx} {meshSettings.domain.maxy} {meshSettings.domain.minz})
    ({meshSettings.domain.minx} {meshSettings.domain.miny} {meshSettings.domain.maxz})
    ({meshSettings.domain.maxx} {meshSettings.domain.miny} {meshSettings.domain.maxz})
    ({meshSettings.domain.maxx} {meshSettings.domain.maxy} {meshSettings.domain.maxz})
    ({meshSettings.domain.minx} {meshSettings.domain.maxy} {meshSettings.domain.maxz})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({meshSettings.domain.nx} {meshSettings.domain.ny} {meshSettings.domain.nz}) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
"""
    for name, patch in meshSettings.patches.items():
        blockMeshDict += f"""\n    {name}
{{
        type {patch.type};
        faces
        (
                ({patch.faces[0]} {patch.faces[1]} {patch.faces[2]} {patch.faces[3]})
        );
}}\n"""

    blockMeshDict += """);
mergePatchPairs
(
);

// ************************************************************************* //
"""

    return blockMeshDict


# Generate blockMeshDict
# read in data to meshSettings from meshSettings.yaml
if __name__ == "__main__":
    meshSettings = MeshSettings.model_validate(
        ampersandPrimitives.yaml_to_dict("examples/basic/meshSettings.yaml"))
    blockMeshDict = create_blockMeshDict(meshSettings)

    # Save to file
    with open("outputs/blockMeshDict", "w") as f:
        f.write(blockMeshDict)

    print("blockMeshDict file created.")
