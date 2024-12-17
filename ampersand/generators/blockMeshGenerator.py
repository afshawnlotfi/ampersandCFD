# %%
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

from ampersand.models.settings import MeshSettings
from ampersand.utils.generation import GenerationUtils


def create_blockMeshDict(mesh_settings: MeshSettings) -> str:
    return f"""{GenerationUtils.createFoamHeader("dictionary", "blockMeshDict")}

// ********* Domain *********
scale {mesh_settings.scale};

vertices
(
    ({mesh_settings.domain.minx} {mesh_settings.domain.miny} {mesh_settings.domain.minz})
    ({mesh_settings.domain.maxx} {mesh_settings.domain.miny} {mesh_settings.domain.minz})
    ({mesh_settings.domain.maxx} {mesh_settings.domain.maxy} {mesh_settings.domain.minz})
    ({mesh_settings.domain.minx} {mesh_settings.domain.maxy} {mesh_settings.domain.minz})
    ({mesh_settings.domain.minx} {mesh_settings.domain.miny} {mesh_settings.domain.maxz})
    ({mesh_settings.domain.maxx} {mesh_settings.domain.miny} {mesh_settings.domain.maxz})
    ({mesh_settings.domain.maxx} {mesh_settings.domain.maxy} {mesh_settings.domain.maxz})
    ({mesh_settings.domain.minx} {mesh_settings.domain.maxy} {mesh_settings.domain.maxz})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({mesh_settings.domain.nx} {mesh_settings.domain.ny} {mesh_settings.domain.nz}) simpleGrading (1 1 1)
);

edges
(
);

boundary
({"\n".join(
        f"""
    {patch_name}
    {{
        type {patch.type};
        faces
        (
            ({patch.faces[0]} {patch.faces[1]} {patch.faces[2]} {patch.faces[3]})
        );
    }}
""" for patch_name, patch in mesh_settings.patches.items())}

);
mergePatchPairs
(
);

// ************************************************************************* //
"""


# Generate blockMeshDict
# read in data to meshSettings from meshSettings.yaml
if __name__ == "__main__":
    from ampersand.primitives import AmpersandUtils
    from ampersand.models.settings import SnappyHexMeshSettings

    meshSettings = SnappyHexMeshSettings.model_validate(AmpersandUtils.yaml_to_dict("examples/basic/meshSettings.yaml"))
    blockMeshDict = create_blockMeshDict(meshSettings)

    # Save to file
    with open("outputs/blockMeshDict", "w") as f:
        f.write(blockMeshDict)

    print("blockMeshDict file created.")

# %%
