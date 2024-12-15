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

from ampersand.models.settings import SearchableBoxGeometry, TriSurfaceMeshGeometry, MeshSettings
from primitives import AmpersandUtils
from ampersand.utils.generation import GenerationUtils


def create_snappyHexMeshDict(meshSettings: MeshSettings):
    """
    Create a snappyHexMeshDict for OpenFOAM.

    Parameters:
    stl_files (list): A list of STL file names.
    refinement_levels (dict): A dictionary where keys are STL file names and values are refinement levels.
    layer_counts (dict): A dictionary where keys are STL file names and values are layer counts.

    Returns:
    str: The content of the snappyHexMeshDict file as a string.
    """
    snappyHexMeshDict = f""
    trueFalse = {True: "true", False: "false"}
    header = GenerationUtils.createFoamHeader(
        className="dictionary", objectName="snappyHexMeshDict")

    steps = f"""
castellatedMesh {meshSettings.snappyHexSteps.castellatedMesh};
snap            {meshSettings.snappyHexSteps.snap};
addLayers       {meshSettings.snappyHexSteps.addLayers};"""

    features = ""
    refinementSurfaces = ""
    added_geo = ""
    # maxRefinementLevel = 1
    # minRefinementRegionLevel = 2
    geometry_str = f"""\ngeometry\n{{"""
    for geometry_name, geometry in meshSettings.geometry.items():
        # For STL surfaces, featureEdges and refinementSurfaces are added
        # maxRefinementLevel = max(maxRefinementLevel, geometry_patch.refineMax)
        if (isinstance(geometry, TriSurfaceMeshGeometry)):
            added_geo = f"""\n
    {geometry_name}
    {{
        type {geometry.type};
        name {geometry_name[:-4]};
        regions
        {{
            {geometry_name[:-4]}
            {{
                name {geometry_name[:-4]};
            }}
        }}
    }}"""
            # Add features and refinement surfaces
            if (geometry.featureEdges):
                features += f"""

        {{
            file \"{geometry_name[:-4]}.eMesh\";
            level {geometry.featureLevel};
        }}"""
            if (geometry.purpose == 'inlet' or geometry.purpose == 'outlet'):
                patchType = 'patch'
                refinementSurfaces += f"""
        {geometry_name[:-4]}
        {{
            level (0 0);
            regions
            {{
                {geometry_name[:-4]}
                {{
                    level ({geometry.refineMin} {geometry.refineMax});
                    patchInfo
                    {{
                        type {patchType};
                    }}
                }}
            }}
        }}"""
            elif (geometry.purpose == 'cellZone'):
                patchType = 'cellZone'
                if geometry.property[1] == True:  # patches will be added
                    refinementSurfaces += f"""
        {geometry_name[:-4]}
        {{
            level (0 0);
            cellZone {geometry_name[:-4]};
            faceZone {geometry_name[:-4]};
            cellZoneInside inside;
            boundary internal;
            faceType boundary;
        }}"""
                else:  # no patches. Just cellZone
                    refinementSurfaces += f"""
        {geometry_name[:-4]}
        {{
            level (0 0);
            cellZone {geometry_name[:-4]};
            faceZone {geometry_name[:-4]};
            cellZoneInside inside;
            boundary internal;
        }}"""
                # if refinementSurface or region, do not add here
            elif (geometry.purpose == 'refinementRegion' or geometry.purpose == 'refinementSurface'):
                pass

            elif (geometry.purpose == 'baffle'):
                patchType = 'wall'
                refinementSurfaces += f"""
        {geometry_name[:-4]}
        {{
            level (0 0);
            regions
            {{
                {geometry_name[:-4]}
                {{
                    faceType baffles;
                    faceZone {geometry_name[:-4]};
                    level ({geometry.refineMin} {geometry.refineMax});
                }}
            }}
        }}"""

            else:
                patchType = 'wall'
                refinementSurfaces += f"""
        {geometry_name[:-4]}
        {{
            level (0 0);
            regions
            {{
                {geometry_name[:-4]}
                {{
                    level ({geometry.refineMin} {geometry.refineMax});
                    patchInfo
                    {{
                        type {patchType};
                    }}
                }}
            }}
        }}"""

        # For searchable boxes, min and max are added
        elif (isinstance(geometry, SearchableBoxGeometry)):
            added_geo = f"""
    {geometry_name}
    {{
        type {geometry.type};
        min ({geometry.min[0]} {geometry.min[1]} {geometry.min[2]});
        max ({geometry.max[0]} {geometry.max[1]} {geometry.max[2]});
    }}"""
        geometry_str += added_geo
    geometry_str += f"""

}}"""

    refinementRegions = f""
    for geometry_name, geometry in meshSettings.geometry.items():
        if (isinstance(geometry, SearchableBoxGeometry)):
            refinementRegions += f"""
        {geometry_name}
        {{
            mode inside;
            levels ((1E15 {geometry.refineMax}));
        }}"""
        elif (isinstance(geometry, TriSurfaceMeshGeometry)):
            if (geometry.purpose == 'refinementSurface'):
                refinementRegions += f"""
        {geometry_name[:-4]}
        {{
            mode distance;
            levels ((1E-4 {geometry.property}));
        }}"""
            elif (geometry.purpose == 'refinementRegion'):
                refinementRegions += f"""
        {geometry_name[:-4]}
        {{
            mode inside;
            levels ((1E15 {geometry.property}));
        }}"""
            elif (geometry.purpose == 'cellZone'):
                refinementRegions += f"""
        {geometry_name[:-4]}
        {{
            mode inside;
            levels ((1E15 {geometry.property[0]}));
        }}"""

        else:
            pass

    castellatedMeshControls = f"""\ncastellatedMeshControls
{{
    maxLocalCells {meshSettings.castellatedMeshControls.maxLocalCells};
    maxGlobalCells {meshSettings.castellatedMeshControls.maxGlobalCells};
    minRefinementCells {meshSettings.castellatedMeshControls.minRefinementCells};
    maxLoadUnbalance {meshSettings.castellatedMeshControls.maxLoadUnbalance};
    nCellsBetweenLevels {meshSettings.castellatedMeshControls.nCellsBetweenLevels};
    features
    (
        {features}
    );
    refinementSurfaces
    {{
        {refinementSurfaces}
    }}
    resolveFeatureAngle {meshSettings.castellatedMeshControls.resolveFeatureAngle};
    refinementRegions
    {{
        {refinementRegions}
    }};
    locationInMesh ({meshSettings.castellatedMeshControls.locationInMesh[0]} {meshSettings.castellatedMeshControls.locationInMesh[1]} {meshSettings.castellatedMeshControls.locationInMesh[2]});
    allowFreeStandingZoneFaces {meshSettings.castellatedMeshControls.allowFreeStandingZoneFaces};
}}"""

    snapControls = f"""\nsnapControls
{{
    nSmoothPatch {meshSettings.snapControls.nSmoothPatch};
    tolerance {meshSettings.snapControls.tolerance};
    nSolveIter {meshSettings.snapControls.nSolveIter};
    nRelaxIter {meshSettings.snapControls.nRelaxIter};
    nFeatureSnapIter {meshSettings.snapControls.nFeatureSnapIter};
    implicitFeatureSnap {meshSettings.snapControls.implicitFeatureSnap};
    explicitFeatureSnap {meshSettings.snapControls.explicitFeatureSnap};
    multiRegionFeatureSnap {meshSettings.snapControls.multiRegionFeatureSnap};
}}"""
    layerControls = f"""\naddLayersControls
{{
    relativeSizes {meshSettings.addLayersControls.relativeSizes};
    layers
    {{"""
    for geometry_name, geometry in meshSettings.geometry.items():
        if (isinstance(geometry, TriSurfaceMeshGeometry)):
            if (geometry.purpose == 'wall'):  # If the surface is a wall, add layers
                layerControls += f"""
            "{geometry_name[:-4]}.*"
            {{
                nSurfaceLayers {geometry.nLayers};
            }}"""
            elif (geometry.purpose == 'baffle'):  # If the surface is a baffle, add layers
                layerControls += f"""
            "{geometry_name[:-4]}.*"
            {{
                nSurfaceLayers {1};
            }}"""
            elif (geometry.purpose == 'cellZone'):
                layerControls += f"""
            "{geometry_name[:-4]}.*"
            {{
                nSurfaceLayers {1};
            }}"""
            else:
                pass
    layerControls += f"""
    }};
    expansionRatio {meshSettings.addLayersControls.expansionRatio};
    finalLayerThickness {meshSettings.addLayersControls.finalLayerThickness};
    //firstLayerThickness {meshSettings.addLayersControls.firstLayerThickness};
    minThickness {meshSettings.addLayersControls.minThickness};
    nGrow {meshSettings.addLayersControls.nGrow};
    featureAngle {meshSettings.addLayersControls.featureAngle};
    slipFeatureAngle {meshSettings.addLayersControls.slipFeatureAngle};
    nRelaxIter {meshSettings.addLayersControls.nRelaxIter};
    nSmoothSurfaceNormals {meshSettings.addLayersControls.nSmoothSurfaceNormals};
    nSmoothNormals {meshSettings.addLayersControls.nSmoothNormals};
    nSmoothThickness {meshSettings.addLayersControls.nSmoothThickness};
    maxFaceThicknessRatio {meshSettings.addLayersControls.maxFaceThicknessRatio};
    maxThicknessToMedialRatio {meshSettings.addLayersControls.maxThicknessToMedialRatio};
    minMedianAxisAngle {meshSettings.addLayersControls.minMedianAxisAngle};
    minMedialAxisAngle {meshSettings.addLayersControls.minMedianAxisAngle};
    nBufferCellsNoExtrude {meshSettings.addLayersControls.nBufferCellsNoExtrude};
    nLayerIter {meshSettings.addLayersControls.nLayerIter};
}}"""
    meshQualityControls = f"""\nmeshQualityControls
{{
    maxNonOrtho {meshSettings.meshQualityControls.maxNonOrtho};
    maxBoundarySkewness {meshSettings.meshQualityControls.maxBoundarySkewness};
    maxInternalSkewness {meshSettings.meshQualityControls.maxInternalSkewness};
    maxConcave {meshSettings.meshQualityControls.maxConcave};
    minVol {meshSettings.meshQualityControls.minVol};
    minTetQuality {meshSettings.meshQualityControls.minTetQuality};
    minArea {meshSettings.meshQualityControls.minArea};
    minTwist {meshSettings.meshQualityControls.minTwist};
    minDeterminant {meshSettings.meshQualityControls.minDeterminant};
    minFaceWeight {meshSettings.meshQualityControls.minFaceWeight};
    minVolRatio {meshSettings.meshQualityControls.minVolRatio};
    minTriangleTwist {meshSettings.meshQualityControls.minTriangleTwist};
    nSmoothScale {meshSettings.meshQualityControls.nSmoothScale};
    errorReduction {meshSettings.meshQualityControls.errorReduction};
}}"""
    debug = f"""
writeFlags
(
    scalarLevels
    layerSets
    layerFields     // write volScalarField for layer coverage
);
debug {meshSettings.debug};
mergeTolerance {meshSettings.mergeTolerance};"""
    snappyHexMeshDict += header+steps+geometry_str+castellatedMeshControls + \
        snapControls+layerControls+meshQualityControls+debug
    return snappyHexMeshDict


# Example usage
if __name__ == "__main__":
    meshSettings = MeshSettings.model_validate(
        AmpersandUtils.yaml_to_dict("examples/basic/meshSettings.yaml")
    )

    snappy_hex_mesh_dict_content = create_snappyHexMeshDict(meshSettings)
    with open("outputs/snappyHexMeshDict", "w") as f:
        f.write(snappy_hex_mesh_dict_content)
    print("snappyHexMeshDict file created.")
