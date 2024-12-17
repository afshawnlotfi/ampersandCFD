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

from ampersand.models.settings import SearchableBoxGeometry, TriSurfaceMeshGeometry, SnappyHexMeshSettings
from primitives import AmpersandUtils
from ampersand.utils.generation import GenerationUtils


def create_snappyHexMeshDict(meshSettings: SnappyHexMeshSettings):
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
    meshSettings = SnappyHexMeshSettings.model_validate(
        AmpersandUtils.yaml_to_dict("examples/basic/meshSettings.yaml")
    )

    snappy_hex_mesh_dict_content = create_snappyHexMeshDict(meshSettings)
    with open("outputs/snappyHexMeshDict", "w") as f:
        f.write(snappy_hex_mesh_dict_content)
    print("snappyHexMeshDict file created.")


from typing import Dict, Any, List, Union

def create_snappy_hex_mesh_dict_python(mesh_settings):
    """
    Convert SnappyHexMeshSettings to a Python dictionary representation.

    Parameters:
    mesh_settings (SnappyHexMeshSettings): Mesh generation settings object

    Returns:
    Dict[str, Any]: A Python dictionary representing the snappyHexMeshDict configuration
    """
    # Initialize the dictionary to store mesh configuration
    snappy_hex_mesh_dict = {}

    # Add mesh generation steps
    snappy_hex_mesh_dict['steps'] = {
        'castellatedMesh': mesh_settings.snappyHexSteps.castellatedMesh,
        'snap': mesh_settings.snappyHexSteps.snap,
        'addLayers': mesh_settings.snappyHexSteps.addLayers
    }

    # Process geometry
    snappy_hex_mesh_dict['geometry'] = {}
    features = []
    refinement_surfaces = {}
    refinement_regions = {}

    for geometry_name, geometry in mesh_settings.geometry.items():
        # Process TriSurfaceMesh geometries
        if hasattr(geometry, 'type') and geometry.type == 'triSurfaceMesh':
            # Prepare geometry details
            geo_details = {
                'type': geometry.type,
                'name': geometry_name[:-4],
                'purpose': geometry.purpose
            }
            
            # Add to main geometry dictionary
            snappy_hex_mesh_dict['geometry'][geometry_name] = geo_details

            # Process feature edges
            if getattr(geometry, 'featureEdges', False):
                features.append({
                    'file': f"{geometry_name[:-4]}.eMesh",
                    'level': geometry.featureLevel
                })

            # Process refinement surfaces
            refinement_surface_config = {
                'level': [0, 0],
                'refineMin': geometry.refineMin if hasattr(geometry, 'refineMin') else 0,
                'refineMax': geometry.refineMax if hasattr(geometry, 'refineMax') else 0
            }

            # Add special handling for different purposes
            if geometry.purpose == 'inlet' or geometry.purpose == 'outlet':
                refinement_surface_config['patchType'] = 'patch'
            elif geometry.purpose == 'baffle':
                refinement_surface_config['patchType'] = 'baffles'
            elif geometry.purpose in ['wall', None]:
                refinement_surface_config['patchType'] = 'wall'
            
            refinement_surfaces[geometry_name[:-4]] = refinement_surface_config

            # Process refinement regions
            if geometry.purpose in ['refinementSurface', 'refinementRegion', 'cellZone']:
                region_config = {
                    'mode': 'inside' if geometry.purpose in ['refinementRegion', 'cellZone'] else 'distance',
                    'level': geometry.property if hasattr(geometry, 'property') else 0
                }
                refinement_regions[geometry_name[:-4]] = region_config

        # Process SearchableBox geometries
        elif hasattr(geometry, 'type') and geometry.type == 'searchableBox':
            box_details = {
                'type': geometry.type,
                'min': geometry.min,
                'max': geometry.max
            }
            snappy_hex_mesh_dict['geometry'][geometry_name] = box_details

            # Add refinement region for searchable boxes
            refinement_regions[geometry_name] = {
                'mode': 'inside',
                'level': geometry.refineMax if hasattr(geometry, 'refineMax') else 0
            }

    # Castellated Mesh Controls
    snappy_hex_mesh_dict['castellatedMeshControls'] = {
        key: getattr(mesh_settings.castellatedMeshControls, key) 
        for key in [
            'maxLocalCells', 'maxGlobalCells', 'minRefinementCells', 
            'maxLoadUnbalance', 'nCellsBetweenLevels', 'resolveFeatureAngle',
            'locationInMesh', 'allowFreeStandingZoneFaces'
        ]
    }

    # Add extracted collections
    snappy_hex_mesh_dict['features'] = features
    snappy_hex_mesh_dict['refinementSurfaces'] = refinement_surfaces
    snappy_hex_mesh_dict['refinementRegions'] = refinement_regions

    # Snap Controls
    snappy_hex_mesh_dict['snapControls'] = {
        key: getattr(mesh_settings.snapControls, key)
        for key in [
            'nSmoothPatch', 'tolerance', 'nSolveIter', 'nRelaxIter', 
            'nFeatureSnapIter', 'implicitFeatureSnap', 'explicitFeatureSnap', 
            'multiRegionFeatureSnap'
        ]
    }

    # Add Layers Controls
    snappy_hex_mesh_dict['addLayersControls'] = {
        key: getattr(mesh_settings.addLayersControls, key)
        for key in [
            'relativeSizes', 'expansionRatio', 'finalLayerThickness', 
            'minThickness', 'nGrow', 'featureAngle', 'slipFeatureAngle', 
            'nRelaxIter', 'nSmoothSurfaceNormals', 'nSmoothNormals', 
            'nSmoothThickness', 'maxFaceThicknessRatio', 
            'maxThicknessToMedialRatio', 'minMedianAxisAngle', 
            'nBufferCellsNoExtrude', 'nLayerIter'
        ]
    }

    # Add layer specifications for specific geometries
    layer_specs = {}
    for geometry_name, geometry in mesh_settings.geometry.items():
        if hasattr(geometry, 'purpose') and geometry.purpose in ['wall', 'baffle', 'cellZone']:
            layer_specs[f"{geometry_name[:-4]}.*"] = {
                'nSurfaceLayers': getattr(geometry, 'nLayers', 1)
            }
    snappy_hex_mesh_dict['addLayersControls']['layers'] = layer_specs

    # Mesh Quality Controls
    snappy_hex_mesh_dict['meshQualityControls'] = {
        key: getattr(mesh_settings.meshQualityControls, key)
        for key in [
            'maxNonOrtho', 'maxBoundarySkewness', 'maxInternalSkewness', 
            'maxConcave', 'minVol', 'minTetQuality', 'minArea', 'minTwist', 
            'minDeterminant', 'minFaceWeight', 'minVolRatio', 
            'minTriangleTwist', 'nSmoothScale', 'errorReduction'
        ]
    }

    # Debug and merge tolerance
    snappy_hex_mesh_dict['debug'] = mesh_settings.debug
    snappy_hex_mesh_dict['mergeTolerance'] = mesh_settings.mergeTolerance

    return snappy_hex_mesh_dict

# Example usage
if __name__ == "__main__":
    from ampersand.models.settings import SnappyHexMeshSettings
    from primitives import AmpersandUtils

    # Load mesh settings from a YAML file
    mesh_settings = SnappyHexMeshSettings.model_validate(
        AmpersandUtils.yaml_to_dict("examples/basic/meshSettings.yaml")
    )

    # Convert to Python dictionary
    python_dict = create_snappy_hex_mesh_dict_python(mesh_settings)
    
    # Optionally, print or save the dictionary
    import json
    print(json.dumps(python_dict, indent=2))