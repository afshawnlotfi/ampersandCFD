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
from pathlib import Path
from typing import Union
import vtk
import numpy as np
import math
from src.thirdparty.stlToOpenFOAM import find_inside_point, is_point_inside, read_stl_file
from src.thirdparty.stlToOpenFOAM import extract_curvature_data, compute_curvature
from src.primitives import AmpersandIO


class StlAnalysis:
    @staticmethod
    def roundl(x):
        return float(np.around(x, decimals=1))

    # to calculate the domain size for blockMeshDict

    @staticmethod
    def calc_domain_size(stlBoundingBox, sizeFactor=1, onGround=False,
                         internalFlow=False, halfModel=False):
        stlMinX, stlMaxX, stlMinY, stlMaxY, stlMinZ, stlMaxZ = stlBoundingBox
        # this part is for external flow
        bbX = stlMaxX - stlMinX
        bbY = stlMaxY - stlMinY
        bbZ = stlMaxZ - stlMinZ
        characLength = max(bbX, bbY, bbZ)
        minX = stlMinX - 3.0*characLength*sizeFactor
        maxX = stlMaxX + 9.0*characLength*sizeFactor
        minY = stlMinY - 2.0*characLength*sizeFactor
        maxY = stlMaxY + 2.0*characLength*sizeFactor
        minZ = stlMinZ - 2.0*characLength*sizeFactor
        maxZ = stlMaxZ + 2.0*characLength*sizeFactor

        if (internalFlow):
            minX = stlMinX - 0.1*bbX*sizeFactor
            maxX = stlMaxX + 0.1*bbX*sizeFactor
            minY = stlMinY - 0.1*bbY*sizeFactor
            maxY = stlMaxY + 0.1*bbY*sizeFactor
            minZ = stlMinZ - 0.1*bbZ*sizeFactor
            maxZ = stlMaxZ + 0.1*bbZ*sizeFactor

        if onGround:  # the the body is touching the ground
            minZ = stlMinZ
            maxZ = stlMaxZ + 4.0*characLength*sizeFactor
        if halfModel:
            maxY = (maxY+minY)/2.
        domain_size = (minX, maxX, minY, maxY, minZ, maxZ)
        return domain_size

    # to calculate the max length of STL
    @staticmethod
    def getMaxSTLDim(stlBoundingBox):
        stlMinX, stlMaxX, stlMinY, stlMaxY, stlMinZ, stlMaxZ = stlBoundingBox
        bbX = stlMaxX - stlMinX
        bbY = stlMaxY - stlMinY
        bbZ = stlMaxZ - stlMinZ
        return max(bbX, bbY, bbZ)

    # to calculate the min size of stl
    @staticmethod
    def getMinSTLDim(stlBoundingBox):
        stlMinX, stlMaxX, stlMinY, stlMaxY, stlMinZ, stlMaxZ = stlBoundingBox
        bbX = stlMaxX - stlMinX
        bbY = stlMaxY - stlMinY
        bbZ = stlMaxZ - stlMinZ
        return min(bbX, bbY, bbZ)

    # to calculate the refinement box for snappyHexMeshDict
    @staticmethod
    def getRefinementBox(stlBoundingBox):
        stlMinX, stlMaxX, stlMinY, stlMaxY, stlMinZ, stlMaxZ = stlBoundingBox
        bbX = stlMaxX - stlMinX
        bbY = stlMaxY - stlMinY
        bbZ = stlMaxZ - stlMinZ
        boxMinX = stlMinX - 0.7*bbX
        boxMaxX = stlMaxX + 15*bbX
        boxMinY = stlMinY - 1.0*bbY
        boxMaxY = stlMaxY + 1.0*bbY
        boxMinZ = stlMinZ - 1.0*bbZ
        boxMaxZ = stlMaxZ + 1.0*bbZ
        return (boxMinX, boxMaxX, boxMinY, boxMaxY, boxMinZ, boxMaxZ)

    @staticmethod
    def getRefinementBoxClose(stlBoundingBox):
        stlMinX, stlMaxX, stlMinY, stlMaxY, stlMinZ, stlMaxZ = stlBoundingBox
        bbX = stlMaxX - stlMinX
        bbY = stlMaxY - stlMinY
        bbZ = stlMaxZ - stlMinZ
        boxMinX = stlMinX - 0.2*bbX
        boxMaxX = stlMaxX + 3.0*bbX
        boxMinY = stlMinY - 0.45*bbY
        boxMaxY = stlMaxY + 0.45*bbY
        boxMinZ = stlMinZ - 0.45*bbZ
        boxMaxZ = stlMaxZ + 0.45*bbZ
        return (boxMinX, boxMaxX, boxMinY, boxMaxY, boxMinZ, boxMaxZ)

    # to add refinement box to mesh settings
    @staticmethod
    def addRefinementBoxToMesh(meshSettings, stl_path, boxName='refinementBox', refLevel=2, internalFlow=False):
        if (internalFlow):
            return meshSettings
        stlBoundingBox = StlAnalysis.compute_bounding_box(stl_path)
        box = StlAnalysis.getRefinementBox(stlBoundingBox)
        meshSettings['geometry'].append({'name': boxName, 'type': 'searchableBox', 'purpose': 'refinement',
                                         'min': [box[0], box[2], box[4]], 'max': [box[1], box[3], box[5]],
                                         'refineMax': refLevel-1})

        fineBox = StlAnalysis.getRefinementBoxClose(stlBoundingBox)
        meshSettings['geometry'].append({'name': 'fineBox', 'type': 'searchableBox', 'purpose': 'refinement',
                                         'min': [fineBox[0], fineBox[2], fineBox[4]], 'max': [fineBox[1], fineBox[3], fineBox[5]],
                                         'refineMax': refLevel})

        return meshSettings

    # refinement box for the ground for external automotive flows
    @staticmethod
    def addGroundRefinementBoxToMesh(meshSettings, stl_path, refLevel=2):
        # if(internalFlow):
        #    return meshSettings
        boxName = 'groundBox'
        stlBoundingBox = StlAnalysis.compute_bounding_box(stl_path)
        xmin, xmax, ymin, ymax, zmin, zmax = stlBoundingBox
        z = meshSettings['domain']['minz']
        z_delta = 0.2*(zmax-zmin)
        box = [-1000.0, 1000., -1000, 1000, z-z_delta, z+z_delta]
        meshSettings['geometry'].append({'name': boxName, 'type': 'searchableBox', 'purpose': 'refinement',
                                         'min': [box[0], box[2], box[4]], 'max': [box[1], box[3], box[5]],
                                         'refineMax': refLevel})
        return meshSettings

    # to calculate nearest wall thickness for a target yPlus value
    @staticmethod
    def calc_y(nu=1e-6, rho=1000., L=1.0, u=1.0, target_yPlus=200):
        # rho = fluid_properties['rho']
        # nu = fluid_properties['nu']
        Re = u*L/nu
        Cf = 0.0592*Re**(-1./5.)
        tau = 0.5*rho*Cf*u**2.
        uStar = np.sqrt(tau/rho)
        y = target_yPlus*nu/uStar
        return y

    # to calculate yPlus value for a given first layer thickness
    @staticmethod
    def calc_yPlus(nu=1e-6, L=1.0, u=1.0, y=0.001):
        Re = u*L/nu
        Cf = 0.0592*Re**(-1./5.)
        tau = 0.5*Cf*u**2.
        uStar = np.sqrt(tau)
        yPlus = uStar*y/nu
        return yPlus

    # calculate nearest cell size for a given expansion ratio and layer count
    @staticmethod
    def calc_cell_size(y_=0.001, nLayers=5, expRatio=1.2, thicknessRatio=0.3):
        max_y = y_*expRatio**(nLayers)
        return max_y/thicknessRatio

    @staticmethod
    def calc_refinement_levels(max_cell_size=0.1, target_cell_size=0.001):
        size_ratio = max_cell_size / target_cell_size
        n = np.log(size_ratio)/np.log(2.)
        # print(n)
        return int(np.ceil(n))

    @staticmethod
    def calc_nx_ny_nz(domain_size, target_cell_size):
        (minX, maxX, minY, maxY, minZ, maxZ) = domain_size
        nx = (maxX-minX)/target_cell_size
        ny = (maxY-minY)/target_cell_size
        nz = (maxZ-minZ)/target_cell_size
        nx, ny, nz = int(math.ceil(nx)), int(math.ceil(ny)), int(math.ceil(nz))
        # it is better to have even number of cells
        if nx // 2 != 0:
            nx += 1
        if ny // 2 != 0:
            ny += 1
        if nz // 2 != 0:
            nz += 1
        return (nx, ny, nz)

    # Function to read STL file and compute bounding box
    @staticmethod
    def compute_bounding_box(stl_file_path):
        # Check if the file exists
        if not os.path.exists(stl_file_path):
            raise FileNotFoundError(
                f"File not found: {stl_file_path}. Make sure the file exists.")
        # Create a reader for the STL file
        reader = vtk.vtkSTLReader()
        reader.SetFileName(stl_file_path)
        reader.Update()

        # Get the output data from the reader
        poly_data = reader.GetOutput()

        # Calculate the bounding box
        bounds = poly_data.GetBounds()
        # xmin, xmax, ymin, ymax, zmin, zmax = bounds
        # Optionally, return the bounding box as a tuple
        return bounds

    # this is the wrapper function to check if a point is inside the mesh
    @staticmethod
    def is_point_inside(stl_file_path, point):
        # Check if the file exists
        if not os.path.exists(stl_file_path):
            raise FileNotFoundError(
                f"File not found: {stl_file_path}. Make sure the file exists.")
        # Create a reader for the STL file
        reader = vtk.vtkSTLReader()
        reader.SetFileName(stl_file_path)
        reader.Update()

        # Get the output data from the reader
        poly_data = reader.GetOutput()
        # Calculate the bounding box
        bounds = poly_data.GetBounds()
        # Check if the point is inside the bounding box
        xmin, xmax, ymin, ymax, zmin, zmax = bounds
        if point[0] < xmin or point[0] > xmax:
            return False
        if point[1] < ymin or point[1] > ymax:
            return False
        if point[2] < zmin or point[2] > zmax:
            return False
        # Check if the point is inside the mesh
        return is_point_inside(poly_data, point)

    @staticmethod
    def read_stl(stl_file_path):
        # Check if the file exists
        if not os.path.exists(stl_file_path):
            raise FileNotFoundError(
                f"File not found: {stl_file_path}. Make sure the file exists.")
        # Create a reader for the STL file
        reader = vtk.vtkSTLReader()
        reader.SetFileName(stl_file_path)
        reader.Update()

        # Get the output data from the reader
        poly_data = reader.GetOutput()
        return poly_data

    @staticmethod
    def calc_nLayer(yFirst=0.001, targetCellSize=0.1, expRatio=1.2):
        n = np.log(targetCellSize*0.4/yFirst)/np.log(expRatio)
        return int(np.ceil(n))

    @staticmethod
    def calc_delta(U=1.0, nu=1e-6, L=1.0):
        Re = U*L/nu
        delta = 0.37*L/Re**(0.2)
        return delta

    @staticmethod
    # calculates N layers and final layer thickness
    # yFirst: first layer thickness
    # delta: boundary layer thickness
    # expRatio: expansion ratio
    def calc_layers(yFirst=0.001, delta=0.01, expRatio=1.2):
        currentThickness = yFirst*2.0  # initial thickness. Twice the yPlus value
        currentDelta = 0
        N = 0
        for i in range(1, 50):
            currentThickness = currentThickness*expRatio**(i)
            currentDelta = currentDelta + currentThickness
            if (currentDelta > delta):
                N = i
                break
        finalLayerThickness = currentThickness

        return N, finalLayerThickness

    @staticmethod
    def calc_layers_from_cell_size(yFirst=0.001, targetCellSize=0.1, expRatio=1.2):

        firstLayerThickness = yFirst*2.0
        finalLayerThickness = targetCellSize*0.35
        nLayers = int(np.log(finalLayerThickness /
                      firstLayerThickness)/np.log(expRatio))
        nLayers = max(1, nLayers)
        return nLayers, finalLayerThickness

    # this function calculates the smallest curvature of the mesh
    # This function calls stlToOpenFOAM functions to read the mesh and calculate curvature

    @staticmethod
    def calc_smallest_curvature(stlFile):
        mesh = read_stl_file(stlFile)
        curved_mesh = compute_curvature(mesh, curvature_type='mean')
        curvature_values = extract_curvature_data(curved_mesh)
        print(f"Curvature values: {curvature_values}")
        min_curvature = np.min(curvature_values)
        return min_curvature

    # to calculate the mesh settings for blockMeshDict and snappyHexMeshDict

    @staticmethod
    def calc_mesh_settings(stlBoundingBox, nu=1e-6, rho=1000., U=1.0, maxCellSize=0.5, sizeFactor=1.0,
                           expansion_ratio=1.5, onGround=False, internalFlow=False, refinement=1,
                           nLayers=5, halfModel=False, thicknessRatio=0.3, GUI=False):
        maxSTLLength = StlAnalysis.getMaxSTLDim(stlBoundingBox)
        minSTLLength = StlAnalysis.getMinSTLDim(stlBoundingBox)
        if (maxCellSize < 0.001):
            maxCellSize = maxSTLLength/4.
        domain_size = StlAnalysis.calc_domain_size(stlBoundingBox=stlBoundingBox, sizeFactor=sizeFactor,
                                                   onGround=onGround, internalFlow=internalFlow, halfModel=halfModel)
        if (refinement == 0):
            if (internalFlow):
                if maxSTLLength/minSTLLength > 10:  # if the geometry is very slender
                    backgroundCellSize = min(maxSTLLength/50., maxCellSize)
                else:
                    backgroundCellSize = min(minSTLLength/8., maxCellSize)
            else:
                # this is the size of largest blockMesh cells
                backgroundCellSize = min(minSTLLength/3., maxCellSize)
            target_yPlus = 70
            # nLayers = 2
            refLevel = 2
        elif (refinement == 1):
            if (internalFlow):
                if maxSTLLength/minSTLLength > 10:  # if the geometry is very slender
                    backgroundCellSize = min(maxSTLLength/70., maxCellSize)
                else:
                    backgroundCellSize = min(minSTLLength/12., maxCellSize)
            else:
                backgroundCellSize = min(minSTLLength/5., maxCellSize)
            target_yPlus = 50
            # nLayers = 5
            refLevel = 4
        elif (refinement == 2):
            if (internalFlow):
                if maxSTLLength/minSTLLength > 10:  # if the geometry is very slender
                    backgroundCellSize = min(maxSTLLength/90., maxCellSize)
                else:
                    backgroundCellSize = min(minSTLLength/16., maxCellSize)
            else:
                backgroundCellSize = min(minSTLLength/7., maxCellSize)
            target_yPlus = 30
            # nLayers = 7
            refLevel = 6
        else:  # medium settings for default
            if (internalFlow):
                backgroundCellSize = min(maxSTLLength/12., maxCellSize)
            else:
                backgroundCellSize = min(maxSTLLength/8., maxCellSize)
            # nLayers = 5
            target_yPlus = 70
            refLevel = 4

        nx, ny, nz = StlAnalysis.calc_nx_ny_nz(domain_size, backgroundCellSize)
        backgroundCellSize = (domain_size[1]-domain_size[0])/nx
        L = maxSTLLength  # this is the characteristic length to be used in Re calculations
        # this is the thickness of closest cell
        target_y = StlAnalysis.calc_y(nu, rho, L, U, target_yPlus=target_yPlus)
        delta = StlAnalysis.calc_delta(U, nu, L)
        # nLayers,finalLayerThickness = stlAnalysis.calc_layers(yFirst=target_y,delta=delta,expRatio=expansion_ratio)
        if (refinement == 0):
            refLevel = max(2, refLevel)
        elif (refinement == 1):
            refLevel = max(4, refLevel)
        elif (refinement == 2):
            refLevel = max(6, refLevel)
        else:
            refLevel = max(3, refLevel)
        targetCellSize = backgroundCellSize/2.**refLevel
        nLayers, finalLayerThickness = StlAnalysis.calc_layers_from_cell_size(
            yFirst=target_y, targetCellSize=targetCellSize, expRatio=expansion_ratio)
        nLayers = max(1, nLayers)

        # adjust refinement levels based on coarse, medium, fine settings

        adjustedNearWallThickness = finalLayerThickness / \
            expansion_ratio**(nLayers-1)
        adjustedYPlus = StlAnalysis.calc_yPlus(
            nu, L, U, adjustedNearWallThickness/2.)

        # minVolumeSize = backgroundCellSize**3/(8.**refLevel*20.)
        # print the summary of results
        AmpersandIO.printMessage(
            "\n-----------------Mesh Settings-----------------")
        AmpersandIO.printMessage(f"Domain size: x({domain_size[0]:6.3f}~{domain_size[1]:6.3f}) y({domain_size[2]:6.3f}~{
                                 domain_size[3]:6.3f}) z({domain_size[4]:6.3f}~{domain_size[5]:6.3f})")
        AmpersandIO.printMessage(f"Nx Ny Nz: {nx},{ny},{nz}")
        AmpersandIO.printMessage(f"Max cell size: {backgroundCellSize}")
        AmpersandIO.printMessage(f"Min cell size: {targetCellSize}")
        AmpersandIO.printMessage(f"Refinement Level:{refLevel}")

        AmpersandIO.printMessage(
            "\n-----------------Turbulence-----------------")
        AmpersandIO.printMessage(f"Target yPlus:{target_yPlus}")
        AmpersandIO.printMessage(f'Reynolds number:{U*L/nu}')
        AmpersandIO.printMessage(f"Boundary layer thickness:{delta}")
        AmpersandIO.printMessage(f"First layer thickness:{adjustedNearWallThickness}")
        AmpersandIO.printMessage(f"Final layer thickness:{finalLayerThickness}")
        AmpersandIO.printMessage(f"YPlus:{adjustedYPlus}")
        AmpersandIO.printMessage(f"Number of layers:{nLayers}")
        return domain_size, nx, ny, nz, refLevel, finalLayerThickness, nLayers

    @staticmethod
    def set_layer_thickness(meshSettings, thickness=0.01):
        meshSettings['addLayersControls']['finalLayerThickness'] = thickness
        minThickness = max(0.0001, thickness/100.)
        meshSettings['addLayersControls']['minThickness'] = minThickness
        return meshSettings

    @staticmethod
    def set_min_vol(meshSettings, minVol=1e-15):
        meshSettings['meshQualityControls']['minVol'] = 1e-15  # minVol/100.
        return meshSettings

    # to set mesh settings for blockMeshDict and snappyHexMeshDict
    @staticmethod
    def set_mesh_settings(meshSettings, domain_size, nx, ny, nz, refLevel, featureLevel=1, nLayers=None):
        meshSettings['domain'] = {'minx': domain_size[0], 'maxx': domain_size[1], 'miny': domain_size[2],
                                  'maxy': domain_size[3], 'minz': domain_size[4], 'maxz': domain_size[5], 'nx': nx, 'ny': ny, 'nz': nz}
        # meshSettings['domain']['nx'] = nx
        # meshSettings['domain']['ny'] = ny
        # meshSettings['domain']['nz'] = nz

        refMin = max(1, refLevel)
        refMax = max(2, refLevel)
        for geometry in meshSettings['geometry']:
            if geometry['type'] == 'triSurfaceMesh':
                geometry['refineMin'] = refMin
                geometry['refineMax'] = refMax
                # geometry['featureEdges'] = 'true'
                geometry['featureLevel'] = featureLevel
                if nLayers is not None:
                    geometry['nLayers'] = nLayers
        return meshSettings

    @staticmethod
    def calc_center_of_mass(mesh: vtk.vtkPolyData):
        center_of_mass_filter = vtk.vtkCenterOfMass()
        center_of_mass_filter.SetInputData(mesh)
        center_of_mass_filter.Update()
        center_of_mass = center_of_mass_filter.GetCenter()
        return center_of_mass

    @staticmethod
    def analyze_stl(stl_file_path):
        mesh = StlAnalysis.read_stl(stl_file_path)
        bounds = StlAnalysis.compute_bounding_box(stl_file_path)
        stlMinX, stlMaxX, stlMinY, stlMaxY, stlMinZ, stlMaxZ = bounds
        outsideX = stlMaxX + 0.05*(stlMaxX-stlMinX)
        outsideY = stlMinY*0.95  # (stlMaxY - stlMinY)/2.
        outsideZ = (stlMaxZ - stlMinZ)/2.
        outsidePoint = (outsideX, outsideY, outsideZ)
        center_of_mass = StlAnalysis.calc_center_of_mass(mesh)
        insidePoint = find_inside_point(
            mesh, center_of_mass, min_bounds=None, max_bounds=None)
        return center_of_mass, insidePoint, outsidePoint

    @staticmethod
    def set_mesh_location(meshSettings, stl_file_path, internalFlow=False):
        center_of_mass, insidePoint, outsidePoint = StlAnalysis.analyze_stl(stl_file_path)
        if internalFlow:
            meshSettings['castellatedMeshControls']['locationInMesh'] = [
                insidePoint[0], insidePoint[1], insidePoint[2]]
        else:
            meshSettings['castellatedMeshControls']['locationInMesh'] = [
                outsidePoint[0], outsidePoint[1], outsidePoint[2]]
        return meshSettings

    @staticmethod 
    def set_stl_solid_name(stl_file: Union[str, Path]) -> int:
        """
        Sets/updates the solid name in an STL file to match the filename.
        
        Args:
            stl_file: Path to STL file as string or Path object
            
        Returns:
            0 on success, -1 on failure
            
        Raises:
            FileNotFoundError: If STL file does not exist
            PermissionError: If lacking write permissions
            UnicodeDecodeError: If file has invalid encoding
        """
        try:
            stl_path = Path(stl_file)
            AmpersandIO.printMessage(f"Setting solid name for {stl_path}")
            
            if not stl_path.exists():
                raise FileNotFoundError(f"STL file not found: {stl_path}")

            if not stl_path.is_file():
                raise ValueError(f"Path is not a file: {stl_path}")

            # Extract solid name from filename without extension
            solid_name = stl_path.stem

            # Read and process file contents
            lines = stl_path.read_text().splitlines()
            new_lines = []
            for line in lines:
                if 'endsolid' in line.lower():
                    line = f"endsolid {solid_name}"
                elif line.lower().lstrip().startswith('solid'):
                    line = f"solid {solid_name}"
                new_lines.append(line + '\n')
            
            # Write back to file
            stl_path.write_text(''.join(new_lines))
            return 0

        except (FileNotFoundError, PermissionError, UnicodeDecodeError) as e:
            AmpersandIO.printMessage(f"Error processing STL file: {e}")
            return -1


if __name__ == "__main__":
    stl_file = "stl/ahmed.stl"
    minCurv = StlAnalysis.calc_smallest_curvature(stl_file)
