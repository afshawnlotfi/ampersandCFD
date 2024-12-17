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
from typing import Literal, Sequence
import vtk
import numpy as np
import math
from ampersand.models.settings import AddLayersControls, BCPatch, BoundingBox, MeshSettings, SnappyHexMeshSettings, SearchableBoxGeometry
from thirdparty.stlToOpenFOAM import find_inside_point, find_outside_point, is_point_inside, read_stl_file
from thirdparty.stlToOpenFOAM import extract_curvature_data, compute_curvature
from primitives import AmpersandIO


RefinementAmount = Literal["coarse", "medium", "fine"]

# @staticmethod
# def find_outside_point(mesh):
#     bounds = STLAnalysis.compute_bounding_box(stl_file_path)
#     outsideX = bounds.maxx + 0.05*(bounds.maxx-bounds.minx)
#     outsideY = bounds.miny*0.95  # (stlMaxY - stlMinY)/2.
#     outsideZ = (bounds.maxz - bounds.minz)/2.
#     outsidePoint = (outsideX, outsideY, outsideZ)
#     return np.array(outsidePoint)

class STLAnalysis:

    # to calculate the domain size for blockMeshDict
    @staticmethod
    def calc_domain_bbox(stl_bbox: BoundingBox, size_factor=1.0, on_ground=False,
                         internal_flow=False, half_model=False):
        max_length = stl_bbox.max_length
        max_size_factor = max_length*size_factor

        if (internal_flow):
            domain_bbox = stl_bbox.scale_dimensions(-0.1*size_factor, 0.1*size_factor, -0.1*size_factor, 0.1*size_factor, -0.1*size_factor, 0.1*size_factor)
        else:
            domain_bbox = stl_bbox.scale_dimensions(-3*max_size_factor, 9*max_size_factor, -2*max_size_factor, 2*max_size_factor, -2*max_size_factor, 2*max_size_factor)


        if on_ground:  # the the body is touching the ground
            domain_bbox.minz = stl_bbox.minz
            domain_bbox.maxz = stl_bbox.maxz + 4.0*max_size_factor
        
        if half_model:
            domain_bbox.maxy = (domain_bbox.maxy+domain_bbox.miny)/2.
        return domain_bbox


    # to calculate the refinement box for snappyHexMeshDict
    @staticmethod
    def getRefinementBox(stl_bbox: BoundingBox):
        return stl_bbox.scale_dimensions(-0.7, 15.0, -1.0, 1.0, -1.0, 1.0)


    @staticmethod
    def getRefinementBoxClose(stl_bbox: BoundingBox):
        return stl_bbox.scale_dimensions(-0.2, 3.0, -0.45, 0.45, -0.45, 0.45)

    # to calculate nearest wall thickness for a target yPlus value
    @staticmethod
    def calc_y_first(nu=1e-6, rho=1000., L=1.0, u=1.0, target_y_plus=200):
        Re = u*L/nu
        Cf = 0.0592*Re**(-1./5.)
        tau = 0.5*rho*Cf*u**2.
        uStar = np.sqrt(tau/rho)
        y = target_y_plus*nu/uStar
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
    def calc_cell_size(y_=0.001, n_layers=5, exp_ratio=1.2, thickness_ratio=0.3):
        max_y = y_*exp_ratio**(n_layers)
        return max_y/thickness_ratio

    @staticmethod
    def calc_refinement_levels(max_cell_size=0.1, target_cell_size=0.001):
        size_ratio = max_cell_size / target_cell_size
        n = np.log(size_ratio)/np.log(2.)
        # print(n)
        return int(np.ceil(n))

    @staticmethod
    def calc_nx_ny_nz(domain_bbox: BoundingBox, target_cell_size: float):
        nx = (domain_bbox.maxx-domain_bbox.minx)/target_cell_size
        ny = (domain_bbox.maxy-domain_bbox.miny)/target_cell_size
        nz = (domain_bbox.maxz-domain_bbox.minz)/target_cell_size
        nx, ny, nz = int(math.ceil(nx)), int(math.ceil(ny)), int(math.ceil(nz))
        # it is better to have even number of cells
        if nx % 2:  # if nx is odd
            nx += 1
        if ny % 2:  # if ny is odd
            ny += 1
        if nz % 2:  # if nz is odd
            nz += 1
        return (nx, ny, nz)

    # Function to read STL file and compute bounding box
    @staticmethod
    def compute_bounding_box(mesh):
        # Calculate the bounding box
        bounds = mesh.GetBounds()
        # xmin, xmax, ymin, ymax, zmin, zmax = bounds
        # Optionally, return the bounding box as a tuple
        return BoundingBox(
            minx=bounds[0],
            maxx=bounds[1],
            miny=bounds[2],
            maxy=bounds[3],
            minz=bounds[4],
            maxz=bounds[5]
        )

    # this is the wrapper function to check if a point is inside the mesh
    @staticmethod
    def is_point_inside(stl_file_path: str, point: Sequence[float]):
        # Check if the file exists
        if not os.path.exists(stl_file_path):
            raise FileNotFoundError(
                f"File not found: {stl_file_path}. Make sure the file exists.")
        # Create a reader for the STL file
        reader = vtk.vtkSTLReader()  # type: ignore
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
    def calc_num_layers(final_layer_thickness: float, yFirst=0.001, expansion_ratio=1.2):
        firstLayerThickness = yFirst*2.0
        return max(1, int(np.log(final_layer_thickness / firstLayerThickness)/np.log(expansion_ratio)))

    # this function calculates the smallest curvature of the mesh
    # This function calls stlToOpenFOAM functions to read the mesh and calculate curvature
    @staticmethod
    def calc_smallest_curvature(mesh):
        curved_mesh = compute_curvature(mesh, curvature_type='mean')
        curvature_values = extract_curvature_data(curved_mesh)
        print(f"Curvature values: {curvature_values}")
        min_curvature = np.min(curvature_values)
        return min_curvature

    # to calculate the mesh settings for blockMeshDict and snappyHexMeshDict
    @staticmethod
    def calc_background_cell_size(refinement_amount: RefinementAmount, domain_bbox: BoundingBox, maxCellSize: float, internalFlow: bool):
        max_length = domain_bbox.max_length
        min_length = domain_bbox.min_length
        
        if (refinement_amount == "coarse"):
            if (internalFlow):
                if max_length/min_length > 10:  # if the geometry is very slender
                    return min(max_length/50., maxCellSize)
                else:
                    return min(min_length/8., maxCellSize)
            else:
                # this is the size of largest blockMesh cells
                return min(min_length/3., maxCellSize)
        elif (refinement_amount == "medium"):
            if (internalFlow):
                if max_length/min_length > 10:  # if the geometry is very slender
                    return min(max_length/70., maxCellSize)
                else:
                    return min(min_length/12., maxCellSize)
            else:
                return min(min_length/5., maxCellSize)
        else:
            if (internalFlow):
                if max_length/min_length > 10:  # if the geometry is very slender
                    return min(max_length/90., maxCellSize)
                else:
                    return min(min_length/16., maxCellSize)
            else:
                return min(min_length/7., maxCellSize)

    @staticmethod
    def calc_ref_level(refinement_amount: RefinementAmount):

        ref_levels_options: dict[RefinementAmount, int] = {
            "coarse": 2,
            "medium": 4,
            "fine": 6,
        }

        if (refinement_amount == "coarse"):
            return max(2, ref_levels_options[refinement_amount])
        elif (refinement_amount == "medium"):
            return max(4, ref_levels_options[refinement_amount])
        elif (refinement_amount == "fine"):
            return max(6, ref_levels_options[refinement_amount])


    @staticmethod
    def calc_mesh_settings(
        stl_bbox: BoundingBox, 
        nu=1e-6, 
        rho=1000.0, 
        U=1.0, 
        max_cell_size=0.5, 
        size_factor=1.0,
        expansion_ratio=1.5, 
        on_ground=False, 
        internal_flow=False, 
        refinement_index=1, 
        half_model=False
    ):
        max_bbox_length = stl_bbox.max_length

        if (max_cell_size < 0.001):
            max_cell_size = max_bbox_length/4.
        domain_bbox = STLAnalysis.calc_domain_bbox(stl_bbox, size_factor,on_ground, internal_flow, half_model)

        refinement_options: list[RefinementAmount] = ["coarse", "medium", "fine"]  # type: ignore
        refinement_amount: RefinementAmount = refinement_options[refinement_index]


        background_cell_size = STLAnalysis.calc_background_cell_size(refinement_amount, domain_bbox, max_cell_size, internal_flow)
        nx, ny, nz = STLAnalysis.calc_nx_ny_nz(domain_bbox, background_cell_size)
        
        x_background_cell_size = (domain_bbox.maxx-domain_bbox.minx)/nx
        L = max_bbox_length  # this is the characteristic length to be used in Re calculations
        

        # Calculate Y Plus
        target_y_plus_options: dict[RefinementAmount, int] = {
            "coarse": 70,
            "medium": 50,
            "fine": 30,
        }
        target_y_plus = target_y_plus_options[refinement_amount]
        
        # this is the thickness of closest cell
        target_y_first = STLAnalysis.calc_y_first(nu, rho, L, U, target_y_plus)

        ref_level =  STLAnalysis.calc_ref_level(refinement_amount)
        target_cell_size = x_background_cell_size/2.**ref_level
        final_layer_thickness = target_cell_size*0.35


        num_layers = STLAnalysis.calc_num_layers(final_layer_thickness, target_y_first, expansion_ratio)

        # adjust refinement levels based on coarse, medium, fine settings
        # adjustedNearWallThickness = final_layer_thickness / expansion_ratio**(nLayers-1)
        # adjustedYPlus = STLAnalysis.calc_yPlus(nu, L, U, adjustedNearWallThickness/2.)

        featureLevel = max(ref_level,1)
        refMin = max(1, ref_level)
        refMax = max(2, ref_level)
        

        # return domain_bbox, nx, ny, nz, ref_level, final_layer_thickness, nLayers



    @staticmethod
    def calc_layer_controls(thickness=0.01):
        return AddLayersControls(
            finalLayerThickness=thickness,
            minThickness=max(0.0001, thickness / 100.),
        )


  # # to add refinement box to mesh settings
    # @staticmethod
    # def addRefinementBoxToMesh(meshSettings: SnappyHexMeshSettings, stl_path: str, boxName='refinementBox', refLevel=2, internalFlow=False):
    #     if internalFlow:
    #         return meshSettings
    #     stlBoundingBox = STLAnalysis.compute_bounding_box(stl_path)
    #     box = STLAnalysis.getRefinementBox(stlBoundingBox)
    #     meshSettings.geometry.append(
    #         SearchableBoxGeometry(
    #             type='searchableBox', 
    #             name=boxName
    #             purpose='refinement', 
    #             property=({'min': [box[0], box[2], box[4]], 'max': [box[1], box[3], box[5]], 'refineMax':  refLevel-1})
    #         )
    #     )

    #     fineBox = STLAnalysis.getRefinementBoxClose(stlBoundingBox)
    #     meshSettings.geometry['fineBox'] = refineMax': refLevel(
    #         purpose='refinement', property=({'min': [fineBox[0], fineBox[2], fineBox[4]], 'max': [fineBox[1], fineBox[3], fineBox[5]], 'refineMax': refLevel})
    #     )

    #     return meshSettings

    # # TODO: this is assuming the z is the the directrion of up
    # # refinement box for the ground for external automotive flows
    # @staticmethod
    # def addGroundRefinementBoxToMesh(meshSettings: SnappyHexMeshSettings, stl_path: str, refLevel=2):
    #     boxName = 'groundBox'
    #     stlBoundingBox = STLAnalysis.compute_bounding_box(stl_path)
    #     z = meshSettings.domain.minz
    #     z_delta = 0.2 * (stlBoundingBox.maxz - stlBoundingBox.minz)
    #     box = [-1000.0, 1000.0, -1000.0, 1000.0, z - z_delta, z + z_delta]
    #     meshSettings.geometry[boxName] = SearchableBoxGeometry(
    #         purpose='refinement',
    #         refineMax= refLevel,
    #         min= [[box[0], box[2], box[4]],
    #         max= [box[1], box[3], box[5]]
    #     )
    #     return meshSettings


    @staticmethod
    def calc_center_of_mass(mesh):
        center_of_mass_filter = vtk.vtkCenterOfMass()  # type: ignore
        center_of_mass_filter.SetInputData(mesh)
        center_of_mass_filter.Update()
        return center_of_mass_filter.GetCenter()

    @staticmethod
    def calc_location_in_mesh(mesh, internalFlow=False):
        center_of_mass = STLAnalysis.calc_center_of_mass(mesh)

        # TODO: make sure this is the same as above function
        if internalFlow:
            return find_inside_point(mesh, center_of_mass, min_bounds=None, max_bounds=None)
        else:
            return tuple(find_outside_point(mesh, center_of_mass, min_bounds=None, max_bounds=None))

    # @staticmethod
    # def read_stl(stl_file_path: str):
    #     # Check if the file exists
    #     if not os.path.exists(stl_file_path):
    #         raise FileNotFoundError(
    #             f"File not found: {stl_file_path}. Make sure the file exists.")
    #     # Create a reader for the STL file
    #     reader = vtk.vtkSTLReader()  # type: ignore
    #     reader.SetFileName(stl_file_path)
    #     reader.Update()

    #     # Get the output data from the reader
    #     poly_data = reader.GetOutput()
    #     return poly_data



    # @staticmethod
    # def set_stl_solid_name(stl_file='input.stl'):
    #     AmpersandIO.printMessage(f"Setting solid name for {stl_file}")
    #     # if the file does not exist, return -1
    #     if not os.path.exists(stl_file):
    #         print(f"File not found: {stl_file}")
    #         return -1
    #     # if exists, extract file name by removing the directory path
    #     new_lines = []
    #     new_stl_file = stl_file[:-4] + ".stl"
    #     solid_name = os.path.basename(stl_file)[:-4]

    #     # open the file
    #     try:
    #         with open(stl_file, 'r') as f:
    #             lines = f.readlines()
    #     except FileNotFoundError:
    #         print(f"File not found: {stl_file}")
    #         return -1
    #     # find the solid name
    #     for line in lines:
    #         if 'endsolid' in line:
    #             # replace the solid name using above solid_name
    #             line = f"endsolid {solid_name}\n"
    #         # replace the endsolid name
    #         elif 'solid' in line:
    #             line = f"solid {solid_name}\n"
    #         else:
    #             pass
    #         new_lines.append(line)
    #     # print(f"Solid name: {solid_name}")
    #     # write the new lines to the file
    #     try:
    #         with open(new_stl_file, 'w') as f:
    #             f.writelines(new_lines)
    #     except FileNotFoundError:
    #         print(f"File not found: {new_stl_file}")
    #         return -1
    #     return 0




#   # minVolumeSize = backgroundCellSize**3/(8.**refLevel*20.)
#         # print the summary of results
#         AmpersandIO.printMessage(
#             "\n-----------------Mesh Settings-----------------", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(f"Domain size: x({domain_bbox.minx:6.3f}~{domain_bbox.maxx:6.3f}) y({domain_bbox.miny:6.3f}~{
#                                  domain_bbox.maxy:6.3f}) z({domain_bbox.minz:6.3f}~{domain_bbox.maxz:6.3f})", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(f"Nx Ny Nz: {nx},{ny},{
#                                  nz}", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(
#             f"Max cell size: {x_background_cell_size}", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(
#             f"Min cell size: {targetCellSize}", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(
#             f"Refinement Level:{refLevel}", GUIMode=GUI, window=window)

#         AmpersandIO.printMessage(
#             "\n-----------------Turbulence-----------------", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(
#             f"Target yPlus:{target_yPlus}", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(
#             f'Reynolds number:{U*L/nu}', GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(f"Boundary layer thickness:{
#                                  delta}", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(f"First layer thickness:{
#                                  adjustedNearWallThickness}", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(f"Final layer thickness:{
#                                  finalLayerThickness}", GUIMode=GUI, window=window)
#         AmpersandIO.printMessage(
#             f"YPlus:{adjustedYPlus}", GUIMode=GUI, window=window)

#         AmpersandIO.printMessage(
#             f"Number of layers:{nLayers}", GUIMode=GUI, window=window)


    # @staticmethod
    # def set_layer_thickness(meshSettings, thickness=0.01):
    #     meshSettings.addLayersControls.finalLayerThickness = thickness
    #     minThickness = max(0.0001, thickness/100.)
    #     meshSettings.addLayersControls.minThickness = minThickness
    #     return meshSettings

    # @staticmethod
    # def set_min_vol(meshSettings, minVol=1e-15):
    #     meshSettings.meshQualityControls.minVol = 1e-15  # minVol/100.
    #     return meshSettings


if __name__ == "__main__":
    minCurv = STLAnalysis.calc_smallest_curvature("stl/ahmed.stl")
    print(minCurv)
