import vtk
import SimpleITK as sitk
import numpy as np
from scipy.spatial import distance
from nipype.interfaces.freesurfer import MRIsConvert
import os
import sys


def read_poly_data(filename):
    # Check which PolyData reader should be used
    if ".vtk" in filename:
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()
    elif ".vtp" in filename:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()
    else:
        print("ERROR: Failed to read in polydata")
        return sys.exit(os.EX_IOERR)


def ras_to_lps(point):
    surf_x, surf_y, surf_z = point
    point = (-surf_x, -surf_y, surf_z)  # must flip y axis to convert from VTK to ITK
    return point


# Find the label of a given vtk point from a label map
def vtk_point_to_label(point, labelmap):
    point = ras_to_lps(point)
    index = labelmap.TransformPhysicalPointToIndex(point)
    x = int(index[0])
    y = int(index[1])
    z = int(index[2])
    return labelmap.GetPixel(x, y, z)


def build_kd_tree(mesh):
    kd_tree = vtk.vtkKdTreePointLocator()
    kd_tree.SetDataSet(mesh)
    kd_tree.BuildLocator()
    return kd_tree


def convert_fs_surface(in_surf, out_surf, to_scanner=True):
    if os.path.isfile(os.path.abspath(out_surf)):
        return os.path.abspath(out_surf)
    mris_convert = MRIsConvert()
    mris_convert.inputs.in_file = in_surf
    mris_convert.inputs.out_file = os.path.abspath(out_surf)
    mris_convert.inputs.to_scanner = to_scanner
    result = mris_convert.run()
    return result.outputs.converted


def get_vtk_file_name(fs_file_name):
    fs_dir, fs_basename = os.path.split(fs_file_name)
    return os.path.join(fs_dir, fs_basename.replace(".", "_") + ".vtk")


def fs_to_vtk(fs_surface):
    output_file = get_vtk_file_name(fs_surface)
    return convert_fs_surface(fs_surface, output_file)


def get_surf(surf_dir, hemisphere, surf):
    return os.path.join(surf_dir, "{0}.{1}".format(hemisphere, surf))


def get_white(surf_dir, hemisphere):
    return get_surf(surf_dir, hemisphere, "white")


def get_pial(surf_dir, hemisphere):
    return get_surf(surf_dir, hemisphere, "pial")


def get_white_and_pial_fs_files(surf_dir, hemisphere):
    fs_white = get_white(surf_dir, hemisphere)
    fs_pial = get_pial(surf_dir, hemisphere)
    return fs_white, fs_pial


def get_white_and_pial_vtk_files(surf_dir, hemisphere):
    fs_white, fs_pial = get_white_and_pial_fs_files(surf_dir, hemisphere)
    return fs_to_vtk(fs_white), fs_to_vtk(fs_pial)


def get_white_and_pial(surf_dir, hemisphere):
    vtk_white, vtk_pial = get_white_and_pial_vtk_files(surf_dir, hemisphere)
    white = read_poly_data(vtk_white)
    pial = read_poly_data(vtk_pial)
    return white, pial


def compute_thickness(wmP, kdTreegm, kdTreewm):
    # Find the closest point to the gray matter surface point
    gmIndex = kdTreegm.FindClosestPoint(wmP)
    gmP = kdTreegm.GetDataSet().GetPoint(gmIndex)
    # compute the distance
    # distance from wm point to gm point
    dst1 = distance.euclidean(wmP, gmP)
    wmIndex = kdTreewm.FindClosestPoint(gmP)
    wmP2 = kdTreegm.GetDataSet().GetPoint(wmIndex)
    # distnace from gm to closest wm point
    dst2 = distance.euclidean(gmP, wmP2)
    # average the two distances
    thickness = (dst1 + dst2)/float(2)
    return thickness


def create_thickness_array():
    thicknesses = vtk.vtkFloatArray()
    thicknesses.SetName("thickness")
    return thicknesses


def calculate_distance(white, pial):
    # setup KdTrees for each surface
    # this will help in finding the closest points
    kd_tree_white = build_kd_tree(white)
    kd_tree_pial = build_kd_tree(pial)

    white_points = white.GetPoints()
    white_count = white.GetNumberOfPoints()
    white_point_data = white.GetPointData()
    thicknesses = create_thickness_array()

    for i in range(0, white_count):
        white_matter_point = white_points.GetPoint(i)

        # compute the thickness
        thickness = compute_thickness(white_matter_point, kd_tree_pial, kd_tree_white)
        thicknesses.InsertNextValue(thickness)

    white_point_data.AddArray(thicknesses)
    return white


def get_surf_dir(subjects_dir, subject_id):
    return os.path.join(subjects_dir, subject_id, "surf")


def write_vtk_file(polydata, file_name):
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(file_name)
    writer.SetInputData(polydata)
    writer.Update()
    return os.path.abspath(writer.GetFileName())


def get_thickness_file(subjects_dir, subject_id, hemisphere):
    surf_dir = get_surf_dir(subjects_dir, subject_id)
    white, pial = get_white_and_pial(surf_dir, hemisphere)
    thickness = calculate_distance(white, pial)
    return write_vtk_file(thickness, os.path.join(surf_dir, "{0}_thickness.vtk".format(hemisphere)))


def get_thickness_files_for_both_hemispheres(subjects_dir, subject_id):
    lh_thickness = get_thickness_file(subjects_dir, subject_id, 'lh')
    rh_thickness = get_thickness_file(subjects_dir, subject_id, 'rh')
    return lh_thickness, rh_thickness


def masked_thickness_values(thickness_file, mask_image_file, array_index=None):
    thickness = read_poly_data(thickness_file)
    mask = sitk.ReadImage(mask_image_file)

    inside_mask_values = list()
    outside_mask_values = list()

    thickness_point_data = thickness.GetPointData()
    if not array_index:
        # set the array index to the last array added to the poly data
        array_index = thickness_point_data.GetNumberOfArrays() - 1
    thickness_values = thickness.GetPointData().GetArray(array_index)

    for point_index in range(thickness.GetNumberOfPoints()):
        point = thickness.GetPoint(point_index)
        mask_value = vtk_point_to_label(point, mask)
        thickness_value = thickness_values.GetValue(point_index)
        if mask_value == 1:
            inside_mask_values.append(thickness_value)
        else:
            outside_mask_values.append(thickness_value)
    return inside_mask_values, outside_mask_values


def calculate_stats(values):
    if values:
        values_array = np.array(values)
        return dict(mean=values_array.mean(), std=values_array.std(), min=values_array.min(), max=values_array.max())
    else:
        return dict(mean=None, std=None, min=None, max=None)


def masked_thickness_stats(thickness_file, mask_image_file):
    inside_mask_values, outside_mask_values = masked_thickness_values(thickness_file, mask_image_file)
    stats = dict()
    stats['inside'] = calculate_stats(inside_mask_values)
    stats['outside'] = calculate_stats(outside_mask_values)
    return stats


def get_thickness_stats_for_both_hemispheres(subjects_dir, subject_id, mask_file):
    stats = dict()
    lh_thickness, rh_thickness = get_thickness_files_for_both_hemispheres(subjects_dir, subject_id)
    stats['lh'] = masked_thickness_stats(lh_thickness, mask_file)
    stats['rh'] = masked_thickness_stats(rh_thickness, mask_file)
    return stats


def main():
    os.environ['PATH'] += ":/Shared/sinapse/sharedopt/apps/freesurfer/Darwin/x86_64/6.0-beta/20150915/bin/"
    mask_file = "/Shared/sinapse/CACHE/20160712_AtrophySimulation_Results/2559/58661/simulation_1/atrophy_regions.nii.gz"
    subj_dir = "/Shared/sinapse/CACHE/20160713_AtrophySimulation_BAW_base_Results/PHD_024/2559_58661/79/"
    print(get_thickness_stats_for_both_hemispheres(subj_dir, "FreeSurfer", mask_file))
    print("done")

if __name__ == "__main__":
    main()
