from nipype.interfaces.base import (BaseInterface, BaseInterfaceInputSpec,
                                    traits, File, TraitedSpec, InputMultiPath,
                                    CommandLineInputSpec, CommandLine, isdefined)
from scipy.spatial import distance
import os
import csv
import vtk
import numpy as np
from os.path import abspath, isfile
import SimpleITK as sitk
from freesurfer_utils import create_label_watershed


# method to parse the labels xml info
def parse_labels_xml(xml_file):
    labels_dict = dict()
    with open(xml_file, 'r') as xml_reader:
        for line in xml_reader:
            if '<Label>' == line.strip():
                name = xml_reader.next().strip().replace('<Name>', '').replace('</Name>', '')
                code = xml_reader.next().strip().replace('<Number>', '').replace('</Number>', '')
                hemi = xml_reader.next().strip().replace('<Hemisphere>', '').replace('</Hemisphere>', '')
                location = xml_reader.next().strip().replace('<Location>', '').replace('</Location>', '')
                labels_dict[code] = dict(name=name, hemisphere=hemi, location=location)
    return labels_dict


def parse_lookup_table(lookup_table_file):
    """Parses a lookup table to determine regions.
    This allows the hemisphere splitting to adapt with updated lookup tables."""
    labels_dict = dict()
    with open(lookup_table_file, 'r') as lookup_table:
        for line in lookup_table:

            # parse line for label code
            row = line.split(' ')
            for i in range(row.count('')):
                row.remove('')
            code = row[0]

            # continue if the code is a number
            if code.isalnum():
                name = row[1]

                # determine hemisphere
                if 'Left' in name or 'lh' in name:
                    hemisphere = 'lh'
                elif 'Right' in name or 'rh' in name:
                    hemisphere = 'rh'
                else:
                    hemisphere = 'N/A'

                # determine location
                # set location to None. Then update it depending on the name.
                location = None

                if 'wm' in name:
                    location = 'wm'
                elif 'ctx' in name or 'gyrus' in name:
                    location = 'gm'
                elif 'CC' in name:
                    location = "cc"
                elif 'Ventricle' in name:
                    location = "ventricle"

                cerebellum_names = ['Cbm', 'Cerebellum', 'Cerebellum', 'Cerebellar', '4th-Ventricle', 'Brain-Stem',
                                    'VentralDC']
                subcortical_names = ['Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Hippocampus', 'Amygdala',
                                     'Accumbens', 'Inf-Lat-Vent']

                for designated_name, list_of_locations in [('cerebellum', cerebellum_names),
                                                           ('subcortical', subcortical_names)]:
                    for location_name in list_of_locations:
                        if location_name in name:
                            location = designated_name

                if not location:
                    location = "UNKNOWN"

                labels_dict[code] = dict(name=name, hemisphere=hemisphere, location=location)

    return labels_dict


def parse_atlas_info(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == '.txt':
        return parse_lookup_table(in_file)
    elif ext == '.xml':
        return parse_labels_xml(in_file)
    else:
        print("Could not parse {0}".format(in_file))


class WMMaskingInputSpec(BaseInterfaceInputSpec):
    atlas_file = File(
        exists=True, mandatory=True,
        desc='Label map used to divide hemispheres')
    posterior_files = traits.Dict(mandatory=True, desc='Posterior probability files')
    brainlabels_file = File(exists=True, mandatory=True,
                            desc='BRAINSABC brain labels')
    atlas_info = File(exists=True, mandatory=True,
                      desc='input label information in xml format')
    dilation = traits.Int(
        default_value=0, desc="""
        Parameter to adjust the dilation of the boundary mask (default=0)
        A negative value will erode the boundary mask.
        """, use_default=True)
    csf_threshold = traits.Float(
        default_value=0.9, desc="""
        Posterior probabilities above this threshold will be considered CSF
        """, use_default=True)
    hncma_file = File(exists=True, desc="HNCMA atlas is used to define ventricles.")


class WMMaskingOutputSpec(TraitedSpec):
    lh_boundary = File(
        desc="""
        Binary mask setting hard boundaries for the outer cerebral cortex
        surfaces for the left hemisphere
        """)
    rh_boundary = File(
        desc="""
        Binary mask setting hard boundaries for the outer cerebral cortex
        surfaces for the right hemisphere""")
    lh_wm = File(
        desc="Binary mask of the white matter for the left hemisphere")
    rh_wm = File(
        desc="Binary mask of the white matter for the right hemisphere")


class WMMasking(BaseInterface):

    """
    Takes in a brainslabel map from BRAINSTools AutoWorkup as well
    as a csf posterior probability map and a label map and outputs
    the wm mask to be used by LOGISMOS-B.
    """

    input_spec = WMMaskingInputSpec
    output_spec = WMMaskingOutputSpec

    def _run_interface(self, runtime):
        atlas_file = self.inputs.atlas_file
        csf_file = self.inputs.posterior_files['CSF']
        brainlabels_file = self.inputs.brainlabels_file
        atlas_info = self.inputs.atlas_info

        # White Matter Masking Script

        # This script creates a white matter mask and a boundary mask for
        # LOGISMOS-B

        import SimpleITK as sitk
        import os

        # Helpful methods
        # method to find largest connected component
        def largest_connected_component(image, minSize=1000):
            return sitk.RelabelComponent(sitk.ConnectedComponent(image), minSize) == 1

        # method to fill holes in the mask
        def fill_mask_holes(image, minSize=0):
            negConnected = largest_connected_component(1 - image, minSize)
            return 1 - negConnected

        def fillLateralVentricle(ventricle, boundary, ventricleDilation=1, dilationFactor=1):
            ventTemp = largest_connected_component(ventricle)
            # Fill the ventricle region
            for i in range(ventricleDilation):
                ventTemp = sitk.DilateObjectMorphology(ventTemp, dilationFactor) * boundary
                ventTemp = largest_connected_component(ventTemp)
            return ventTemp

        # Read input images
        # Read images
        malf_image = sitk.ReadImage(atlas_file, sitk.sitkUInt32)
        csf_posteriors_image = sitk.ReadImage(csf_file, sitk.sitkFloat64)
        brainlabelsImage = sitk.ReadImage(brainlabels_file, sitk.sitkUInt32)

        # Create label dictionary, hemisphere, and cerebellum masks
        RightTemplate = malf_image < 0
        LeftTemplate = malf_image < 0
        CerebellumMask = malf_image < 0
        LeftInsulaGM = malf_image < 0
        RightInsulaGM = malf_image < 0
        # Define regions not to subtract from hemisphere mask
        # preserve the ventricle boundary mask
        preservedRegions = malf_image < 0

        # Brain labels subcortical regions
        subcortical_regions = (brainlabelsImage == 23) + (brainlabelsImage == 24) + \
            (brainlabelsImage == 25) + (brainlabelsImage == 21)

        atlas_dict = parse_atlas_info(atlas_info)
        for code in atlas_dict.iterkeys():
            location = atlas_dict[code]['location']
            hemi = atlas_dict[code]['hemisphere']
            name = atlas_dict[code]['name']
            if location == 'subcortical':
                subcortical_regions = subcortical_regions + (malf_image == code)
            if location == 'cerebellum':
                CerebellumMask = CerebellumMask + (malf_image == code)
            elif hemi == 'lh':
                LeftTemplate = LeftTemplate + (malf_image == code)
                if location == 'ventricle':
                    left_ventricle_label_mask = malf_image == code
                elif 'insula' in name and location == 'gm':
                    LeftInsulaGM = LeftInsulaGM + (malf_image == code)
            elif hemi == 'rh':
                RightTemplate = RightTemplate + (malf_image == code)
                if location == 'ventricle':
                    right_ventricle_label_mask = malf_image == code
                elif 'insula' in name and location == 'gm':
                    RightInsulaGM = RightInsulaGM + (malf_image == code)
            elif location == 'ventricle' or location == 'cc':
                preservedRegions = preservedRegions + (malf_image == code)

        sitk.WriteImage(CerebellumMask, "Cerebellum.nii.gz")

        filled_right_template = fill_mask_holes(RightTemplate, 1000)
        filled_left_template = fill_mask_holes(LeftTemplate, 1000)

        def create_hemisphere_splits(right_template, left_template):
            template = (right_template > 0) + (left_template > 0)*2
            splits = create_label_watershed(template)
            right_split = splits == 1
            left_split = splits == 2
            return right_split, left_split

        right_hemisphere, left_hemisphere = create_hemisphere_splits(filled_right_template, filled_left_template)

        # Create left and right hemisphere WM masks

        # Define the latereral ventricles
        # Extract the right lateral and inferior ventricle from the label map
        if isdefined(self.inputs.hncma_file):
            hncma_left_ventricle_code = 4
            hncma_right_ventricle_code = 43
            hncma_atlas = sitk.ReadImage(self.inputs.hncma_file)
            right_ventricle_label_mask = (right_ventricle_label_mask + (hncma_atlas == hncma_right_ventricle_code)) > 0
            left_ventricle_label_mask = (left_ventricle_label_mask + (hncma_atlas == hncma_left_ventricle_code)) > 0

        right_ventricle_boundary = right_hemisphere
        right_ventricle_final = fillLateralVentricle(right_ventricle_label_mask, right_ventricle_boundary)

        # Extract the left lateral and inferior ventricle from the label map
        left_ventricle_boundary = left_hemisphere
        left_ventricle_final = fillLateralVentricle(left_ventricle_label_mask, left_ventricle_boundary)

        # Add subcortical regions and lateral ventricles to the WM mask
        white_matter = brainlabelsImage == 1
        complete_white_matter = (white_matter + subcortical_regions + right_ventricle_final + left_ventricle_final) > 0
        white_matter_final = largest_connected_component(fill_mask_holes(complete_white_matter, 1000))

        # Regions not included in white matter
        left_white_matter_template = left_hemisphere * (CerebellumMask == 0)
        right_white_matter_template = right_hemisphere * (CerebellumMask == 0)

        # Split the hemispheres
        # left hemisphere white matter mask
        LeftWhiteMatter = largest_connected_component(white_matter_final * left_white_matter_template)
        # right hemisphere white matter mask
        RightWhiteMatter = largest_connected_component(white_matter_final * right_white_matter_template)

        # Create right and left boundary masks
        # dilate preserved regions
        preservedDilation = 1
        preservedRegions = sitk.DilateObjectMorphology(preservedRegions, preservedDilation)
        # add the whtie matter masks to the preserved regions
        preservedRegions = preservedRegions + RightWhiteMatter + LeftWhiteMatter > 0

        # Define CSF regions
        # Threshold for CSF
        Thresh = self.inputs.csf_threshold
        CSFMask = sitk.BinaryThreshold(csf_posteriors_image, Thresh)
        CSFMask = CSFMask * (1 - preservedRegions)

        # Remove mask around cerebellum and brain stem
        cerebellumDilation = 1
        cerebellumMaskDilated = sitk.DilateObjectMorphology(CerebellumMask, cerebellumDilation)
        leftBoundaryMask = left_hemisphere * (1 - cerebellumMaskDilated)
        leftBoundaryMask = largest_connected_component(leftBoundaryMask)
        rightBoundaryMask = right_hemisphere * (1 - cerebellumMaskDilated)
        rightBoundaryMask = largest_connected_component(rightBoundaryMask)

        # Dilate or erode boundary masks
        boundaryDilation = self.inputs.dilation
        if boundaryDilation < 0:
            # Erode brainlabels
            # LOGISMOS-B dilates the input labels, so this can be offset
            # by eroding the brainlabels mask
            brainlabelserosion = 1
            brainlabelsmask = sitk.ErodeObjectMorphology(brainlabelsImage > 0, brainlabelserosion)
            brainlabelsImage = sitk.Cast(brainlabelsmask, sitk.sitkUInt32) * brainlabelsImage
        elif boundaryDilation > 0:
            leftBoundaryMask = sitk.DilateObjectMorphology(leftBoundaryMask, boundaryDilation)
            rightBoundaryMask = sitk.DilateObjectMorphology(rightBoundaryMask, boundaryDilation)

        # Remove CSF from hemisphere masks
        leftBoundaryMask = leftBoundaryMask * (1 - CSFMask)
        leftBoundaryMask = largest_connected_component(leftBoundaryMask)
        rightBoundaryMask = rightBoundaryMask * (1 - CSFMask)
        rightBoundaryMask = largest_connected_component(rightBoundaryMask)

        # Convert mask to brains label map
        leftBrainLabels = sitk.Cast(leftBoundaryMask, sitk.sitkUInt32) * brainlabelsImage
        rightBrainLabels = sitk.Cast(rightBoundaryMask, sitk.sitkUInt32) * brainlabelsImage
        sitk.WriteImage(leftBrainLabels, self._list_outputs()['lh_boundary'])
        sitk.WriteImage(rightBrainLabels, self._list_outputs()['rh_boundary'])
        sitk.WriteImage(LeftWhiteMatter, self._list_outputs()['lh_wm'])
        sitk.WriteImage(RightWhiteMatter, self._list_outputs()['rh_wm'])

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['lh_boundary'] = abspath('LeftBrainLabels.nii.gz')
        outputs['rh_boundary'] = abspath('RightBrainLabels.nii.gz')
        outputs['lh_wm'] = abspath('LeftWhiteMatter.nii.gz')
        outputs['rh_wm'] = abspath('RightWhiteMatter.nii.gz')

        return outputs


# Define a useful function for reading PolyData vtk files
# Read a PolyData file and output a vtk PolyData object
def readPolyData(filename):
    "Reading polydata: " + filename
    # Check which PolyData reader should be used
    if ".vtk" in filename:
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()
    else:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()


# Find the label of a given vtk point from a label map
def vtkPoint_to_label(point, labelmap):
    surfx, surfy, surfz = point
    point = (-surfx, -surfy, surfz) # must flip y axis to convert from VTK to ITK
    index = labelmap.TransformPhysicalPointToIndex(point)
    x = int(index[0])
    y = int(index[1])
    z = int(index[2])
    label = int(labelmap.GetPixel(x, y, z))
    return label


# As advised on the ITK mailing list, label dilation can be implemented via
# distance transforms and watershed transforms. This algorithm is illustrated
# in SimpleITK python code below (courtesy of Bradely Lowekamp)
def multilabel_dilation(img, radius=1, kernel=sitk.BinaryDilateImageFilter.Ball):
    distImg = sitk.SignedMaurerDistanceMap(img != 0,
                                           insideIsPositive=False,
                                           squaredDistance=False,
                                           useImageSpacing=False)
    dilatImg = sitk.BinaryDilate(img != 0, radius, kernel)
    wsImg = sitk.MorphologicalWatershedFromMarkers(distImg, img)
    return sitk.Cast(dilatImg, sitk.sitkUInt64) * wsImg


class CreateGMLabelMapInputSpec(BaseInterfaceInputSpec):
    atlas_file = File(
        exists=True, mandatory=True,
        desc='Label map used to define gray matter regions')
    atlas_info = File(exists=True, mandatory=True,
                      desc='input label information in xml format')


class CreateGMLabelMapOutputSpec(TraitedSpec):
    out_file = File(desc="gray matter label map")


class CreateGMLabelMap(BaseInterface):

    """
    Selects the gray matter labels and then dilates them
    """

    input_spec = CreateGMLabelMapInputSpec
    output_spec = CreateGMLabelMapOutputSpec

    def _run_interface(self, runtime):
        atlas_file = self.inputs.atlas_file
        atlas_dict = parse_labels_xml(self.inputs.atlas_info)
        atlas_img = sitk.Cast(sitk.ReadImage(atlas_file), sitk.sitkUInt64)
        gm_mask = atlas_img < 0
        for code in atlas_dict.iterkeys():
            location = atlas_dict[code]['location']
            if location == 'gm':
                gm_mask = gm_mask + (atlas_img == code)
        gm_mask = gm_mask > 0
        gm_labels = sitk.Cast(gm_mask, sitk.sitkUInt64) * atlas_img
        out_img = multilabel_dilation(sitk.Cast(gm_labels, sitk.sitkUInt64))
        sitk.WriteImage(out_img, self._list_outputs()['out_file'])

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath('gm_labels.nii.gz')
        return outputs


class ComputeDistanceInputSpec(BaseInterfaceInputSpec):
    wm_file = File(exists=True, desc='vtk polydata mesh surface', mandatory=True)
    gm_file = File(exists=True, desc='vtk polydata mesh surface', mandatory=True)
    labels_file = File(exists=True, desc='label image file', mandatory=True)
    hemisphere = traits.Enum('lh', 'rh', desc='hemisphere being processed', mandatory=True)
    atlas_info = File(exists=True, mandatory=False,
                      desc='input label information in xml format')


class ComputeDistanceOutputSpec(TraitedSpec):
    out_file = File(desc="vtk polydata mesh surface with distance scalars")


class ComputeDistance(BaseInterface):

    """
    Nipype wrappers for a comparing 2 surfaces
    Compute the surface to surface distance between 2 using similar to FreeSurfer
    """

    input_spec = ComputeDistanceInputSpec
    output_spec = ComputeDistanceOutputSpec

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(
            '{0}_ctx_results.csv'.format(self.inputs.hemisphere))
        return outputs

    def _run_interface(self, runtime):
        labelmap = sitk.ReadImage(self.inputs.labels_file)
        uncleanwm = readPolyData(self.inputs.wm_file)
        uncleangm = readPolyData(self.inputs.gm_file)
        if isdefined(self.inputs.atlas_info):
            atlas_dict = parse_labels_xml(self.inputs.atlas_info)
        # Clean the data
        cleanwm = vtk.vtkCleanPolyData()
        cleanwm.SetInputData(uncleanwm)
        cleangm = vtk.vtkCleanPolyData()
        cleangm.SetInputData(uncleangm)
        cleanwm.Update()
        cleangm.Update()
        wmsurf = cleanwm.GetOutput()
        gmsurf = cleangm.GetOutput()
        # setup KdTrees for each surface
        # this will help in finding the closest points
        kdTreewm = vtk.vtkKdTreePointLocator()
        kdTreewm.SetDataSet(wmsurf)
        kdTreewm.BuildLocator()
        kdTreegm = vtk.vtkKdTreePointLocator()
        kdTreegm.SetDataSet(gmsurf)
        kdTreegm.BuildLocator()
        measurements = dict()
        wmPD = wmsurf.GetPointData()
        wmPoints = wmsurf.GetPoints()
        wmCount = wmPD.GetNumberOfTuples()
        for i in range(0, wmCount):
            wmP = wmPoints.GetPoint(i)
            # Find the closest point to the gray matter surface point
            gmIndex = kdTreegm.FindClosestPoint(wmP)
            gmP = kdTreegm.GetDataSet().GetPoint(gmIndex)
            # Get the gray matter label from the label map
            gmlabel = vtkPoint_to_label(gmP, labelmap)
            if gmlabel != 0:
                label = str(gmlabel)
            else:
                # if the gray matter label is not defined try the wm label
                wmlabel = vtkPoint_to_label(wmP, labelmap)
                if wmlabel != 0:
                    label = str(wmlabel)
                else:
                    # label is not known
                    label = 'UNKNOWN'
            # compute the distance
            # distance from wm point to gm point
            dst1 = distance.euclidean(wmP, gmP)
            wmIndex = kdTreewm.FindClosestPoint(gmP)
            wmP2 = kdTreegm.GetDataSet().GetPoint(wmIndex)
            # distnace from gm to closest wm point
            dst2 = distance.euclidean(gmP, wmP2)
            # average the two distances
            thickness = (dst1 + dst2) / 2
            if not measurements.has_key(label):
                # first point in a labeled region
                measurements[label] = [thickness]
            else:
                measurements[label].append(thickness)

        mu = ["mean"]
        median = ["median"]
        std = ["std"]
        count = ["points"]
        minimum = ["min"]
        maximum = ["max"]
        labels = ["label"]
        for key in measurements.iterkeys():
            labels.append(key)
            data = np.array(measurements[key])
            mu.append(np.mean(data))
            median.append(np.median(data))
            std.append(np.std(data))
            count.append(len(data))
            minimum.append(np.min(data))
            maximum.append(np.max(data))

        out_csv = self._list_outputs()['out_file']
        with open(out_csv, 'w') as CSV_file:
            writer = csv.writer(CSV_file)
            writer.writerows([labels, mu, std, count, minimum, maximum])

        return runtime


class LOGISMOSBInputSpec(CommandLineInputSpec):
    t1_file = File(exists=True, desc='T1 scan output by BAW', argstr='--inputT1 %s', mandatory=True)
    t2_file = File(exists=True, genfile=True, desc='T2 scan output by BAW', argstr='--inputT2 %s', mandatory=False)
    mesh_file = File(exists=True, desc='final mesh of the white matter surface (must have a genus equal to 0)',
                     argstr='-m %s', mandatory=True)
    wm_file = File(exists=True, desc='final binary image of the white matter surface (must have a genus equal to 0)',
                   argstr='-b %s', mandatory=True)
    atlas_file = File(exists=True, desc='hcnma atlas to define brain regions. If different atlas is used, thick ' +
                                        'regions must be defined',
                      argstr='-z %s', mandatory=False)
    brainlabels_file = File(exists=True, desc='skullstripped brainlabels file', argstr='--inputABCLabels %s',
                            mandatory=True)
    smoothnessConstraint = traits.Int(desc='smoothness constraint',
                                      argstr='--smoothnessConstraint %d', mandatory=True)
    nColumns = traits.Int(desc="number of vertices", argstr="--nColumns %d", Mandatory=False)
    columnChoice = traits.String(desc="some parameter", argstr="--columnChoice %s", Mandatory=False)
    columnHeight = traits.Int(desc="column height", argstr="--columnHeight %d", Mandatory=False)
    nodeSpacing = traits.Float(desc="node spacing", argstr="--nodeSpacing %.2f", Mandatory=False)
    w = traits.Float(desc="w", argstr="-w %.2f", Mandatory=False)
    a = traits.Float(desc="a", argstr="-a %.2f", Mandatory=False)
    nPropagate = traits.Int(desc="number of propagations", argstr="--nPropagate %d", Mandatory=False)
    basename = traits.String(desc="basename for output files", argstr="--outputBase %s", Mandatory=True)
    thick_regions = traits.List(traits.Int(), argstr="-r %s", mandatory=False, sep=',',
                                desc="List of regions in the atlas file to that will be thicker")
    useHNCMALabels = traits.Bool(argstr="--useHNCMALabels", desc="Uses HCNMA label map to define thick regions")
    wm_proba_file = File(exist=True, argstr='--wmProbaMap %s', desc="White matter pobability map.")
    gm_proba_file = File(exist=True, argstr='--gmProbaMap %s', desc="Gray matter pobability map.")


class LOGISMOSBOutputSpec(TraitedSpec):
    gmsurface_file = File(desc="path/name of GM surface file")
    wmsurface_file = File(desc="path/name of WM surface file")
    profile_file = File(desc="output profile file")


class LOGISMOSB(CommandLine):
    _cmd = 'LOGISMOS-B'
    input_spec = LOGISMOSBInputSpec
    output_spec = LOGISMOSBOutputSpec

    def _gen_filename(self, name):
        if name == "t2_file":
            return self.inputs.t1_file
        return None

    def _format_arg(self, name, spec, value):
        if name == "t2_file" and not os.path.isfile(value):
            print "Using T1 as T2 file"
            value = self.inputs.t1_file
        return super(LOGISMOSB, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['gmsurface_file'] = os.path.abspath(self.inputs.basename + "_GMresult.vtp")
        outputs['wmsurface_file'] = os.path.abspath(self.inputs.basename + "_WMresult.vtp")
        outputs['profile_file'] = os.path.abspath(self.inputs.basename + "_profile.vtk")
        return outputs


class BSGInputSpec(CommandLineInputSpec):
    in_file = File(exists=True, mandatory=True,
                   desc="binary ITK image mask file",
                   argstr="--inputImageFile %s")

    out_file = File(desc="output vtk polydata surface mesh file",
                    argstr="--outputSurface %s", mandatory=True)

    smoothSurface = traits.Bool(desc="smooth the surface (default: off)", argstr="--smoothSurface")

    numIterations = traits.Int(desc="number of iterations to smooth the surface (default: 5)",
                               argstr="--numIterations %d")

    surfaceValue = traits.Float(desc="The iso-surface value for the resulting surface (default: 0.5)",
                                argstr="--surfaceValue %.2f")

    decimateSurface = traits.Bool(desc="decimate the surface (default: off)",
                                  argstr="--decimateSurface")

    numberOfElements = traits.Int(desc="Number of faces desired after decimation (default: 70000)",
                                  argstr="--numberOfElements %d")

    relaxationFactor = traits.Float(dec="The Relaxation Factor Used in Smoothing (default: 0.1)",
                                    argstr="--relaxationFactor %.2f")


class BSGOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="output vtk mesh surface file")


class BRAINSSurfaceGeneration(CommandLine):
    _cmd = 'BRAINSSurfaceGeneration'
    input_spec = BSGInputSpec
    output_spec = BSGOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs


class Genus0InputSpec(CommandLineInputSpec):
    in_file = File(exists=True,
                   argstr="--inputVolume %s",
                   desc="Input the image volume to be topologically corrected",
                   mandatory=True)
    out_mask = File(argstr="--outputVolume %s",
                    desc="Topologically corrected image volume output (ignored if computeSurface is set)")
    out_mesh = File(argstr="--vtkOutput %s",
                    desc="File to write a VTK triangluated mesh to (ignored if computeSurface is not set)")
    stl = File(argstr="--stl %s",
               desc="File to write an STL triangluated mesh to.")
    cutLoops = traits.Bool(argstr="--cutLoops",
                           desc="Cut loops instead of patching holes (default: OFF)")
    connectedComponent = traits.Bool(argstr="--connectedComponent",
                                     desc="Extract largest connected component before processing (default: OFF)")
    connectivity = traits.Int(argstr="--connectivity %d",
                              desc="Controls the discrete connectivity model (18|6 default: 18). 18-connectivity only allows for the output to be a vtk surface (computeSurface must be OFF).")
    computeSurface = traits.Bool(argstr="--computeSurface",
                                 desc="Compute VTK surface instead of corrected image volume (default: OFF)")
    extractFinalConnectedComponent = traits.Bool(argstr="--extractFinalConnectedComponent",
                                                 desc="Extracts the largest connected component after the processing. (default: 0)")
    biggestComponent = traits.Bool(argstr="--biggestComponent",
                                   desc="Extract largest component of the triangulated result. (The volume result needs to be followed by an extraction of the largest connected component if desired; use 'extractFinalConnectedComponent'.) (default: 0)")
    returnParameterFile = File(argstr="--returnparameterfile %s",
                               desc="Filename in which to write simple return parameters (int, float, int-vector, etc.) as opposed to bulk return parameters (image, geometry, transform, measurement, table).")
    processInformationAddress = File(argstr="--processinformationaddress %s",
                                     desc="Address of a structure to store process information (progress, abort, etc.). (default: 0)")
    xml = traits.Bool(argstr="--xml",
                      desc="Produce xml description of command line arguments (default: 0)")
    echo = traits.Bool(argstr="--echo",
                       desc="Echo the command line arguments (default: 0)")
    commandHelp = traits.Bool(argstr="--help",
                              desc="Displays the parameters to run this command")

class Genus0OutputSpec(TraitedSpec):
    out_file = File(desc="white matter binary mask image (.nii.gz) or a vtk mesh")


class GenusZeroImageFilter(CommandLine):
    _cmd = 'GenusZeroImageFilter'
    input_spec = Genus0InputSpec
    output_spec = Genus0OutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if self.inputs.computeSurface:
            outputs['out_file'] = os.path.abspath(self.inputs.out_mesh)
        else:
            outputs['out_file'] = os.path.abspath(self.inputs.out_mask)
        return outputs
