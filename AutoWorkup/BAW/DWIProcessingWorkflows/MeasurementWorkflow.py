"""
MeasurementWorkflow.py
============================
Description:
    The purpose of this is to...
    
Author:

Usage:

"""

import os
from functools import reduce

import SimpleITK as sitk
import nipype
import nipype.interfaces.io as nio  # Data i/oS
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces import ants
from nipype.interfaces.base import (
    CommandLine,
    CommandLineInputSpec,
    TraitedSpec,
    File,
    Directory,
)
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.semtools import *
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
from past.utils import old_div


def create_measurement_workflow(WFname, LABELS_CONFIG_FILE):
    """
    This Function takes in...

    :param WFname:
    :param LABELS_CONFIG_FILE:
    :return:
    """
    # \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    ###### UTILITY FUNCTIONS #######
    # This function returns a label map that only covers the FOV of the input DWI scan
    def create_dwi_label_map(T2LabelMapVolume, DWIBrainMask):
        """
        This Function takes in...

        :param T2LabelMapVolume:
        :param DWIBrainMask:
        :return:
        """
        import os
        import SimpleITK as sitk

        T2LabelMapVolume = sitk.ReadImage(
            T2LabelMapVolume, sitk.sitkUInt16.encode("ascii", "replace")
        )  # FreeSurfer labelmap needs uint-16
        DWIBrainMask = sitk.ReadImage(DWIBrainMask)
        # 1- Dilate input DWI mask
        dilateFilter = sitk.BinaryDilateImageFilter()
        dilateFilter.SetKernelRadius(1)
        dilated_mask = dilateFilter.Execute(DWIBrainMask)
        # 2- Resample dilated mask to the space of T2LabelMap (1x1x1)
        # Use Linear interpolation + thresholding
        resFilt = sitk.ResampleImageFilter()
        resFilt.SetReferenceImage(T2LabelMapVolume)
        resFilt.SetOutputPixelType(sitk.sitkFloat32)
        resFilt.SetInterpolator(sitk.sitkLinear)
        resampled_dilated_mask = resFilt.Execute(dilated_mask)
        # Thresholding by 0
        threshFilt = sitk.BinaryThresholdImageFilter()
        thresh_resampled_dilated_mask = threshFilt.Execute(
            resampled_dilated_mask, 0.0001, 1.0, 1, 0
        )
        # 3- Cast the thresholded image to uInt-16
        castFilt = sitk.CastImageFilter()
        castFilt.SetOutputPixelType(sitk.sitkUInt16)
        casted_thresh_resampled_dilated_mask = castFilt.Execute(
            thresh_resampled_dilated_mask
        )
        # 4- Multiply this binary mask to the T2 labelmap volume
        mulFilt = sitk.MultiplyImageFilter()
        DWILabelMapVolume = mulFilt.Execute(
            casted_thresh_resampled_dilated_mask, T2LabelMapVolume
        )
        # write the output label map
        outputVolume = os.path.realpath("DWILabelMapVolume.nrrd")
        sitk.WriteImage(DWILabelMapVolume, outputVolume)
        return outputVolume

    def make_resampled_in_file_list(
        FAImage,
        MDImage,
        RDImage,
        FrobeniusNormImage,
        Lambda1Image,
        Lambda2Image,
        Lambda3Image,
    ):
        """
        This Function takes in...

        :param FAImage:
        :param MDImage:
        :param RDImage:
        :param FrobeniusNormImage:
        :param Lambda1Image:
        :param Lambda2Image:
        :param Lambda3Image:
        :return:
        """
        RISsList = [
            FAImage,
            MDImage,
            RDImage,
            FrobeniusNormImage,
            Lambda1Image,
            Lambda2Image,
            Lambda3Image,
        ]
        return RISsList

    # This functions computes statistics of each input RIS volume over all input labels
    # and writes the results as a CSV file
    def compute_statistics(
        inputVolume, T2LabelMapVolume, DWILabelMapVolume, labelCodesFile
    ):
        """
        This Function takes in...

        :param inputVolume:
        :param T2LabelMapVolume:
        :param DWILabelMapVolume:
        :param labelCodesFile:
        :return:
        """
        import os
        import SimpleITK as sitk
        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        #### Util Funcs ####
        def create_labels_dictionary(labelCodesFile):
            """
            This Function takes in...

            :param labelCodesFile:
            :return:
            """
            import csv

            labelsDictionary = OrderedDict()
            with open(labelCodesFile) as lf:
                reader = csv.reader(lf, delimiter=",")
                for line in reader:
                    if line[0][0] == "#":
                        continue
                    else:
                        labelsDictionary[line[0]] = line[1]
            return labelsDictionary

        def compute_voxel_volume(inputVolume):
            """
            This Function takes in...

            :param inputVolume:
            :return:
            """
            import operator

            return reduce(operator.mul, inputVolume.GetSpacing())

        def return_statistics_list(
            labelID, voxelVolume, resampledRISVolume, DWILabelMap, T2LabelMap
        ):
            """
            This Function takes in...

            :param labelID:
            :param voxelVolume:
            :param resampledRISVolume:
            :param DWILabelMap:
            :param T2LabelMap:
            :return:
            """
            from past.utils import old_div
            from collections import (
                OrderedDict,
            )  # Need OrderedDict internally to ensure consistent ordering

            statFilter = sitk.LabelStatisticsImageFilter()
            # RIS stats over input label ID
            statFilter.Execute(resampledRISVolume, DWILabelMap)
            mean = statFilter.GetMean(labelID)
            std = statFilter.GetSigma(labelID)
            maximum = statFilter.GetMaximum(labelID)
            minimum = statFilter.GetMinimum(labelID)
            median = statFilter.GetMedian(labelID)
            effectiveVolume = statFilter.GetCount(labelID) * voxelVolume
            # compute total volume of input label ID in the non-cropped labelmap (T2LabelMap)
            statFilter.Execute(resampledRISVolume, T2LabelMap)
            totalVolume = statFilter.GetCount(labelID) * voxelVolume
            # if effectiveVolume is 0, the label is missed in dwi scan, or it doesn't exists in current labelmaps
            # in both cases we need zero confidence coefficient for that.
            if effectiveVolume == 0:
                confidence_coeficient = 0
                maximum = 0
                minimum = 0
            else:
                if totalVolume == 0:
                    raise ValueError(
                        "Label {0} is not found in T2 labels map, but exists in DWI labels map!".format(
                            labelID
                        )
                    )
                confidence_coeficient = old_div(effectiveVolume, totalVolume)
            # Now create statistics list
            statsList = [
                format(mean, ".4f"),
                format(std, ".4f"),
                format(maximum, ".4f"),
                format(minimum, ".4f"),
                format(median, ".4f"),
                effectiveVolume,
                totalVolume,
                format(confidence_coeficient, ".3f"),
            ]
            return statsList, totalVolume

        def write_label_statistics(filename, statisticsDictionary):
            """
            This Function takes in...

            :param filename:
            :param statisticsDictionary:
            :return:
            """
            import csv

            with open(filename, "wb") as lf:
                headerdata = [
                    [
                        "#Label",
                        "mean",
                        "std",
                        "max",
                        "min",
                        "median",
                        "effective_volume",
                        "total_volume",
                        "confidence_coeficient",
                    ]
                ]
                wr = csv.writer(lf, delimiter=",")
                wr.writerows(headerdata)
                for key, value in sorted(statisticsDictionary.items()):
                    wr.writerows([[key] + value])

        #### #### #### ####
        resampledRISVolume = sitk.ReadImage(inputVolume)
        T2LabelMap = sitk.ReadImage(T2LabelMapVolume)
        DWILabelMap = sitk.ReadImage(DWILabelMapVolume)
        labelsDictionary = create_labels_dictionary(labelCodesFile)
        statisticsDictionary = OrderedDict()
        voxelVolume = compute_voxel_volume(resampledRISVolume)
        for key in labelsDictionary:
            labelID = int(key)
            [statisticsList, total_volume] = return_statistics_list(
                labelID, voxelVolume, resampledRISVolume, DWILabelMap, T2LabelMap
            )
            if total_volume != 0:
                statisticsDictionary[labelsDictionary[key]] = statisticsList
        # Create output file name
        inputBaseName = os.path.basename(inputVolume)
        inputName = os.path.splitext(inputBaseName)[0]
        RISName = inputName.split("_", 1)[0]
        CSVStatisticsFile = os.path.realpath(RISName + "_statistics.csv")
        write_label_statistics(CSVStatisticsFile, statisticsDictionary)
        return CSVStatisticsFile

    # This function helps to pick desirable output from the output list
    def pick_from_file(inputlist, item):
        """
        This Function takes in...

        :param inputList:
        :param item:
        :return:
        """
        return inputlist[item]

    def resample_ris_volumes(referenceVolume, inputVolume):
        """
        This Function takes in...

        :param referenceVolume:
        :param inputVolume:
        :return: outputVolume
        """
        import os
        import SimpleITK as sitk

        refVolume = sitk.ReadImage(referenceVolume)
        RISVolume = sitk.ReadImage(inputVolume)
        # 0- because of numberical precision error, we have small negative values
        # in rotationary invariant scalars e.g. -1e-13,
        # so set voxel value x to 0 if -1e-5<x<1e-5
        binaryThreshFilt = sitk.BinaryThresholdImageFilter()
        maskNearZero = binaryThreshFilt.Execute(RISVolume, -1e-5, 1e-5, 0, 1)
        RISVolume = RISVolume * sitk.Cast(maskNearZero, sitk.sitkFloat64)
        # Linear interpolation
        resFilt = sitk.ResampleImageFilter()
        resFilt.SetReferenceImage(refVolume)
        resFilt.SetOutputPixelType(RISVolume.GetPixelIDValue())
        resFilt.SetInterpolator(sitk.sitkLinear)
        RIS_resampled = resFilt.Execute(RISVolume)
        """
        # sqrt voxel-wise + cubic BSpline + square voxel-wise
        # 1- voxel-wise square root of input volume
        sqrtFilt = sitk.SqrtImageFilter()
        RIS_sqrt = sqrtFilt.Execute(RISVolume)
        # 2-resample squared image using cubic BSpline
        resFilt = sitk.ResampleImageFilter()
        resFilt.SetReferenceImage(refVolume)
        resFilt.SetInterpolator(sitk.sitkBSpline)
        RIS_sqrt_res = resFilt.Execute(RIS_sqrt)
        # 3- square the resampled RIS volume voxel-wise
        squarFilt = sitk.SquareImageFilter()
        RIS_resampled = squarFilt.Execute(RIS_sqrt_res)
        """
        # Create output file name
        inputBaseName = os.path.basename(inputVolume)
        RISName = os.path.splitext(inputBaseName)[0]
        outputVolume = os.path.realpath(RISName + "_res.nrrd")
        sitk.WriteImage(RIS_resampled, outputVolume)
        assert os.path.isfile(outputVolume), (
            "Resampled RIS file is not found: %s" % outputVolume
        )
        return outputVolume

    #################################
    # \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    MeasurementWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=[
                "T2LabelMapVolume",
                "DWIBrainMask",
                "LabelsConfigFile",
                "FAImage",
                "MDImage",
                "RDImage",
                "FrobeniusNormImage",
                "Lambda1Image",
                "Lambda2Image",
                "Lambda3Image",
            ]
        ),
        name="inputsSpec",
    )
    inputsSpec.inputs.LabelsConfigFile = LABELS_CONFIG_FILE

    outputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=[
                "FA_stats",
                "MD_stats",
                "RD_stats",
                "FrobeniusNorm_stats",
                "Lambda1_stats",
                "Lambda2_stats",
                "Lambda3_stats",
            ]
        ),
        name="outputsSpec",
    )

    # Step1: Create the labelmap volume for DWI scan
    CreateDWILabelMapNode = pe.Node(
        interface=Function(
            function=create_dwi_label_map,
            input_names=["T2LabelMapVolume", "DWIBrainMask"],
            output_names=["DWILabelMapVolume"],
        ),
        name="create_dwi_label_map",
    )
    MeasurementWF.connect(
        inputsSpec, "T2LabelMapVolume", CreateDWILabelMapNode, "T2LabelMapVolume"
    )
    MeasurementWF.connect(
        inputsSpec, "DWIBrainMask", CreateDWILabelMapNode, "DWIBrainMask"
    )

    # Now we have two labelmap volumes (both have 1x1x1 voxel lattice):
    # (1) T2LabelMap: Used to compute total_volume for each label
    # (2) DWILabelMap: It is probably cropped and missed some part of labels,
    #                  and is used to compute all stats like [mean,std,max,min,median,effective_volume].

    # Step2: Resample each RIS to T2LabelmapVolume voxel lattice
    MakeResamplerInFilesListNode = pe.Node(
        interface=Function(
            function=make_resampled_in_file_list,
            input_names=[
                "FAImage",
                "MDImage",
                "RDImage",
                "FrobeniusNormImage",
                "Lambda1Image",
                "Lambda2Image",
                "Lambda3Image",
            ],
            output_names=["RISsList"],
        ),
        name="MakeResamplerInFilesListNode",
    )
    MeasurementWF.connect(
        [
            (
                inputsSpec,
                MakeResamplerInFilesListNode,
                [
                    ("FAImage", "FAImage"),
                    ("MDImage", "MDImage"),
                    ("RDImage", "RDImage"),
                    ("FrobeniusNormImage", "FrobeniusNormImage"),
                    ("Lambda1Image", "Lambda1Image"),
                    ("Lambda2Image", "Lambda2Image"),
                    ("Lambda3Image", "Lambda3Image"),
                ],
            )
        ]
    )
    # To resample RIS volumes we should consider that the output of resampling
    # should not have any negative intensity value becuase negative values have no
    # meaning in rotationally invariant scalar measures.
    # There are 3 options:
    # 1- Use linear interpolation (commented out part using BRAINSResample)
    # 2- Use Gaussian interpolation
    # 3- Use cubic BSpline interpolation
    # Third option is chosen here, but with some considerations:
    # "voxel-wise squared root of intensity values" +
    # cubic BSpline interpolation +
    # "voxel-wise square of intesity values"
    ResampleRISsNode = pe.MapNode(
        interface=Function(
            function=resample_ris_volumes,
            input_names=["referenceVolume", "inputVolume"],
            output_names=["outputVolume"],
        ),
        name="ResampleRISs",
        iterfield=["inputVolume"],
    )
    MeasurementWF.connect(
        inputsSpec, "T2LabelMapVolume", ResampleRISsNode, "referenceVolume"
    )
    MeasurementWF.connect(
        MakeResamplerInFilesListNode, "RISsList", ResampleRISsNode, "inputVolume"
    )
    """
    ResampleRISsNode = pe.MapNode(interface=BRAINSResample(), name="ResampleRISs",
                                  iterfield=['inputVolume', 'outputVolume'])
    ResampleRISsNode.inputs.interpolationMode = 'Linear'
    ResampleRISsNode.inputs.pixelType = 'float'
    ResampleRISsNode.inputs.outputVolume = ['FA_res.nrrd','MD_res.nrrd','RD_res.nrrd','frobenius_norm_res.nrrd',
                                            'lambda1_res.nrrd','lambda2_res.nrrd','lambda3_res.nrrd']
    MeasurementWF.connect(inputsSpec,'T2LabelMapVolume',ResampleRISsNode,'referenceVolume')
    MeasurementWF.connect(MakeResamplerInFilesListNode,'RISsList',ResampleRISsNode,'inputVolume')
    """
    # Step3: Computes statistics of each resampled RIS over all input labels
    # and writes the results as a CSV file (a csv file for each RIS)
    ComputeStatisticsNode = pe.MapNode(
        interface=Function(
            function=compute_statistics,
            input_names=[
                "inputVolume",
                "T2LabelMapVolume",
                "DWILabelMapVolume",
                "labelCodesFile",
            ],
            output_names=["CSVStatisticsFile"],
        ),
        name="compute_statistics",
        iterfield=["inputVolume"],
    )
    MeasurementWF.connect(
        ResampleRISsNode, "outputVolume", ComputeStatisticsNode, "inputVolume"
    )
    MeasurementWF.connect(
        inputsSpec, "T2LabelMapVolume", ComputeStatisticsNode, "T2LabelMapVolume"
    )
    MeasurementWF.connect(
        CreateDWILabelMapNode,
        "DWILabelMapVolume",
        ComputeStatisticsNode,
        "DWILabelMapVolume",
    )
    MeasurementWF.connect(
        inputsSpec, "LabelsConfigFile", ComputeStatisticsNode, "labelCodesFile"
    )
    MeasurementWF.connect(
        ComputeStatisticsNode,
        ("CSVStatisticsFile", pick_from_file, 0),
        outputsSpec,
        "FA_stats",
    )
    MeasurementWF.connect(
        ComputeStatisticsNode,
        ("CSVStatisticsFile", pick_from_file, 1),
        outputsSpec,
        "MD_stats",
    )
    MeasurementWF.connect(
        ComputeStatisticsNode,
        ("CSVStatisticsFile", pick_from_file, 2),
        outputsSpec,
        "RD_stats",
    )
    MeasurementWF.connect(
        ComputeStatisticsNode,
        ("CSVStatisticsFile", pick_from_file, 3),
        outputsSpec,
        "FrobeniusNorm_stats",
    )
    MeasurementWF.connect(
        ComputeStatisticsNode,
        ("CSVStatisticsFile", pick_from_file, 4),
        outputsSpec,
        "Lambda1_stats",
    )
    MeasurementWF.connect(
        ComputeStatisticsNode,
        ("CSVStatisticsFile", pick_from_file, 5),
        outputsSpec,
        "Lambda2_stats",
    )
    MeasurementWF.connect(
        ComputeStatisticsNode,
        ("CSVStatisticsFile", pick_from_file, 6),
        outputsSpec,
        "Lambda3_stats",
    )

    return MeasurementWF
