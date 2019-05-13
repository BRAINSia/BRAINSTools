## \author Ali Ghayoor
##
## This workflow corrects the succeptiblity artifacts of the input DWI data,
## and aligns that to corresponding t2 image space
## Important outputs: CorrectedDWI_in_T2Space, DWIBrainMask

"""
CorrectionWorkflow.py
============================
Description:
    The purpose of this is to...
    
Author:
Usage:
"""

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

from utilities.misc import common_ants_registration_settings


def create_correction_workflow(WFname):
    """
    This Function takes in...

    :param WFname:
    :return: CorrectionWF
    """
    ###### UTILITY FUNCTIONS #######
    # \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
    # remove the skull from the T2 volume
    def extract_brain_from_head(RawScan, BrainLabels):
        """
        This function will remove the skull from the T2 volume

        :param Rawscan:
        :param BrainLabels:
        :return:
        """
        import os
        import SimpleITK as sitk

        # Remove skull from the head scan
        assert os.path.exists(RawScan), "File not found: %s" % RawScan
        assert os.path.exists(BrainLabels), "File not found: %s" % BrainLabels
        headImage = sitk.ReadImage(RawScan)
        labelsMap = sitk.ReadImage(BrainLabels)
        label_mask = labelsMap > 0
        brainImage = sitk.Cast(headImage, sitk.sitkInt16) * sitk.Cast(
            label_mask, sitk.sitkInt16
        )
        outputVolume = os.path.realpath("T2Stripped.nrrd")
        sitk.WriteImage(brainImage, outputVolume)
        return outputVolume

    def make_resampled_in_file_list(inputT2, inputLabelMap):
        """
        This function..

        :param inputT2:
        :param inputLabelMap:
        :return:
        """
        imagesList = [inputT2, inputLabelMap]
        return imagesList

    # This function helps to pick desirable output from the output list
    def pick_from_file(inlist, item):
        """
        This function will

        :param inlist:
        :param item:
        :return:
        """
        return inlist[item]

    # Create registration mask for ANTs from resampled label map image
    def create_ants_registration_mask(brainMask):
        """
        This function will

        :param brainmask:
        :return:
        """
        import os
        import SimpleITK as sitk

        assert os.path.exists(brainMask), "File not found: %s" % brainMask
        labelsMap = sitk.ReadImage(brainMask)
        label_mask = labelsMap > 0
        # dilate the label mask
        dilateFilter = sitk.BinaryDilateImageFilter()
        dilateFilter.SetKernelRadius(12)
        dilated_mask = dilateFilter.Execute(label_mask)
        regMask = dilated_mask
        registrationMask = os.path.realpath("registrationMask.nrrd")
        sitk.WriteImage(regMask, registrationMask)
        return registrationMask

    # Save direction cosine for the input volume
    def save_direction_cosine_to_matrix(inputVolume):
        """
        This function will return the direction cosine for the input volume

        :param inputVolume:
        :return:
        """
        import os
        import SimpleITK as sitk

        assert os.path.exists(inputVolume), "File not found: %s" % inputVolume
        t2 = sitk.ReadImage(inputVolume)
        directionCosine = t2.GetDirection()
        return directionCosine

    def make_force_dc_file_list(inputB0, inputT2, inputLabelMap):
        """
        This function will

        :param inputB0:
        :param inputT2:
        :param inputLabelMap:
        :return:
        """
        import os

        assert os.path.exists(inputB0), "File not found: %s" % inputB0
        assert os.path.exists(inputT2), "File not found: %s" % inputT2
        assert os.path.exists(inputLabelMap), "File not found: %s" % inputLabelMap
        imagesList = [inputB0, inputT2, inputLabelMap]
        return imagesList

    # Force DC to ID
    def force_dc_to_id(inputVolume):
        """
        This function will force DC to ID

        :param inputVolume:
        :return:
        """
        import os
        import SimpleITK as sitk

        inImage = sitk.ReadImage(inputVolume)
        inImage.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
        outputVolume = os.path.realpath("IDDC_" + os.path.basename(inputVolume))
        sitk.WriteImage(inImage, outputVolume)
        return outputVolume

    def restore_dc_from_saved_matrix(inputVolume, inputDirectionCosine):
        """
        This function will

        :param inputVolume:
        :param inputDirectionCosine:
        :return:
        """
        import os
        import SimpleITK as sitk

        inImage = sitk.ReadImage(inputVolume)
        inImage.SetDirection(inputDirectionCosine)
        outputVolume = os.path.realpath("CorrectedDWI.nrrd")
        sitk.WriteImage(inImage, outputVolume)
        return outputVolume

    def get_rigid_transform_inverse(inputTransform):
        """
        This function will

        :param inputTransform:
        :return:
        """
        import os
        import SimpleITK as sitk

        inputTx = sitk.ReadTransform(inputTransform)
        versorRigidTx = sitk.VersorRigid3DTransform()
        versorRigidTx.SetFixedParameters(inputTx.GetFixedParameters())
        versorRigidTx.SetParameters(inputTx.GetParameters())
        invTx = versorRigidTx.GetInverse()
        inverseTransform = os.path.realpath(
            "Inverse_" + os.path.basename(inputTransform)
        )
        sitk.WriteTransform(invTx, inverseTransform)
        return inverseTransform

    #################################
    # \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    CorrectionWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(
        interface=IdentityInterface(fields=["T2Volume", "DWIVolume", "LabelMapVolume"]),
        name="inputsSpec",
    )

    outputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=["CorrectedDWI", "CorrectedDWI_in_T2Space", "DWIBrainMask"]
        ),
        name="outputsSpec",
    )

    # Step0: remove the skull from the T2 volume
    ExtractBRAINFromHeadNode = pe.Node(
        interface=Function(
            function=extract_brain_from_head,
            input_names=["RawScan", "BrainLabels"],
            output_names=["outputVolume"],
        ),
        name="extract_brain_from_head",
    )

    CorrectionWF.connect(inputsSpec, "T2Volume", ExtractBRAINFromHeadNode, "RawScan")
    CorrectionWF.connect(
        inputsSpec, "LabelMapVolume", ExtractBRAINFromHeadNode, "BrainLabels"
    )

    # Step1: extract B0 from DWI volume
    EXTRACT_B0 = pe.Node(interface=extractNrrdVectorIndex(), name="EXTRACT_B0")
    EXTRACT_B0.inputs.vectorIndex = 0
    EXTRACT_B0.inputs.outputVolume = "B0_Image.nrrd"
    CorrectionWF.connect(inputsSpec, "DWIVolume", EXTRACT_B0, "inputVolume")

    # Step2: Register T2 to B0 space using BRAINSFit
    BFit_T2toB0 = pe.Node(interface=BRAINSFit(), name="BFit_T2toB0")
    BFit_T2toB0.inputs.costMetric = "MMI"
    BFit_T2toB0.inputs.numberOfSamples = 100000
    BFit_T2toB0.inputs.numberOfIterations = [1500]
    BFit_T2toB0.inputs.numberOfHistogramBins = 50
    BFit_T2toB0.inputs.maximumStepLength = 0.2
    BFit_T2toB0.inputs.minimumStepLength = [0.00005]
    BFit_T2toB0.inputs.useRigid = True
    BFit_T2toB0.inputs.useAffine = True
    BFit_T2toB0.inputs.maskInferiorCutOffFromCenter = 65
    BFit_T2toB0.inputs.maskProcessingMode = "ROIAUTO"
    BFit_T2toB0.inputs.ROIAutoDilateSize = 13
    BFit_T2toB0.inputs.backgroundFillValue = 0.0
    BFit_T2toB0.inputs.initializeTransformMode = "useCenterOfHeadAlign"
    BFit_T2toB0.inputs.strippedOutputTransform = "T2ToB0_RigidTransform.h5"
    BFit_T2toB0.inputs.writeOutputTransformInFloat = True
    CorrectionWF.connect(EXTRACT_B0, "outputVolume", BFit_T2toB0, "fixedVolume")
    CorrectionWF.connect(
        ExtractBRAINFromHeadNode, "outputVolume", BFit_T2toB0, "movingVolume"
    )

    # Step3: Use T_rigid to "resample" T2 and label map images to B0 image space
    MakeResamplerInFilesListNode = pe.Node(
        Function(
            function=make_resampled_in_file_list,
            input_names=["inputT2", "inputLabelMap"],
            output_names=["imagesList"],
        ),
        name="MakeResamplerInFilesListNode",
    )
    CorrectionWF.connect(
        [
            (
                ExtractBRAINFromHeadNode,
                MakeResamplerInFilesListNode,
                [("outputVolume", "inputT2")],
            ),
            (
                inputsSpec,
                MakeResamplerInFilesListNode,
                [("LabelMapVolume", "inputLabelMap")],
            ),
        ]
    )

    ResampleToB0Space = pe.MapNode(
        interface=BRAINSResample(),
        name="ResampleToB0Space",
        iterfield=["inputVolume", "pixelType", "outputVolume"],
    )
    ResampleToB0Space.inputs.interpolationMode = "Linear"
    ResampleToB0Space.inputs.outputVolume = ["T2toB0.nrrd", "BRAINMaskToB0.nrrd"]
    ResampleToB0Space.inputs.pixelType = ["ushort", "binary"]
    CorrectionWF.connect(
        BFit_T2toB0, "strippedOutputTransform", ResampleToB0Space, "warpTransform"
    )
    CorrectionWF.connect(
        EXTRACT_B0, "outputVolume", ResampleToB0Space, "referenceVolume"
    )
    CorrectionWF.connect(
        MakeResamplerInFilesListNode, "imagesList", ResampleToB0Space, "inputVolume"
    )

    # Step4: Create registration mask from resampled label map image
    CreateRegistrationMask = pe.Node(
        interface=Function(
            function=create_ants_registration_mask,
            input_names=["brainMask"],
            output_names=["registrationMask"],
        ),
        name="create_ants_registration_mask",
    )
    CorrectionWF.connect(
        ResampleToB0Space,
        ("outputVolume", pick_from_file, 1),
        CreateRegistrationMask,
        "brainMask",
    )

    # Step5: Save direction cosine for the resampled T2 image
    SaveDirectionCosineToMatrixNode = pe.Node(
        interface=Function(
            function=save_direction_cosine_to_matrix,
            input_names=["inputVolume"],
            output_names=["directionCosine"],
        ),
        name="save_direction_cosine_to_matrix",
    )
    CorrectionWF.connect(
        ResampleToB0Space,
        ("outputVolume", pick_from_file, 0),
        SaveDirectionCosineToMatrixNode,
        "inputVolume",
    )

    # Step6: Force DC to ID
    MakeForceDCFilesListNode = pe.Node(
        Function(
            function=make_force_dc_file_list,
            input_names=["inputB0", "inputT2", "inputLabelMap"],
            output_names=["imagesList"],
        ),
        name="MakeForceDCFilesListNode",
    )
    CorrectionWF.connect(
        [
            (EXTRACT_B0, MakeForceDCFilesListNode, [("outputVolume", "inputB0")]),
            (
                ResampleToB0Space,
                MakeForceDCFilesListNode,
                [(("outputVolume", pick_from_file, 0), "inputT2")],
            ),
            (
                CreateRegistrationMask,
                MakeForceDCFilesListNode,
                [("registrationMask", "inputLabelMap")],
            ),
        ]
    )

    ForceDCtoIDNode = pe.MapNode(
        interface=Function(
            function=force_dc_to_id,
            input_names=["inputVolume"],
            output_names=["outputVolume"],
        ),
        name="force_dc_to_id",
        iterfield=["inputVolume"],
    )
    CorrectionWF.connect(
        MakeForceDCFilesListNode, "imagesList", ForceDCtoIDNode, "inputVolume"
    )

    # Step7: Run antsRegistration in one direction
    antsReg_B0ToTransformedT2 = pe.Node(
        interface=ants.Registration(), name="antsReg_B0ToTransformedT2"
    )
    antsReg_B0ToTransformedT2.inputs.interpolation = "Linear"
    antsReg_B0ToTransformedT2.inputs.dimension = 3
    antsReg_B0ToTransformedT2.inputs.transforms = ["SyN"]
    antsReg_B0ToTransformedT2.inputs.transform_parameters = [(0.25, 3.0, 0.0)]
    antsReg_B0ToTransformedT2.inputs.metric = ["MI"]
    antsReg_B0ToTransformedT2.inputs.sampling_strategy = [None]
    antsReg_B0ToTransformedT2.inputs.sampling_percentage = [1.0]
    antsReg_B0ToTransformedT2.inputs.metric_weight = [1.0]
    antsReg_B0ToTransformedT2.inputs.radius_or_number_of_bins = [32]
    antsReg_B0ToTransformedT2.inputs.number_of_iterations = [[70, 50, 40]]
    antsReg_B0ToTransformedT2.inputs.convergence_threshold = [1e-6]
    antsReg_B0ToTransformedT2.inputs.convergence_window_size = [10]
    antsReg_B0ToTransformedT2.inputs.use_histogram_matching = [True]
    antsReg_B0ToTransformedT2.inputs.shrink_factors = [[3, 2, 1]]
    antsReg_B0ToTransformedT2.inputs.smoothing_sigmas = [[2, 1, 0]]
    antsReg_B0ToTransformedT2.inputs.sigma_units = ["vox"]
    antsReg_B0ToTransformedT2.inputs.use_estimate_learning_rate_once = [False]
    antsReg_B0ToTransformedT2.inputs.write_composite_transform = True
    antsReg_B0ToTransformedT2.inputs.collapse_output_transforms = False
    antsReg_B0ToTransformedT2.inputs.initialize_transforms_per_stage = False
    antsReg_B0ToTransformedT2.inputs.output_transform_prefix = "Tsyn"
    antsReg_B0ToTransformedT2.inputs.winsorize_lower_quantile = 0.01
    antsReg_B0ToTransformedT2.inputs.winsorize_upper_quantile = 0.99
    antsReg_B0ToTransformedT2.inputs.float = True
    antsReg_B0ToTransformedT2.inputs.num_threads = -1
    antsReg_B0ToTransformedT2.inputs.args = "--restrict-deformation 0x1x0"
    CorrectionWF.connect(
        ForceDCtoIDNode,
        ("outputVolume", pick_from_file, 1),
        antsReg_B0ToTransformedT2,
        "fixed_image",
    )
    CorrectionWF.connect(
        ForceDCtoIDNode,
        ("outputVolume", pick_from_file, 2),
        antsReg_B0ToTransformedT2,
        "fixed_image_masks",
    )
    CorrectionWF.connect(
        ForceDCtoIDNode,
        ("outputVolume", pick_from_file, 0),
        antsReg_B0ToTransformedT2,
        "moving_image",
    )

    # Step8: Now, all necessary transforms are acquired. It's a time to
    #        transform input DWI image into T2 image space
    # {DWI} --> force_dc_to_id --> gtractResampleDWIInPlace(using SyN transfrom)
    # --> Restore DirectionCosine From Saved Matrix --> gtractResampleDWIInPlace(inverse of T_rigid from BFit)
    # --> {CorrectedDW_in_T2Space}
    DWI_ForceDCtoIDNode = pe.Node(
        interface=Function(
            function=force_dc_to_id,
            input_names=["inputVolume"],
            output_names=["outputVolume"],
        ),
        name="DWI_ForceDCtoIDNode",
    )
    CorrectionWF.connect(inputsSpec, "DWIVolume", DWI_ForceDCtoIDNode, "inputVolume")

    gtractResampleDWI_SyN = pe.Node(
        interface=gtractResampleDWIInPlace(), name="gtractResampleDWI_SyN"
    )
    CorrectionWF.connect(
        DWI_ForceDCtoIDNode, "outputVolume", gtractResampleDWI_SyN, "inputVolume"
    )
    CorrectionWF.connect(
        antsReg_B0ToTransformedT2,
        "composite_transform",
        gtractResampleDWI_SyN,
        "warpDWITransform",
    )
    CorrectionWF.connect(
        ForceDCtoIDNode,
        ("outputVolume", pick_from_file, 1),
        gtractResampleDWI_SyN,
        "referenceVolume",
    )  # fixed image of antsRegistration
    gtractResampleDWI_SyN.inputs.outputVolume = "IDDC_correctedDWI.nrrd"

    RestoreDCFromSavedMatrixNode = pe.Node(
        interface=Function(
            function=restore_dc_from_saved_matrix,
            input_names=["inputVolume", "inputDirectionCosine"],
            output_names=["outputVolume"],
        ),
        name="restore_dc_from_saved_matrix",
    )
    CorrectionWF.connect(
        gtractResampleDWI_SyN,
        "outputVolume",
        RestoreDCFromSavedMatrixNode,
        "inputVolume",
    )
    CorrectionWF.connect(
        SaveDirectionCosineToMatrixNode,
        "directionCosine",
        RestoreDCFromSavedMatrixNode,
        "inputDirectionCosine",
    )
    CorrectionWF.connect(
        RestoreDCFromSavedMatrixNode, "outputVolume", outputsSpec, "CorrectedDWI"
    )

    GetRigidTransformInverseNode = pe.Node(
        interface=Function(
            function=get_rigid_transform_inverse,
            input_names=["inputTransform"],
            output_names=["inverseTransform"],
        ),
        name="get_rigid_transform_inverse",
    )
    CorrectionWF.connect(
        BFit_T2toB0,
        "strippedOutputTransform",
        GetRigidTransformInverseNode,
        "inputTransform",
    )

    gtractResampleDWIInPlace_Trigid = pe.Node(
        interface=gtractResampleDWIInPlace(), name="gtractResampleDWIInPlace_Trigid"
    )
    CorrectionWF.connect(
        RestoreDCFromSavedMatrixNode,
        "outputVolume",
        gtractResampleDWIInPlace_Trigid,
        "inputVolume",
    )
    CorrectionWF.connect(
        GetRigidTransformInverseNode,
        "inverseTransform",
        gtractResampleDWIInPlace_Trigid,
        "inputTransform",
    )  # Inverse of rigid transform from BFit
    gtractResampleDWIInPlace_Trigid.inputs.outputVolume = (
        "CorrectedDWI_in_T2Space_estimate.nrrd"
    )
    gtractResampleDWIInPlace_Trigid.inputs.outputResampledB0 = (
        "CorrectedDWI_in_T2Space_estimate_B0.nrrd"
    )

    # Setp9: An extra registration step to tune the alignment between the CorrecetedDWI_in_T2Space image and T2 image.
    BFit_TuneRegistration = pe.Node(interface=BRAINSFit(), name="BFit_TuneRegistration")
    BFit_TuneRegistration.inputs.costMetric = "MMI"
    BFit_TuneRegistration.inputs.numberOfSamples = 100000
    BFit_TuneRegistration.inputs.numberOfIterations = [1500]
    BFit_TuneRegistration.inputs.numberOfHistogramBins = 50
    BFit_TuneRegistration.inputs.maximumStepLength = 0.2
    BFit_TuneRegistration.inputs.minimumStepLength = [0.00005]
    BFit_TuneRegistration.inputs.useRigid = True
    BFit_TuneRegistration.inputs.useAffine = True
    BFit_TuneRegistration.inputs.maskInferiorCutOffFromCenter = 65
    BFit_TuneRegistration.inputs.maskProcessingMode = "ROIAUTO"
    BFit_TuneRegistration.inputs.ROIAutoDilateSize = 13
    BFit_TuneRegistration.inputs.backgroundFillValue = 0.0
    BFit_TuneRegistration.inputs.initializeTransformMode = "useCenterOfHeadAlign"
    BFit_TuneRegistration.inputs.strippedOutputTransform = (
        "CorrectedB0inT2Space_to_T2_RigidTransform.h5"
    )
    BFit_TuneRegistration.inputs.writeOutputTransformInFloat = True
    CorrectionWF.connect(
        ExtractBRAINFromHeadNode, "outputVolume", BFit_TuneRegistration, "fixedVolume"
    )  # T2 brain volume
    CorrectionWF.connect(
        gtractResampleDWIInPlace_Trigid,
        "outputResampledB0",
        BFit_TuneRegistration,
        "movingVolume",
    )  # CorrectedB0_in_T2Space

    gtractResampleDWIInPlace_TuneRigidTx = pe.Node(
        interface=gtractResampleDWIInPlace(),
        name="gtractResampleDWIInPlace_TuneRigidTx",
    )
    CorrectionWF.connect(
        gtractResampleDWIInPlace_Trigid,
        "outputVolume",
        gtractResampleDWIInPlace_TuneRigidTx,
        "inputVolume",
    )
    CorrectionWF.connect(
        BFit_TuneRegistration,
        "strippedOutputTransform",
        gtractResampleDWIInPlace_TuneRigidTx,
        "inputTransform",
    )
    gtractResampleDWIInPlace_TuneRigidTx.inputs.outputVolume = (
        "CorrectedDWI_in_T2Space.nrrd"
    )
    gtractResampleDWIInPlace_TuneRigidTx.inputs.outputResampledB0 = (
        "CorrectedDWI_in_T2Space_B0.nrrd"
    )

    # Finally we pass the outputs of the gtractResampleDWIInPlace_TuneRigidTx to the outputsSpec
    CorrectionWF.connect(
        gtractResampleDWIInPlace_TuneRigidTx,
        "outputVolume",
        outputsSpec,
        "CorrectedDWI_in_T2Space",
    )

    # Step10: Create brain mask from the input labelmap
    DWIBRAINMASK = pe.Node(interface=BRAINSResample(), name="DWIBRAINMASK")
    DWIBRAINMASK.inputs.interpolationMode = "Linear"
    DWIBRAINMASK.inputs.outputVolume = "BrainMaskForDWI.nrrd"
    DWIBRAINMASK.inputs.pixelType = "binary"
    CorrectionWF.connect(
        gtractResampleDWIInPlace_TuneRigidTx,
        "outputResampledB0",
        DWIBRAINMASK,
        "referenceVolume",
    )
    CorrectionWF.connect(inputsSpec, "LabelMapVolume", DWIBRAINMASK, "inputVolume")
    CorrectionWF.connect(DWIBRAINMASK, "outputVolume", outputsSpec, "DWIBrainMask")

    return CorrectionWF
