#!/usr/bin/env python
# ################################################################################
## Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
## Language:  Python
##
## Author:  Hans J. Johnson, David Welch
##
##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.
##
#################################################################################


"""
 baseline.py
============================
Description:
    This software is distributed WITHOUT ANY WARRANTY; without even
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     See the above copyright notices for more information.

Author:
    Hans J. Johnson
    David Welch

Usage:

"""

import os
from builtins import str

import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.utils.misc import package_check

# package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
package_check("numpy", "1.3", "tutorial1")
package_check("scipy", "0.7", "tutorial1")
# THIS IS NOT REQUIRED package_check('matplotlib','1.4','turorial1')
package_check("networkx", "1.0", "tutorial1")
package_check("IPython", "0.10", "tutorial1")

from utilities.distributed import modify_qsub_args
from utilities.image_processing import fix_wm_partitioning
from PipeLineFunctionHelpers import accumulate_like_tissue_posteriors
from PipeLineFunctionHelpers import unwrap_posterior_images_from_list_tuples_function

from .WorkupT1T2LandmarkInitialization import create_landmark_initialize_workflow
from .WorkupT1T2TissueClassify import create_tissue_classify_workflow
from .WorkupJointFusion import create_joint_fusion_workflow
from .WorkupAddsonBrainStem import create_brainstem_workflow

from utilities.misc import *

try:
    from nipype.interfaces.semtools import *
except ImportError:
    from AutoWorkup.semtools import *

from nipype.interfaces.semtools.registration.brainsresample import BRAINSResample

# from nipype.interfaces.semtools.filtering.denoising import UnbiasedNonLocalMeans
from nipype.interfaces.ants.segmentation import DenoiseImage
from nipype.interfaces.ants.segmentation import N4BiasFieldCorrection
from nipype.interfaces.semtools.segmentation.specialized import (
    BRAINSCreateLabelMapFromProbabilityMaps,
)


def get_list_element(nestedList, index):
    """
    This function...

    :param nestedList:
    :param index:
    :return:
    """
    return nestedList[index]


def get_all_t1s_length(allT1s):
    """
    This function...

    :param allT1s:
    :return:
    """
    return len(allT1s)


##INFO:  Move to module that can be re-used
##       get_largest_label is copied elsewhere
def create_left_right_wm_hemispheres(
    BRAINLABELSFile,
    HDCMARegisteredVentricleMaskFN,
    LeftHemisphereMaskName,
    RightHemisphereMaskName,
    WM_LeftHemisphereFileName,
    WM_RightHemisphereFileName,
):
    """
    This function..
    :param BRAINLABELSFile:
    :param HDCMARegisteredVentricleMaskFN:
    :param LeftHemisphereMaskName:
    :param RightHemisphereMaskName:
    :param WM_LeftHemisphereFileName:
    :param WM_RightHemisphereFileName:
    :return:
    """
    import SimpleITK as sitk
    import os

    def get_largest_label(inputMask, UseErosionCleaning):
        """
        This function...

        :param inputMask:
        :param UseErosionCleaning:
        :return:
        """
        LargestComponentCode = 1
        if UseErosionCleaning:
            erosionMask = sitk.ErodeObjectMorphology(inputMask, 1)
        else:
            erosionMask = inputMask
        CC = sitk.ConnectedComponent(erosionMask)
        Rlabel = sitk.RelabelComponent(CC)
        largestMask = Rlabel == LargestComponentCode
        if UseErosionCleaning:
            dilateMask = sitk.DilateObjectMorphology(largestMask, 1)
        else:
            dilateMask = largestMask

        return largestMask * dilateMask > 0

    ABCLabelsImage = sitk.Cast(sitk.ReadImage(BRAINLABELSFile), sitk.sitkUInt32)
    # # Remove brain stem and cerebellum
    BS = ABCLabelsImage == 30
    Cerebellum_GM = ABCLabelsImage == 12
    Cerebellum_WM = ABCLabelsImage == 11
    KeepRegion = sitk.Cast((1 - (BS + Cerebellum_GM + Cerebellum_WM)), sitk.sitkUInt32)

    ABCLabelsImage = KeepRegion * ABCLabelsImage

    HDCMARegisteredVentricleLabels = sitk.Cast(
        sitk.ReadImage(HDCMARegisteredVentricleMaskFN), sitk.sitkUInt32
    )
    ABCCSFLabelCode = 4
    HDMCALeftVentricleCode = 4
    HDMCARightVentricleCode = 43
    HDCMAMask = (HDCMARegisteredVentricleLabels == HDMCALeftVentricleCode) + (
        HDCMARegisteredVentricleLabels == HDMCARightVentricleCode
    )
    ExpandVentValue = 5
    HDCMAMask_d5 = sitk.DilateObjectMorphology(HDCMAMask, ExpandVentValue)
    CSFMaskImage = ABCLabelsImage == ABCCSFLabelCode
    VentricleMask = (HDCMAMask_d5 * CSFMaskImage + HDCMAMask) > 0
    VentricleMask_d2 = sitk.DilateObjectMorphology(VentricleMask, 2)
    ABCWMLabelCode = 1
    WMMaskImage = ABCLabelsImage == ABCWMLabelCode

    subcorticalRegions = (
        ABCLabelsImage >= 12
    )  # All subcortical regions are listed greater than equal to values of 12
    WMSubcortFilled = (WMMaskImage + subcorticalRegions) > 0

    WMSubcortFilled_CC = get_largest_label(WMSubcortFilled, False)
    WMSubcortFilled_CC_Ventricles = (WMSubcortFilled_CC + VentricleMask_d2) > 0
    neg_WMSubcortFilled_CC = 1 - WMSubcortFilled_CC
    neg_WMSubcortFilled_CC_bg = get_largest_label(neg_WMSubcortFilled_CC, False)
    neg_WMSubcortFilled_CC_bg_holes = neg_WMSubcortFilled_CC - neg_WMSubcortFilled_CC_bg

    WM_Final = neg_WMSubcortFilled_CC_bg_holes + WMSubcortFilled_CC_Ventricles > 0

    ####################################
    ### START WM

    # Template masks for left and right hemispheres
    Left_template = (
        sitk.Cast(sitk.ReadImage(LeftHemisphereMaskName), sitk.sitkUInt32) > 0
    )
    Right_template = (
        sitk.Cast(sitk.ReadImage(RightHemisphereMaskName), sitk.sitkUInt32) > 0
    )

    # Split into left and right hemispheres
    WM_left = WM_Final * Left_template > 0
    WM_right = WM_Final * Right_template > 0

    WM_Largest_left = get_largest_label(WM_left, False)
    WM_Largest_right = get_largest_label(WM_right, False)

    WM_left_extras = WM_left - WM_Largest_left
    WM_right_extras = WM_right - WM_Largest_right

    WM_Largest_left = get_largest_label(WM_Largest_left + WM_right_extras, False)
    WM_Largest_right = get_largest_label(WM_Largest_right + WM_left_extras, False)

    ## Write todisk
    WM_LeftHemisphereFileName = os.path.abspath(WM_LeftHemisphereFileName)
    sitk.WriteImage(WM_Largest_left, WM_LeftHemisphereFileName)

    WM_RightHemisphereFileName = os.path.abspath(WM_RightHemisphereFileName)
    sitk.WriteImage(WM_Largest_right, WM_RightHemisphereFileName)

    ## TODO Add splitting into hemispheres code here
    return WM_LeftHemisphereFileName, WM_RightHemisphereFileName


def generate_single_session_template_wf(
    projectid,
    subjectid,
    sessionid,
    onlyT1,
    hasPDs,
    hasFLs,
    master_config,
    phase,
    interpMode,
    pipeline_name,
    doDenoise=True,
    badT2=False,
    useEMSP=None,
):
    """
    Run autoworkup on a single sessionid

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.

    :param projectid:
    :param subjectid:
    :param sessionid:
    :param onlyT1:
    :param hasPDs:
    :param hasFLs:
    :param master_config:
    :param pipeline_name:
    :param doDenoise:
    :param badT2:
    :param useEMSP:
    :return:
    """

    # if  not 'landmark' in master_config['components'] or not 'auxlmk' in master_config['components'] or not 'tissue_classify' in master_config['components']:
    #    print "Baseline DataSink requires 'AUXLMK' and/or 'TISSUE_CLASSIFY'!!!"
    #    raise NotImplementedError
    # master_config['components'].append('auxlmk')
    # master_config['components'].append('tissue_classify')

    assert phase in [
        "atlas-based-reference",
        "subject-based-reference",
    ], "Unknown phase! Valid entries: 'atlas-based-reference', 'subject-based-reference'"

    if "tissue_classify" in master_config["components"]:
        assert (
            "landmark" in master_config["components"]
        ), "tissue_classify Requires landmark step!"
    # NOT TRUE if 'landmark' in master_config['components']:
    #    assert 'denoise' in master_config['components'], "landmark Requires denoise step!"

    if "jointfusion_2015_wholebrain" in master_config["components"]:
        assert (
            "warp_atlas_to_subject" in master_config["components"]
        ), "jointfusion_2015_wholebrain requires warp_atlas_to_subject!"

    from workflows.atlasNode import make_atlas_node

    baw201 = pe.Workflow(name=pipeline_name)

    inputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=[
                "atlasLandmarkFilename",
                "atlasWeightFilename",
                "LLSModel",
                "inputTemplateModel",
                "template_t1_denoised_gaussian",
                "atlasDefinition",
                "T1s",
                "T2s",
                "PDs",
                "FLs",
                "OTHERs",
                "EMSP",
                "hncma_atlas",
                "template_rightHemisphere",
                "template_leftHemisphere",
                "template_WMPM2_labels",
                "template_nac_labels",
                "template_ventricles",
                "template_headregion",
            ]
        ),
        run_without_submitting=True,
        name="inputspec",
    )

    outputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=[
                "t1_average",
                "t2_average",
                "pd_average",
                "fl_average",
                "posteriorImages",
                "outputLabels",
                "outputHeadLabels",
                "atlasToSubjectTransform",
                "atlasToSubjectInverseTransform",
                "atlasToSubjectRegistrationState",
                "BCD_ACPC_T1_CROPPED",
                "outputLandmarksInACPCAlignedSpace",
                "outputLandmarksInInputSpace",
                "output_tx",
                "LMIatlasToSubject_tx",
                "writeBranded2DImage",
                "brainStemMask",
                "UpdatedPosteriorsList",  # Longitudinal
            ]
        ),
        run_without_submitting=True,
        name="outputspec",
    )

    dsName = "{0}_ds_{1}_lmks".format(phase, sessionid)
    DataSinkLandmarks = pe.Node(name=dsName, interface=nio.DataSink())
    DataSinkLandmarks.overwrite = master_config["ds_overwrite"]
    DataSinkLandmarks.inputs.container = "{0}/{1}/{2}".format(
        projectid, subjectid, sessionid
    )
    DataSinkLandmarks.inputs.base_directory = master_config["resultdir"]
    del dsName

    atlas_static_directory = master_config["atlascache"]
    if master_config["workflow_phase"] == "atlas-based-reference":
        PostACPCAlignToAtlas = False
        atlas_warped_directory = master_config["atlascache"]
        atlasABCNode_XML = make_atlas_node(
            atlas_warped_directory,
            "BABCXMLAtlas_{0}".format(sessionid),
            ["W_BRAINSABCSupport"],
        )
        baw201.connect(
            atlasABCNode_XML,
            "ExtendedAtlasDefinition_xml",
            inputsSpec,
            "atlasDefinition",
        )

        atlasABCNode_W = make_atlas_node(
            atlas_warped_directory,
            "BABCAtlas_W{0}".format(sessionid),
            ["W_BRAINSABCSupport", "W_LabelMapsSupport"],
        )
        baw201.connect(
            [
                (
                    atlasABCNode_W,
                    inputsSpec,
                    [
                        ("hncma_atlas", "hncma_atlas"),
                        ("template_headregion", "template_headregion"),
                        ("template_leftHemisphere", "template_leftHemisphere"),
                        ("template_rightHemisphere", "template_rightHemisphere"),
                        ("template_WMPM2_labels", "template_WMPM2_labels"),
                        ("template_nac_labels", "template_nac_labels"),
                        ("template_ventricles", "template_ventricles"),
                    ],
                )
            ]
        )
        ## These landmarks are only relevant for the atlas-based-reference case
        atlasBCDNode_W = make_atlas_node(
            atlas_warped_directory, "BBCDAtlas_W{0}".format(sessionid), ["W_BCDSupport"]
        )
        baw201.connect(
            [
                (
                    atlasBCDNode_W,
                    inputsSpec,
                    [
                        (
                            "template_t1_denoised_gaussian",
                            "template_t1_denoised_gaussian",
                        ),
                        ("template_landmarks_50Lmks_fcsv", "atlasLandmarkFilename"),
                    ],
                )
            ]
        )

    elif master_config["workflow_phase"] == "subject-based-reference":
        PostACPCAlignToAtlas = True  # Use this subjects atlas image to align landmarks
        print((master_config["previousresult"]))
        atlas_warped_directory = os.path.join(
            master_config["previousresult"], subjectid, "Atlas"
        )

        template_DG = pe.Node(
            interface=nio.DataGrabber(
                infields=["subject"],
                outfields=[
                    "outAtlasXMLFullPath",
                    "hncma_atlas",
                    "template_leftHemisphere",
                    "template_rightHemisphere",
                    "template_WMPM2_labels",
                    "template_nac_labels",
                    "template_ventricles",
                    "template_t1_denoised_gaussian",
                    "template_landmarks_50Lmks_fcsv",
                    "template_headregion",
                ],
            ),
            name="Template_DG",
        )
        template_DG.inputs.base_directory = master_config["previousresult"]
        template_DG.inputs.subject = subjectid
        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        template_DG.inputs.field_template = OrderedDict(
            {
                "outAtlasXMLFullPath": "%s/Atlas/AtlasDefinition_%s.xml",
                "hncma_atlas": "%s/Atlas/AVG_hncma_atlas.nii.gz",
                "template_leftHemisphere": "%s/Atlas/AVG_template_leftHemisphere.nii.gz",
                "template_rightHemisphere": "%s/Atlas/AVG_template_rightHemisphere.nii.gz",
                "template_WMPM2_labels": "%s/Atlas/AVG_template_WMPM2_labels.nii.gz",
                "template_nac_labels": "%s/Atlas/AVG_template_nac_labels.nii.gz",
                "template_ventricles": "%s/Atlas/AVG_template_ventricles.nii.gz",
                "template_t1_denoised_gaussian": "%s/Atlas/AVG_T1.nii.gz",
                "template_landmarks_50Lmks_fcsv": "%s/Atlas/AVG_LMKS.fcsv",
                "template_headregion": "%s/Atlas/AVG_template_headregion.nii.gz",
            }
        )
        template_DG.inputs.template_args = OrderedDict(
            {
                "outAtlasXMLFullPath": [["subject", "subject"]],
                "hncma_atlas": [["subject"]],
                "template_leftHemisphere": [["subject"]],
                "template_rightHemisphere": [["subject"]],
                "template_WMPM2_labels": [["subject"]],
                "template_nac_labels": [["subject"]],
                "template_ventricles": [["subject"]],
                "template_t1_denoised_gaussian": [["subject"]],
                "template_landmarks_50Lmks_fcsv": [["subject"]],
                "template_headregion": [["subject"]],
            }
        )
        template_DG.inputs.template = "*"
        template_DG.inputs.sort_filelist = True
        template_DG.inputs.raise_on_empty = True

        baw201.connect(
            template_DG, "outAtlasXMLFullPath", inputsSpec, "atlasDefinition"
        )
        baw201.connect(
            [
                (
                    template_DG,
                    inputsSpec,
                    [
                        ## Already connected ('template_t1_denoised_gaussian','template_t1_denoised_gaussian'),
                        ("hncma_atlas", "hncma_atlas"),
                        ("template_leftHemisphere", "template_leftHemisphere"),
                        ("template_rightHemisphere", "template_rightHemisphere"),
                        ("template_WMPM2_labels", "template_WMPM2_labels"),
                        ("template_nac_labels", "template_nac_labels"),
                        ("template_ventricles", "template_ventricles"),
                    ],
                )
            ]
        )
        ## These landmarks are only relevant for the atlas-based-reference case
        baw201.connect(
            [
                (
                    template_DG,
                    inputsSpec,
                    [
                        (
                            "template_t1_denoised_gaussian",
                            "template_t1_denoised_gaussian",
                        ),
                        ("template_landmarks_50Lmks_fcsv", "atlasLandmarkFilename"),
                        ("template_headregion", "template_headregion"),
                    ],
                )
            ]
        )

    else:
        assert 0 == 1, "Invalid workflow type specified for singleSession"

    atlasBCDNode_S = make_atlas_node(
        atlas_static_directory, "BBCDAtlas_S{0}".format(sessionid), ["S_BCDSupport"]
    )
    baw201.connect(
        [
            (
                atlasBCDNode_S,
                inputsSpec,
                [
                    ("template_weights_50Lmks_wts", "atlasWeightFilename"),
                    ("LLSModel_50Lmks_h5", "LLSModel"),
                    ("T1_50Lmks_mdl", "inputTemplateModel"),
                ],
            )
        ]
    )

    if doDenoise:
        print("\ndenoise image filter\n")
        makeDenoiseInImageList = pe.Node(
            Function(
                function=make_out_from_file,
                input_names=[
                    "T1List",
                    "T2List",
                    "PDList",
                    "FLList",
                    "OTHERList",
                    "postfix",
                    "postfixBFC",
                    "postfixUnwrapped",
                    "PrimaryT1",
                    "ListOutType",
                ],
                output_names=[
                    "inImageList",
                    "outImageList",
                    "outBFCImageList",
                    "outUnwrappedImageList",
                    "imageTypeList",
                ],
            ),
            run_without_submitting=True,
            name="99_makeDenoiseInImageList",
        )
        baw201.connect(inputsSpec, "T1s", makeDenoiseInImageList, "T1List")
        baw201.connect(inputsSpec, "T2s", makeDenoiseInImageList, "T2List")
        baw201.connect(inputsSpec, "PDs", makeDenoiseInImageList, "PDList")
        baw201.connect(inputsSpec, "FLs", makeDenoiseInImageList, "FLList")
        baw201.connect(inputsSpec, "OTHERs", makeDenoiseInImageList, "OTHERList")
        makeDenoiseInImageList.inputs.ListOutType = False
        makeDenoiseInImageList.inputs.postfix = "_ants_denoised.nii.gz"
        makeDenoiseInImageList.inputs.postfixBFC = "_N4_BFC.nii.gz"
        makeDenoiseInImageList.inputs.postfixUnwrapped = "_unwrapped.nii.gz"
        makeDenoiseInImageList.inputs.PrimaryT1 = None  # an emptyList HACK

        print("\nDenoise:\n")
        DenoiseInputImgs = pe.MapNode(
            interface=DenoiseImage(),
            name="denoiseInputImgs",
            iterfield=["input_image", "output_image"],
        )
        DenoiseInputImgs.plugin_args = {
            "qsub_args": modify_qsub_args(master_config["queue"], 2, 4, 8),
            "overwrite": True,
        }
        DenoiseInputImgs.inputs.num_threads = -1
        DenoiseInputImgs.synchronize = True
        DenoiseInputImgs.inputs.dimension = 3

        # Rician has a bug in it as of 2016-02-08 DenoiseInputImgs.inputs.noise_model= 'Rician'
        # Rician bug fixed by Nick Tustison 2016-02-15
        DenoiseInputImgs.inputs.noise_model = "Rician"
        # DenoiseInputImgs.inputs.save_noise=True # we do need this until NIPYPE is fixed
        DenoiseInputImgs.inputs.save_noise = (
            False  # we don't need the noise image for BAW
        )
        DenoiseInputImgs.inputs.shrink_factor = 1  # default
        baw201.connect(
            [
                (
                    makeDenoiseInImageList,
                    DenoiseInputImgs,
                    [("inImageList", "input_image")],
                ),
                (
                    makeDenoiseInImageList,
                    DenoiseInputImgs,
                    [("outImageList", "output_image")],
                ),
            ]
        )

        print("\nN4BiasFieldCorrection:\n")
        N4BFC = pe.MapNode(
            interface=N4BiasFieldCorrection(),
            name="N4BFC",
            iterfield=["input_image", "output_image"],
        )
        N4BFC.plugin_args = {
            "qsub_args": modify_qsub_args(master_config["queue"], 4, 8, 8),
            "overwrite": True,
        }
        N4BFC.inputs.num_threads = -1
        N4BFC.inputs.dimension = 3
        N4BFC.inputs.bspline_fitting_distance = 200
        N4BFC.inputs.shrink_factor = 2
        N4BFC.inputs.n_iterations = [100, 100, 100, 75]
        N4BFC.inputs.convergence_threshold = 0.0000000001

        baw201.connect(
            [
                (DenoiseInputImgs, N4BFC, [("output_image", "input_image")]),
                (makeDenoiseInImageList, N4BFC, [("outBFCImageList", "output_image")]),
            ]
        )

        print("\nMerge all T1 and T2 List\n")
        makePreprocessingOutList = pe.Node(
            Function(
                function=generate_separate_image_type_list,
                input_names=["inFileList", "inTypeList"],
                output_names=["T1s", "T2s", "PDs", "FLs", "OTHERs"],
            ),
            run_without_submitting=False,
            name="99_makePreprocessingOutList",
        )
        baw201.connect(N4BFC, "output_image", makePreprocessingOutList, "inFileList")
        baw201.connect(
            makeDenoiseInImageList,
            "imageTypeList",
            makePreprocessingOutList,
            "inTypeList",
        )

    else:
        makePreprocessingOutList = inputsSpec

    if "landmark" in master_config["components"]:
        DoReverseMapping = False  # Set to true for debugging outputs
        if "auxlmk" in master_config["components"]:
            DoReverseMapping = True
        myLocalLMIWF = create_landmark_initialize_workflow(
            "LandmarkInitialize",
            master_config,
            interpMode,
            PostACPCAlignToAtlas,
            DoReverseMapping,
            useEMSP,
            Debug=False,
        )

        baw201.connect(
            [
                (
                    makePreprocessingOutList,
                    myLocalLMIWF,
                    [(("T1s", get_list_element, 0), "inputspec.inputVolume")],
                ),
                (
                    inputsSpec,
                    myLocalLMIWF,
                    [
                        ("atlasLandmarkFilename", "inputspec.atlasLandmarkFilename"),
                        ("atlasWeightFilename", "inputspec.atlasWeightFilename"),
                        ("LLSModel", "inputspec.LLSModel"),
                        ("inputTemplateModel", "inputspec.inputTemplateModel"),
                        ("template_t1_denoised_gaussian", "inputspec.atlasVolume"),
                        ("EMSP", "inputspec.EMSP"),
                    ],
                ),
                (
                    myLocalLMIWF,
                    outputsSpec,
                    [
                        (
                            "outputspec.outputResampledCroppedVolume",
                            "BCD_ACPC_T1_CROPPED",
                        ),
                        (
                            "outputspec.outputLandmarksInACPCAlignedSpace",
                            "outputLandmarksInACPCAlignedSpace",
                        ),
                        (
                            "outputspec.outputLandmarksInInputSpace",
                            "outputLandmarksInInputSpace",
                        ),
                        ("outputspec.outputTransform", "output_tx"),
                        ("outputspec.atlasToSubjectTransform", "LMIatlasToSubject_tx"),
                        ("outputspec.writeBranded2DImage", "writeBranded2DImage"),
                    ],
                ),
            ]
        )
        baw201.connect(
            [
                (
                    outputsSpec,
                    DataSinkLandmarks,  # INFO: change to myLocalLMIWF -> DataSink
                    [
                        (
                            "outputLandmarksInACPCAlignedSpace",
                            "ACPCAlign.@outputLandmarks_ACPC",
                        ),
                        ("writeBranded2DImage", "ACPCAlign.@writeBranded2DImage"),
                        ("BCD_ACPC_T1_CROPPED", "ACPCAlign.@BCD_ACPC_T1_CROPPED"),
                        (
                            "outputLandmarksInInputSpace",
                            "ACPCAlign.@outputLandmarks_Input",
                        ),
                        ("output_tx", "ACPCAlign.@output_tx"),
                        ("LMIatlasToSubject_tx", "ACPCAlign.@LMIatlasToSubject_tx"),
                    ],
                )
            ]
        )

    if "tissue_classify" in master_config["components"]:
        dsName = "{0}_ds_{1}_tissue".format(phase, sessionid)
        DataSinkTissue = pe.Node(name=dsName, interface=nio.DataSink())
        DataSinkTissue.overwrite = master_config["ds_overwrite"]
        DataSinkTissue.inputs.container = "{0}/{1}/{2}".format(
            projectid, subjectid, sessionid
        )
        DataSinkTissue.inputs.base_directory = master_config["resultdir"]
        del dsName

        useRegistrationMask = master_config["use_registration_masking"]

        myLocalTCWF = create_tissue_classify_workflow(
            "TissueClassify", master_config, interpMode, useRegistrationMask
        )
        baw201.connect(
            [
                (makePreprocessingOutList, myLocalTCWF, [("T1s", "inputspec.T1List")]),
                (makePreprocessingOutList, myLocalTCWF, [("T2s", "inputspec.T2List")]),
                (makePreprocessingOutList, myLocalTCWF, [("PDs", "inputspec.PDList")]),
                (makePreprocessingOutList, myLocalTCWF, [("FLs", "inputspec.FLList")]),
                (
                    makePreprocessingOutList,
                    myLocalTCWF,
                    [("OTHERs", "inputspec.OTHERList")],
                ),
                (
                    inputsSpec,
                    myLocalTCWF,
                    [
                        ("atlasDefinition", "inputspec.atlasDefinition"),
                        ("template_t1_denoised_gaussian", "inputspec.atlasVolume"),
                        ("template_headregion", "inputspec.atlasheadregion"),
                        (("T1s", get_all_t1s_length), "inputspec.T1_count"),
                    ],
                ),
                (
                    myLocalLMIWF,
                    myLocalTCWF,
                    [
                        (
                            "outputspec.outputResampledCroppedVolume",
                            "inputspec.PrimaryT1",
                        ),
                        (
                            "outputspec.atlasToSubjectTransform",
                            "inputspec.atlasToSubjectInitialTransform",
                        ),
                    ],
                ),
                (
                    myLocalTCWF,
                    outputsSpec,
                    [
                        ("outputspec.t1_average", "t1_average"),
                        ("outputspec.t2_average", "t2_average"),
                        ("outputspec.pd_average", "pd_average"),
                        ("outputspec.fl_average", "fl_average"),
                        ("outputspec.posteriorImages", "posteriorImages"),
                        ("outputspec.outputLabels", "outputLabels"),
                        ("outputspec.outputHeadLabels", "outputHeadLabels"),
                        (
                            "outputspec.atlasToSubjectTransform",
                            "atlasToSubjectTransform",
                        ),
                        (
                            "outputspec.atlasToSubjectInverseTransform",
                            "atlasToSubjectInverseTransform",
                        ),
                        (
                            "outputspec.atlasToSubjectRegistrationState",
                            "atlasToSubjectRegistrationState",
                        ),
                    ],
                ),
            ]
        )

        dsName = "{0}_ds_{1}_tissue_t1".format(phase, sessionid)
        DataSinkTissueT1 = pe.Node(name=dsName, interface=nio.DataSink())
        DataSinkTissueT1.overwrite = master_config["ds_overwrite"]
        DataSinkTissueT1.inputs.container = "{0}/{1}/{2}".format(
            projectid, subjectid, sessionid
        )
        DataSinkTissueT1.inputs.base_directory = master_config["resultdir"]
        del dsName
        baw201.connect(
            outputsSpec, "t1_average", DataSinkTissueT1, "TissueClassify.@t1"
        )

        if not onlyT1:
            dsName = "{0}_ds_{1}_tissue_t2".format(phase, sessionid)
            DataSinkTissueT2 = pe.Node(name=dsName, interface=nio.DataSink())
            DataSinkTissueT2.overwrite = master_config["ds_overwrite"]
            DataSinkTissueT2.inputs.container = "{0}/{1}/{2}".format(
                projectid, subjectid, sessionid
            )
            DataSinkTissueT2.inputs.base_directory = master_config["resultdir"]
            # DataSinkTissueT2.inputs.ignore_exception = True
            del dsName
            baw201.connect(
                outputsSpec, "t2_average", DataSinkTissueT2, "TissueClassify.@t2"
            )

        if hasPDs:
            dsName = "{0}_ds_{1}_tissue_pd".format(phase, sessionid)
            DataSinkTissuePD = pe.Node(name=dsName, interface=nio.DataSink())
            DataSinkTissuePD.overwrite = master_config["ds_overwrite"]
            DataSinkTissuePD.inputs.container = "{0}/{1}/{2}".format(
                projectid, subjectid, sessionid
            )
            DataSinkTissuePD.inputs.base_directory = master_config["resultdir"]
            # DataSinkTissuePD.inputs.ignore_exception = True
            del dsName
            baw201.connect(
                outputsSpec, "pd_average", DataSinkTissuePD, "TissueClassify.@pd"
            )

        if hasFLs:
            dsName = "{0}_ds_{1}_tissue_fl".format(phase, sessionid)
            DataSinkTissueFL = pe.Node(name=dsName, interface=nio.DataSink())
            DataSinkTissueFL.overwrite = master_config["ds_overwrite"]
            DataSinkTissueFL.inputs.container = "{0}/{1}/{2}".format(
                projectid, subjectid, sessionid
            )
            DataSinkTissueFL.inputs.base_directory = master_config["resultdir"]
            # DataSinkTissueFL.inputs.ignore_exception = True
            del dsName
            baw201.connect(
                outputsSpec, "fl_average", DataSinkTissueFL, "TissueClassify.@fl"
            )

        currentFixWMPartitioningName = "_".join(
            ["fix_wm_partitioning", str(subjectid), str(sessionid)]
        )
        FixWMNode = pe.Node(
            interface=Function(
                function=fix_wm_partitioning,
                input_names=["brainMask", "PosteriorsList"],
                output_names=[
                    "UpdatedPosteriorsList",
                    "MatchingFGCodeList",
                    "MatchingLabelList",
                    "nonAirRegionMask",
                ],
            ),
            name=currentFixWMPartitioningName,
        )

        baw201.connect(
            [
                (
                    myLocalTCWF,
                    FixWMNode,
                    [
                        ("outputspec.outputLabels", "brainMask"),
                        (
                            (
                                "outputspec.posteriorImages",
                                unwrap_posterior_images_from_list_tuples_function,
                            ),
                            "PosteriorsList",
                        ),
                    ],
                ),
                (
                    FixWMNode,
                    outputsSpec,
                    [("UpdatedPosteriorsList", "UpdatedPosteriorsList")],
                ),
            ]
        )

        currentBRAINSCreateLabelMapName = (
            "BRAINSCreateLabelMapFromProbabilityMaps_"
            + str(subjectid)
            + "_"
            + str(sessionid)
        )
        BRAINSCreateLabelMapNode = pe.Node(
            interface=BRAINSCreateLabelMapFromProbabilityMaps(),
            name=currentBRAINSCreateLabelMapName,
        )

        ## INFO:  Fix the file names
        BRAINSCreateLabelMapNode.inputs.dirtyLabelVolume = "fixed_headlabels_seg.nii.gz"
        BRAINSCreateLabelMapNode.inputs.cleanLabelVolume = (
            "fixed_brainlabels_seg.nii.gz"
        )

        baw201.connect(
            [
                (
                    FixWMNode,
                    BRAINSCreateLabelMapNode,
                    [
                        ("UpdatedPosteriorsList", "inputProbabilityVolume"),
                        ("MatchingFGCodeList", "foregroundPriors"),
                        ("MatchingLabelList", "priorLabelCodes"),
                        ("nonAirRegionMask", "nonAirRegionMask"),
                    ],
                )
            ]
        )
        baw201.connect(
            [
                (
                    BRAINSCreateLabelMapNode,
                    DataSinkTissue,
                    [  # brainstem code below replaces this ('cleanLabelVolume', 'TissueClassify.@outputLabels'),
                        ("dirtyLabelVolume", "TissueClassify.@outputHeadLabels")
                    ],
                )
            ]
        )
        baw201.connect(
            [
                (
                    myLocalTCWF,
                    DataSinkTissue,
                    [
                        (
                            "outputspec.atlasToSubjectTransform",
                            "TissueClassify.@atlas2session_tx",
                        ),
                        (
                            "outputspec.atlasToSubjectInverseTransform",
                            "TissueClassify.@atlas2sessionInverse_tx",
                        ),
                    ],
                )
            ]
        )
        baw201.connect(
            [
                (
                    FixWMNode,
                    DataSinkTissue,
                    [("UpdatedPosteriorsList", "TissueClassify.@posteriors")],
                )
            ]
        )

        currentAccumulateLikeTissuePosteriorsName = (
            "AccumulateLikeTissuePosteriors_" + str(subjectid) + "_" + str(sessionid)
        )
        AccumulateLikeTissuePosteriorsNode = pe.Node(
            interface=Function(
                function=accumulate_like_tissue_posteriors,
                input_names=["posteriorImages"],
                output_names=["AccumulatePriorsList", "AccumulatePriorsNames"],
            ),
            name=currentAccumulateLikeTissuePosteriorsName,
        )

        baw201.connect(
            [
                (
                    FixWMNode,
                    AccumulateLikeTissuePosteriorsNode,
                    [("UpdatedPosteriorsList", "posteriorImages")],
                ),
                (
                    AccumulateLikeTissuePosteriorsNode,
                    DataSinkTissue,
                    [
                        (
                            "AccumulatePriorsList",
                            "ACCUMULATED_POSTERIORS.@AccumulateLikeTissuePosteriorsOutputDir",
                        )
                    ],
                ),
            ]
        )

        """
        brain stem adds on feature
        inputs:
            - landmark (fcsv) file
            - fixed brainlabels seg.nii.gz
        output:
            - complete_brainlabels_seg.nii.gz Segmentation
        """
        myLocalBrainStemWF = create_brainstem_workflow(
            "BrainStem", master_config["queue"], "complete_brainlabels_seg.nii.gz"
        )

        baw201.connect(
            [
                (
                    myLocalLMIWF,
                    myLocalBrainStemWF,
                    [
                        (
                            "outputspec.outputLandmarksInACPCAlignedSpace",
                            "inputspec.inputLandmarkFilename",
                        )
                    ],
                ),
                (
                    BRAINSCreateLabelMapNode,
                    myLocalBrainStemWF,
                    [("cleanLabelVolume", "inputspec.inputTissueLabelFilename")],
                ),
            ]
        )

        baw201.connect(
            myLocalBrainStemWF,
            "outputspec.ouputTissuelLabelFilename",
            DataSinkTissue,
            "TissueClassify.@complete_brainlabels_seg",
        )

    dsName = "{0}_ds_{1}_seg".format(phase, sessionid)
    DataSinkSegmentation = pe.Node(name=dsName, interface=nio.DataSink())
    DataSinkSegmentation.overwrite = master_config["ds_overwrite"]
    DataSinkSegmentation.inputs.container = "{0}/{1}/{2}".format(
        projectid, subjectid, sessionid
    )
    DataSinkSegmentation.inputs.base_directory = master_config["resultdir"]
    del dsName

    ###########################
    assert (
        "segmentation" not in master_config["components"]
    ), "ERROR, segmentation (aka BRAINSCut) no longer supported"

    if "warp_atlas_to_subject" in master_config["components"]:
        ##
        ##~/src/NEP-build/bin/BRAINSResample
        # --warpTransform AtlasToSubjectPreBABC_Composite.h5
        #  --inputVolume  /Shared/sinapse/CACHE/x20141001_KIDTEST_base_CACHE/Atlas/hncma-atlas.nii.gz
        #  --referenceVolume  /Shared/sinapse/CACHE/x20141001_KIDTEST_base_CACHE/singleSession_KID1_KT1/LandmarkInitialize/BROIAuto_cropped/Cropped_BCD_ACPC_Aligned.nii.gz
        # !--outputVolume hncma.nii.gz
        # !--interpolationMode NearestNeighbor
        # !--pixelType short
        ##
        ##

        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        BResample = OrderedDict()
        AtlasLabelMapsToResample = [
            "hncma_atlas",
            "template_WMPM2_labels",
            "template_nac_labels",
        ]

        for atlasImage in AtlasLabelMapsToResample:
            BResample[atlasImage] = pe.Node(
                interface=BRAINSResample(), name="BRAINSResample_" + atlasImage
            )
            BResample[atlasImage].plugin_args = {
                "qsub_args": modify_qsub_args(master_config["queue"], 1, 1, 1),
                "overwrite": True,
            }
            BResample[atlasImage].inputs.pixelType = "short"
            BResample[atlasImage].inputs.interpolationMode = "NearestNeighbor"
            BResample[atlasImage].inputs.outputVolume = atlasImage + ".nii.gz"

            baw201.connect(
                myLocalTCWF,
                "outputspec.t1_average",
                BResample[atlasImage],
                "referenceVolume",
            )
            baw201.connect(inputsSpec, atlasImage, BResample[atlasImage], "inputVolume")
            baw201.connect(
                myLocalTCWF,
                "outputspec.atlasToSubjectTransform",
                BResample[atlasImage],
                "warpTransform",
            )
            baw201.connect(
                BResample[atlasImage],
                "outputVolume",
                DataSinkSegmentation,
                "WarpedAtlas2Subject.@" + atlasImage,
            )

        AtlasBinaryMapsToResample = [
            "template_rightHemisphere",
            "template_leftHemisphere",
            "template_ventricles",
            "template_headregion",
        ]

        for atlasImage in AtlasBinaryMapsToResample:
            BResample[atlasImage] = pe.Node(
                interface=BRAINSResample(), name="BRAINSResample_" + atlasImage
            )
            BResample[atlasImage].plugin_args = {
                "qsub_args": modify_qsub_args(master_config["queue"], 1, 1, 1),
                "overwrite": True,
            }
            BResample[atlasImage].inputs.pixelType = "binary"
            BResample[atlasImage].inputs.interpolationMode = (
                "Linear"  ## Conversion to distance map, so use linear to resample distance map
            )
            BResample[atlasImage].inputs.outputVolume = atlasImage + ".nii.gz"

            baw201.connect(
                myLocalTCWF,
                "outputspec.t1_average",
                BResample[atlasImage],
                "referenceVolume",
            )
            baw201.connect(inputsSpec, atlasImage, BResample[atlasImage], "inputVolume")
            baw201.connect(
                myLocalTCWF,
                "outputspec.atlasToSubjectTransform",
                BResample[atlasImage],
                "warpTransform",
            )
            baw201.connect(
                BResample[atlasImage],
                "outputVolume",
                DataSinkSegmentation,
                "WarpedAtlas2Subject.@" + atlasImage,
            )

        WhiteMatterHemisphereNode = pe.Node(
            interface=Function(
                function=create_left_right_wm_hemispheres,
                input_names=[
                    "BRAINLABELSFile",
                    "HDCMARegisteredVentricleMaskFN",
                    "LeftHemisphereMaskName",
                    "RightHemisphereMaskName",
                    "WM_LeftHemisphereFileName",
                    "WM_RightHemisphereFileName",
                ],
                output_names=[
                    "WM_LeftHemisphereFileName",
                    "WM_RightHemisphereFileName",
                ],
            ),
            name="WhiteMatterHemisphere",
        )
        WhiteMatterHemisphereNode.inputs.WM_LeftHemisphereFileName = (
            "left_hemisphere_wm.nii.gz"
        )
        WhiteMatterHemisphereNode.inputs.WM_RightHemisphereFileName = (
            "right_hemisphere_wm.nii.gz"
        )

        baw201.connect(
            myLocalBrainStemWF,
            "outputspec.ouputTissuelLabelFilename",
            WhiteMatterHemisphereNode,
            "BRAINLABELSFile",
        )
        baw201.connect(
            BResample["hncma_atlas"],
            "outputVolume",
            WhiteMatterHemisphereNode,
            "HDCMARegisteredVentricleMaskFN",
        )
        baw201.connect(
            BResample["template_leftHemisphere"],
            "outputVolume",
            WhiteMatterHemisphereNode,
            "LeftHemisphereMaskName",
        )
        baw201.connect(
            BResample["template_rightHemisphere"],
            "outputVolume",
            WhiteMatterHemisphereNode,
            "RightHemisphereMaskName",
        )

        baw201.connect(
            WhiteMatterHemisphereNode,
            "WM_LeftHemisphereFileName",
            DataSinkSegmentation,
            "WarpedAtlas2Subject.@LeftHemisphereWM",
        )
        baw201.connect(
            WhiteMatterHemisphereNode,
            "WM_RightHemisphereFileName",
            DataSinkSegmentation,
            "WarpedAtlas2Subject.@RightHemisphereWM",
        )

    if (
        "jointfusion_2015_wholebrain" in master_config["components"]
    ):  ## HACK Do JointFusion labeling
        ## HACK FOR NOW SHOULD BE MORE ELEGANT FROM THE .config file
        if badT2:
            onlyT1 = True
        if onlyT1:
            print("T1 only processing in jointFusion")
        else:
            print("Multimodal processing in jointFusion")

        myLocalJointFusion = create_joint_fusion_workflow(
            "JointFusion", onlyT1, master_config
        )
        baw201.connect(
            myLocalTCWF,
            "outputspec.t1_average",
            myLocalJointFusion,
            "inputspec.subj_t1_image",
        )
        baw201.connect(
            myLocalTCWF,
            "outputspec.t2_average",
            myLocalJointFusion,
            "inputspec.subj_t2_image",
        )
        baw201.connect(
            myLocalBrainStemWF,
            "outputspec.ouputTissuelLabelFilename",
            myLocalJointFusion,
            "inputspec.subj_fixed_head_labels",
        )
        baw201.connect(
            myLocalTCWF,
            "outputspec.posteriorImages",
            myLocalJointFusion,
            "inputspec.subj_posteriors",
        )

        baw201.connect(
            BResample["template_leftHemisphere"],
            "outputVolume",
            myLocalJointFusion,
            "inputspec.subj_left_hemisphere",
        )
        baw201.connect(
            myLocalLMIWF,
            "outputspec.outputLandmarksInACPCAlignedSpace",
            myLocalJointFusion,
            "inputspec.subj_lmks",
        )
        baw201.connect(
            atlasBCDNode_S,
            "template_weights_50Lmks_wts",
            myLocalJointFusion,
            "inputspec.atlasWeightFilename",
        )

        inputLabelFileJointFusionnameSpec = pe.Node(
            interface=IdentityInterface(fields=["labelBaseFilename"]),
            run_without_submitting=True,
            name="inputLabelFileJointFusionnameSpec",
        )
        baw201.connect(
            inputLabelFileJointFusionnameSpec,
            "labelBaseFilename",
            myLocalJointFusion,
            "inputspec.labelBaseFilename",
        )

        # baw201.connect(myLocalJointFusion,'outputspec.JointFusion_HDAtlas20_2015_label',DataSinkSegmentation,'JointFusion.@JointFusion_HDAtlas20_2015_label')
        # baw201.connect(myLocalJointFusion,'outputspec.JointFusion_HDAtlas20_2015_CSFVBInjected_label',DataSinkSegmentation,'JointFusion.@JointFusion_HDAtlas20_2015_CSFVBInjected_label')
        baw201.connect(
            myLocalJointFusion,
            "outputspec.JointFusion_HDAtlas20_2015_fs_standard_label",
            DataSinkSegmentation,
            "JointFusion.@JointFusion_HDAtlas20_2015_fs_standard_label",
        )
        baw201.connect(
            myLocalJointFusion,
            "outputspec.JointFusion_HDAtlas20_2015_lobe_label",
            DataSinkSegmentation,
            "JointFusion.@JointFusion_HDAtlas20_2015_lobe_label",
        )
        # baw201.connect(myLocalJointFusion,'outputspec.JointFusion_extended_snapshot',DataSinkSegmentation,'JointFusion.@JointFusion_extended_snapshot')
        baw201.connect(
            myLocalJointFusion,
            "outputspec.JointFusion_HDAtlas20_2015_dustCleaned_label",
            DataSinkSegmentation,
            "JointFusion.@JointFusion_HDAtlas20_2015_dustCleaned_label",
        )

        baw201.connect(
            myLocalJointFusion,
            "outputspec.JointFusion_volumes_csv",
            DataSinkSegmentation,
            "JointFusion.allVol.@JointFusion_volumesCSV",
        )
        baw201.connect(
            myLocalJointFusion,
            "outputspec.JointFusion_volumes_json",
            DataSinkSegmentation,
            "JointFusion.allVol.@JointFusion_volumesJSON",
        )
        baw201.connect(
            myLocalJointFusion,
            "outputspec.JointFusion_lobe_volumes_csv",
            DataSinkSegmentation,
            "JointFusion.lobeVol.@JointFusion_lobe_volumesCSV",
        )
        baw201.connect(
            myLocalJointFusion,
            "outputspec.JointFusion_lobe_volumes_json",
            DataSinkSegmentation,
            "JointFusion.lobeVol.@JointFusion_lobe_volumesJSON",
        )

        if "logismosb" in master_config["components"]:
            print("LOGISMOSB Workflow")
            from logismosb import create_logb_workflow

            # connect LOGISMOSB inputs
            myLocalLOGISMOSBWF = create_logb_workflow(
                master_config=master_config,
                plugin_args={
                    "qsub_args": modify_qsub_args(
                        queue=master_config["queue"],
                        memoryGB=8,
                        minThreads=16,
                        maxThreads=16,
                    ),
                    "overwrite": True,
                },
            )

            # In order to allow comparisons with FreeSurfer, the input files to LOGISMOSB need to be resampled
            resample_logb_inputs = True
            if resample_logb_inputs:

                def select_first_file(file_list):
                    first_file = file_list[0]
                    return first_file

                select_t1_node = pe.Node(
                    Function(["file_list"], ["first_file"], select_first_file),
                    "SelectSingleT14LOGB",
                )
                baw201.connect(inputsSpec, "T1s", select_t1_node, "file_list")

                resample_brain_labels = pe.Node(
                    BRAINSResample(), name="ResampleBrainLabels4LOGB"
                )

                baw201.connect(
                    myLocalBrainStemWF,
                    "outputspec.ouputTissuelLabelFilename",
                    resample_brain_labels,
                    "inputVolume",
                )
                baw201.connect(
                    myLocalLMIWF,
                    "outputspec.outputTransform",
                    resample_brain_labels,
                    "warpTransform",
                )
                baw201.connect(
                    select_t1_node,
                    "first_file",
                    resample_brain_labels,
                    "referenceVolume",
                )
                resample_brain_labels.inputs.outputVolume = (
                    "brain_labels_in_original_space.nii.gz"
                )
                resample_brain_labels.inputs.pixelType = "short"
                resample_brain_labels.inputs.interpolationMode = "NearestNeighbor"
                resample_brain_labels.inputs.inverseTransform = True
                baw201.connect(
                    resample_brain_labels,
                    "outputVolume",
                    myLocalLOGISMOSBWF,
                    "inputspec.brainlabels_file",
                )

                resample_hncma = pe.Node(BRAINSResample(), name="ResampleHNCMA4LOGB")
                baw201.connect(
                    BResample["hncma_atlas"],
                    "outputVolume",
                    resample_hncma,
                    "inputVolume",
                )
                baw201.connect(
                    myLocalLMIWF,
                    "outputspec.outputTransform",
                    resample_hncma,
                    "warpTransform",
                )
                baw201.connect(
                    select_t1_node, "first_file", resample_hncma, "referenceVolume"
                )
                resample_hncma.inputs.outputVolume = "hncma_in_original_space.nii.gz"
                resample_hncma.inputs.pixelType = "short"
                resample_hncma.inputs.interpolationMode = "NearestNeighbor"
                resample_hncma.inputs.inverseTransform = True
                baw201.connect(
                    resample_hncma,
                    "outputVolume",
                    myLocalLOGISMOSBWF,
                    "inputspec.hncma_atlas",
                )

                for struct in ("t1", "t2"):
                    resample_struct = pe.Node(
                        BRAINSResample(), name="Resample{0}".format(struct.upper())
                    )
                    resample_struct.inputs.outputVolume = (
                        "{0}_in_original_space.nii.gz".format(struct)
                    )
                    resample_struct.inputs.pixelType = "short"
                    resample_struct.inputs.interpolationMode = "Linear"
                    resample_struct.inputs.inverseTransform = True
                    baw201.connect(
                        select_t1_node, "first_file", resample_struct, "referenceVolume"
                    )
                    baw201.connect(
                        myLocalLMIWF,
                        "outputspec.outputTransform",
                        resample_struct,
                        "warpTransform",
                    )
                    baw201.connect(
                        myLocalTCWF,
                        "outputspec.{0}_average".format(struct),
                        resample_struct,
                        "inputVolume",
                    )
                    baw201.connect(
                        resample_struct,
                        "outputVolume",
                        myLocalLOGISMOSBWF,
                        "inputspec.{0}_file".format(struct),
                    )

                resample_joint_fusion = pe.Node(BRAINSResample(), "ResampleJointFusion")
                resample_joint_fusion.inputs.outputVolume = "malf_resampled.nii.gz"
                resample_joint_fusion.inputs.pixelType = "short"
                resample_joint_fusion.inputs.interpolationMode = "NearestNeighbor"
                resample_joint_fusion.inputs.inverseTransform = True
                baw201.connect(
                    select_t1_node,
                    "first_file",
                    resample_joint_fusion,
                    "referenceVolume",
                )
                baw201.connect(
                    myLocalLMIWF,
                    "outputspec.outputTransform",
                    resample_joint_fusion,
                    "warpTransform",
                )
                baw201.connect(
                    myLocalJointFusion,
                    "outputspec.JointFusion_HDAtlas20_2015_fs_standard_label",
                    resample_joint_fusion,
                    "inputVolume",
                )
                baw201.connect(
                    resample_joint_fusion,
                    "outputVolume",
                    myLocalLOGISMOSBWF,
                    "inputspec.joint_fusion_file",
                )

                def select_posterior(posterior_files, key="CSF"):
                    out_file = posterior_files[key]
                    return out_file

                def to_dict(string, key="CSF"):
                    out_dict = {key: string}
                    return out_dict

                select_csf = pe.Node(
                    Function(["posterior_files"], ["out_file"], select_posterior),
                    "SelectCSF",
                )
                baw201.connect(
                    myLocalTCWF,
                    "outputspec.posteriorImages",
                    select_csf,
                    "posterior_files",
                )

                resample_CSF = pe.Node(BRAINSResample(), "ResampleCSF")
                resample_CSF.inputs.outputVolume = "CSF_resampled.nii.gz"
                resample_CSF.inputs.pixelType = "float"
                resample_CSF.inputs.interpolationMode = "Linear"
                resample_CSF.inputs.inverseTransform = True
                baw201.connect(
                    select_t1_node, "first_file", resample_CSF, "referenceVolume"
                )
                baw201.connect(
                    myLocalLMIWF,
                    "outputspec.outputTransform",
                    resample_CSF,
                    "warpTransform",
                )
                baw201.connect(select_csf, "out_file", resample_CSF, "inputVolume")

                csf_to_dict = pe.Node(
                    Function(["string"], ["out_dict"], to_dict), "CSFDict"
                )
                baw201.connect(resample_CSF, "outputVolume", csf_to_dict, "string")
                baw201.connect(
                    csf_to_dict,
                    "out_dict",
                    myLocalLOGISMOSBWF,
                    "inputspec.posterior_files",
                )

            else:
                baw201.connect(
                    myLocalTCWF,
                    "outputspec.t1_average",
                    myLocalLOGISMOSBWF,
                    "inputspec.t1_file",
                )
                baw201.connect(
                    myLocalTCWF,
                    "outputspec.t2_average",
                    myLocalLOGISMOSBWF,
                    "inputspec.t2_file",
                )
                baw201.connect(
                    myLocalJointFusion,
                    "outputspec.JointFusion_HDAtlas20_2015_fs_standard_label",
                    myLocalLOGISMOSBWF,
                    "inputspec.joint_fusion_file",
                )
                baw201.connect(
                    BResample["hncma_atlas"],
                    "outputVolume",
                    myLocalLOGISMOSBWF,
                    "inputspec.hncma_atlas",
                )
                baw201.connect(
                    myLocalBrainStemWF,
                    "outputspec.ouputTissuelLabelFilename",
                    myLocalLOGISMOSBWF,
                    "inputspec.brainlabels_file",
                )
                baw201.connect(
                    myLocalTCWF,
                    "outputspec.posteriorImages",
                    myLocalLOGISMOSBWF,
                    "inputspec.posterior_files",
                )

            # connect LOGISMOSB outputs to the data sink
            baw201.connect(
                myLocalLOGISMOSBWF,
                "outputspec.lh_gmsurface_file",
                DataSink,
                "LOGISMOSB.@lh_gm_surf",
            )
            baw201.connect(
                myLocalLOGISMOSBWF,
                "outputspec.lh_wmsurface_file",
                DataSink,
                "LOGISMOSB.@lh_wm_surf",
            )
            baw201.connect(
                myLocalLOGISMOSBWF,
                "outputspec.rh_wmsurface_file",
                DataSink,
                "LOGISMOSB.@rh_wm_surf",
            )
            baw201.connect(
                myLocalLOGISMOSBWF,
                "outputspec.rh_gmsurface_file",
                DataSink,
                "LOGISMOSB.@rh_gm_surf",
            )

    if "fs_nipype" in master_config["components"]:
        from nipype.workflows.smri.freesurfer import create_reconall_workflow

        num_threads = 12
        # HACK to convert subject_dir to supported string type
        old_str = type("")
        subject_dir = old_str(
            os.path.join(master_config["resultdir"], projectid, subjectid, sessionid)
        )

        # HACK: Select first T2 from the list of input T2s
        def select_first_file(file_list):
            first_file = file_list[0]
            return first_file

        select_t2_node = pe.Node(
            Function(["file_list"], ["first_file"], select_first_file), "SelectSingleT2"
        )
        baw201.connect(inputsSpec, "T2s", select_t2_node, "file_list")

        reconall = create_reconall_workflow(
            plugin_args={
                "qsub_args": modify_qsub_args(
                    queue=master_config["queue"],
                    memoryGB=8,
                    minThreads=num_threads,
                    maxThreads=num_threads,
                ),
                "overwrite": True,
            }
        )

        baw201.connect(
            [
                (inputsSpec, reconall, [("T1s", "inputspec.T1_files")]),
                (select_t2_node, reconall, [("first_file", "inputspec.T2_file")]),
            ]
        )
        if not os.path.exists(subject_dir):
            os.makedirs(subject_dir)
        reconall.inputs.inputspec.subjects_dir = subject_dir
        reconall.inputs.inputspec.num_threads = num_threads
        reconall.inputs.inputspec.subject_id = "FreeSurfer"

    return baw201
