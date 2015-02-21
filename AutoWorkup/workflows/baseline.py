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

import os
#"""Import necessary modules from nipype."""
# from nipype.utils.config import config
# config.set('logging', 'log_to_file', 'false')
# config.set_log_dir(os.getcwd())
#--config.set('logging', 'workflow_level', 'DEBUG')
#--config.set('logging', 'interface_level', 'DEBUG')
#--config.set('execution','remove_unnecessary_outputs','false')

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio

from nipype.interfaces.utility import IdentityInterface, Function

from nipype.utils.misc import package_check
# package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
package_check('numpy', '1.3', 'tutorial1')
package_check('scipy', '0.7', 'tutorial1')
# THIS IS NOT REQUIRED package_check('matplotlib','1.4','turorial1')
package_check('networkx', '1.0', 'tutorial1')
package_check('IPython', '0.10', 'tutorial1')

from utilities.distributed import modify_qsub_args
from PipeLineFunctionHelpers import convertToList, FixWMPartitioning, AccumulateLikeTissuePosteriors
from PipeLineFunctionHelpers import UnwrapPosteriorImagesFromDictionaryFunction as flattenDict

from WorkupT1T2LandmarkInitialization import CreateLandmarkInitializeWorkflow
from WorkupT1T2TissueClassify import CreateTissueClassifyWorkflow
from WorkupT1T2MALF import CreateMALFWorkflow
from WorkupAddsonBrainStem import CreateBrainstemWorkflow

from utilities.misc import *

try:
    from SEMTools import *
except ImportError:
    from AutoWorkup.SEMTools import *

from SEMTools.registration.brainsresample import BRAINSResample

from SEMTools.filtering.denoising import UnbiasedNonLocalMeans
from SEMTools.segmentation.specialized import BRAINSCreateLabelMapFromProbabilityMaps


def get_list_element(nestedList, index):
    return nestedList[index]


def DetermineIfSegmentationShouldBeDone(master_config):
    """ This function is in a trival state right now, but
    more complicated rulesets may be necessary in the furture
    to determine when segmentation should be run.
    This is being left so that anticipated future
    changes are easier to implement.
    """
    do_BRAINSCut_Segmentation = False
    if master_config['workflow_phase'] == 'atlas-based-reference':
        if 'segmentation' in master_config['components']:
            do_BRAINSCut_Segmentation = True
    elif master_config['workflow_phase'] == 'subject-based-reference':
        if 'segmentation' in master_config['components']:
            do_BRAINSCut_Segmentation = True
    return do_BRAINSCut_Segmentation


def getAllT1sLength(allT1s):
    return len(allT1s)


##TODO:  Move to module that can be re-used
##       GetLargestLabel is copied elsewhere
def CreateLeftRightWMHemispheres(BRAINLABELSFile,
                                HDCMARegisteredVentricleMaskFN,
                                LeftHemisphereMaskName,
                                RightHemisphereMaskName,

                                WM_LeftHemisphereFileName,
                                WM_RightHemisphereFileName):
    import SimpleITK as sitk
    import os

    def GetLargestLabel(inputMask, UseErosionCleaning):
        LargestComponentCode = 1
        if UseErosionCleaning:
            erosionMask = sitk.ErodeObjectMorphology(inputMask, 1)
        else:
            erosionMask = inputMask
        CC = sitk.ConnectedComponent(erosionMask)
        Rlabel = sitk.RelabelComponent(CC)
        largestMask = ( Rlabel == LargestComponentCode)
        if UseErosionCleaning:
            dilateMask = sitk.DilateObjectMorphology(largestMask, 1)
        else:
            dilateMask = largestMask

        return (largestMask * dilateMask > 0)

    ABCLabelsImage = sitk.Cast(sitk.ReadImage(BRAINLABELSFile), sitk.sitkUInt32)
    # # Remove brain stem and cerebellum
    BS = (ABCLabelsImage == 30)
    Cerebellum_GM = (ABCLabelsImage == 12)
    Cerebellum_WM = (ABCLabelsImage == 11)
    KeepRegion = sitk.Cast((1 - (BS + Cerebellum_GM + Cerebellum_WM)), sitk.sitkUInt32)

    ABCLabelsImage = KeepRegion * ABCLabelsImage

    HDCMARegisteredVentricleLabels = sitk.Cast(sitk.ReadImage(HDCMARegisteredVentricleMaskFN), sitk.sitkUInt32)
    ABCCSFLabelCode = 4
    HDMCALeftVentricleCode = 4
    HDMCARightVentricleCode = 43
    HDCMAMask = ( HDCMARegisteredVentricleLabels == HDMCALeftVentricleCode ) + (
    HDCMARegisteredVentricleLabels == HDMCARightVentricleCode)
    ExpandVentValue = 5
    HDCMAMask_d5 = sitk.DilateObjectMorphology(HDCMAMask, ExpandVentValue)
    CSFMaskImage = (ABCLabelsImage == ABCCSFLabelCode)
    VentricleMask = ( ( HDCMAMask_d5 * CSFMaskImage + HDCMAMask ) > 0 )
    VentricleMask_d2 = sitk.DilateObjectMorphology(VentricleMask, 2)
    ABCWMLabelCode = 1
    WMMaskImage = (ABCLabelsImage == ABCWMLabelCode)

    subcorticalRegions = (
    ABCLabelsImage >= 12 )  # All subcortical regions are listed greater than equal to values of 12
    WMSubcortFilled = ( ( WMMaskImage + subcorticalRegions ) > 0 )

    WMSubcortFilled_CC = GetLargestLabel(WMSubcortFilled, False)
    WMSubcortFilled_CC_Ventricles = ( ( WMSubcortFilled_CC + VentricleMask_d2 ) > 0 )
    neg_WMSubcortFilled_CC = ( 1 - WMSubcortFilled_CC )
    neg_WMSubcortFilled_CC_bg = GetLargestLabel(neg_WMSubcortFilled_CC, False)
    neg_WMSubcortFilled_CC_bg_holes = (neg_WMSubcortFilled_CC - neg_WMSubcortFilled_CC_bg )

    WM_Final = (neg_WMSubcortFilled_CC_bg_holes + WMSubcortFilled_CC_Ventricles > 0 )

    ####################################
    ### START WM

    # Template masks for left and right hemispheres
    Left_template = (sitk.Cast(sitk.ReadImage(LeftHemisphereMaskName), sitk.sitkUInt32) > 0 )
    Right_template = (sitk.Cast(sitk.ReadImage(RightHemisphereMaskName), sitk.sitkUInt32) > 0 )

    # Split into left and right hemispheres
    WM_left = ( WM_Final * Left_template > 0 )
    WM_right = ( WM_Final * Right_template > 0 )

    WM_Largest_left = GetLargestLabel(WM_left, False)
    WM_Largest_right = GetLargestLabel(WM_right, False)

    WM_left_extras = WM_left - WM_Largest_left
    WM_right_extras = WM_right - WM_Largest_right

    WM_Largest_left = GetLargestLabel(WM_Largest_left + WM_right_extras, False)
    WM_Largest_right = GetLargestLabel(WM_Largest_right + WM_left_extras, False)

    ## Write todisk
    WM_LeftHemisphereFileName = os.path.abspath(WM_LeftHemisphereFileName)
    sitk.WriteImage(WM_Largest_left, WM_LeftHemisphereFileName)

    WM_RightHemisphereFileName = os.path.abspath(WM_RightHemisphereFileName)
    sitk.WriteImage(WM_Largest_right, WM_RightHemisphereFileName)

    ## TODO Add splitting into hemispheres code here
    return WM_LeftHemisphereFileName, WM_RightHemisphereFileName



def generate_single_session_template_WF(projectid, subjectid, sessionid, onlyT1, master_config, phase, interpMode,
                                        pipeline_name, doDenoise=True):
    """
    Run autoworkup on a single sessionid

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """

    #if  not 'landmark' in master_config['components'] or not 'auxlmk' in master_config['components'] or not 'tissue_classify' in master_config['components']:
    #    print "Baseline DataSink requires 'AUXLMK' and/or 'TISSUE_CLASSIFY'!!!"
    #    raise NotImplementedError
    # master_config['components'].append('auxlmk')
    # master_config['components'].append('tissue_classify')

    assert phase in ['atlas-based-reference',
                     'subject-based-reference'], "Unknown phase! Valid entries: 'atlas-based-reference', 'subject-based-reference'"

    if 'tissue_classify' in master_config['components']:
        assert ('landmark' in master_config['components'] ), "tissue_classify Requires landmark step!"
    if 'landmark' in master_config['components']:
        assert 'denoise' in master_config['components'], "landmark Requires denoise step!"

    from workflows.atlasNode import MakeAtlasNode

    baw201 = pe.Workflow(name=pipeline_name)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['atlasLandmarkFilename', 'atlasWeightFilename',
                                                             'LLSModel', 'inputTemplateModel', 'template_t1',
                                                             'atlasDefinition', 'T1s', 'T2s', 'PDs', 'FLs', 'OTHERs',
                                                             'hncma_atlas',
                                                             'template_rightHemisphere',
                                                             'template_leftHemisphere',
                                                             'template_WMPM2_labels',
                                                             'template_nac_labels',
                                                             'template_ventricles']),
                         run_without_submitting=True, name='inputspec')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['t1_average', 't2_average', 'pd_average', 'fl_average',
                                                              'posteriorImages', 'outputLabels', 'outputHeadLabels',
                                                              'atlasToSubjectTransform',
                                                              'atlasToSubjectInverseTransform',
                                                              'atlasToSubjectRegistrationState',
                                                              'BCD_ACPC_T1_CROPPED',
                                                              'outputLandmarksInACPCAlignedSpace',
                                                              'outputLandmarksInInputSpace',
                                                              'output_tx', 'LMIatlasToSubject_tx',
                                                              'writeBranded2DImage',
                                                              'brainStemMask',
                                                              'UpdatedPosteriorsList'  # Longitudinal
    ]),
                          run_without_submitting=True, name='outputspec')

    dsName = "{0}_ds_{1}".format(phase, sessionid)
    DataSink = pe.Node(name=dsName, interface=nio.DataSink())
    DataSink.overwrite = master_config['ds_overwrite']
    DataSink.inputs.container = '{0}/{1}/{2}'.format(projectid, subjectid, sessionid)
    DataSink.inputs.base_directory = master_config['resultdir']

    atlas_static_directory = master_config['atlascache']
    if master_config['workflow_phase'] == 'atlas-based-reference':
        atlas_warped_directory = master_config['atlascache']
        atlasABCNode_XML = MakeAtlasNode(atlas_warped_directory, 'BABCXMLAtlas_{0}'.format(sessionid),
                                         ['W_BRAINSABCSupport'])
        baw201.connect(atlasABCNode_XML, 'ExtendedAtlasDefinition_xml', inputsSpec, 'atlasDefinition')

        atlasABCNode_W = MakeAtlasNode(atlas_warped_directory, 'BABCAtlas_W{0}'.format(sessionid),
                                       ['W_BRAINSABCSupport', 'W_LabelMapsSupport'])
        baw201.connect([( atlasABCNode_W, inputsSpec, [
            ('hncma_atlas', 'hncma_atlas'),
            ('template_leftHemisphere', 'template_leftHemisphere'),
            ('template_rightHemisphere', 'template_rightHemisphere'),
            ('template_WMPM2_labels', 'template_WMPM2_labels'),
            ('template_nac_labels', 'template_nac_labels'),
            ('template_ventricles', 'template_ventricles')]
                        )]
        )
        ## These landmarks are only relevant for the atlas-based-reference case
        atlasBCDNode_W = MakeAtlasNode(atlas_warped_directory, 'BBCDAtlas_W{0}'.format(sessionid),
                                       ['W_BCDSupport'])
        baw201.connect([(atlasBCDNode_W, inputsSpec,
                         [('template_t1', 'template_t1'),
                          ('template_landmarks_50Lmks_fcsv', 'atlasLandmarkFilename'),
                         ]),
        ])
        ## Needed for both segmentation and template building prep
        atlasBCUTNode_W = MakeAtlasNode(atlas_warped_directory,
                                        'BBCUTAtlas_W{0}'.format(sessionid), ['W_BRAINSCutSupport'])



    elif master_config['workflow_phase'] == 'subject-based-reference':
        print master_config['previousresult']
        atlas_warped_directory = os.path.join(master_config['previousresult'], subjectid, 'Atlas')

        atlasBCUTNode_W = pe.Node(interface=nio.DataGrabber(infields=['subject'],
                                                            outfields=[
                                                                "l_accumben_ProbabilityMap",
                                                                "r_accumben_ProbabilityMap",
                                                                "l_caudate_ProbabilityMap",
                                                                "r_caudate_ProbabilityMap",
                                                                "l_globus_ProbabilityMap",
                                                                "r_globus_ProbabilityMap",
                                                                "l_hippocampus_ProbabilityMap",
                                                                "r_hippocampus_ProbabilityMap",
                                                                "l_putamen_ProbabilityMap",
                                                                "r_putamen_ProbabilityMap",
                                                                "l_thalamus_ProbabilityMap",
                                                                "r_thalamus_ProbabilityMap",
                                                                "phi",
                                                                "rho",
                                                                "theta"
                                                            ]),
                                  name='PerSubject_atlasBCUTNode_W')
        atlasBCUTNode_W.inputs.base_directory = master_config['previousresult']
        atlasBCUTNode_W.inputs.subject = subjectid
        atlasBCUTNode_W.inputs.field_template = {
            'l_accumben_ProbabilityMap': '%s/Atlas/AVG_l_accumben_ProbabilityMap.nii.gz',
            'r_accumben_ProbabilityMap': '%s/Atlas/AVG_r_accumben_ProbabilityMap.nii.gz',
            'l_caudate_ProbabilityMap': '%s/Atlas/AVG_l_caudate_ProbabilityMap.nii.gz',
            'r_caudate_ProbabilityMap': '%s/Atlas/AVG_r_caudate_ProbabilityMap.nii.gz',
            'l_globus_ProbabilityMap': '%s/Atlas/AVG_l_globus_ProbabilityMap.nii.gz',
            'r_globus_ProbabilityMap': '%s/Atlas/AVG_r_globus_ProbabilityMap.nii.gz',
            'l_hippocampus_ProbabilityMap': '%s/Atlas/AVG_l_hippocampus_ProbabilityMap.nii.gz',
            'r_hippocampus_ProbabilityMap': '%s/Atlas/AVG_r_hippocampus_ProbabilityMap.nii.gz',
            'l_putamen_ProbabilityMap': '%s/Atlas/AVG_l_putamen_ProbabilityMap.nii.gz',
            'r_putamen_ProbabilityMap': '%s/Atlas/AVG_r_putamen_ProbabilityMap.nii.gz',
            'l_thalamus_ProbabilityMap': '%s/Atlas/AVG_l_thalamus_ProbabilityMap.nii.gz',
            'r_thalamus_ProbabilityMap': '%s/Atlas/AVG_r_thalamus_ProbabilityMap.nii.gz',
            'phi': '%s/Atlas/AVG_phi.nii.gz',
            'rho': '%s/Atlas/AVG_rho.nii.gz',
            'theta': '%s/Atlas/AVG_theta.nii.gz'
        }
        atlasBCUTNode_W.inputs.template_args = {
            'l_accumben_ProbabilityMap': [['subject']],
            'r_accumben_ProbabilityMap': [['subject']],
            'l_caudate_ProbabilityMap': [['subject']],
            'r_caudate_ProbabilityMap': [['subject']],
            'l_globus_ProbabilityMap': [['subject']],
            'r_globus_ProbabilityMap': [['subject']],
            'l_hippocampus_ProbabilityMap': [['subject']],
            'r_hippocampus_ProbabilityMap': [['subject']],
            'l_putamen_ProbabilityMap': [['subject']],
            'r_putamen_ProbabilityMap': [['subject']],
            'l_thalamus_ProbabilityMap': [['subject']],
            'r_thalamus_ProbabilityMap': [['subject']],
            'phi': [['subject']],
            'rho': [['subject']],
            'theta': [['subject']]
        }
        atlasBCUTNode_W.inputs.template = '*'
        atlasBCUTNode_W.inputs.sort_filelist = True
        atlasBCUTNode_W.inputs.raise_on_empty = True

        template_DG = pe.Node(interface=nio.DataGrabber(infields=['subject'],
                                                        outfields=['outAtlasXMLFullPath',
                                                                   'hncma_atlas',
                                                                   'template_leftHemisphere',
                                                                   'template_rightHemisphere',
                                                                   'template_WMPM2_labels',
                                                                   'template_nac_labels',
                                                                   'template_ventricles',
                                                                   'template_t1',
                                                                   'template_landmarks_50Lmks_fcsv'
                                                        ]),
                              name='Template_DG')
        template_DG.inputs.base_directory = master_config['previousresult']
        template_DG.inputs.subject = subjectid
        template_DG.inputs.field_template = {'outAtlasXMLFullPath': '%s/Atlas/AtlasDefinition_%s.xml',
                                             'hncma_atlas': '%s/Atlas/AVG_hncma_atlas.nii.gz',
                                             'template_leftHemisphere': '%s/Atlas/AVG_template_leftHemisphere.nii.gz',
                                             'template_rightHemisphere': '%s/Atlas/AVG_template_rightHemisphere.nii.gz',
                                             'template_WMPM2_labels': '%s/Atlas/AVG_template_WMPM2_labels.nii.gz',
                                             'template_nac_labels': '%s/Atlas/AVG_template_nac_labels.nii.gz',
                                             'template_ventricles': '%s/Atlas/AVG_template_ventricles.nii.gz',
                                             'template_t1': '%s/Atlas/AVG_T1.nii.gz',
                                             'template_landmarks_50Lmks_fcsv': '%s/Atlas/AVG_LMKS.fcsv',
        }
        template_DG.inputs.template_args = {'outAtlasXMLFullPath': [['subject', 'subject']],
                                            'hncma_atlas': [['subject']],
                                            'template_leftHemisphere': [['subject']],
                                            'template_rightHemisphere': [['subject']],
                                            'template_WMPM2_labels': [['subject']],
                                            'template_nac_labels': [['subject']],
                                            'template_ventricles': [['subject']],
                                            'template_t1': [['subject']],
                                            'template_landmarks_50Lmks_fcsv': [['subject']]
        }
        template_DG.inputs.template = '*'
        template_DG.inputs.sort_filelist = True
        template_DG.inputs.raise_on_empty = True

        baw201.connect(template_DG, 'outAtlasXMLFullPath', inputsSpec, 'atlasDefinition')
        baw201.connect([(template_DG, inputsSpec, [
            ## Already connected ('template_t1','template_t1'),
            ('hncma_atlas', 'hncma_atlas'),
            ('template_leftHemisphere', 'template_leftHemisphere'),
            ('template_rightHemisphere', 'template_rightHemisphere'),
            ('template_WMPM2_labels', 'template_WMPM2_labels'),
            ('template_nac_labels', 'template_nac_labels'),
            ('template_ventricles', 'template_ventricles')]
                        )]
        )
        ## These landmarks are only relevant for the atlas-based-reference case
        baw201.connect([(template_DG, inputsSpec,
                         [('template_t1', 'template_t1'),
                          ('template_landmarks_50Lmks_fcsv', 'atlasLandmarkFilename'),
                         ]),
        ])

    else:
        assert 0 == 1, "Invalid workflow type specified for singleSession"

    atlasBCDNode_S = MakeAtlasNode(atlas_static_directory, 'BBCDAtlas_S{0}'.format(sessionid),
                                   ['S_BCDSupport'])
    baw201.connect([(atlasBCDNode_S, inputsSpec,
                     [('template_weights_50Lmks_wts', 'atlasWeightFilename'),
                      ('LLSModel_50Lmks_h5', 'LLSModel'),
                      ('T1_50Lmks_mdl', 'inputTemplateModel')
                     ]),
    ])

    if doDenoise:
        print("\ndenoise image filter\n")
        makeDenoiseInImageList = pe.Node(Function(function=MakeOutFileList,
                                                  input_names=['T1List', 'T2List', 'PDList', 'FLList',
                                                               'OtherList', 'postfix', 'PrimaryT1'],
                                                  output_names=['inImageList', 'outImageList', 'imageTypeList']),
                                         run_without_submitting=True, name="99_makeDenoiseInImageList")
        baw201.connect(inputsSpec, 'T1s', makeDenoiseInImageList, 'T1List')
        baw201.connect(inputsSpec, 'T2s', makeDenoiseInImageList, 'T2List')
        baw201.connect(inputsSpec, 'PDs', makeDenoiseInImageList, 'PDList')
        makeDenoiseInImageList.inputs.FLList = []  # an emptyList HACK
        makeDenoiseInImageList.inputs.PrimaryT1 = None  # an emptyList HACK
        makeDenoiseInImageList.inputs.postfix = "_UNM_denoised.nii.gz"
        # HACK baw201.connect( inputsSpec, 'FLList', makeDenoiseInImageList, 'FLList' )
        baw201.connect(inputsSpec, 'OTHERs', makeDenoiseInImageList, 'OtherList')

        print("\nDenoise:\n")
        DenoiseInputImgs = pe.MapNode(interface=UnbiasedNonLocalMeans(),
                                      name='denoiseInputImgs',
                                      iterfield=['inputVolume',
                                                 'outputVolume'])
        DenoiseInputImgs.inputs.rc = [1, 1, 1]
        DenoiseInputImgs.inputs.rs = [4, 4, 4]
        DenoiseInputImgs.plugin_args = {'qsub_args': modify_qsub_args(master_config['queue'], .2, 1, 1),
                                        'overwrite': True}
        baw201.connect([(makeDenoiseInImageList, DenoiseInputImgs, [('inImageList', 'inputVolume')]),
                        (makeDenoiseInImageList, DenoiseInputImgs, [('outImageList', 'outputVolume')])
        ])
        print("\nMerge all T1 and T2 List\n")
        makePreprocessingOutList = pe.Node(Function(function=GenerateSeparateImageTypeList,
                                                    input_names=['inFileList', 'inTypeList'],
                                                    output_names=['T1s', 'T2s', 'PDs', 'FLs', 'OtherList']),
                                           run_without_submitting=True, name="99_makePreprocessingOutList")
        baw201.connect(DenoiseInputImgs, 'outputVolume', makePreprocessingOutList, 'inFileList')
        baw201.connect(makeDenoiseInImageList, 'imageTypeList', makePreprocessingOutList, 'inTypeList')

    else:
        makePreprocessingOutList = inputsSpec

    if 'landmark' in master_config['components']:
        DoReverseMapping = False  # Set to true for debugging outputs
        if 'auxlmk' in master_config['components']:
            DoReverseMapping = True
        myLocalLMIWF = CreateLandmarkInitializeWorkflow("LandmarkInitialize", interpMode, DoReverseMapping)

        baw201.connect([(makePreprocessingOutList, myLocalLMIWF,
                         [(('T1s', get_list_element, 0), 'inputspec.inputVolume' )]),
                        (inputsSpec, myLocalLMIWF,
                         [('atlasLandmarkFilename', 'inputspec.atlasLandmarkFilename'),
                          ('atlasWeightFilename', 'inputspec.atlasWeightFilename'),
                          ('LLSModel', 'inputspec.LLSModel'),
                          ('inputTemplateModel', 'inputspec.inputTemplateModel'),
                          ('template_t1', 'inputspec.atlasVolume')]),
                        (myLocalLMIWF, outputsSpec,
                         [('outputspec.outputResampledCroppedVolume', 'BCD_ACPC_T1_CROPPED'),
                          ('outputspec.outputLandmarksInACPCAlignedSpace',
                           'outputLandmarksInACPCAlignedSpace'),
                          ('outputspec.outputLandmarksInInputSpace',
                           'outputLandmarksInInputSpace'),
                          ('outputspec.outputTransform', 'output_tx'),
                          ('outputspec.atlasToSubjectTransform', 'LMIatlasToSubject_tx'),
                          ('outputspec.writeBranded2DImage', 'writeBranded2DImage')])
        ])
        baw201.connect([(outputsSpec, DataSink,  # TODO: change to myLocalLMIWF -> DataSink
                         [('outputLandmarksInACPCAlignedSpace', 'ACPCAlign.@outputLandmarks_ACPC'),
                          ('writeBranded2DImage', 'ACPCAlign.@writeBranded2DImage'),
                          ('BCD_ACPC_T1_CROPPED', 'ACPCAlign.@BCD_ACPC_T1_CROPPED'),
                          ('outputLandmarksInInputSpace', 'ACPCAlign.@outputLandmarks_Input'),
                          ('output_tx', 'ACPCAlign.@output_tx'),
                          ('LMIatlasToSubject_tx', 'ACPCAlign.@LMIatlasToSubject_tx'), ]
                        )
        ]
        )

    if 'tissue_classify' in master_config['components']:
        useRegistrationMask = master_config['use_registration_masking']

        myLocalTCWF = CreateTissueClassifyWorkflow("TissueClassify", master_config, interpMode,useRegistrationMask)
        baw201.connect([(makePreprocessingOutList, myLocalTCWF, [('T1s', 'inputspec.T1List')]),
                        (makePreprocessingOutList, myLocalTCWF, [('T2s', 'inputspec.T2List')]),
                        (inputsSpec, myLocalTCWF, [('atlasDefinition', 'inputspec.atlasDefinition'),
                                                   ('template_t1', 'inputspec.atlasVolume'),
                                                   (('T1s', getAllT1sLength), 'inputspec.T1_count'),
                                                   ('PDs', 'inputspec.PDList'),
                                                   ('FLs', 'inputspec.FLList'),
                                                   ('OTHERs', 'inputspec.OtherList')
                        ]),
                        (myLocalLMIWF, myLocalTCWF, [('outputspec.outputResampledCroppedVolume', 'inputspec.PrimaryT1'),
                                                     ('outputspec.atlasToSubjectTransform',
                                                      'inputspec.atlasToSubjectInitialTransform')]),
                        (myLocalTCWF, outputsSpec, [('outputspec.t1_average', 't1_average'),
                                                    ('outputspec.t2_average', 't2_average'),
                                                    ('outputspec.pd_average', 'pd_average'),
                                                    ('outputspec.fl_average', 'fl_average'),
                                                    ('outputspec.posteriorImages', 'posteriorImages'),
                                                    ('outputspec.outputLabels', 'outputLabels'),
                                                    ('outputspec.outputHeadLabels', 'outputHeadLabels'),
                                                    ('outputspec.atlasToSubjectTransform', 'atlasToSubjectTransform'),
                                                    ('outputspec.atlasToSubjectInverseTransform',
                                                     'atlasToSubjectInverseTransform'),
                                                    ('outputspec.atlasToSubjectRegistrationState',
                                                     'atlasToSubjectRegistrationState')
                        ]),
        ])

        baw201.connect([(outputsSpec, DataSink,  # TODO: change to myLocalTCWF -> DataSink
                         [(('t1_average', convertToList), 'TissueClassify.@t1'),
                          (('t2_average', convertToList), 'TissueClassify.@t2'),
                          (('pd_average', convertToList), 'TissueClassify.@pd'),
                          (('fl_average', convertToList), 'TissueClassify.@fl')])
        ])

        currentFixWMPartitioningName = "_".join(['FixWMPartitioning', str(subjectid), str(sessionid)])
        FixWMNode = pe.Node(interface=Function(function=FixWMPartitioning,
                                               input_names=['brainMask', 'PosteriorsList'],
                                               output_names=['UpdatedPosteriorsList', 'MatchingFGCodeList',
                                                             'MatchingLabelList', 'nonAirRegionMask']),
                            name=currentFixWMPartitioningName)

        baw201.connect([(myLocalTCWF, FixWMNode, [('outputspec.outputLabels', 'brainMask'),
                                                  (('outputspec.posteriorImages', flattenDict), 'PosteriorsList')]),
                        (FixWMNode, outputsSpec, [('UpdatedPosteriorsList', 'UpdatedPosteriorsList')]),
        ])

        currentBRAINSCreateLabelMapName = 'BRAINSCreateLabelMapFromProbabilityMaps_' + str(subjectid) + "_" + str(
            sessionid)
        BRAINSCreateLabelMapNode = pe.Node(interface=BRAINSCreateLabelMapFromProbabilityMaps(),
                                           name=currentBRAINSCreateLabelMapName)

        ## TODO:  Fix the file names
        BRAINSCreateLabelMapNode.inputs.dirtyLabelVolume = 'fixed_headlabels_seg.nii.gz'
        BRAINSCreateLabelMapNode.inputs.cleanLabelVolume = 'fixed_brainlabels_seg.nii.gz'

        baw201.connect([(FixWMNode, BRAINSCreateLabelMapNode, [('UpdatedPosteriorsList', 'inputProbabilityVolume'),
                                                               ('MatchingFGCodeList', 'foregroundPriors'),
                                                               ('MatchingLabelList', 'priorLabelCodes'),
                                                               ('nonAirRegionMask', 'nonAirRegionMask')]),
                        (BRAINSCreateLabelMapNode, DataSink,
                         [  # brainstem code below replaces this ('cleanLabelVolume', 'TissueClassify.@outputLabels'),
                            ('dirtyLabelVolume', 'TissueClassify.@outputHeadLabels')]),
                        (myLocalTCWF, DataSink, [('outputspec.atlasToSubjectTransform',
                                                  'TissueClassify.@atlas2session_tx'),
                                                 ('outputspec.atlasToSubjectInverseTransform',
                                                  'TissueClassify.@atlas2sessionInverse_tx')]),
                        (FixWMNode, DataSink, [('UpdatedPosteriorsList', 'TissueClassify.@posteriors')]),
        ])

        currentAccumulateLikeTissuePosteriorsName = 'AccumulateLikeTissuePosteriors_' + str(subjectid) + "_" + str(
            sessionid)
        AccumulateLikeTissuePosteriorsNode = pe.Node(interface=Function(function=AccumulateLikeTissuePosteriors,
                                                                        input_names=['posteriorImages'],
                                                                        output_names=['AccumulatePriorsList',
                                                                                      'AccumulatePriorsNames']),
                                                     name=currentAccumulateLikeTissuePosteriorsName)

        baw201.connect([(FixWMNode, AccumulateLikeTissuePosteriorsNode, [('UpdatedPosteriorsList', 'posteriorImages')]),
                        (AccumulateLikeTissuePosteriorsNode, DataSink, [('AccumulatePriorsList',
                                                                         'ACCUMULATED_POSTERIORS.@AccumulateLikeTissuePosteriorsOutputDir')])])

        """
        brain stem adds on feature
        inputs:
            - landmark (fcsv) file
            - fixed brainlabels seg.nii.gz
        output:
            - complete_brainlabels_seg.nii.gz Segmentation
        """
        myLocalBrainStemWF = CreateBrainstemWorkflow("BrainStem",
                                                     master_config['queue'],
                                                     "complete_brainlabels_seg.nii.gz")

        baw201.connect([(myLocalLMIWF, myLocalBrainStemWF, [('outputspec.outputLandmarksInACPCAlignedSpace',
                                                             'inputspec.inputLandmarkFilename')]),
                        (BRAINSCreateLabelMapNode, myLocalBrainStemWF, [('cleanLabelVolume',
                                                                         'inputspec.inputTissueLabelFilename')])
        ])

        baw201.connect(myLocalBrainStemWF, 'outputspec.ouputTissuelLabelFilename', DataSink,
                       'TissueClassify.@complete_brainlabels_seg')


    ###########################
    do_BRAINSCut_Segmentation = DetermineIfSegmentationShouldBeDone(master_config)
    if do_BRAINSCut_Segmentation:
        from workflows.segmentation import segmentation
        from workflows.WorkupT1T2BRAINSCut import GenerateWFName

        sname = 'segmentation'
        segWF = segmentation(projectid, subjectid, sessionid, master_config, onlyT1, pipeline_name=sname)

        baw201.connect([(inputsSpec, segWF,
                         [
                             ('template_t1', 'inputspec.template_t1')
                         ])
        ])
        baw201.connect([(atlasBCUTNode_W, segWF,
                         [
                             ('rho', 'inputspec.rho'),
                             ('phi', 'inputspec.phi'),
                             ('theta', 'inputspec.theta'),
                             ('l_caudate_ProbabilityMap', 'inputspec.l_caudate_ProbabilityMap'),
                             ('r_caudate_ProbabilityMap', 'inputspec.r_caudate_ProbabilityMap'),
                             ('l_hippocampus_ProbabilityMap', 'inputspec.l_hippocampus_ProbabilityMap'),
                             ('r_hippocampus_ProbabilityMap', 'inputspec.r_hippocampus_ProbabilityMap'),
                             ('l_putamen_ProbabilityMap', 'inputspec.l_putamen_ProbabilityMap'),
                             ('r_putamen_ProbabilityMap', 'inputspec.r_putamen_ProbabilityMap'),
                             ('l_thalamus_ProbabilityMap', 'inputspec.l_thalamus_ProbabilityMap'),
                             ('r_thalamus_ProbabilityMap', 'inputspec.r_thalamus_ProbabilityMap'),
                             ('l_accumben_ProbabilityMap', 'inputspec.l_accumben_ProbabilityMap'),
                             ('r_accumben_ProbabilityMap', 'inputspec.r_accumben_ProbabilityMap'),
                             ('l_globus_ProbabilityMap', 'inputspec.l_globus_ProbabilityMap'),
                             ('r_globus_ProbabilityMap', 'inputspec.r_globus_ProbabilityMap')
                         ]
                        )])

        atlasBCUTNode_S = MakeAtlasNode(atlas_static_directory,
                                        'BBCUTAtlas_S{0}'.format(sessionid), ['S_BRAINSCutSupport'])
        baw201.connect(atlasBCUTNode_S, 'trainModelFile_txtD0060NT0060_gz',
                       segWF, 'inputspec.trainModelFile_txtD0060NT0060_gz')

        ## baw201_outputspec = baw201.get_node('outputspec')
        baw201.connect([(myLocalTCWF, segWF, [('outputspec.t1_average', 'inputspec.t1_average'),
                                              ('outputspec.atlasToSubjectRegistrationState',
                                               'inputspec.atlasToSubjectRegistrationState'),
                                              ('outputspec.outputLabels', 'inputspec.inputLabels'),
                                              ('outputspec.posteriorImages', 'inputspec.posteriorImages'),
                                              ('outputspec.outputHeadLabels', 'inputspec.inputHeadLabels')
        ] ),
                        (myLocalLMIWF, segWF, [('outputspec.atlasToSubjectTransform', 'inputspec.LMIatlasToSubject_tx')
                        ] ),
                        (FixWMNode, segWF, [('UpdatedPosteriorsList', 'inputspec.UpdatedPosteriorsList')
                        ] ),
        ])
        if not onlyT1:
            baw201.connect([(myLocalTCWF, segWF, [('outputspec.t2_average', 'inputspec.t2_average')])])

    if 'warp_atlas_to_subject' in master_config['components']:
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

        ## TODO : SHOULD USE BRAINSCut transform that was refined even further!

        BResample = dict()
        AtlasLabelMapsToResample = [
            'hncma_atlas',
            'template_WMPM2_labels',
            'template_nac_labels',
        ]


        for atlasImage in AtlasLabelMapsToResample:
            BResample[atlasImage] = pe.Node(interface=BRAINSResample(), name="BRAINSResample_" + atlasImage)
            BResample[atlasImage].plugin_args = {'qsub_args': modify_qsub_args(master_config['queue'], 1, 1, 1),
                                                 'overwrite': True}
            BResample[atlasImage].inputs.pixelType = 'short'
            BResample[atlasImage].inputs.interpolationMode = 'NearestNeighbor'
            BResample[atlasImage].inputs.outputVolume = atlasImage + ".nii.gz"

            baw201.connect(myLocalTCWF, 'outputspec.t1_average', BResample[atlasImage], 'referenceVolume')
            baw201.connect(inputsSpec, atlasImage, BResample[atlasImage], 'inputVolume')
            baw201.connect(myLocalTCWF, 'outputspec.atlasToSubjectTransform',
                           BResample[atlasImage], 'warpTransform')
            baw201.connect(BResample[atlasImage], 'outputVolume', DataSink, 'WarpedAtlas2Subject.@' + atlasImage)

        AtlasBinaryMapsToResample = [
            'template_rightHemisphere',
            'template_leftHemisphere',
            'template_ventricles']

        for atlasImage in AtlasBinaryMapsToResample:
            BResample[atlasImage] = pe.Node(interface=BRAINSResample(), name="BRAINSResample_" + atlasImage)
            BResample[atlasImage].plugin_args = {'qsub_args': modify_qsub_args(master_config['queue'], 1, 1, 1),
                                                 'overwrite': True}
            BResample[atlasImage].inputs.pixelType = 'binary'
            BResample[
                atlasImage].inputs.interpolationMode = 'Linear'  ## Conversion to distance map, so use linear to resample distance map
            BResample[atlasImage].inputs.outputVolume = atlasImage + ".nii.gz"

            baw201.connect(myLocalTCWF, 'outputspec.t1_average', BResample[atlasImage], 'referenceVolume')
            baw201.connect(inputsSpec, atlasImage, BResample[atlasImage], 'inputVolume')
            baw201.connect(myLocalTCWF, 'outputspec.atlasToSubjectTransform', BResample[atlasImage], 'warpTransform')
            baw201.connect(BResample[atlasImage], 'outputVolume', DataSink, 'WarpedAtlas2Subject.@' + atlasImage)

        BRAINSCutAtlasImages = [
            'rho',
            'phi',
            'theta',
            'l_caudate_ProbabilityMap',
            'r_caudate_ProbabilityMap',
            'l_hippocampus_ProbabilityMap',
            'r_hippocampus_ProbabilityMap',
            'l_putamen_ProbabilityMap',
            'r_putamen_ProbabilityMap',
            'l_thalamus_ProbabilityMap',
            'r_thalamus_ProbabilityMap',
            'l_accumben_ProbabilityMap',
            'r_accumben_ProbabilityMap',
            'l_globus_ProbabilityMap',
            'r_globus_ProbabilityMap'
        ]
        for atlasImage in BRAINSCutAtlasImages:
            BResample[atlasImage] = pe.Node(interface=BRAINSResample(), name="BCUTBRAINSResample_" + atlasImage)
            BResample[atlasImage].plugin_args = {'qsub_args': modify_qsub_args(master_config['queue'], 1, 1, 1),
                                                 'overwrite': True}
            BResample[atlasImage].inputs.pixelType = 'float'
            BResample[
                atlasImage].inputs.interpolationMode = 'Linear'  ## Conversion to distance map, so use linear to resample distance map
            BResample[atlasImage].inputs.outputVolume = atlasImage + ".nii.gz"

            baw201.connect(myLocalTCWF, 'outputspec.t1_average', BResample[atlasImage], 'referenceVolume')
            baw201.connect(atlasBCUTNode_W, atlasImage, BResample[atlasImage], 'inputVolume')
            baw201.connect(myLocalTCWF, 'outputspec.atlasToSubjectTransform', BResample[atlasImage], 'warpTransform')
            baw201.connect(BResample[atlasImage], 'outputVolume', DataSink, 'WarpedAtlas2Subject.@' + atlasImage)


        WhiteMatterHemisphereNode = pe.Node(interface=Function(function=CreateLeftRightWMHemispheres,
                                                                        input_names=['BRAINLABELSFile',
                                                                                     'HDCMARegisteredVentricleMaskFN',
                                                                                     'LeftHemisphereMaskName',
                                                                                     'RightHemisphereMaskName',
                                                                                     'WM_LeftHemisphereFileName',
                                                                                     'WM_RightHemisphereFileName'],
                                                                        output_names=['WM_LeftHemisphereFileName',
                                                                                      'WM_RightHemisphereFileName']),
                                                     name="WhiteMatterHemisphere")
        WhiteMatterHemisphereNode.inputs.WM_LeftHemisphereFileName ="left_hemisphere_wm.nii.gz"
        WhiteMatterHemisphereNode.inputs.WM_RightHemisphereFileName ="right_hemisphere_wm.nii.gz"

        baw201.connect(myLocalBrainStemWF,'outputspec.ouputTissuelLabelFilename',WhiteMatterHemisphereNode,'BRAINLABELSFile')
        baw201.connect(BResample['hncma_atlas'],'outputVolume',WhiteMatterHemisphereNode,'HDCMARegisteredVentricleMaskFN')
        baw201.connect(BResample['template_leftHemisphere'],'outputVolume',WhiteMatterHemisphereNode,'LeftHemisphereMaskName')
        baw201.connect(BResample['template_rightHemisphere'],'outputVolume',WhiteMatterHemisphereNode,'RightHemisphereMaskName')

        baw201.connect(WhiteMatterHemisphereNode,'WM_LeftHemisphereFileName',DataSink,'WarpedAtlas2Subject.@LeftHemisphereWM')
        baw201.connect(WhiteMatterHemisphereNode,'WM_RightHemisphereFileName',DataSink,'WarpedAtlas2Subject.@RightHemisphereWM')

    if 'malf_2012_neuro' in master_config['components']:  ## HACK Do MALF labeling
        good_subjects = [
            '1001', '1004', '1005','1011',
            '1012', '1018', '1019', '1102',
            '1103', '1104', '1120', '1129',
            '1009', '1010', '1013', '1014',
            '1036', '1109', '1117', '1122']

        ## HACK FOR NOW SHOULD BE MORE ELEGANT FROM THE .config file
        BASE_DATA_GRABBER_DIR='/Shared/johnsonhj/HDNI/Neuromorphometrics/20141116_Neuromorphometrics_base_Results/Neuromorphometrics/2012Subscription'

        myLocalMALF = CreateMALFWorkflow("MALF", master_config,good_subjects,BASE_DATA_GRABBER_DIR)
        baw201.connect(myLocalTCWF,'outputspec.t1_average',myLocalMALF,'inputspec.subj_t1_image')
        baw201.connect(myLocalBrainStemWF, 'outputspec.ouputTissuelLabelFilename',myLocalMALF,'inputspec.subj_fixed_head_labels')

        baw201.connect(BResample['template_leftHemisphere'],'outputVolume',myLocalMALF,'inputspec.subj_left_hemisphere')
        baw201.connect(myLocalLMIWF, 'outputspec.outputLandmarksInACPCAlignedSpace' ,myLocalMALF,'inputspec.subj_lmks')
        baw201.connect(atlasBCDNode_S,'template_weights_50Lmks_wts',myLocalMALF,'inputspec.atlasWeightFilename')
        baw201.connect(myLocalMALF,'outputspec.MALF_neuro2012_labelmap',DataSink,'TissueClassify.@MALF_neuro2012_labelmap')
        baw201.connect(myLocalMALF,'outputspec.MALF_fswm_extended_neuro2012_labelmap',DataSink,'TissueClassify.@MALF_fswm_extended_neuro2012_labelmap')
        baw201.connect(myLocalMALF,'outputspec.MALF_fswm_standard_neuro2012_labelmap',DataSink,'TissueClassify.@MALF_fswm_standard_neuro2012_labelmap')
        baw201.connect(myLocalMALF,'outputspec.MALF_extended_snapshot',DataSink,'TissueClassify.@MALF_extended_snapshot')

    return baw201
