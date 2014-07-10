#!/usr/bin/env python
#################################################################################
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

import sys
import string
#"""Import necessary modules from nipype."""
# from nipype.utils.config import config
# config.set('logging', 'log_to_file', 'false')
# config.set_log_dir(os.getcwd())
#--config.set('logging', 'workflow_level', 'DEBUG')
#--config.set('logging', 'interface_level', 'DEBUG')
#--config.set('execution','remove_unnecessary_outputs','false')

from nipype.utils.misc import package_check
# package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
package_check('numpy', '1.3', 'tutorial1')
package_check('scipy', '0.7', 'tutorial1')
package_check('networkx', '1.0', 'tutorial1')
package_check('IPython', '0.10', 'tutorial1')

try:
    from SEMTools import *
except ImportError:
    from AutoWorkup.SEMTools import *

def get_list_element(nestedList, index):
    return nestedList[index]

def getAllT1sLength(allT1s):
    return len(allT1s)

def baseline_workflow(projectid, subjectid, sessionid, master_config, phase='baseline', interpMode='Linear', pipeline_name=''):
    """
    Run autoworkup on a single session

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """

    if not 'auxlmk' in master_config['components'] or not 'tissue_classify' in master_config['components']:
        print "Baseline DataSink requires at least 'AUXLMK' or 'TISSUE_CLASSIFY'"
        master_config['components'].append('auxlmk')
        master_config['components'].append('tissue_classify')
    assert phase in ['baseline', 'longitudinal'], "Unknown phase! Valid entries: 'baseline', 'longitudinal'"

    from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, Directory
    from nipype.interfaces.base import traits, isdefined, BaseInterface
    from nipype.interfaces.utility import Split, Rename, IdentityInterface, Function
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from baseline import get_list_element, getAllT1sLength  # Can we replace with len()?
    from utilities.misc import GenerateWFName
    from PipeLineFunctionHelpers import convertToList, FixWMPartitioning
    from PipeLineFunctionHelpers import UnwrapPosteriorImagesFromDictionaryFunction as flattenDict


    baw201 = pe.Workflow(name=pipeline_name)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['atlasLandmarkFilename', 'atlasWeightFilename',
                                                             'LLSModel', 'inputTemplateModel', 'template_t1',
                                                             'atlasDefinition', 'T1s', 'T2s', 'PDs', 'FLs', 'OTHERs']),
                         run_without_submitting=True, name='inputspec')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['t1_average', 't2_average', 'pd_average', 'fl_average',
                                                              'posteriorImages', 'outputLabels', 'outputHeadLabels',
                                                              'tc_atlas2session_tx',
                                                              'tc_atlas2sessionInverse_tx',
                                                              'BCD_ACPC_T1_CROPPED',
                                                              'outputLandmarksInACPCAlignedSpace',
                                                              'outputLandmarksInInputSpace',
                                                              'output_tx', 'LMIatlasToSubject_tx',
                                                              'writeBranded2DImage',
                                                              'UpdatedPosteriorsList'  # Longitudinal
                                                              ]),
                          run_without_submitting=True, name='outputspec')

    from WorkupT1T2LandmarkInitialization import CreateLandmarkInitializeWorkflow
    DoReverseMapping = False   # Set to true for debugging outputs
    if 'auxlmk' in master_config['components']:
        DoReverseMapping = True
    myLocalLMIWF = CreateLandmarkInitializeWorkflow("LandmarkInitialize", interpMode, DoReverseMapping)

    baw201.connect([(inputsSpec, myLocalLMIWF,
                               [(('T1s', get_list_element, 0), 'inputspec.inputVolume'),
                                ('atlasLandmarkFilename', 'inputspec.atlasLandmarkFilename'),
                                ('atlasWeightFilename', 'inputspec.atlasWeightFilename'),
                                ('LLSModel', 'inputspec.LLSModel'),
                                ('inputTemplateModel', 'inputspec.inputTemplateModel'),
                                ('template_t1', 'inputspec.atlasVolume')]),
                    (myLocalLMIWF, outputsSpec, [('outputspec.outputResampledCroppedVolume',
                                                  'BCD_ACPC_T1_CROPPED'),
                                                 ('outputspec.outputLandmarksInACPCAlignedSpace',
                                                  'outputLandmarksInACPCAlignedSpace'),
                                                 ('outputspec.outputLandmarksInInputSpace',
                                                  'outputLandmarksInInputSpace'),
                                                 ('outputspec.outputTransform', 'output_tx'),
                                                 ('outputspec.atlasToSubjectTransform',
                                                  'LMIatlasToSubject_tx'),
                                                 ('outputspec.writeBranded2DImage', 'writeBranded2DImage')])
                        ])

    if 'tissue_classify' in master_config['components']:
        from WorkupT1T2TissueClassify import CreateTissueClassifyWorkflow
        myLocalTCWF = CreateTissueClassifyWorkflow("TissueClassify", master_config['queue'], master_config['long_q'], interpMode)
        import os
        baw201.connect([(inputsSpec, myLocalTCWF, [('atlasDefinition', 'inputspec.atlasDefinition'),
                                                   ('T1s', 'inputspec.T1List'),
                                                   (('T1s', getAllT1sLength), 'inputspec.T1_count'),
                                                   ('T2s', 'inputspec.T2List'),
                                                   ('PDs', 'inputspec.PDList'),
                                                   ('FLs', 'inputspec.FLList'),
                                                   ('OTHERs', 'inputspec.OtherList')]),
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
                                                    ('outputspec.atlasToSubjectTransform', 'tc_atlas2session_tx'),
                                                    ('outputspec.atlasToSubjectInverseTransform',
                                                     'tc_atlas2sessionInverse_tx')]),
                       ])

    dsName = "{0}_ds_{1}".format(phase, sessionid)
    DataSink = pe.Node(name=dsName, interface=nio.DataSink())
    DataSink.overwrite = master_config['ds_overwrite']
    DataSink.inputs.container = '{0}/{1}/{2}'.format(projectid, subjectid, sessionid)
    DataSink.inputs.base_directory = master_config['resultdir']

    if phase == 'baseline':
        baw201.connect([(outputsSpec, DataSink,  # TODO: change to myLocalTCWF -> DataSink
                                   [(('t1_average', convertToList), 'Baseline.@t1_average'),
                                    (('t2_average', convertToList), 'Baseline.@t2_average'),
                                    (('pd_average', convertToList), 'Baseline.@pd_average'),
                                    (('fl_average', convertToList), 'Baseline.@fl_average'),
                                    (('outputLabels', convertToList), 'Baseline.@labels'),
                                    (('posteriorImages', flattenDict),'TissueClassify')]),
                                 ])
    elif phase == 'longitudinal':
        from PipeLineFunctionHelpers import AccumulateLikeTissuePosteriors
        baw201.connect([(outputsSpec, DataSink, # TODO: change to myLocalTCWF -> DataSink
                                   [(('t1_average', convertToList), 'Longitudinal.@t1_average'),
                                    (('t2_average', convertToList), 'Longitudinal.@t2_average'),
                                    (('pd_average', convertToList), 'Longitudinal.@pd_average'),
                                    (('fl_average', convertToList), 'Longitudinal.@fl_average'),
                                    (('outputLabels', convertToList), 'Longitudinal.@labels'),
                                    (('posteriorImages', flattenDict),'TissueClassify')
                                   ]
                        ),
                       ]
                      )
    else:
        raise NotImplementedError("Missing valid pipeline stage! Options: 'baseline', 'longitudinal'")
    baw201.connect([(outputsSpec, DataSink, # TODO: change to myLocalLMIWF -> DataSink
                                [('outputLandmarksInACPCAlignedSpace', 'ACPCAlign.@outputLandmarks_ACPC'),
                                 ('writeBranded2DImage', 'ACPCAlign.@writeBranded2DImage'),
                                 ('BCD_ACPC_T1_CROPPED', 'ACPCAlign.@BCD_ACPC_T1_CROPPED'),
                                 ('outputLandmarksInInputSpace', 'ACPCAlign.@outputLandmarks_Input'),
                                 ('output_tx', 'ACPCAlign.@output_tx'),
                                 ('LMIatlasToSubject_tx', 'ACPCAlign.@LMIatlasToSubject_tx'),]
                    )
                   ]
                  )

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

    currentBRAINSCreateLabelMapName = 'BRAINSCreateLabelMapFromProbabilityMaps_' + str(subjectid) + "_" + str(sessionid)
    BRAINSCreateLabelMapNode = pe.Node(interface=BRAINSCreateLabelMapFromProbabilityMaps(),
                                       name=currentBRAINSCreateLabelMapName)
    ## TODO:  Fix the file names
    BRAINSCreateLabelMapNode.inputs.dirtyLabelVolume = 'fixed_headlabels_seg.nii.gz'
    BRAINSCreateLabelMapNode.inputs.cleanLabelVolume = 'fixed_brainlabels_seg.nii.gz'

    baw201.connect([(FixWMNode, BRAINSCreateLabelMapNode, [('UpdatedPosteriorsList','inputProbabilityVolume'),
                                                           ('MatchingFGCodeList', 'foregroundPriors'),
                                                           ('MatchingLabelList', 'priorLabelCodes'),
                                                           ('nonAirRegionMask', 'nonAirRegionMask')]),
                    (BRAINSCreateLabelMapNode, DataSink, [('cleanLabelVolume', 'TissueClassify.@outputLabels'),
                                                          ('dirtyLabelVolume',
                                                           'TissueClassify.@outputHeadLabels')]),
                    (myLocalTCWF, DataSink, [('outputspec.atlasToSubjectTransform',
                                              'TissueClassify.@atlas2session_tx'),
                                             ('outputspec.atlasToSubjectInverseTransform',
                                              'TissueClassify.@atlas2sessionInverse_tx')]),
                    (FixWMNode, DataSink, [('UpdatedPosteriorsList', 'TissueClassify.@posteriors')]),
                   ])

    currentAccumulateLikeTissuePosteriorsName = 'AccumulateLikeTissuePosteriors_' + str(subjectid) + "_" + str(sessionid)
    AccumulateLikeTissuePosteriorsNode = pe.Node(interface=Function(function=AccumulateLikeTissuePosteriors,
                                                                    input_names=['posteriorImages'],
                                                                    output_names=['AccumulatePriorsList',
                                                                                  'AccumulatePriorsNames']),
                                                 name=currentAccumulateLikeTissuePosteriorsName)

    baw201.connect([(FixWMNode, AccumulateLikeTissuePosteriorsNode, [('UpdatedPosteriorsList', 'posteriorImages')]),
                    (AccumulateLikeTissuePosteriorsNode, DataSink, [('AccumulatePriorsList',
                                                                     'ACCUMULATED_POSTERIORS.@AccumulateLikeTissuePosteriorsOutputDir')])])
    return baw201
