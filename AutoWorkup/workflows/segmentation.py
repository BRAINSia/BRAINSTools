#!/usr/bin/env python
#################################################################################
## Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
## Language:  Python
##
## Author:  David Welch
##
##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.
##
#################################################################################

#"""Import necessary modules from nipype."""
# from nipype.utils.config import config
# config.set('logging', 'log_to_file', 'false')
# config.set_log_dir(os.getcwd())
#--config.set('logging', 'workflow_level', 'DEBUG')
#--config.set('logging', 'interface_level', 'DEBUG')
#--config.set('execution','remove_unnecessary_outputs','true')

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from nipype.utils.misc import package_check
# package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
package_check('numpy', '1.3', 'tutorial1')
package_check('scipy', '0.7', 'tutorial1')
package_check('networkx', '1.0', 'tutorial1')
package_check('IPython', '0.10', 'tutorial1')

from utilities.misc import CommonANTsRegistrationSettings


def segmentation(projectid, subjectid, sessionid, master_config, onlyT1=True, pipeline_name=''):
    import os.path
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from nipype.interfaces import ants
    from nipype.interfaces.utility import IdentityInterface, Function, Merge
    # Set universal pipeline options
    from nipype import config
    config.update_config(master_config)

    from PipeLineFunctionHelpers import ClipT1ImageWithBrainMask
    from .WorkupT1T2BRAINSCut import CreateBRAINSCutWorkflow
    from utilities.distributed import modify_qsub_args
    from nipype.interfaces.semtools import BRAINSSnapShotWriter

    #CLUSTER_QUEUE=master_config['queue']
    CLUSTER_QUEUE_LONG=master_config['long_q']
    baw200 = pe.Workflow(name=pipeline_name)

    # HACK: print for debugging
    for key, itme in list(master_config.items()):
        print("-" * 30)
        print(key, ":", itme)
    print("-" * 30)
    #END HACK

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['t1_average',
                                                             't2_average',
                                                             'template_t1',
                                                             'hncma_atlas',
                                                             'LMIatlasToSubject_tx',
                                                             'inputLabels',
                                                             'inputHeadLabels',
                                                             'posteriorImages',
                                                             'UpdatedPosteriorsList',
                                                             'atlasToSubjectRegistrationState',
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
                                                             'r_globus_ProbabilityMap',
                                                             'trainModelFile_txtD0060NT0060_gz',
    ]),
                         run_without_submitting=True, name='inputspec')

    # outputsSpec = pe.Node(interface=IdentityInterface(fields=[...]),
    #                       run_without_submitting=True, name='outputspec')

    currentClipT1ImageWithBrainMaskName = 'ClipT1ImageWithBrainMask_' + str(subjectid) + "_" + str(sessionid)
    ClipT1ImageWithBrainMaskNode = pe.Node(interface=Function(function=ClipT1ImageWithBrainMask,
                                                              input_names=['t1_image', 'brain_labels',
                                                                           'clipped_file_name'],
                                                              output_names=['clipped_file']),
                                            name=currentClipT1ImageWithBrainMaskName)
    ClipT1ImageWithBrainMaskNode.inputs.clipped_file_name = 'clipped_from_BABC_labels_t1.nii.gz'

    baw200.connect([(inputsSpec, ClipT1ImageWithBrainMaskNode, [('t1_average', 't1_image'),
                                                                ('inputLabels', 'brain_labels')])])

    currentA2SantsRegistrationPostABCSyN = 'A2SantsRegistrationPostABCSyN_' + str(subjectid) + "_" + str(sessionid)
    ## TODO: It would be great to update the BRAINSABC atlasToSubjectTransform at this point, but
    ##       That requires more testing, and fixes to ANTS to properly collapse transforms.
    ##       For now we are simply creating a dummy node to pass through


    A2SantsRegistrationPostABCSyN = pe.Node(interface=ants.Registration(), name=currentA2SantsRegistrationPostABCSyN)

    many_cpu_ANTsSyN_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE_LONG,8,8,16), 'overwrite': True}
    A2SantsRegistrationPostABCSyN.plugin_args = many_cpu_ANTsSyN_options_dictionary
    CommonANTsRegistrationSettings(
                  antsRegistrationNode=A2SantsRegistrationPostABCSyN,
                  registrationTypeDescription="A2SantsRegistrationPostABCSyN",
                  output_transform_prefix='AtlasToSubjectPostBABC_SyN',
                  output_warped_image='atlas2subjectPostBABC.nii.gz',
                  output_inverse_warped_image='subject2atlasPostBABC.nii.gz',
                  save_state='SavedInternalSyNStatePostBABC.h5',
                  invert_initial_moving_transform=None)

    ## TODO: Try multi-modal registration here
    baw200.connect([(inputsSpec, A2SantsRegistrationPostABCSyN, [('atlasToSubjectRegistrationState', 'restore_state'),
                                                                  ('t1_average', 'fixed_image'),
                                                                  ('template_t1', 'moving_image')])
                   ])


    myLocalSegWF = CreateBRAINSCutWorkflow(projectid,
                                           subjectid,
                                           sessionid,
                                           master_config['queue'],
                                           master_config['long_q'],
                                           "Segmentation",
                                           onlyT1)
    MergeStage2AverageImagesName = "99_mergeAvergeStage2Images_" + str(sessionid)
    MergeStage2AverageImages = pe.Node(interface=Merge(2), run_without_submitting=True,
                                       name=MergeStage2AverageImagesName)

    baw200.connect([(inputsSpec, myLocalSegWF, [('t1_average', 'inputspec.T1Volume'),
                                                ('template_t1', 'inputspec.template_t1'),
                                                ('posteriorImages', "inputspec.posteriorDictionary"),
                                                ('inputLabels', 'inputspec.RegistrationROI'),]),
                    (inputsSpec, MergeStage2AverageImages, [('t1_average', 'in1')]),
                    (A2SantsRegistrationPostABCSyN, myLocalSegWF, [('composite_transform',
                                                                     'inputspec.atlasToSubjectTransform')])
                   ])

    baw200.connect([(inputsSpec, myLocalSegWF,
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
                         ('r_globus_ProbabilityMap', 'inputspec.r_globus_ProbabilityMap'),
                         ('trainModelFile_txtD0060NT0060_gz', 'inputspec.trainModelFile_txtD0060NT0060_gz')
                     ]
                    )]
    )

    if not onlyT1:
        baw200.connect([(inputsSpec, myLocalSegWF, [('t2_average', 'inputspec.T2Volume')]),
                        (inputsSpec, MergeStage2AverageImages, [('t2_average', 'in2')])])
        file_count = 15  # Count of files to merge into MergeSessionSubjectToAtlas
    else:
        file_count = 14  # Count of files to merge into MergeSessionSubjectToAtlas


    ## NOTE: Element 0 of AccumulatePriorsList is the accumulated GM tissue
    # baw200.connect([(AccumulateLikeTissuePosteriorsNode, myLocalSegWF,
    #               [(('AccumulatePriorsList', getListIndex, 0), "inputspec.TotalGM")]),
    #               ])

    ### Now define where the final organized outputs should go.
    DataSink = pe.Node(nio.DataSink(), name="CleanedDenoisedSegmentation_DS_" + str(subjectid) + "_" + str(sessionid))
    DataSink.overwrite = master_config['ds_overwrite']
    DataSink.inputs.base_directory = master_config['resultdir']
    # DataSink.inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'BRAINSCut')
    # DataSink.inputs.regexp_substitutions = GenerateBRAINSCutImagesOutputPattern(projectid, subjectid, sessionid)
    DataSink.inputs.substitutions = [('Segmentations', os.path.join(projectid, subjectid, sessionid, 'CleanedDenoisedRFSegmentations')),
                                     ('subjectANNLabel_', ''),
                                     ('ANNContinuousPrediction', ''),
                                     ('subject.nii.gz', '.nii.gz'),
                                     ('_seg.nii.gz', '_seg.nii.gz'),
                                     ('.nii.gz', '_seg.nii.gz'),
                                     ('_seg_seg', '_seg')]

    baw200.connect([(myLocalSegWF, DataSink, [('outputspec.outputBinaryLeftCaudate', 'Segmentations.@LeftCaudate'),
                                              ('outputspec.outputBinaryRightCaudate', 'Segmentations.@RightCaudate'),
                                              ('outputspec.outputBinaryLeftHippocampus', 'Segmentations.@LeftHippocampus'),
                                              ('outputspec.outputBinaryRightHippocampus', 'Segmentations.@RightHippocampus'),
                                              ('outputspec.outputBinaryLeftPutamen', 'Segmentations.@LeftPutamen'),
                                              ('outputspec.outputBinaryRightPutamen', 'Segmentations.@RightPutamen'),
                                              ('outputspec.outputBinaryLeftThalamus', 'Segmentations.@LeftThalamus'),
                                              ('outputspec.outputBinaryRightThalamus', 'Segmentations.@RightThalamus'),
                                              ('outputspec.outputBinaryLeftAccumben', 'Segmentations.@LeftAccumben'),
                                              ('outputspec.outputBinaryRightAccumben', 'Segmentations.@RightAccumben'),
                                              ('outputspec.outputBinaryLeftGlobus', 'Segmentations.@LeftGlobus'),
                                              ('outputspec.outputBinaryRightGlobus', 'Segmentations.@RightGlobus'),
                                              ('outputspec.outputLabelImageName', 'Segmentations.@LabelImageName'),
                                              ('outputspec.outputCSVFileName', 'Segmentations.@CSVFileName')]),
                    # (myLocalSegWF, DataSink, [('outputspec.cleaned_labels', 'Segmentations.@cleaned_labels')])
                   ])


    MergeStage2BinaryVolumesName = "99_MergeStage2BinaryVolumes_" + str(sessionid)
    MergeStage2BinaryVolumes = pe.Node(interface=Merge(12), run_without_submitting=True,
                                       name=MergeStage2BinaryVolumesName)

    baw200.connect([(myLocalSegWF, MergeStage2BinaryVolumes, [('outputspec.outputBinaryLeftAccumben', 'in1'),
                                                              ('outputspec.outputBinaryLeftCaudate', 'in2'),
                                                              ('outputspec.outputBinaryLeftPutamen', 'in3'),
                                                              ('outputspec.outputBinaryLeftGlobus', 'in4'),
                                                              ('outputspec.outputBinaryLeftThalamus', 'in5'),
                                                              ('outputspec.outputBinaryLeftHippocampus', 'in6'),
                                                              ('outputspec.outputBinaryRightAccumben', 'in7'),
                                                              ('outputspec.outputBinaryRightCaudate', 'in8'),
                                                              ('outputspec.outputBinaryRightPutamen', 'in9'),
                                                              ('outputspec.outputBinaryRightGlobus', 'in10'),
                                                              ('outputspec.outputBinaryRightThalamus', 'in11'),
                                                              ('outputspec.outputBinaryRightHippocampus', 'in12')])
                   ])

    ## SnapShotWriter for Segmented result checking:
    SnapShotWriterNodeName = "SnapShotWriter_" + str(sessionid)
    SnapShotWriter = pe.Node(interface=BRAINSSnapShotWriter(), name=SnapShotWriterNodeName)

    SnapShotWriter.inputs.outputFilename = 'snapShot' + str(sessionid) + '.png'  # output specification
    SnapShotWriter.inputs.inputPlaneDirection = [2, 1, 1, 1, 1, 0, 0]
    SnapShotWriter.inputs.inputSliceToExtractInPhysicalPoint = [-3, -7, -3, 5, 7, 22, -22]

    baw200.connect([(MergeStage2AverageImages, SnapShotWriter, [('out', 'inputVolumes')]),
                    (MergeStage2BinaryVolumes, SnapShotWriter, [('out', 'inputBinaryVolumes')]),
                    (SnapShotWriter, DataSink, [('outputFilename', 'Segmentations.@outputSnapShot')])
                    ])

    # currentAntsLabelWarpToSubject = 'AntsLabelWarpToSubject' + str(subjectid) + "_" + str(sessionid)
    # AntsLabelWarpToSubject = pe.Node(interface=ants.ApplyTransforms(), name=currentAntsLabelWarpToSubject)
    #
    # AntsLabelWarpToSubject.inputs.num_threads = -1
    # AntsLabelWarpToSubject.inputs.dimension = 3
    # AntsLabelWarpToSubject.inputs.output_image = 'warped_hncma_atlas_seg.nii.gz'
    # AntsLabelWarpToSubject.inputs.interpolation = "MultiLabel"
    #
    # baw200.connect([(A2SantsRegistrationPostABCSyN, AntsLabelWarpToSubject, [('composite_transform', 'transforms')]),
    #                 (inputsSpec, AntsLabelWarpToSubject, [('t1_average', 'reference_image'),
    #                                                       ('hncma_atlas', 'input_image')])
    #                 ])
    # #####
    # ### Now define where the final organized outputs should go.
    # AntsLabelWarpedToSubject_DSName = "AntsLabelWarpedToSubject_DS_" + str(sessionid)
    # AntsLabelWarpedToSubject_DS = pe.Node(nio.DataSink(), name=AntsLabelWarpedToSubject_DSName)
    # AntsLabelWarpedToSubject_DS.overwrite = master_config['ds_overwrite']
    # AntsLabelWarpedToSubject_DS.inputs.base_directory = master_config['resultdir']
    # AntsLabelWarpedToSubject_DS.inputs.substitutions = [('AntsLabelWarpedToSubject', os.path.join(projectid, subjectid, sessionid, 'AntsLabelWarpedToSubject'))]
    #
    # baw200.connect([(AntsLabelWarpToSubject, AntsLabelWarpedToSubject_DS, [('output_image', 'AntsLabelWarpedToSubject')])])

    MergeSessionSubjectToAtlasName = "99_MergeSessionSubjectToAtlas_" + str(sessionid)
    MergeSessionSubjectToAtlas = pe.Node(interface=Merge(file_count), run_without_submitting=True,
                                         name=MergeSessionSubjectToAtlasName)

    baw200.connect([(myLocalSegWF, MergeSessionSubjectToAtlas, [('outputspec.outputBinaryLeftAccumben', 'in1'),
                                                                ('outputspec.outputBinaryLeftCaudate', 'in2'),
                                                                ('outputspec.outputBinaryLeftPutamen', 'in3'),
                                                                ('outputspec.outputBinaryLeftGlobus', 'in4'),
                                                                ('outputspec.outputBinaryLeftThalamus', 'in5'),
                                                                ('outputspec.outputBinaryLeftHippocampus', 'in6'),
                                                                ('outputspec.outputBinaryRightAccumben', 'in7'),
                                                                ('outputspec.outputBinaryRightCaudate', 'in8'),
                                                                ('outputspec.outputBinaryRightPutamen', 'in9'),
                                                                ('outputspec.outputBinaryRightGlobus', 'in10'),
                                                                ('outputspec.outputBinaryRightThalamus', 'in11'),
                                                                ('outputspec.outputBinaryRightHippocampus', 'in12')]),
                    # (FixWMPartitioningNode, MergeSessionSubjectToAtlas, [('UpdatedPosteriorsList', 'in13')]),
                    (inputsSpec, MergeSessionSubjectToAtlas, [('UpdatedPosteriorsList', 'in13')]),
                    (inputsSpec, MergeSessionSubjectToAtlas, [('t1_average', 'in14')])
                    ])

    if not onlyT1:
        assert file_count == 15
        baw200.connect([(inputsSpec, MergeSessionSubjectToAtlas, [('t2_average', 'in15')])])

    LinearSubjectToAtlasANTsApplyTransformsName = 'LinearSubjectToAtlasANTsApplyTransforms_' + str(sessionid)
    LinearSubjectToAtlasANTsApplyTransforms = pe.MapNode(interface=ants.ApplyTransforms(), iterfield=['input_image'],
                                                         name=LinearSubjectToAtlasANTsApplyTransformsName)
    LinearSubjectToAtlasANTsApplyTransforms.inputs.num_threads = -1
    LinearSubjectToAtlasANTsApplyTransforms.inputs.interpolation = 'Linear'

    baw200.connect([(A2SantsRegistrationPostABCSyN, LinearSubjectToAtlasANTsApplyTransforms, [('inverse_composite_transform',
                                                                                              'transforms')]),
                    (inputsSpec, LinearSubjectToAtlasANTsApplyTransforms, [('template_t1', 'reference_image')]),
                    (MergeSessionSubjectToAtlas, LinearSubjectToAtlasANTsApplyTransforms, [('out', 'input_image')])
                    ])

    MergeMultiLabelSessionSubjectToAtlasName = "99_MergeMultiLabelSessionSubjectToAtlas_" + str(sessionid)
    MergeMultiLabelSessionSubjectToAtlas = pe.Node(interface=Merge(2), run_without_submitting=True,
                                                   name=MergeMultiLabelSessionSubjectToAtlasName)

    baw200.connect([(inputsSpec, MergeMultiLabelSessionSubjectToAtlas, [('inputLabels', 'in1'),
                                                                        ('inputHeadLabels', 'in2')])
                   ])

    ### This is taking this sessions RF label map back into NAC atlas space.
    #{
    MultiLabelSubjectToAtlasANTsApplyTransformsName = 'MultiLabelSubjectToAtlasANTsApplyTransforms_' + str(sessionid) + '_map'
    MultiLabelSubjectToAtlasANTsApplyTransforms = pe.MapNode(interface=ants.ApplyTransforms(), iterfield=['input_image'],
                                                             name=MultiLabelSubjectToAtlasANTsApplyTransformsName)
    MultiLabelSubjectToAtlasANTsApplyTransforms.inputs.num_threads = -1
    MultiLabelSubjectToAtlasANTsApplyTransforms.inputs.interpolation = 'MultiLabel'

    baw200.connect([(A2SantsRegistrationPostABCSyN, MultiLabelSubjectToAtlasANTsApplyTransforms,
                     [('inverse_composite_transform', 'transforms')]),
                      (inputsSpec, MultiLabelSubjectToAtlasANTsApplyTransforms, [('template_t1', 'reference_image')]),
                      (MergeMultiLabelSessionSubjectToAtlas, MultiLabelSubjectToAtlasANTsApplyTransforms,
                       [('out', 'input_image')])
                   ])
    #}
    ### Now we must take the sessions to THIS SUBJECTS personalized atlas.
    #{
    #}

    ### Now define where the final organized outputs should go.
    Subj2Atlas_DSName = "SubjectToAtlas_DS_" + str(sessionid)
    Subj2Atlas_DS = pe.Node(nio.DataSink(), name=Subj2Atlas_DSName)
    Subj2Atlas_DS.overwrite = master_config['ds_overwrite']
    Subj2Atlas_DS.inputs.base_directory = master_config['resultdir']
    Subj2Atlas_DS.inputs.regexp_substitutions = [(r'_LinearSubjectToAtlasANTsApplyTransforms_[^/]*',
                                                  r'' + sessionid + '/')]

    baw200.connect([(LinearSubjectToAtlasANTsApplyTransforms, Subj2Atlas_DS,
                     [('output_image', 'SubjectToAtlasWarped.@linear_output_images')])])

    Subj2AtlasTransforms_DSName = "SubjectToAtlasTransforms_DS_" + str(sessionid)
    Subj2AtlasTransforms_DS = pe.Node(nio.DataSink(), name=Subj2AtlasTransforms_DSName)
    Subj2AtlasTransforms_DS.overwrite = master_config['ds_overwrite']
    Subj2AtlasTransforms_DS.inputs.base_directory = master_config['resultdir']
    Subj2AtlasTransforms_DS.inputs.regexp_substitutions = [(r'SubjectToAtlasWarped',
                                                            r'SubjectToAtlasWarped/' + sessionid + '/')]

    baw200.connect([(A2SantsRegistrationPostABCSyN, Subj2AtlasTransforms_DS,
                     [('composite_transform', 'SubjectToAtlasWarped.@composite_transform'),
                      ('inverse_composite_transform', 'SubjectToAtlasWarped.@inverse_composite_transform')])])
    # baw200.connect([(MultiLabelSubjectToAtlasANTsApplyTransforms, Subj2Atlas_DS, [('output_image', 'SubjectToAtlasWarped.@multilabel_output_images')])])

    if master_config['plugin_name'].startswith('SGE'):  # for some nodes, the qsub call needs to be modified on the cluster
        A2SantsRegistrationPostABCSyN.plugin_args = {'template': master_config['plugin_args']['template'], 'overwrite': True,
                                                      'qsub_args': modify_qsub_args(master_config['queue'], 8, 8, 24)}
        SnapShotWriter.plugin_args = {'template': master_config['plugin_args']['template'], 'overwrite': True,
                                      'qsub_args': modify_qsub_args(master_config['queue'], 1, 1, 1 )}
        LinearSubjectToAtlasANTsApplyTransforms.plugin_args = {'template': master_config['plugin_args']['template'],
                                                               'overwrite': True,
                                                               'qsub_args': modify_qsub_args(master_config['queue'], 1, 1, 1 )}
        MultiLabelSubjectToAtlasANTsApplyTransforms.plugin_args = {'template': master_config['plugin_args']['template'],
                                                                   'overwrite': True,
                                                                   'qsub_args': modify_qsub_args(master_config['queue'], 1, 1, 1 )}

    return baw200
