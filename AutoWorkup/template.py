#! /usr/bin/env python
"""
template.py
=========
This program is used to generate the subject- and session-specific workflows for BRAINSTool processing

Usage:
  template.py subject [--rewrite-datasinks] [--wfrun PLUGIN] [--use-sentinal] [--use-shuffle] [--dotfilename PFILE] --pe ENV --ExperimentConfig FILE SUBJECTS...
  template.py population [--rewrite-datasinks] [--wfrun PLUGIN] [--use-sentinal] [--use-shuffle] [--dotfilename PFILE] --pe ENV --ExperimentConfig FILE -p POPFILE SUBJECTS...
  template.py -v | --version
  template.py -h | --help

Arguments:
  SUBJECTS              List of subject IDs to process

Options:
  -h, --help            Show this help and exit
  -v, --version         Print the version and exit
  --dotfilename=PFILE   Turn on printing pipeline to file PFILE
  --rewrite-datasinks   Turn on the Nipype option to overwrite all files in the 'results' directory
  --use-sentinal        Use the t1_average file as a marker to determine if session needs to be run
  --use-shuffle         Shuffle the subjects randomly to minimize multiple runs collision probability
  --pe=ENV              The processing environment to use from configuration file
  --wfrun=PLUGIN        The name of the workflow plugin option (default: 'local')
  --ExperimentConfig=FILE   The configuration file
  -p POPFILE            File containing dictionary 'groups' in Python syntax

Examples:
  $ template.py population --pe OSX --ExperimentConfig my_baw.config -p groupfile.txt all
  $ template.py subject --wfrun helium_all.q --pe OSX --ExperimentConfig my_baw.config 1058 1059
  $ template.py subject --rewrite-datasinks --pe OSX --ExperimentConfig my_baw.config 2001

"""
import glob
import os
import sys
import traceback

from baw_exp import OpenSubjectDatabase

def get_processed(resultdir):
    import glob
    # resultdir/subject_dir/Atlas/AVG_T1.nii.gz
    sential_file_pattern = "*/Atlas/AVG_template_rightHemisphere.nii.gz"
    processedPaths = glob.glob( os.path.join(resultdir, sential_file_pattern) )
    processed = [ os.path.basename(os.path.dirname(os.path.dirname(s))) for s in processedPaths ]
    print("SKIPPING COMPLETED: {0}".format(processed))
    return processed

def generate_group_dictionary(input_subjects, cache, resultdir, prefix, dbfile, useSentinal, shuffle=False, groupfile=None):
    import random
    _temp = OpenSubjectDatabase(cache, ['all'], prefix, dbfile)
    if "all" in input_subjects:
        input_subjects =  _temp.getAllSubjects();
    if useSentinal:
        print("="*80)
        print("Using Sentinal Files to Limit Jobs Run")
        _all_subjects = set(input_subjects )
        processed = set( get_processed( resultdir ) )
        subjects = list( _all_subjects - processed) #NOTE - in set operation notation removes values
    else:
        subjects = input_subjects
    if shuffle:
        random.shuffle(subjects)  # randomly shuffle to get max cluster efficiency

    if groupfile is not None:
        group_dictionary = _temp.getSessionsFromGroup(groupfile)
    else:
        group_dictionary = dict()
        for subject in subjects:
            group_dictionary[subject]=_temp.getSessionsFromSubject(subject)
    return subjects, group_dictionary

## Merge the different groups together
def MergeByExtendListElements(t1s, t2s, pds, fls, labels, posteriors, passive_intensities, passive_masks):
    """
    *** NOTE:  All input lists MUST have the same number of elements (even if they are null) ***

    output = [{'T1':        os.path.join(mydatadir, '01_T1_half.nii.gz'),
               'INV_T1':    os.path.join(mydatadir, '01_T1_inv_half.nii.gz'),
               'LABEL_MAP': os.path.join(mydatadir, '01_T1_inv_half.nii.gz')
              },
              {'T1':        os.path.join(mydatadir, '02_T1_half.nii.gz'),
               'INV_T1':    os.path.join(mydatadir, '02_T1_inv_half.nii.gz'),
               'LABEL_MAP': os.path.join(mydatadir, '02_T1_inv_half.nii.gz')
              },
              {'T1':        os.path.join(mydatadir, '03_T1_half.nii.gz'),
               'INV_T1':    os.path.join(mydatadir, '03_T1_inv_half.nii.gz'),
               'LABEL_MAP': os.path.join(mydatadir, '03_T1_inv_half.nii.gz')
              }
             ]
    #          SUBJECT_01                    SUBJECT_02                        SUBJECT_03
    labels = ['brain_label_seg.nii.gz',      'brain_label_seg.nii.gz',          ...      ]
    pds    = [None,                          None,                              ...      ]
    t1s    = ['t1_average_BRAINSABC.nii.gz', 't1_average_BRAINSABC.nii.gz',     ...      ]
    t2s    = ['t2_average_BRAINSABC.nii.gz', 't2_average_BRAINSABC.nii.gz',     ...      ]

    """
    # print "t1s", t1s
    # print "t2s", t2s
    # print "pds", pds
    # print "fls", fls
    # print "labels", labels
    # print "$$$$$$$$$$$$$$$$$$$$$$$"
    # print "posteriors", posteriors
    ListOfImagesDictionaries = [dict() for i in t1s]  # Initial list with empty dictionaries
    ## HACK:  Need to make it so that AVG_AIR.nii.gz has a background value of 1
    registrationImageTypes = ['T1']  # ['T1','T2'] someday.
    DefaultContinuousInterpolationType = 'Linear'  # or 'LanczosWindowedSinc' ('Linear' for speed)
    interpolationMapping = {'T1': DefaultContinuousInterpolationType,
                            'T2': DefaultContinuousInterpolationType,
                            'PD': DefaultContinuousInterpolationType,
                            'FL': DefaultContinuousInterpolationType,
                            'BRAINMASK': 'MultiLabel'
                            }
    for index in range(len(t1s)):
        if t1s[index] is not None:
            ListOfImagesDictionaries[index]['T1'] = t1s[index]
        if isinstance(t2s, list) and t2s[index] is not None:
            ListOfImagesDictionaries[index]['T2'] = t2s[index]
        if isinstance(pds, list) and pds[index] is not None:
            ListOfImagesDictionaries[index]['PD'] = pds[index]
        if isinstance(fls, list) and fls[index] is not None:
            ListOfImagesDictionaries[index]['FL'] = fls[index]
        if labels[index] is not None:
            ListOfImagesDictionaries[index]['BRAINMASK'] = labels[index]
        print ListOfImagesDictionaries[index]
        for key, value in posteriors.items():
            # print "key:", key, " -> value:", value
            ListOfImagesDictionaries[index][key] = value[index]
            interpolationMapping[key] = DefaultContinuousInterpolationType
        for key, value in passive_intensities.items():
            # print "key:", key, " -> value:", value
            ListOfImagesDictionaries[index][key] = value[index]
            interpolationMapping[key] = DefaultContinuousInterpolationType
        for key, value in passive_masks.items():
            # print "key:", key, " -> value:", value
            ListOfImagesDictionaries[index][key] = value[index]
            interpolationMapping[key] = 'MultiLabel'

    # print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    #print "ListOfImagesDictionaries", ListOfImagesDictionaries
    #print "registrationImageTypes", registrationImageTypes
    #print "interpolationMapping", interpolationMapping
    return ListOfImagesDictionaries, registrationImageTypes, interpolationMapping


def xml_filename(group):
    return 'AtlasDefinition_{0}.xml'.format(group)

def get_dict_value(dictonary, key):
    return dictionary[key]

def _template_runner(argv, environment, experiment, pipeline_options, cluster):
    print "Getting subjects from database..."
    # subjects = argv["--subjects"].split(',')
    subjects, group_dictionary = generate_group_dictionary(argv['SUBJECTS'],
            experiment['cachedir'],
            experiment['resultdir'],
            environment['prefix'],
            experiment['dbfile'],
            argv['--use-sentinal'], argv['--use-shuffle']
            ) # Build database before parallel section
    useSentinal = argv['--use-sentinal']
    for group_key, group_value in group_dictionary.keys():
        print("Processing atlas generation for {0}".format(group_key))
        print("="*80)
        print "Copying Atlas directory and determining appropriate Nipype options..."
        pipeline_options = nipype_options(argv, pipeline_options, cluster, experiment, environment)  # Generate Nipype options
        print "Dispatching jobs to the system..."
        ######
        ###### Now start workflow construction # FIXME: def construct_workflow():
        ######
        # Set universal pipeline options
        nipype_config.update_config(pipeline_options)

        ready_for_template_building = True
        for session in group_value:
            path_test = os.path.join(experiment['previousresult'],'*/*/{0}/TissueClassify/t1_average_BRAINSABC.nii.gz'.format(session))
            t1_file_result = glob.glob(path_test)
            try:
                assert len(t1_file_result) == 1, "Only one t1_average_BRAINSABC.nii.gz per session"
            except AssertionError:
                print("Incorrect number of t1 images found for data grabber\n\tat {0}".format(path_test))
                ready_for_template_building = False
        if not ready_for_template_building:
            print("TEMPORARY SKIPPING:  Not ready to process {0}".format(group_key))
            continue

        base_output_directory = os.path.join(pipeline_options['logging']['log_directory'], group_key)
        if argv['--groupfile'] is None:
            template = pe.Workflow(name='SubjectAtlas_Template_' + group_key)
        else:
            template = pe.Workflow(name='PopulationAtlas_Template_' + group_key)
        template.base_dir = base_output_directory

        subjectNode = pe.Node(interface=IdentityInterface(fields=['subject']), run_without_submitting=True, name='99_subjectIterator')
        subjectNode.inputs.subject = group_key

        sessionsExtractorNode = pe.Node(Function(function=get_dict_value,
                                                 input_names=['dictionary','key'],
                                                 output_names=['sessions']),
                                        run_without_submitting=True,
                                        name="99_sessionsExtractor")
        sessionsExtractorNode.inputs.dictionary = group_dictionary



        baselineOptionalDG = pe.MapNode(nio.DataGrabber(infields=['session'],
                                                        outfields=['t2_average', 'pd_average', 'fl_average'],
                                                        run_without_submitting=True),
                                        run_without_submitting=True,
                                        iterfield=['session'],
                                        name='BaselineOptional_DG')
        baselineOptionalDG.inputs.base_directory = experiment['previousresult']
        baselineOptionalDG.inputs.sort_filelist = True
        baselineOptionalDG.inputs.raise_on_empty = False
        baselineOptionalDG.inputs.template = '*'
        baselineOptionalDG.inputs.field_template = {'t2_average':'*/*/%s/TissueClassify/t2_average_BRAINSABC.nii.gz',
                                                    'pd_average':'*/*/%s/TissueClassify/pd_average_BRAINSABC.nii.gz',
                                                    'fl_average':'*/*/%s/TissueClassify/fl_average_BRAINSABC.nii.gz'}
        baselineOptionalDG.inputs.template_args  = {'t2_average':[['session']],
                                                    'pd_average':[['session']],
                                                    'fl_average':[['session']]}


        baselineRequiredDG = pe.MapNode(nio.DataGrabber(infields=['session'],
                                                        outfields=['t1_average', 'brainMaskLabels', 'posteriorImages',
                                                                   'passive_intensities','passive_masks', 'BCD_ACPC_Landmarks_fcsv'],
                                                        run_without_submitting=True),
                                        run_without_submitting=True,
                                        iterfield=['session'], name='Baseline_DG')
        baselineRequiredDG.inputs.base_directory = experiment['previousresult']
        baselineRequiredDG.inputs.sort_filelist = True
        baselineRequiredDG.inputs.raise_on_empty = True
        baselineRequiredDG.inputs.template = '*'
        posterior_files = ['AIR', 'BASAL', 'CRBLGM', 'CRBLWM', 'CSF', 'GLOBUS', 'HIPPOCAMPUS',
                           'NOTCSF', 'NOTGM', 'NOTVB', 'NOTWM', 'SURFGM', 'THALAMUS', 'VB', 'WM']
        passive_intensities_files = [
            'rho.nii.gz',
            'phi.nii.gz',
            'theta.nii.gz',
            'l_thalamus_ProbabilityMap.nii.gz',
            'r_accumben_ProbabilityMap.nii.gz',
            'l_globus_ProbabilityMap.nii.gz',
            'l_accumben_ProbabilityMap.nii.gz',
            'l_caudate_ProbabilityMap.nii.gz',
            'l_putamen_ProbabilityMap.nii.gz',
            'r_thalamus_ProbabilityMap.nii.gz',
            'r_putamen_ProbabilityMap.nii.gz',
            'r_caudate_ProbabilityMap.nii.gz',
            'r_hippocampus_ProbabilityMap.nii.gz',
            'r_globus_ProbabilityMap.nii.gz',
            'l_hippocampus_ProbabilityMap.nii.gz'
            ]
        passive_mask_files = [
            'template_WMPM2_labels.nii.gz',
            'hncma_atlas.nii.gz',
            'template_nac_labels.nii.gz',
            'template_leftHemisphere.nii.gz',
            'template_rightHemisphere.nii.gz',
            'template_ventricles.nii.gz'
            ]
        baselineRequiredDG.inputs.field_template = {'t1_average':             '*/*/%s/TissueClassify/t1_average_BRAINSABC.nii.gz',
                                                    'brainMaskLabels':        '*/*/%s/TissueClassify/complete_brainlabels_seg.nii.gz',
                                                    'BCD_ACPC_Landmarks_fcsv':'*/*/%s/ACPCAlign/BCD_ACPC_Landmarks.fcsv',
                                                    'posteriorImages':        '*/*/%s/TissueClassify/POSTERIOR_%s.nii.gz',
                                                    'passive_intensities':    '*/*/%s/WarpedAtlas2Subject/%s',
                                                    'passive_masks':          '*/*/%s/WarpedAtlas2Subject/%s'}
        baselineRequiredDG.inputs.template_args  = {'t1_average':[['session']],
                                                    'brainMaskLabels':[['session']],
                                                    'BCD_ACPC_Landmarks_fcsv':[['session']],
                                                    'posteriorImages':[['session', posterior_files]],
                                                    'passive_intensities':[['session', passive_intensities_files]],
                                                    'passive_masks':[['session', passive_mask_files]]}

        MergeByExtendListElementsNode = pe.Node(Function(function=MergeByExtendListElements,
                                                         input_names=['t1s', 't2s', 'pds', 'fls', 'labels', 'posteriors',
                                                                      'passive_intensities', 'passive_masks'],
                                                         output_names=['ListOfImagesDictionaries', 'registrationImageTypes',
                                                                       'interpolationMapping']),
                                                run_without_submitting=True,
                                                name="99_MergeByExtendListElements")

        template.connect([(subjectNode, sessionsExtractorNode, [('subject','key')]),
                          (sessionsExtractorNode, baselineRequiredDG, [('sessions', 'session')]),
                          (sessionsExtractorNode, baselineOptionalDG, [('sessions', 'session')]),
                          (baselineRequiredDG, MergeByExtendListElementsNode,
                                    [('t1_average', 't1s'),
                                     ('brainMaskLabels', 'labels'),
                                     (('posteriorImages',
                                        ConvertSessionsListOfPosteriorListToDictionaryOfSessionLists), 'posteriors')
                                     ]),
                          (baselineOptionalDG, MergeByExtendListElementsNode,
                                    [
                                     ('t2_average', 't2s'),
                                     ('pd_average', 'pds'),
                                     ('fl_average', 'fls')
                                     ]),
                          (baselineRequiredDG, MergeByExtendListElementsNode,
                                     [
                                      (('passive_intensities',
                                        ConvertSessionsListOfPosteriorListToDictionaryOfSessionLists), 'passive_intensities')
                                     ]),
                          (baselineRequiredDG, MergeByExtendListElementsNode,
                                     [
                                     (('passive_masks',
                                        ConvertSessionsListOfPosteriorListToDictionaryOfSessionLists), 'passive_masks')
                                     ])
                        ])

        myInitAvgWF = pe.Node(interface=ants.AverageImages(), name='Atlas_antsSimpleAverage')  # was 'Phase1_antsSimpleAverage'
        myInitAvgWF.inputs.dimension = 3
        myInitAvgWF.inputs.normalize = True
        template.connect(baselineRequiredDG, 't1_average', myInitAvgWF, "images")
        ####################################################################################################
        # TEMPLATE_BUILD_RUN_MODE = 'MULTI_IMAGE'
        # if numSessions == 1:
        #     TEMPLATE_BUILD_RUN_MODE = 'SINGLE_IMAGE'
        ####################################################################################################
        buildTemplateIteration1 = BAWantsRegistrationTemplateBuildSingleIterationWF('iteration01')
        # buildTemplateIteration2 = buildTemplateIteration1.clone(name='buildTemplateIteration2')
        buildTemplateIteration2 = BAWantsRegistrationTemplateBuildSingleIterationWF('Iteration02')

        CreateAtlasXMLAndCleanedDeformedAveragesNode = pe.Node(interface=Function(function=CreateAtlasXMLAndCleanedDeformedAverages,
                                                              input_names=['t1_image', 'deformed_list', 'AtlasTemplate', 'outDefinition'],
                                                              output_names=['outAtlasFullPath', 'clean_deformed_list']),
                                           # This is a lot of work, so submit it run_without_submitting=True,
                                           run_without_submitting=True,  # HACK:  THIS NODE REALLY SHOULD RUN ON THE CLUSTER!
                                           name='99_CreateAtlasXMLAndCleanedDeformedAverages')

        if pipeline_options['plugin_name'].startswith('SGE'):  # for some nodes, the qsub call needs to be modified on the cluster

            CreateAtlasXMLAndCleanedDeformedAveragesNode.plugin_args = {'template': pipeline_options['plugin_args']['template'],
                                                    'qsub_args': modify_qsub_args(cluster['queue'], 1, 1, 1),
                                                    'overwrite': True}
            for bt in [buildTemplateIteration1, buildTemplateIteration2]:
                BeginANTS = bt.get_node("BeginANTS")
                BeginANTS.plugin_args = {'template': pipeline_options['plugin_args']['template'], 'overwrite': True,
                                         'qsub_args': modify_qsub_args(cluster['queue'], 4, 2, 4)}
                wimtdeformed = bt.get_node("wimtdeformed")
                wimtdeformed.plugin_args = {'template': pipeline_options['plugin_args']['template'], 'overwrite': True,
                                            'qsub_args': modify_qsub_args(cluster['queue'], 2, 2, 2)}

                #AvgAffineTransform = bt.get_node("AvgAffineTransform")
                #AvgAffineTransform.plugin_args = {'template': pipeline_options['plugin_args']['template'], 'overwrite': True,
                #                                  'qsub_args': modify_qsub_args(cluster['queue'], 2, 1, 1)}

                wimtPassivedeformed = bt.get_node("wimtPassivedeformed")
                wimtPassivedeformed.plugin_args = {'template': pipeline_options['plugin_args']['template'], 'overwrite': True,
                                                    'qsub_args': modify_qsub_args(cluster['queue'], 2, 2, 4)}

        # Running off previous baseline experiment
        NACCommonAtlas = MakeAtlasNode(experiment['atlascache'], 'NACCommonAtlas_{0}'.format('subject'),
                ['S_BRAINSABCSupport'] ) ## HACK : replace 'subject' with subject id once this is a loop rather than an iterable.
        template.connect([(myInitAvgWF, buildTemplateIteration1, [('output_average_image', 'inputspec.fixed_image')]),
                          (MergeByExtendListElementsNode, buildTemplateIteration1, [('ListOfImagesDictionaries', 'inputspec.ListOfImagesDictionaries'),
                                                                                    ('registrationImageTypes', 'inputspec.registrationImageTypes'),
                                                                                    ('interpolationMapping','inputspec.interpolationMapping')]),
                          (buildTemplateIteration1, buildTemplateIteration2, [('outputspec.template', 'inputspec.fixed_image')]),
                          (MergeByExtendListElementsNode, buildTemplateIteration2, [('ListOfImagesDictionaries', 'inputspec.ListOfImagesDictionaries'),
                                                                                    ('registrationImageTypes','inputspec.registrationImageTypes'),
                                                                                    ('interpolationMapping', 'inputspec.interpolationMapping')]),
                          (subjectNode, CreateAtlasXMLAndCleanedDeformedAveragesNode, [(('subject', xml_filename), 'outDefinition')]),
                          (NACCommonAtlas, CreateAtlasXMLAndCleanedDeformedAveragesNode, [('ExtendedAtlasDefinition_xml_in', 'AtlasTemplate')]),
                          (buildTemplateIteration2, CreateAtlasXMLAndCleanedDeformedAveragesNode, [('outputspec.template', 't1_image'),
                                                                               ('outputspec.passive_deformed_templates', 'deformed_list')]),
                          ])


        ## Genearate an average lmks file.
        myAverageLmk = pe.Node(interface = GenerateAverageLmkFile(), name="myAverageLmk" )
        myAverageLmk.inputs.outputLandmarkFile = "AVG_LMKS.fcsv"
        template.connect(baselineRequiredDG,'BCD_ACPC_Landmarks_fcsv',myAverageLmk,'inputLandmarkFiles')

        # Create DataSinks
        SubjectAtlas_DataSink = pe.Node(nio.DataSink(), name="Subject_DS")
        SubjectAtlas_DataSink.overwrite = pipeline_options['ds_overwrite']
        SubjectAtlas_DataSink.inputs.base_directory = experiment['resultdir']

        template.connect([(subjectNode, SubjectAtlas_DataSink, [('subject', 'container')]),
                          (CreateAtlasXMLAndCleanedDeformedAveragesNode, SubjectAtlas_DataSink, [('outAtlasFullPath', 'Atlas.@definitions')]),
                          (CreateAtlasXMLAndCleanedDeformedAveragesNode, SubjectAtlas_DataSink, [('clean_deformed_list', 'Atlas.@passive_deformed_templates')]),

                          (subjectNode, SubjectAtlas_DataSink, [(('subject', outputPattern), 'regexp_substitutions')]),
                          (buildTemplateIteration2, SubjectAtlas_DataSink, [('outputspec.template', 'Atlas.@template')]),
                          (myAverageLmk,SubjectAtlas_DataSink,[('outputLandmarkFile','Atlas.@outputLandmarkFile')]),
                         ])

        dotfilename = argv['--dotfilename']
        if dotfilename is not None:
            print("WARNING: Printing workflow, but not running pipeline")
            print_workflow(template, plugin=pipeline_options['plugin_name'], dotfilename=dotfilename)
        else:
            run_workflow(template, plugin=pipeline_options['plugin_name'], plugin_args=pipeline_options['plugin_args'])

if __name__ == '__main__':
    import sys
    from AutoWorkup import setup_environment

    from docopt import docopt

    argv = docopt(__doc__, version='1.1')
    print argv
    valid_phases = ['subject-template-generation', 'population-template-generation']
    if argv['subject']:
        this_phase = valid_phases[0]
    elif argv['population']:
        this_phase = valid_phases[1]
    else:
        raise RuntimeError("Unknown workphase attempted!")
    print("Running {0} workphase".format(this_phase))
    print '=' * 100
    argv['--workphase'] = this_phase
    environment, experiment, pipeline, cluster = setup_environment(argv)
    from nipype import config as nipype_config
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from nipype.interfaces.utility import IdentityInterface, Function
    import nipype.interfaces.ants as ants

    from PipeLineFunctionHelpers import ConvertSessionsListOfPosteriorListToDictionaryOfSessionLists
    from workflows.atlasNode import MakeAtlasNode, CreateAtlasXMLAndCleanedDeformedAverages
    from utilities.misc import GenerateSubjectOutputPattern as outputPattern
    from utilities.distributed import modify_qsub_args
    from workflows.utils import run_workflow, print_workflow
    from BAWantsRegistrationBuildTemplate import BAWantsRegistrationTemplateBuildSingleIterationWF
    from utilities.configFileParser import nipype_options
    from SEMTools.testing.generateaveragelmkfile import GenerateAverageLmkFile
    exit = _template_runner(argv, environment, experiment, pipeline, cluster)
    sys.exit(exit)
