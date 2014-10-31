#! /usr/bin/env python
"""
template.py
=========
This program is used to generate the subject- and session-specific workflows for BRAINSTool processing

Usage:
  template.py [--rewrite-datasinks] [--wfrun PLUGIN] [--use-sentinal] [--dotfilename PFILE] --workphase WORKPHASE --pe ENV --ExperimentConfig FILE SUBJECTS...
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
  --pe=ENV              The processing environment to use from configuration file
  --wfrun=PLUGIN        The name of the workflow plugin option (default: 'local')
  --workphase WORKPHASE The type of processing to be done only VALID is ['subject-template-generation']
  --ExperimentConfig=FILE   The configuration file

Examples:
  $ template.py --pe OSX --ExperimentConfig my_baw.config all
  $ template.py --wfrun helium_all.q --pe OSX --ExperimentConfig my_baw.config 1058 1059
  $ template.py --rewrite-datasinks --pe OSX --ExperimentConfig my_baw.config 2001

"""
import os
import sys
import traceback

from baw_exp import OpenSubjectDatabase

def get_processed_subjects( resultdir ):
    import glob
    # resultdir/subject_dir/Atlas/AVG_T1.nii.gz
    sential_file_pattern = "*/Atlas/AVG_T1.nii.gz"
    processedSubjectsPaths = glob.glob( os.path.join(resultdir, sential_file_pattern) )
    processedSubjects = [ os.path.basename(os.path.dirname(os.path.dirname(s))) for s in processedSubjectsPaths ]
    return processedSubjects

def get_subjects_sessions_dictionary(subjects, cache, resultdir, prefix, dbfile, useSentinal, shuffle=False):
    import random
    _temp = OpenSubjectDatabase(cache, ['all'], prefix, dbfile)
    if useSentinal:
        print("="*80)
        print("Using Sentinal Files to Limit Jobs Run")
        _all_subjects = set( _temp.getAllSubjects() )
        _processed_subjects = set( get_processed_subjects( resultdir ) )
        subjects = list( _all_subjects - _processed_subjects ) #NOTE - in set operation notation removes values
    elif "all" in subjects:
        subjects = _temp.getAllSubjects()
    if shuffle:
        random.shuffle(subjects)  # randomly shuffle to get max
    subject_sessions_dictionary = dict()
    for subject in subjects:
        subject_sessions_dictionary[subject]=_temp.getSessionsFromSubject(subject)
    return subjects,subject_sessions_dictionary

def MergeByExtendListElements(t1s, t2s, pds, fls, labels, posteriors):
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
    labels = ['brain_label_seg.nii.gz', 'brain_label_seg.nii.gz']
    pds = [None, None]
    t1s = ['t1_average_BRAINSABC.nii.gz', 't1_average_BRAINSABC.nii.gz']
    t2s = ['t2_average_BRAINSABC.nii.gz', 't2_average_BRAINSABC.nii.gz']

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
    for list_index in range(len(t1s)):
        if t1s[list_index] is not None:
            ListOfImagesDictionaries[list_index]['T1'] = t1s[list_index]
        if isinstance(t2s, list) and t2s[list_index] is not None:
            ListOfImagesDictionaries[list_index]['T2'] = t2s[list_index]
        if isinstance(pds, list) and pds[list_index] is not None:
            ListOfImagesDictionaries[list_index]['PD'] = pds[list_index]
        if isinstance(fls, list) and fls[list_index] is not None:
            ListOfImagesDictionaries[list_index]['FL'] = fls[list_index]
        if labels[list_index] is not None:
            ListOfImagesDictionaries[list_index]['BRAINMASK'] = labels[list_index]
        print ListOfImagesDictionaries[list_index]
        for key, value in posteriors.items():
            # print "key:", key, " -> value:", value
            ListOfImagesDictionaries[list_index][key] = value[list_index]

    # print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    # print "ListOfImagesDictionaries", ListOfImagesDictionaries
    # print "registrationImageTypes", registrationImageTypes
    # print "interpolationMapping", interpolationMapping
    return ListOfImagesDictionaries, registrationImageTypes, interpolationMapping


def xml_filename(subject):
    return 'AtlasDefinition_{0}.xml'.format(subject)

def getSessionsFromSubjectDictionary(subject_session_dictionary,subject):
    return subject_session_dictionary[subject]


def _template_runner(argv, environment, experiment, pipeline, cluster):
    print "Getting subjects from database..."
    # subjects = argv["--subjects"].split(',')
    subjects, subjects_sessions_dictionary = get_subjects_sessions_dictionary(argv['SUBJECTS'],
            experiment['cachedir'],
            experiment['resultdir'],
            environment['prefix'],
            experiment['dbfile'],
            argv['--use-sentinal']
            ) # Build database before parallel section
    print("Processing atlas generation for these subjects")
    print(subjects)
    print("="*80)
    print "Copying Atlas directory and determining appropriate Nipype options..."
    pipeline = nipype_options(argv, pipeline, cluster, experiment, environment)  # Generate Nipype options
    print "Dispatching jobs to the system..."

    ######
    ###### Now start workflow construction
    ######
    # Set universal pipeline options
    nipype_config.update_config(pipeline)

    template = pe.Workflow(name='SubjectAtlas_Template')
    template.base_dir = pipeline['logging']['log_directory']

    subjectIterator = pe.Node(interface=IdentityInterface(fields=['subject']), run_without_submitting=True, name='99_subjectIterator')
    subjectIterator.iterables = ('subject', subjects)

    sessionsExtractorNode = pe.Node(Function(function=getSessionsFromSubjectDictionary,
                                                      input_names=['subject_session_dictionary','subject'],
                                                      output_names=['sessions']),
                                   run_without_submitting=True, name="99_sessionsExtractor")
    sessionsExtractorNode.inputs.subject_session_dictionary = subjects_sessions_dictionary

    baselineDG = pe.MapNode(nio.DataGrabber(infields=['subject','session'],
                                            outfields=['t1_average', 't2_average', 'pd_average',
                                                       'fl_average', 'brainMaskLabels',
                                                       'posteriorImages']),
                            iterfield=['session'], name='Baseline_DG')

    baselineDG.inputs.base_directory = experiment['previousresult']
    baselineDG.inputs.sort_filelist = True
    baselineDG.inputs.raise_on_empty = False
    baselineDG.inputs.template = '*'
    posterior_files = ['AIR', 'BASAL', 'CRBLGM', 'CRBLWM', 'CSF', 'GLOBUS', 'HIPPOCAMPUS',
                       'NOTCSF', 'NOTGM', 'NOTVB', 'NOTWM', 'SURFGM', 'THALAMUS', 'VB', 'WM']
    baselineDG.inputs.field_template = {'t1_average':'*/%s/%s/TissueClassify/t1_average_BRAINSABC.nii.gz',
                                        't2_average':'*/%s/%s/TissueClassify/t2_average_BRAINSABC.nii.gz',
                                        'pd_average':'*/%s/%s/TissueClassify/pd_average_BRAINSABC.nii.gz',
                                        'fl_average':'*/%s/%s/TissueClassify/fl_average_BRAINSABC.nii.gz',
                                   'brainMaskLabels':'*/%s/%s/TissueClassify/fixed_brainlabels_seg.nii.gz',
                                   'posteriorImages':'*/%s/%s/TissueClassify/POSTERIOR_%s.nii.gz'
                                   }
    baselineDG.inputs.template_args  = {'t1_average':[['subject','session']],
                                        't2_average':[['subject','session']],
                                        'pd_average':[['subject','session']],
                                        'fl_average':[['subject','session']],
                                   'brainMaskLabels':[['subject','session']],
                                   'posteriorImages':[['subject','session', posterior_files]]
                                   }

    MergeByExtendListElementsNode = pe.Node(Function(function=MergeByExtendListElements,
                                                     input_names=['t1s', 't2s',
                                                                  'pds', 'fls',
                                                                  'labels', 'posteriors'],
                                                     output_names=['ListOfImagesDictionaries', 'registrationImageTypes',
                                                                   'interpolationMapping']),
                                            run_without_submitting=True, name="99_MergeByExtendListElements")

    template.connect([(subjectIterator, baselineDG, [('subject', 'subject')]),
                      (subjectIterator, sessionsExtractorNode, [('subject','subject')]),
                      (sessionsExtractorNode, baselineDG, [('sessions', 'session')]),
                      (baselineDG, MergeByExtendListElementsNode, [('t1_average', 't1s'),
                                                                   ('t2_average', 't2s'),
                                                                   ('pd_average', 'pds'),
                                                                   ('fl_average', 'fls'),
                                                                   ('brainMaskLabels', 'labels'),
                                                                   (('posteriorImages', ConvertSessionsListOfPosteriorListToDictionaryOfSessionLists), 'posteriors')])
                    ])

    myInitAvgWF = pe.Node(interface=ants.AverageImages(), name='Atlas_antsSimpleAverage')  # was 'Phase1_antsSimpleAverage'
    myInitAvgWF.inputs.dimension = 3
    myInitAvgWF.inputs.normalize = True
    template.connect(baselineDG, 't1_average', myInitAvgWF, "images")
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

    if pipeline['plugin_name'].startswith('SGE'):  # for some nodes, the qsub call needs to be modified on the cluster

        CreateAtlasXMLAndCleanedDeformedAveragesNode.plugin_args = {'template': pipeline['plugin_args']['template'],
                                                'qsub_args': modify_qsub_args(cluster['queue'], 1, 1, 1),
                                                'overwrite': True}
        for bt in [buildTemplateIteration1, buildTemplateIteration2]:
            ##################################################
            # *** Hans, is this TODO already addressed? ***  #
            # ---->  # TODO:  Change these parameters  <---- #
            ##################################################
            BeginANTS = bt.get_node("BeginANTS")
            BeginANTS.plugin_args = {'template': pipeline['plugin_args']['template'], 'overwrite': True,
                                     'qsub_args': modify_qsub_args(cluster['queue'], 8, 8, 24)}
            wimtdeformed = bt.get_node("wimtdeformed")
            wimtdeformed.plugin_args = {'template': pipeline['plugin_args']['template'], 'overwrite': True,
                                        'qsub_args': modify_qsub_args(cluster['queue'], 2, 2, 2)}
            AvgAffineTransform = bt.get_node("AvgAffineTransform")
            AvgAffineTransform.plugin_args = {'template': pipeline['plugin_args']['template'], 'overwrite': True,
                                              'qsub_args': modify_qsub_args(cluster['queue'], 2, 1, 1)}
            wimtPassivedeformed = bt.get_node("wimtPassivedeformed")
            wimtPassivedeformed.plugin_args = {'template': pipeline['plugin_args']['template'], 'overwrite': True,
                                                'qsub_args': modify_qsub_args(cluster['queue'], 2, 2, 2)}

    # Running off previous baseline experiment
    NACCommonAtlas = MakeAtlasNode(experiment['atlascache'], 'NACCommonAtlas_{0}'.format('subject'), 'TemplateBuildSupport') ## HACK : replace 'subject' with subject id once this is a loop rather than an iterable.
    template.connect([(myInitAvgWF, buildTemplateIteration1, [('output_average_image', 'inputspec.fixed_image')]),
                      (MergeByExtendListElementsNode, buildTemplateIteration1, [('ListOfImagesDictionaries', 'inputspec.ListOfImagesDictionaries'),
                                                                                ('registrationImageTypes', 'inputspec.registrationImageTypes'),
                                                                                ('interpolationMapping','inputspec.interpolationMapping')]),
                      (buildTemplateIteration1, buildTemplateIteration2, [('outputspec.template', 'inputspec.fixed_image')]),
                      (MergeByExtendListElementsNode, buildTemplateIteration2, [('ListOfImagesDictionaries', 'inputspec.ListOfImagesDictionaries'),
                                                                                ('registrationImageTypes','inputspec.registrationImageTypes'),
                                                                                ('interpolationMapping', 'inputspec.interpolationMapping')]),
                      (subjectIterator, CreateAtlasXMLAndCleanedDeformedAveragesNode, [(('subject', xml_filename), 'outDefinition')]),
                      (NACCommonAtlas, CreateAtlasXMLAndCleanedDeformedAveragesNode, [('ExtendedAtlasDefinition_xml_in', 'AtlasTemplate')]),
                      (buildTemplateIteration2, CreateAtlasXMLAndCleanedDeformedAveragesNode, [('outputspec.template', 't1_image'),
                                                                           ('outputspec.passive_deformed_templates', 'deformed_list')]),
                      ])

    # Create DataSinks
    SubjectAtlas_DataSink = pe.Node(nio.DataSink(), name="Subject_DS")
    SubjectAtlas_DataSink.overwrite = pipeline['ds_overwrite']
    SubjectAtlas_DataSink.inputs.base_directory = experiment['resultdir']

    template.connect([(subjectIterator, SubjectAtlas_DataSink, [('subject', 'container')]),
                      (CreateAtlasXMLAndCleanedDeformedAveragesNode, SubjectAtlas_DataSink, [('outAtlasFullPath', 'Atlas.@definitions')]),
                      (CreateAtlasXMLAndCleanedDeformedAveragesNode, SubjectAtlas_DataSink, [('clean_deformed_list', 'Atlas.@passive_deformed_templates')]),

                      (subjectIterator, SubjectAtlas_DataSink, [(('subject', outputPattern), 'regexp_substitutions')]),
                      (buildTemplateIteration2, SubjectAtlas_DataSink, [('outputspec.template', 'Atlas.@template')]),
                     ])

    dotfilename = argv['--dotfilename']
    if dotfilename is not None:
        print("WARNING: Printing workflow, but not running pipeline")
        print_workflow(template, plugin=pipeline['plugin_name'], dotfilename=dotfilename)
    else:
        run_workflow(template, plugin=pipeline['plugin_name'], plugin_args=pipeline['plugin_args'])

if __name__ == '__main__':
    import sys
    from AutoWorkup import setup

    from docopt import docopt

    argv = docopt(__doc__, version='1.1')
    print argv
    print '=' * 100
    configs = setup(argv)
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

    exit = _template_runner(argv, *configs)
    sys.exit(exit)
