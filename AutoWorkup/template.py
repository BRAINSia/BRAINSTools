#! /usr/bin/env python
"""
template.py
=========
This program is used to generate the subject- and session-specific workflows for BRAINSTool processing

Usage:
  template.py [--rewrite-datasinks] [--wfrun PLUGIN] [--dotfilename PFILE] --subjects ID --pe ENV --ExperimentConfig FILE
  template.py -v | --version
  template.py -h | --help

Arguments:


Options:
  -h, --help            Show this help and exit
  -v, --version         Print the version and exit
  --dotfilename PFILE   Turn on printing pipeline to file PFILE
  --rewrite-datasinks   Turn on the Nipype option to overwrite all files in the 'results' directory
  --pe ENV              The processing environment to use from configuration file
  --subjects ID1,ID2    The comma separated list of subject IDs to process
  --wfrun PLUGIN        The name of the workflow plugin option (default: 'local')
  --ExperimentConfig FILE   The configuration file

Examples:
  $ template.py --subjects 1058,1059 --pe OSX --ExperimentConfig my_baw.config
  $ template.py --wfrun helium_all.q --subjects 1058 --pe OSX --ExperimentConfig my_baw.config
  $ template.py --rewrite-datasinks --subjects 1058 --pe OSX --ExperimentConfig my_baw.config

"""
import os
import sys
import traceback

import nipype.config
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
from nipype.interfaces.utility import IdentityInterface, Function
import nipype.interfaces.ants as ants

from template import MergeByExtendListElements, xml_filename
from PipeLineFunctionHelpers import mapPosteriorList
from workflows.atlasNode import GetAtlasNode, MakeNewAtlasTemplate
from utilities.misc import GenerateSubjectOutputPattern as outputPattern
from utilities.distributed import modify_qsub_args

def MergeByExtendListElements(t1s, t2s, pds, fls, labels, posteriors):
    """
    *** NOTE:  ALl input lists MUST have the same number of elements (even if they are null) ***

    output = [{'T1': os.path.join(mydatadir, '01_T1_half.nii.gz'),
                               'INV_T1': os.path.join(mydatadir, '01_T1_inv_half.nii.gz'),
                               'LABEL_MAP': os.path.join(mydatadir, '01_T1_inv_half.nii.gz')
                              },
                              {'T1': os.path.join(mydatadir, '02_T1_half.nii.gz'),
                               'INV_T1': os.path.join(mydatadir, '02_T1_inv_half.nii.gz'),
                               'LABEL_MAP': os.path.join(mydatadir, '02_T1_inv_half.nii.gz')
                              },
                              {'T1': os.path.join(mydatadir, '03_T1_half.nii.gz'),
                               'INV_T1': os.path.join(mydatadir, '03_T1_inv_half.nii.gz'),
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


def main(*args, dotfilename=None):
    subjects, environment, experiment, pipeline, cluster = args
    # Set universal pipeline options
    nipype.config.update_config(pipeline)
    assert nipype.config.get('execution', 'plugin') == pipeline['execution']['plugin']

    template = pe.Workflow(name='SubjectAtlas_Template')
    template.base_dir = pipeline['logging']['log_directory']

    if 'previouscache' in experiment:
        # Running off previous baseline experiment
        BAtlas = GetAtlasNode(experiment['previouscache'], 'BAtlas')
    else:
        # Running after previous baseline experiment
        raise "ERROR, Never use the NAC atlas for template building"
        BAtlas = GetAtlasNode(os.path.dirname(experiment['atlascache']), 'BAtlas')
    inputspec = pe.Node(interface=IdentityInterface(fields=['subject']), name='inputspec')
    inputspec.iterables = ('subject', subjects)

    baselineDG = pe.Node(nio.DataGrabber(infields=['subject'], outfields=['t1_average', 't2_average', 'pd_average',
                                                                            'fl_average', 'outputLabels', 'posteriorImages']),
                         name='Baseline_DG')
    if 'previousresult' in experiment:
        baselineDG.inputs.base_directory = experiment['previousresult']
    else:
        baselineDG.inputs.base_directory = experiment['resultdir']
    baselineDG.inputs.sort_filelist = True
    baselineDG.inputs.raise_on_empty = False
    baselineDG.inputs.template = '*/%s/*/Baseline/%s.nii.gz'
    baselineDG.inputs.template_args['t1_average'] = [['subject', 't1_average_BRAINSABC']]
    baselineDG.inputs.template_args['t2_average'] = [['subject', 't2_average_BRAINSABC']]
    baselineDG.inputs.template_args['pd_average'] = [['subject', 'pd_average_BRAINSABC']]
    baselineDG.inputs.template_args['fl_average'] = [['subject', 'fl_average_BRAINSABC']]
    baselineDG.inputs.template_args['outputLabels'] = [['subject', 'brain_label_seg']]
    baselineDG.inputs.field_template = {'posteriorImages':'*/%s/*/TissueClassify/POSTERIOR_%s.nii.gz'}
    posterior_files = ['AIR', 'BASAL', 'CRBLGM', 'CRBLWM', 'CSF', 'GLOBUS', 'HIPPOCAMPUS', 'NOTCSF', 'NOTGM', 'NOTVB', 'NOTWM',
                       'SURFGM', 'THALAMUS', 'VB', 'WM']
    baselineDG.inputs.template_args['posteriorImages'] = [['subject', posterior_files]]

    MergeByExtendListElementsNode = pe.Node(Function(function=MergeByExtendListElements,
                                                     input_names=['t1s', 't2s',
                                                                  'pds', 'fls',
                                                                  'labels', 'posteriors'],
                                                     output_names=['ListOfImagesDictionaries', 'registrationImageTypes',
                                                                   'interpolationMapping']),
                                            run_without_submitting=True, name="99_MergeByExtendListElements")
    from PipeLineFunctionHelpers import WrapPosteriorImagesFromDictionaryFunction as wrapfunc
    template.connect([(inputspec, baselineDG, [('subject', 'subject')]),
                      (baselineDG, MergeByExtendListElementsNode, [('t1_average', 't1s'),
                                                                   ('t2_average', 't2s'),
                                                                   ('pd_average', 'pds'),
                                                                   ('fl_average', 'fls'),
                                                                   ('outputLabels', 'labels'),
                                                                   (('posteriorImages', wrapfunc), 'posteriors')])
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
    from BAWantsRegistrationBuildTemplate import BAWantsRegistrationTemplateBuildSingleIterationWF as registrationWF
    buildTemplateIteration1 = registrationWF('iteration01')
    # buildTemplateIteration2 = buildTemplateIteration1.clone(name='buildTemplateIteration2')
    buildTemplateIteration2 = registrationWF('Iteration02')

    MakeNewAtlasTemplateNode = pe.Node(interface=Function(function=MakeNewAtlasTemplate,
                                                          input_names=['t1_image', 'deformed_list', 'AtlasTemplate', 'outDefinition'],
                                                          output_names=['outAtlasFullPath', 'clean_deformed_list']),
                                       # This is a lot of work, so submit it run_without_submitting=True,
                                       run_without_submitting=True,  # HACK:  THIS NODE REALLY SHOULD RUN ON THE CLUSTER!
                                       name='99_MakeNewAtlasTemplate')

    if pipeline['execution']['plugin'].startswith('SGE'):  # for some nodes, the qsub call needs to be modified on the cluster

        MakeNewAtlasTemplateNode.plugin_args = {'template': pipeline['plugin_args']['template'],
                                                'qsub_args': modify_qsub_args(cluster['queue'], '1000M', 1, 1),
                                                'overwrite': True}
        for bt in [buildTemplateIteration1, buildTemplateIteration2]:
            ##################################################
            # *** Hans, is this TODO already addressed? ***  #
            # ---->  # TODO:  Change these parameters  <---- #
            ##################################################
            BeginANTS = bt.get_node("BeginANTS")
            BeginANTS.plugin_args = {'template': pipeline['plugin_args']['template'], 'overwrite': True,
                                     'qsub_args': modify_qsub_args(cluster['queue'], '9000M', 4, hard=False)}
            wimtdeformed = bt.get_node("wimtdeformed")
            wimtdeformed.plugin_args = {'template': pipeline['plugin_args']['template'], 'overwrite': True,
                                        'qsub_args': modify_qsub_args(cluster['queue'], '2000M', 1, 2)}
            AvgAffineTransform = bt.get_node("AvgAffineTransform")
            AvgAffineTransform.plugin_args = {'template': pipeline['plugin_args']['template'], 'overwrite': True,
                                              'qsub_args': modify_qsub_args(cluster['queue'], '2000M', 1)}
            wimtPassivedeformed = bt.get_node("wimtPassivedeformed")
            wimtPassivedeformed.plugin_args = {'template': pipeline['plugin_args']['template'], 'overwrite': True,
                                                'qsub_args': modify_qsub_args(cluster['queue'], '2000M', 1, 2)}

    template.connect([(myInitAvgWF, buildTemplateIteration1, [('output_average_image', 'inputspec.fixed_image')]),
                      (MergeByExtendListElementsNode, buildTemplateIteration1, [('ListOfImagesDictionaries', 'inputspec.ListOfImagesDictionaries'),
                                                                                ('registrationImageTypes', 'inputspec.registrationImageTypes'),
                                                                                ('interpolationMapping','inputspec.interpolationMapping')]),
                      (buildTemplateIteration1, buildTemplateIteration2, [('outputspec.template', 'inputspec.fixed_image')]),
                      (MergeByExtendListElementsNode, buildTemplateIteration2, [('ListOfImagesDictionaries', 'inputspec.ListOfImagesDictionaries'),
                                                                                ('registrationImageTypes','inputspec.registrationImageTypes'),
                                                                                ('interpolationMapping', 'inputspec.interpolationMapping')]),
                      (inputspec, MakeNewAtlasTemplateNode, [(('subject', xml_filename), 'outDefinition')]),
                      (BAtlas, MakeNewAtlasTemplateNode, [('ExtendedAtlasDefinition_xml_in', 'AtlasTemplate')]),
                      (buildTemplateIteration2, MakeNewAtlasTemplateNode, [('outputspec.template', 't1_image'),
                                                                           ('outputspec.passive_deformed_templates', 'deformed_list')]),
                      ])

    # Create DataSinks
    Atlas_DataSink = pe.Node(nio.DataSink(), name="Atlas_DS")
    Atlas_DataSink.overwrite = pipeline['ds_overwrite']
    Atlas_DataSink.inputs.base_directory = experiment['resultdir']

    Subject_DataSink = pe.Node(nio.DataSink(), name="Subject_DS")
    Subject_DataSink.overwrite = pipeline['ds_overwrite']
    Subject_DataSink.inputs.base_directory = experiment['resultdir']

    template.connect([(inputspec, Atlas_DataSink, [('subject', 'container')]),
                      (buildTemplateIteration1, Atlas_DataSink, [('outputspec.template', 'Atlas.iteration1')]),  # Unnecessary
                      (MakeNewAtlasTemplateNode, Atlas_DataSink, [('outAtlasFullPath', 'Atlas.definitions')]),
                      (BAtlas, Atlas_DataSink, [('template_landmarks_50Lmks_fcsv', 'Atlas.20111119_BCD.@fcsv'),
                                                ('template_weights_50Lmks_wts', 'Atlas.20111119_BCD.@wts'),
                                                ('LLSModel_50Lmks_hdf5', 'Atlas.20111119_BCD.@hdf5'),
                                                ('T1_50Lmks_mdl', 'Atlas.20111119_BCD.@mdl')]),
                      (inputspec, Subject_DataSink, [(('subject', outputPattern), 'regexp_substitutions')]),
                      (buildTemplateIteration2, Subject_DataSink, [('outputspec.template', 'ANTSTemplate.@template')]),
                      (MakeNewAtlasTemplateNode, Subject_DataSink, [('clean_deformed_list', 'ANTSTemplate.@passive_deformed_templates')]),
                     ])

    from workflows.utils import run_workflow, print_workflow
    if dotfilename is not None:
        return print_workflow(template, plugin=pipeline['execution']['plugin'], dotfilename=dotfilename)
    return run_workflow(template, plugin=pipeline['execution']['plugin'], plugin_args=pipeline['plugin_args'])

def _template_runner(argv, environment, experiment, pipeline, cluster):
    from utilities.configFileParser import nipype_options
    from utilities.misc import add_dict
    from AutoWorkup import get_subjects

    print "Getting subjects from database..."
    subjects = get_subjects(argv, experiment['cachedir'], environment['prefix'], experiment['dbfile']) # Build database before parallel section
    print "Copying Atlas directory and determining appropriate Nipype options..."
    pipeline = nipype_options(argv, pipeline, cluster, experiment, environment)  # Generate Nipype options
    print "Dispatching jobs to the system..."
    return main(subjects, environment, experiment, pipeline, cluster, argv['--dotfilename'])

if __name__ == '__main__':
    import sys
    from AutoWorkup import setup

    from docopt import docopt

    argv = docopt(__doc__, version='1.1')
    print argv
    print '=' * 100
    configs = setup(argv)
    exit = _template_runner(argv, *configs)
    sys.exit(exit)
