#! /usr/bin/env python
"""
singleSession.py
=========
This program is used to generate the subject- and session-specific workflows for BRAINSTool processing

Usage:
  singleSession.py [--rewrite-datasinks] [--wfrun PLUGIN] --workphase WORKPHASE --pe ENV --ExperimentConfig FILE SESSIONS...
  singleSession.py -v | --version
  singleSession.py -h | --help

Arguments:
  SESSIONS              List of sessions to process. Specifying 'all' processes every session in
                        the database (specified in the --ExperimentConfig FILE)

Options:
  -h, --help            Show this help and exit
  -v, --version         Print the version and exit
  --rewrite-datasinks   Turn on the Nipype option to overwrite all files in the 'results' directory
  --pe ENV              The processing environment to use from configuration file
  --wfrun PLUGIN        The name of the workflow plugin option (default: 'local')
  --workphase WORKPHASE The type of processing to be done [atlas-based-reference|subject-based-reference]
  --ExperimentConfig FILE   The configuration file


Examples:
  $ singleSession.py --pe OSX --ExperimentConfig my_baw.config all
  $ singleSession.py --wfrun SGEGraph --pe OSX --ExperimentConfig my_baw.config 00001 00002
  $ singleSession.py --rewrite-datasinks --pe OSX --ExperimentConfig my_baw.config 00003

"""

import os

def create_singleSession(dataDict, master_config, interpMode, pipeline_name):
    """
    create singleSession workflow on a single session

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """
    assert 'tissue_classify' in master_config['components'] or \
      'auxlmk' in master_config['components'] or \
      'segmentation' in master_config['components']

    from nipype import config, logging
    config.update_config(master_config)  # Set universal pipeline options
    logging.update_logging(config)

    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, Directory, traits, isdefined, BaseInterface
    from nipype.interfaces.utility import Split, Rename, IdentityInterface, Function

    from workflows.baseline import generate_single_session_template_WF
    from PipeLineFunctionHelpers import convertToList
    from utilities.misc import GenerateSubjectOutputPattern as outputPattern
    from utilities.misc import GenerateWFName
    from workflows.utils import run_workflow, print_workflow
    from workflows.atlasNode import MakeAtlasNode

    project = dataDict['project']
    subject = dataDict['subject']
    session = dataDict['session']

    pname = "{0}_{1}_{2}".format(master_config['workflow_phase'], subject, session)
    sessionWorkflow = generate_single_session_template_WF(project, subject, session, master_config,
                             phase=master_config['workflow_phase'],
                             interpMode=interpMode,
                             pipeline_name=pipeline_name)
    sessionWorkflow.base_dir = master_config['cachedir']

    SSinputsSpecPtr = sessionWorkflow.get_node('inputspec')
    SSinputsSpecPtr.inputs.T1s = dataDict['T1s']
    SSinputsSpecPtr.inputs.T2s = dataDict['T2s']
    SSinputsSpecPtr.inputs.PDs = dataDict['PDs']
    SSinputsSpecPtr.inputs.FLs = dataDict['FLs']
    SSinputsSpecPtr.inputs.OTHERs = dataDict['OTs']
    atlasBCDNode = MakeAtlasNode(master_config['atlascache'], 'BBCDAtlas_{0}'.format(session),['BCDSupport'])
    sessionWorkflow.connect([(atlasBCDNode, SSinputsSpecPtr, [('template_t1', 'template_t1'),
                                                      ('template_landmarks_50Lmks_fcsv', 'atlasLandmarkFilename'),
                                                      ('template_weights_50Lmks_wts', 'atlasWeightFilename'),
                                                      ('LLSModel_50Lmks_hdf5', 'LLSModel'),
                                                      ('T1_50Lmks_mdl', 'inputTemplateModel')]),
                            ])
    if  master_config['workflow_phase'] == 'atlas-based-reference':
        atlasABCNode = MakeAtlasNode(master_config['atlascache'], 'BABCAtlas_{0}'.format(session),['BRAINSABCSupport'])  # TODO: input atlas csv

        sessionWorkflow.connect([(atlasABCNode, SSinputsSpecPtr,
                                    [ ('ExtendedAtlasDefinition_xml', 'atlasDefinition') ]
                                 ),
                               ])
    elif master_config['workflow_phase'] == 'subject-based-reference':
        print master_config['previousresult']
        template_DG = pe.Node(interface=nio.DataGrabber(infields=['subject'],
                                                        outfields=['outAtlasXMLFullPath']),
                              name='Template_DG')
        template_DG.inputs.base_directory = master_config['previousresult']
        template_DG.inputs.subject = subject
        template_DG.inputs.template = '%s/Atlas/AtlasDefinition_%s.xml'
        template_DG.inputs.template_args['outAtlasXMLFullPath'] = [['subject','subject']]
        template_DG.inputs.sort_filelist = True
        template_DG.inputs.raise_on_empty = True

        sessionWorkflow.connect([(template_DG, SSinputsSpecPtr,
                                   [('outAtlasXMLFullPath', 'atlasDefinition')
                                   ]),
                                 ])
    else:
        assert 0 == 1 , "Invalid workflow type specified for singleSession"

    if 'segmentation' in master_config['components']:
        from workflows.segmentation import segmentation
        from workflows.WorkupT1T2BRAINSCut import GenerateWFName
        try:
            bCutInputName = ".".join([GenerateWFName(project, subject, session, 'Segmentation'), 'inputspec'])
        except:
            print project, subject, session
            raise
        sname = 'segmentation'
        onlyT1 = not(len(dataDict['T2s']) > 0)
        atlasBCUTNode = MakeAtlasNode(master_config['atlascache'], 'BBCUTAtlas_{0}'.format(session),['BRAINSCutSupport'])
        segWF = segmentation(project, subject, session, master_config, onlyT1, pipeline_name=sname)
        sessionWorkflow.connect([(atlasBCUTNode, segWF,
                                [('hncma-atlas', 'inputspec.hncma-atlas'),
                                 ('template_t1', 'inputspec.template_t1'),
                                 ('template_t1', bCutInputName + '.template_t1'),
                                 ('rho', bCutInputName + '.rho'),
                                 ('phi', bCutInputName + '.phi'),
                                 ('theta', bCutInputName + '.theta'),
                                 ('l_caudate_ProbabilityMap', bCutInputName + '.l_caudate_ProbabilityMap'),
                                 ('r_caudate_ProbabilityMap', bCutInputName + '.r_caudate_ProbabilityMap'),
                                 ('l_hippocampus_ProbabilityMap', bCutInputName + '.l_hippocampus_ProbabilityMap'),
                                 ('r_hippocampus_ProbabilityMap', bCutInputName + '.r_hippocampus_ProbabilityMap'),
                                 ('l_putamen_ProbabilityMap', bCutInputName + '.l_putamen_ProbabilityMap'),
                                 ('r_putamen_ProbabilityMap', bCutInputName + '.r_putamen_ProbabilityMap'),
                                 ('l_thalamus_ProbabilityMap', bCutInputName + '.l_thalamus_ProbabilityMap'),
                                 ('r_thalamus_ProbabilityMap', bCutInputName + '.r_thalamus_ProbabilityMap'),
                                 ('l_accumben_ProbabilityMap', bCutInputName + '.l_accumben_ProbabilityMap'),
                                 ('r_accumben_ProbabilityMap', bCutInputName + '.r_accumben_ProbabilityMap'),
                                 ('l_globus_ProbabilityMap', bCutInputName + '.l_globus_ProbabilityMap'),
                                 ('r_globus_ProbabilityMap', bCutInputName + '.r_globus_ProbabilityMap'),
                                 ('trainModelFile_txtD0060NT0060_gz', bCutInputName + '.trainModelFile_txtD0060NT0060_gz')])])
        outputSpec = sessionWorkflow.get_node('outputspec')
        sessionWorkflow.connect([(outputSpec, segWF, [('t1_average', 'inputspec.t1_average'),
                                             ('LMIatlasToSubject_tx', 'inputspec.LMIatlasToSubject_tx'),
                                             ('outputLabels', 'inputspec.inputLabels'),
                                             ('posteriorImages', 'inputspec.posteriorImages'),
                                             ('tc_atlas2sessionInverse_tx',
                                              'inputspec.TissueClassifyatlasToSubjectInverseTransform'),
                                             ('UpdatedPosteriorsList', 'inputspec.UpdatedPosteriorsList'),
                                             ('outputHeadLabels', 'inputspec.inputHeadLabels')])
                                ])
        if not onlyT1:
            sessionWorkflow.connect([(outputSpec, segWF, [('t1_average', 'inputspec.t2_average')])])

    return sessionWorkflow

def createAndRun(sessions, environment, experiment, pipeline, cluster):
    from baw_exp import OpenSubjectDatabase
    from utilities.misc import add_dict
    from workflows.utils import run_workflow, print_workflow
    master_config = {}
    for configDict in [environment, experiment, pipeline, cluster]:
        master_config = add_dict(master_config, configDict)
    database = OpenSubjectDatabase(experiment['cachedir'], ['all'], environment['prefix'], experiment['dbfile'])
    database.open_connection()
    try:
        all_sessions = database.getAllSessions()
        if not set(sessions) <= set(all_sessions) and 'all' not in sessions:
            missing = set(sessions) - set(all_sessions)
            assert len(missing) == 0, "Requested sessions are missing from the database: {0}".format(missing)
        elif 'all' in sessions:
            sessions = set(all_sessions)
        else:
            sessions = set(sessions)
        print "!="*40
        print("Doing sessions {0}".format(sessions))
        print "!="*40
        for session in sessions:
            _dict = {}
            _dict['session'] = session
            _dict['project'] = database.getProjFromSession(session)
            _dict['subject'] = database.getSubjFromSession(session)
            _dict['T1s'] = database.getFilenamesByScantype(session, ['T1-15', 'T1-30'])
            _dict['T2s'] = database.getFilenamesByScantype(session, ['T2-15', 'T2-30'])
            _dict['PDs'] = database.getFilenamesByScantype(session, ['PD-15', 'PD-30'])
            _dict['FLs'] = database.getFilenamesByScantype(session, ['FL-15', 'FL-30'])
            _dict['OTs'] = database.getFilenamesByScantype(session, ['OTHER-15', 'OTHER-30'])
            sentinal_file = os.path.join(
                                        master_config['resultdir'],
                                        _dict['project'],
                                        _dict['subject'],
                                        _dict['session'],
                                        "TissueClassify",
                                        "t1_average_BRAINSABC.nii.gz"
                                        )

            if os.path.exists( sentinal_file ):
                print("SKIPPING: {0} exists".format( sentinal_file) )
            else:
                workflow = create_singleSession(_dict, master_config, 'Linear', 'singleSession_{0}_{1}'.format(_dict['subject'], _dict['session']))
                print("Starting session {0}".format(session))
                # HACK Hard-coded to SGEGraph, but --wfrun is ignored completely
                run_workflow(workflow,plugin=master_config['plugin_name'], plugin_args=master_config['plugin_args'])
    except:
        raise
    finally:
        try:
            database.close_connection()
        except:
            pass

def _main(environment, experiment, pipeline, cluster, **kwds):
    from utilities.configFileParser import nipype_options
    from utilities.misc import add_dict

    print "Copying Atlas directory and determining appropriate Nipype options..."
    pipeline = nipype_options(kwds, pipeline, cluster, experiment, environment)  # Generate Nipype options
    print "Getting session(s) from database..."
    createAndRun(kwds['SESSIONS'], environment, experiment, pipeline, cluster)
    return 0

if __name__ == '__main__':
    import sys
    from docopt import docopt
    from AutoWorkup import setup

    argv = docopt(__doc__, version='1.1')
    print argv
    print '=' * 100
    configs = setup(argv)
    exit = _main(*configs, **argv)
    sys.exit(exit)
