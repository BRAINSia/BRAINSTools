#! /usr/bin/env python
"""
singleSession.py
=========
This program is used to generate the subject- and session-specific workflows for BRAINSTool processing

Usage:
  singleSession.py [--rewrite-datasinks] [--wfrun PLUGIN] [--use-sentinal] [--dry-run] --workphase WORKPHASE --pe ENV --ExperimentConfig FILE SESSIONS...
  singleSession.py -v | --version
  singleSession.py -h | --help

Arguments:
  SESSIONS              List of sessions to process. Specifying 'all' processes every session in
                        the database (specified in the --ExperimentConfig FILE)

Options:
  -h, --help            Show this help and exit
  -v, --version         Print the version and exit
  --rewrite-datasinks   Turn on the Nipype option to overwrite all files in the 'results' directory
  --use-sentinal        Use the t1_average file as a marker to determine if session needs to be run
  --dry-run             Do not submit jobs, but print diagnostics about which jobs would be run
  --pe ENV              The processing environment to use from configuration file
  --wfrun PLUGIN        The name of the workflow plugin option (default: 'local')
  --workphase WORKPHASE The type of processing to be done [atlas-based-reference|subject-based-reference]
  --ExperimentConfig FILE   The configuration file


Examples:
  $ singleSession.py --pe OSX --ExperimentConfig my_baw.config all
  $ singleSession.py --use-sentinal --wfrun SGEGraph --pe OSX --ExperimentConfig my_baw.config 00001 00002
  $ singleSession.py --use-sentinal --dry-run --wfrun SGEGraph --pe OSX --ExperimentConfig my_baw.config 00001 00002
  $ singleSession.py --rewrite-datasinks --pe OSX --ExperimentConfig my_baw.config 00003

"""


def _DetermineIfSegmentationShouldBeDone(master_config):
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
    return  do_BRAINSCut_Segmentation

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
        'denoise' in master_config['components'] or \
        'landmark' in master_config['components'] or \
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
    from workflows.atlasNode import MakeAtlasNode

    project = dataDict['project']
    subject = dataDict['subject']
    session = dataDict['session']

    blackListFileName = dataDict['T1s'][0] + '_noDenoise'
    isBlackList=os.path.isfile( blackListFileName )

    pname = "{0}_{1}_{2}".format(master_config['workflow_phase'], subject, session)
    sessionWorkflow = generate_single_session_template_WF(project, subject, session, master_config,
                                                          phase=master_config['workflow_phase'],
                                                          interpMode=interpMode,
                                                          pipeline_name=pipeline_name,
                                                          doDenoise=(not isBlackList))
    sessionWorkflow.base_dir = master_config['cachedir']

    sessionWorkflow_inputsspec = sessionWorkflow.get_node('inputspec')
    sessionWorkflow_inputsspec.inputs.T1s = dataDict['T1s']
    sessionWorkflow_inputsspec.inputs.T2s = dataDict['T2s']
    sessionWorkflow_inputsspec.inputs.PDs = dataDict['PDs']
    sessionWorkflow_inputsspec.inputs.FLs = dataDict['FLs']
    sessionWorkflow_inputsspec.inputs.OTHERs = dataDict['OTs']
    atlasBCDNode = MakeAtlasNode(master_config['atlascache'], 'BBCDAtlas_{0}'.format(session), ['BCDSupport'])
    sessionWorkflow.connect([(atlasBCDNode, sessionWorkflow_inputsspec,
                              [('template_t1', 'template_t1'),
                               ('template_landmarks_50Lmks_fcsv','atlasLandmarkFilename'),
                               ('template_weights_50Lmks_wts', 'atlasWeightFilename'),
                               ('LLSModel_50Lmks_h5', 'LLSModel'),
                               ('T1_50Lmks_mdl', 'inputTemplateModel')
                               ]),
                             ])

    if master_config['workflow_phase'] == 'atlas-based-reference':
        # TODO: input atlas csv
        atlasABCNode = MakeAtlasNode(master_config['atlascache'], 'BABCAtlas_{0}'.format(session), ['BRAINSABCSupport','LabelMapsSupport'])
        sessionWorkflow.connect(atlasABCNode,'ExtendedAtlasDefinition_xml',sessionWorkflow_inputsspec,'atlasDefinition')
        sessionWorkflow.connect([( atlasABCNode,sessionWorkflow_inputsspec, [
                                                ('hncma_atlas','hncma_atlas'),
                                                ('template_leftHemisphere','template_leftHemisphere'),
                                                ('template_rightHemisphere','template_rightHemisphere'),
                                                ('template_WMPM2_labels','template_WMPM2_labels'),
                                                ('template_nac_labels','template_nac_labels'),
                                                ('template_ventricles','template_ventricles') ]
                                 )]
                                )

    elif master_config['workflow_phase'] == 'subject-based-reference':
        print master_config['previousresult']
        template_DG = pe.Node(interface=nio.DataGrabber(infields=['subject'],
                                                        outfields=['outAtlasXMLFullPath']),
                              name='Template_DG')
        template_DG.inputs.base_directory = master_config['previousresult']
        template_DG.inputs.subject = subject
        template_DG.inputs.template = '%s/Atlas/AtlasDefinition_%s.xml'
        template_DG.inputs.template_args['outAtlasXMLFullPath'] = [['subject', 'subject']]
        template_DG.inputs.sort_filelist = True
        template_DG.inputs.raise_on_empty = True

        sessionWorkflow.connect([(template_DG, sessionWorkflow_inputsspec,
                                  [('outAtlasXMLFullPath', 'atlasDefinition')
                                   ]),
                                 ])
    else:
        assert 0 == 1, "Invalid workflow type specified for singleSession"

    do_BRAINSCut_Segmentation = _DetermineIfSegmentationShouldBeDone(master_config)
    if do_BRAINSCut_Segmentation:
        from workflows.segmentation import segmentation
        from workflows.WorkupT1T2BRAINSCut import GenerateWFName
        try:
            bCutInputName = ".".join([GenerateWFName(project, subject, session, 'Segmentation'), 'inputspec'])
        except:
            print project, subject, session
            raise
        sname = 'segmentation'
        onlyT1 = not(len(dataDict['T2s']) > 0)
        atlasBCUTNode = MakeAtlasNode(master_config['atlascache'],
                                      'BBCUTAtlas_{0}'.format(session), ['BRAINSCutSupport'])
        segWF = segmentation(project, subject, session, master_config, onlyT1, pipeline_name=sname)
        ##TODO: sessionWorkflow.connect(WsegWF)
        sessionWorkflow.connect([(atlasBCUTNode, segWF,
                                [
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
                                                      ('atlasToSubjectRegistrationState','inputspec.atlasToSubjectRegistrationState'),
                                                      ('outputLabels', 'inputspec.inputLabels'),
                                                      ('posteriorImages', 'inputspec.posteriorImages'),
                                                      ('UpdatedPosteriorsList', 'inputspec.UpdatedPosteriorsList'),
                                                      ('outputHeadLabels', 'inputspec.inputHeadLabels')

                                 ])
                                 ])
        if not onlyT1:
            sessionWorkflow.connect([(outputSpec, segWF, [('t2_average', 'inputspec.t2_average')])])

    return sessionWorkflow


def createAndRun(sessions, environment, experiment, pipeline, cluster, useSentinal, dryRun):
    from utilities.misc import add_dict
    from workflows.utils import run_workflow, print_workflow
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
        print "!=" * 40
        print("Doing sessions {0}".format(sessions))
        print "!=" * 40
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
            sentinal_file_basedir = os.path.join(
                    master_config['resultdir'],
                    _dict['project'],
                    _dict['subject'],
                    _dict['session']
            )

            sentinal_file = os.path.join( sentinal_file_basedir ) ## NO SENTINAL FILE
            if 'denoise' in master_config['components']:
                pass

            ## Use t1 average sentinal file if  specified.
            if 'landmark' in master_config['components']:
                sentinal_file = os.path.join(
                    sentinal_file_basedir,
                    "ACPCAlign",
                    "landmarkInitializer_atlas_to_subject_transform.h5"
                )

            if 'tissue_classify' in master_config['components']:
                sentinal_file = os.path.join(
                    sentinal_file_basedir,
                    "TissueClassify",
                    "t1_average_BRAINSABC.nii.gz"
                )

            if 'warp_atlas_to_subject' in master_config['components']:
                sentinal_file = os.path.join(
                    sentinal_file_basedir,
                    "WarpedAtlas2Subject",
                    "template_rightHemisphere.nii.gz"
                )

            ## Use different sentinal file if segmentation specified.
            do_BRAINSCut_Segmentation = _DetermineIfSegmentationShouldBeDone(master_config)
            if do_BRAINSCut_Segmentation:
                sentinal_file = os.path.join(
                    sentinal_file_basedir,
                    "CleanedDenoisedRFSegmentations",
                    "allLabels_seg.nii.gz"
                )
            print "#" * 50
            print sentinal_file + " exists? " + str(os.path.exists(sentinal_file))
            print "-" * 50
            if useSentinal and os.path.exists(sentinal_file):
                print("SKIPPING: {0} exists".format(sentinal_file))
            else:
                print("PROCESSING INCOMPLETE:  {0} does not exists".format(sentinal_file))
                if dryRun == False:
                    workflow = create_singleSession(_dict, master_config, 'Linear',
                                                'singleSession_{0}_{1}'.format(_dict['subject'], _dict['session']))
                    print("Starting session {0}".format(session))
                    # HACK Hard-coded to SGEGraph, but --wfrun is ignored completely
                    run_workflow(workflow, plugin=master_config['plugin_name'], plugin_args=master_config['plugin_args'])
                else:
                    print("EXITING WITHOUT WORK DUE TO dryRun flag")
    except:
        raise
    finally:
        try:
            database.close_connection()
        except:
            pass


def _SingleSession_main(environment, experiment, pipeline, cluster, **kwds):

    from utilities.configFileParser import nipype_options


    print "Copying Atlas directory and determining appropriate Nipype options..."
    pipeline = nipype_options(kwds, pipeline, cluster, experiment, environment)  # Generate Nipype options
    print "Getting session(s) from database..."
    createAndRun(kwds['SESSIONS'], environment, experiment, pipeline, cluster, useSentinal=kwds['--use-sentinal'],dryRun=kwds['--dry-run'])
    return 0


######################################
# Set up the environment, process command line options, and start processing
#
if __name__ == '__main__':
    import sys
    import os

    from docopt import docopt
    from AutoWorkup import setup_environment

    argv = docopt(__doc__, version='1.1')
    print argv
    print '=' * 100
    environment, experiment, pipeline, cluster = setup_environment(argv)


    exit = _SingleSession_main(environment, experiment, pipeline, cluster, **argv)
    sys.exit(exit)
