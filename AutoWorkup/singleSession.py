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


def _create_singleSession(dataDict, master_config, interpMode, pipeline_name):
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
           'segmentation' in master_config['components'] or \
           'malf_2012_neuro' in master_config['components']

    from nipype import config, logging

    config.update_config(master_config)  # Set universal pipeline options
    logging.update_logging(config)

    from workflows.baseline import generate_single_session_template_WF

    project = dataDict['project']
    subject = dataDict['subject']
    session = dataDict['session']

    blackListFileName = dataDict['T1s'][0] + '_noDenoise'
    isBlackList = os.path.isfile(blackListFileName)

    pname = "{0}_{1}_{2}".format(master_config['workflow_phase'], subject, session)
    onlyT1 = not (len(dataDict['T2s']) > 0)
    sessionWorkflow = generate_single_session_template_WF(project, subject, session, onlyT1, master_config,
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
    return sessionWorkflow


def createAndRun(sessions, environment, experiment, pipeline, cluster, useSentinal, dryRun):
    from baw_exp import OpenSubjectDatabase
    from utilities.misc import add_dict

    from workflows.utils import run_workflow

    master_config = {}
    for configDict in [environment, experiment, pipeline, cluster]:
        master_config = add_dict(master_config, configDict)
    database = OpenSubjectDatabase(experiment['cachedir'], ['all'], environment['prefix'], experiment['dbfile'])
    database.open_connection()
    try:
        all_sessions = database.getAllSessions()
        if not set(sessions) <= set(all_sessions) and 'all' not in sessions:
            missing = set(sessions) - set(all_sessions)
            assert len(missing) == 0, "Requested sessions are missing from the database: {0}\n\n{1}".format(missing,all_sessions)
        elif 'all' in sessions:
            sessions = set(all_sessions)
        else:
            sessions = set(sessions)
        print "!=" * 40
        print("Doing sessions {0}".format(sessions))
        print "!=" * 40
        for session in sessions:
            _dict = {}
            subject = database.getSubjFromSession(session)
            _dict['session'] = session
            _dict['project'] = database.getProjFromSession(session)
            _dict['subject'] = subject
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

            sentinal_file_list = list()
            sentinal_file_list.append(os.path.join(sentinal_file_basedir))
            if 'denoise' in master_config['components']:
                # # NO SENTINAL FILE
                pass

            # # Use t1 average sentinal file if  specified.
            if 'landmark' in master_config['components']:
                sentinal_file_list.append(os.path.join(
                    sentinal_file_basedir,
                    "ACPCAlign",
                    "landmarkInitializer_atlas_to_subject_transform.h5"
                ))

            if 'tissue_classify' in master_config['components']:
                for tc_file in ["complete_brainlabels_seg.nii.gz", "t1_average_BRAINSABC.nii.gz"]:
                    sentinal_file_list.append(os.path.join(
                        sentinal_file_basedir,
                        "TissueClassify",
                        tc_file
                    ))

            if 'warp_atlas_to_subject' in master_config['components']:
                sentinal_file_list.append(os.path.join(
                    sentinal_file_basedir,
                    "WarpedAtlas2Subject",
                    "rho.nii.gz"
                ))
                sentinal_file_list.append(os.path.join(
                    sentinal_file_basedir,
                    "WarpedAtlas2Subject",
                    "left_hemisphere_wm.nii.gz"
                ))
            if 'malf_2012_neuro' in master_config['components']:
                sentinal_file_list.append(os.path.join(
                    sentinal_file_basedir,
                    "TissueClassify",
                    "neuro2012_20fusion_merge_seg.nii.gz"
                ))
                sentinal_file_list.append(os.path.join(
                    sentinal_file_basedir,
                    "TissueClassify",
                    "fswm_extended_neuro2012_20_merge_seg.nii.gz"
                ))
                sentinal_file_list.append(os.path.join(
                    sentinal_file_basedir,
                    "TissueClassify",
                    "fswm_standard_neuro2012_20_merge_seg.nii.gz"
                ))
                sentinal_file_list.append(os.path.join(
                    sentinal_file_basedir,
                    "TissueClassify",
                    "fswm_extended_neuro2012_labelmap.png"
                ))


            if master_config['workflow_phase'] == 'atlas-based-reference':
                atlasDirectory = os.path.join(master_config['atlascache'], 'spatialImages', 'rho.nii.gz')
            else:
                atlasDirectory = os.path.join(master_config['previousresult'], subject, 'Atlas', 'AVG_rho.nii.gz')

            if os.path.exists(atlasDirectory):
                print "LOOKING FOR DIRECTORY {0}".format(atlasDirectory)
            else:
                print "MISSING REQUIRED ATLAS INPUT {0}".format(atlasDirectory)
                print("SKIPPING: {0} prerequisites missing".format(session))
                continue

            ## Use different sentinal file if segmentation specified.
            from workflows.baseline import DetermineIfSegmentationShouldBeDone

            do_BRAINSCut_Segmentation = DetermineIfSegmentationShouldBeDone(master_config)
            if do_BRAINSCut_Segmentation:
                sentinal_file_list.append(os.path.join(
                    sentinal_file_basedir,
                    "CleanedDenoisedRFSegmentations",
                    "allLabels_seg.nii.gz"
                ))

            def allPathsExists(list_of_paths):
                is_missing = False
                for ff in list_of_paths:
                    if not os.path.exists(ff):
                        is_missing = True
                        print("MISSING: {0}".format(ff))
                return not is_missing

            if useSentinal and allPathsExists(sentinal_file_list):
                print("SKIPPING: {0} exists".format(sentinal_file_list))
            else:
                print("PROCESSING INCOMPLETE: at least 1 required file does not exists")
                if dryRun == False:
                    workflow = _create_singleSession(_dict, master_config, 'Linear',
                                                     'singleSession_{0}_{1}'.format(_dict['subject'], _dict['session']))
                    print("Starting session {0}".format(session))
                    # HACK Hard-coded to SGEGraph, but --wfrun is ignored completely
                    run_workflow(workflow, plugin=master_config['plugin_name'],
                                 plugin_args=master_config['plugin_args'])
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
    createAndRun(kwds['SESSIONS'], environment, experiment, pipeline, cluster, useSentinal=kwds['--use-sentinal'],
                 dryRun=kwds['--dry-run'])
    return 0


# #####################################
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
