#!/usr/bin/env python
#################################################################################
## Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
## Language:  Python
##
## Author:  Hans J. Johnson
##
##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.
##
#################################################################################
from __future__ import print_function
from __future__ import absolute_import
from past.builtins import execfile
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from builtins import object
import os
import re
import sys
import traceback

import multiprocessing
import time
##############################################################################


def OpenSubjectDatabase(ExperimentBaseDirectoryCache, single_subject, mountPrefix, subject_data_file):
    import os.path
    import SessionDB
    subjectDatabaseFile = os.path.join(ExperimentBaseDirectoryCache, 'InternalWorkflowSubjectDB.db')
    ## TODO:  Only make DB if db is older than subject_data_file.
    if (not os.path.exists(subjectDatabaseFile)) or \
      (os.path.getmtime(subjectDatabaseFile) < os.path.getmtime(subject_data_file)):
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
        ExperimentDatabase.MakeNewDB(subject_data_file, mountPrefix)
    else:
        print("Single_subject {0}: Using cached database, {1}".format(single_subject,subjectDatabaseFile))
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
    #print "ENTIRE DB for {_subjid}: ".format(_subjid=ExperimentDatabase.getSubjectFilter())
    #print "^^^^^^^^^^^^^"
    #for row in ExperimentDatabase.getEverything():
    #    print row
    #print "^^^^^^^^^^^^^"
    return ExperimentDatabase

def DoSingleSubjectProcessing(sp_args):

    CACHE_ATLASPATH, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG,QSTAT_IMMEDIATE_EXE,QSTAT_CACHED_EXE, \
          ExperimentBaseDirectoryCache, ExperimentBaseDirectoryResults, subject_data_file, \
          GLOBAL_DATA_SINK_REWRITE, JOB_SCRIPT, WORKFLOW_COMPONENTS, \
          input_arguments, mountPrefix,start_time,subjectid, PreviousExperimentName = sp_args

    while time.time() < start_time :
        time.sleep(start_time-time.time()+1)
        print("Delaying start for {0}".format(subjectid))

    list_with_one_subject = [ subjectid ]
    ExperimentDatabase = OpenSubjectDatabase(ExperimentBaseDirectoryCache, list_with_one_subject, mountPrefix,
                                             subject_data_file)

    import WorkupT1T2  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    if not PreviousExperimentName is None:
        print("Running based on previous experiment results...")
        baw200 = WorkupT1T2.WorkupT1T2(subjectid, mountPrefix,
                                       os.path.join(ExperimentBaseDirectoryCache, str(subjectid)),
                                       ExperimentBaseDirectoryResults,
                                       ExperimentDatabase,
                                       CACHE_ATLASPATH,
                                       GLOBAL_DATA_SINK_REWRITE,
                                       WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS, CLUSTER_QUEUE=CLUSTER_QUEUE,
                                       CLUSTER_QUEUE_LONG=CLUSTER_QUEUE_LONG, SGE_JOB_SCRIPT=JOB_SCRIPT,
                                       PreviousBaseDirectoryResults=PreviousExperimentName)
    else:
        baw200 = WorkupT1T2.WorkupT1T2(subjectid, mountPrefix,
                                       os.path.join(ExperimentBaseDirectoryCache, str(subjectid)),
                                       ExperimentBaseDirectoryResults,
                                       ExperimentDatabase,
                                       CACHE_ATLASPATH,
                                       GLOBAL_DATA_SINK_REWRITE,
                                       WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS, CLUSTER_QUEUE=CLUSTER_QUEUE,
                                       CLUSTER_QUEUE_LONG=CLUSTER_QUEUE_LONG, SGE_JOB_SCRIPT=JOB_SCRIPT)
    print("Start Processing")
    SGEFlavor = 'SGE'
    try:
        if input_arguments.wfrun == 'helium_all.q':
            try:
                #baw200.write_graph()
                pass
            except:
                pass
            baw200.run(plugin=SGEFlavor,
                       plugin_args=dict(template=JOB_SCRIPT,
                                        qsub_args=modify_qsub_args(CLUSTER_QUEUE,2,1,1),
                                        qstatProgramPath=QSTAT_IMMEDIATE_EXE,
                                        qstatCachedProgramPath=QSTAT_CACHED_EXE))
        elif input_arguments.wfrun == 'helium_all.q_graph':
            try:
                #baw200.write_graph()
                pass
            except:
                pass
            SGEFlavor = 'SGEGraph'  # Use the SGEGraph processing
            baw200.run(plugin=SGEFlavor,
                       plugin_args=dict(template=JOB_SCRIPT,
                                        qsub_args=modify_qsub_args(CLUSTER_QUEUE,2,1,1),
                                        qstatProgramPath=QSTAT_IMMEDIATE_EXE,
                                        qstatCachedProgramPath=QSTAT_CACHED_EXE))
        elif input_arguments.wfrun == 'ipl_OSX':
            try:
                baw200.write_graph()
            except:
                pass
            print("Running On ipl_OSX")
            baw200.run(plugin=SGEFlavor,
                       plugin_args=dict(template=JOB_SCRIPT,
                                        qsub_args=modify_qsub_args(CLUSTER_QUEUE,2,1,1),
                                        qstatProgramPath=QSTAT_IMMEDIATE_EXE,
                                        qstatCachedProgramPath=QSTAT_CACHED_EXE))
        elif input_arguments.wfrun == 'local_4':
            try:
                baw200.write_graph()
            except:
                pass
            print("Running with 4 parallel processes on local machine")
            baw200.run(plugin='MultiProc', plugin_args={'n_procs': 4})
        elif input_arguments.wfrun == 'local_12':
            try:
                baw200.write_graph()
            except:
                pass
            print("Running with 12 parallel processes on local machine")
            baw200.run(plugin='MultiProc', plugin_args={'n_procs': 12})
        elif input_arguments.wfrun == 'ds_runner':
            class ds_runner(object):
                def run(self, graph, **kwargs):
                    for node in graph.nodes():
                        if '_ds' in node.name.lower():
                            node.run()

            baw200.run(plugin=ds_runner())
        elif input_arguments.wfrun == 'local':
            try:
                baw200.write_graph()
            except:
                pass
            print("Running sequentially on local machine")
            # baw200.run(updatehash=True)
            baw200.run()
        else:
            print("You must specify the run environment type. [helium_all.q,helium_all.q_graph,ipl_OSX,local_4,local_12,local]")
            print(input_arguments.wfrun)
            sys.exit(-1)
    except:
        print("ERROR: EXCEPTION CAUGHT IN RUNNING SUBJECT {0}".format(subjectid))
        traceback.print_exc(file=sys.stdout)
        return False
    return True


def MasterProcessingController(argv=None):
    import argparse
    import configparser
    import csv
    import string

    if argv == None:
        argv = sys.argv

    # Create and parse input arguments
    parser = argparse.ArgumentParser(description='Runs a mini version of BRAINSAutoWorkup')
    group = parser.add_argument_group('Required')
    group.add_argument('-pe', action="store", dest='processingEnvironment', required=True,
                       help='The name of the processing environment to use from the config file')
    group.add_argument('-wfrun', action="store", dest='wfrun', required=True,
                       help='The name of the workflow running plugin to use')
    group.add_argument('-subject', action="store", dest='subject', required=True,
                       help='The name of the subject to process')
    group.add_argument('-ExperimentConfig', action="store", dest='ExperimentConfig', required=True,
                       help='The path to the file that describes the entire experiment')
    parser.add_argument('-rewrite_datasinks', action='store_true', default=False,
                        help='Use if the datasinks should be forced rerun.\nDefault: value in configuration file')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    config = configparser.ConfigParser(allow_no_value=True)
    config.read(args.ExperimentConfig)

    # Pipeline-specific information
    GLOBAL_DATA_SINK_REWRITE = setDataSinkRewriteValue(args.rewrite_datasinks, config.getboolean('NIPYPE', 'GLOBAL_DATA_SINK_REWRITE'))
    experiment = get_experiment_settings(config)
    # Platform specific information
    environment = get_environment_settings(config)
    if environment['cluster']:
        cluster = get_cluster_settings(config)
    sys.path = environment('PYTHONPATH')
    os.environ['PATH'] = ':'.join(environment['PATH'])
    # Virtualenv
    if not environment['virtualenv_dir'] is None:
        print("Loading virtualenv_dir...")
        execfile(environment['virtualenv_dir'], dict(__file__=environment['virtualenv_dir']))
    ###### Now ensure that all the required packages can be read in from this custom path
    #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # print sys.path
    ## Check to ensure that SimpleITK can be found
    import SimpleITK as sitk
    from nipype import config  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    config.enable_debug_mode()  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    #config.enable_provenance()

    ##############################################################################
    from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
    from nipype.interfaces.base import traits, isdefined, BaseInterface
    from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
    import nipype.interfaces.io as nio   # Data i/o
    import nipype.pipeline.engine as pe  # pypeline engine
    from nipype.interfaces.freesurfer import ReconAll

    from nipype.utils.misc import package_check
    # package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
    package_check('numpy', '1.3', 'tutorial1')
    package_check('scipy', '0.7', 'tutorial1')
    package_check('networkx', '1.0', 'tutorial1')
    package_check('IPython', '0.10', 'tutorial1')

    try:
        verify_empty_freesurfer_env()
    except EnvironmentError:
        raise

    # Define platform specific output write paths
    if not os.path.exists(experiment['output_cache']):
        os.makedirs(experiment['output_cache'])
    if not os.path.exists(experiment['output_results']):
        os.makedirs(experiment['output_results'])
    if 'input_results' in list(experiment.keys()):
        assert os.path.exists(experiment['input_results']), "The previous experiment directory does not exist: {0}".format(experiment['input_results'])

    #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    #    Define platform specific output write paths
    mountPrefix = expConfig.get(input_arguments.processingEnvironment, 'MOUNTPREFIX')
    BASEOUTPUTDIR = expConfig.get(input_arguments.processingEnvironment, 'BASEOUTPUTDIR')
    ExperimentBaseDirectoryPrefix = os.path.realpath(os.path.join(BASEOUTPUTDIR, ExperimentName))
    ExperimentBaseDirectoryCache = ExperimentBaseDirectoryPrefix + "_CACHE"
    ExperimentBaseDirectoryResults = ExperimentBaseDirectoryPrefix + "_Results"
    if not os.path.exists(ExperimentBaseDirectoryCache):
        os.makedirs(ExperimentBaseDirectoryCache)
    if not os.path.exists(ExperimentBaseDirectoryResults):
        os.makedirs(ExperimentBaseDirectoryResults)
    if not PreviousExperimentName is None:
        PreviousBaseDirectoryPrefix = os.path.realpath(os.path.join(BASEOUTPUTDIR, PreviousExperimentName))
        PreviousBaseDirectoryResults = PreviousBaseDirectoryPrefix + "_Results"
        assert os.path.exists(PreviousBaseDirectoryResults), "The previous experiment directory does not exist: {0}".format(PreviousBaseDirectoryResults)
    else:
        PreviousBaseDirectoryResults = None
    #    Define workup common reference data sets
    #    The ATLAS needs to be copied to the ExperimentBaseDirectoryPrefix
    #    The ATLAS pathing must stay constant
    ATLASPATH = expConfig.get(input_arguments.processingEnvironment, 'ATLASPATH')
    if not os.path.exists(ATLASPATH):
        print("ERROR:  Invalid Path for Atlas: {0}".format(ATLASPATH))
        sys.exit(-1)
    CACHE_ATLASPATH = os.path.realpath(os.path.join(ExperimentBaseDirectoryCache, 'Atlas'))
    from distutils.dir_util import copy_tree
    if not os.path.exists(CACHE_ATLASPATH):
        print("Copying a reference of the atlas to the experiment cache directory:\n    from: {0}\n    to: {1}".format(ATLASPATH, CACHE_ATLASPATH))
        copy_tree(ATLASPATH, CACHE_ATLASPATH, preserve_mode=1, preserve_times=1)
        ## Now generate the xml file with the correct pathing
        file_replace(os.path.join(ATLASPATH, 'ExtendedAtlasDefinition.xml.in'), os.path.join(CACHE_ATLASPATH, 'ExtendedAtlasDefinition.xml'), "@ATLAS_INSTALL_DIRECTORY@", CACHE_ATLASPATH)
    else:
        print("Atlas already exists in experiment cache directory: {0}".format(CACHE_ATLASPATH))

    ## Set custom environmental variables so that subproceses work properly (i.e. for FreeSurfer)
    CUSTOM_ENVIRONMENT = eval(environment['misc'])
    # print CUSTOM_ENVIRONMENT
    for key, value in list(CUSTOM_ENVIRONMENT.items()):
        # print "SETTING: ", key, value
        os.putenv(key, value)
        os.environ[key] = value
    # print os.environ
    # sys.exit(-1)

    WORKFLOW_COMPONENTS = experiment['components']
    if 'FREESURFER' in WORKFLOW_COMPONENTS:
        check_freesurfer_environment()

    cluster = setup_cpu(args.wfrun, config)  # None unless wfrun is 'helium*' or 'ipl_OSX', then dict()

    print("Configuring Pipeline")
    ## Ensure that entire db is built and cached before parallel section starts.
    _ignoreme = OpenSubjectDatabase(experiment['output_cache'], [ "all" ], environment['prefix'], environment['subject_data_file'])
    to_do_subjects = args.subject.split(',')
    if to_do_subjects[0] == "all":
        to_do_subjects=_ignoreme.getAllSubjects()
    _ignoreme = None

    ## Create the shell wrapper script for ensuring that all jobs running on remote hosts from SGE
    #  have the same environment as the job submission host.

    JOB_SCRIPT = get_global_sge_script(sys.path, os.environ['PATH'], CUSTOM_ENVIRONMENT, MODULES)
    print(JOB_SCRIPT)

    # Randomly shuffle to_do_subjects to get max
    import random
    random.shuffle(to_do_subjects)

    ## Make a list of all the arguments to be processed
    sp_args_list = list()
    start_time=time.time()
    subj_index = 1
    for subjectid in to_do_subjects:
        delay = 2.5*subj_index
        subj_index += 1
        print("START DELAY: {0}".format(delay))
        sp_args=(CACHE_ATLASPATH, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG,QSTAT_IMMEDIATE_EXE,QSTAT_CACHED_EXE,
                                  experiment['output_cache'], experiment['output_results'], environment['subject_data_file'],
                                  GLOBAL_DATA_SINK_REWRITE, JOB_SCRIPT, WORKFLOW_COMPONENTS, args,
                                  mountPrefix, start_time+delay, subjectid, PreviousBaseDirectoryResult)
        sp_args_list.append(sp_args)
    if 'local' in args.wfrun:
        print("RUNNING WITHOUT POOL BUILDING")
        for sp_args in sp_args_list:
            DoSingleSubjectProcessing(sp_args)
    else:
        ## Make a pool of workers to submit simultaneously
        from multiprocessing import Pool
        myPool = Pool(processes=64,maxtasksperchild=1)
        all_results=myPool.map_async(DoSingleSubjectProcessing,sp_args_list).get(1e100)

        for indx in range(0,len(sp_args_list)):
            if all_results[indx] == False:
                    print("FAILED for {0}".format(sp_args_list[indx][-1]))

    print("THIS RUN OF BAW FOR SUBJS {0} HAS COMPLETED".format(to_do_subjects))
    return 0

if __name__ == "__main__":
    import sys
    main_status = MasterProcessingController(None)
    sys.exit(main_status)
