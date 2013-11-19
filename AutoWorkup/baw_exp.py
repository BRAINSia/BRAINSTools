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
import os
import re
import sys
import traceback

import multiprocessing
import time
##############################################################################


def get_global_sge_script(pythonPathsList, binPathsList, customEnvironment={}):
    """This is a wrapper script for running commands on an SGE cluster
so that all the python modules and commands are pathed properly"""

    custEnvString = ""
    for key, value in customEnvironment.items():
        custEnvString += "export " + key + "=" + value + "\n"

    PYTHONPATH = ":".join(pythonPathsList)
    BASE_BUILDS = ":".join(binPathsList)
    GLOBAL_SGE_SCRIPT = """#!/bin/bash
echo "STARTED at: $(date +'%F-%T')"
echo "Ran on: $(hostname)"
export PATH={BINPATH}
export PYTHONPATH={PYTHONPATH}

echo "========= CUSTOM ENVIORNMENT SETTINGS =========="
echo "export PYTHONPATH={PYTHONPATH}"
echo "export PATH={BINPATH}"
echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

echo "With custom environment:"
echo {CUSTENV}
{CUSTENV}
## NOTE:  nipype inserts the actual commands that need running below this section.
""".format(PYTHONPATH=PYTHONPATH, BINPATH=BASE_BUILDS, CUSTENV=custEnvString)
    return GLOBAL_SGE_SCRIPT

# From http://stackoverflow.com/questions/1597649/replace-strings-in-files-by-python


def file_replace(fname, out_fname, pat, s_after):
    if fname == out_fname:
        print "ERROR: input and output file names can not match"
        sys.exit(-1)

    # first, see if the pattern is even in the file.
    with open(fname) as f:
        if not any(re.search(pat, line) for line in f):
            print "ERROR:  substitution pattern not found in reference file"
            sys.exit(-1)

    # pattern is in the file, so perform replace operation.
    with open(fname) as f:
        out = open(out_fname, "w")
        for line in f:
            out.write(re.sub(pat, s_after, line))
        out.close()


def setDataSinkRewriteValue(cli, cfg):
    """
    Return the behavior for boolean GLOBAL_DATA_SINK_REWRITE

    If the flag '--rewrite_datasinks' is set on the command line, pipeline will force rerun of the
    pipeline datasinks.  If not, then the configuration file entry 'GLOBAL_DATA_SINK_REWRITE' under the
    'PIPELINE' heading will control this behavior

    >>> import baw_exp as baw
    >>> baw.setDataSinkRewriteValue(True, True); baw.GLOBAL_DATA_SINK_REWRITE
    *** Ignoring datasinks for pipeline rewriting ***
    False
    >>> baw.setDataSinkRewriteValue(False, True); baw.GLOBAL_DATA_SINK_REWRITE
    *** Ignoring datasinks for pipeline rewriting ***
    False
    >>> baw.setDataSinkRewriteValue(True, False); baw.GLOBAL_DATA_SINK_REWRITE
    *** Ignoring datasinks for pipeline rewriting ***
    False
    >>> baw.setDataSinkRewriteValue(False, False); baw.GLOBAL_DATA_SINK_REWRITE
    True

    :param cli: command line value
    :type cli: bool
    :param cfg: configuration file value
    :type cfg: bool

    Sets the variable `GLOBAL_DATA_SINK_REWRITE` constant flag used in :mod:`WorkupT1T2()`
    """
    assert isinstance(cli, bool) and isinstance(cfg, bool), "Inputs are not boolean: {0}, {1}".format(cli, cfg)
    if cli or cfg:
        print "*** Force datasinks rewriting for pipeline rewriting ***, commandline= {0}, configfile= {1}".format(cli, cfg)  # TODO: Use logging
        GLOBAL_DATA_SINK_REWRITE = True
    else:
        print "*** Default datasink behavior for pipeline ***, commandline= {0}, configfile= {1}".format(cli, cfg)  # TODO: Use logging
        GLOBAL_DATA_SINK_REWRITE = False
    return GLOBAL_DATA_SINK_REWRITE


def OpenSubjectDatabase(ExperimentBaseDirectoryCache, single_subject, mountPrefix, subject_data_file):
    import SessionDB

    subjectDatabaseFile = os.path.join(ExperimentBaseDirectoryCache, 'InternalWorkflowSubjectDB.db')
    ## TODO:  Only make DB if db is older than subject_data_file.
    if (not os.path.exists(subjectDatabaseFile)) or (
            os.path.getmtime(subjectDatabaseFile) < os.path.getmtime(subject_data_file)):
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
        ExperimentDatabase.MakeNewDB(subject_data_file, mountPrefix)
        ExperimentDatabase = None
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
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

    CACHE_ATLASPATH, CACHE_BCDMODELPATH, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG,QSTAT_IMMEDIATE_EXE,QSTAT_CACHED_EXE, \
          ExperimentBaseDirectoryCache, ExperimentBaseDirectoryResults, subject_data_file, \
          GLOBAL_DATA_SINK_REWRITE, JOB_SCRIPT, WORKFLOW_COMPONENTS, \
          input_arguments, mountPrefix,start_time,subjectid = sp_args

    while time.time() < start_time :
        time.sleep(start_time-time.time()+1)
        print "Delaying start for {0}".format(subjectid)

    list_with_one_subject = [ subjectid ]
    ExperimentDatabase = OpenSubjectDatabase(ExperimentBaseDirectoryCache, list_with_one_subject, mountPrefix,
                                             subject_data_file)

    if input_arguments.doshort:
        import ShortWorkupT1T2
        baw200 = ShortWorkupT1T2.ShortWorkupT1T2(subjectid, mountPrefix,
                                                 os.path.join(ExperimentBaseDirectoryCache, str(subjectid)),
                                                 ExperimentBaseDirectoryResults,
                                                 ExperimentDatabase,
                                                 CACHE_ATLASPATH,
                                                 CACHE_BCDMODELPATH,
                                                 GLOBAL_DATA_SINK_REWRITE,
                                                 WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS, CLUSTER_QUEUE=CLUSTER_QUEUE,
                                                 CLUSTER_QUEUE_LONG=CLUSTER_QUEUE_LONG)
    else:
        import WorkupT1T2  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
        baw200 = WorkupT1T2.WorkupT1T2(subjectid, mountPrefix,
                                       os.path.join(ExperimentBaseDirectoryCache, str(subjectid)),
                                       ExperimentBaseDirectoryResults,
                                       ExperimentDatabase,
                                       CACHE_ATLASPATH,
                                       CACHE_BCDMODELPATH,
                                       GLOBAL_DATA_SINK_REWRITE,
                                       WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS, CLUSTER_QUEUE=CLUSTER_QUEUE,
                                       CLUSTER_QUEUE_LONG=CLUSTER_QUEUE_LONG, SGE_JOB_SCRIPT=JOB_SCRIPT)
    print "Start Processing"
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
                                        qsub_args="-S /bin/bash -cwd -pe smp 1-12 -l h_vmem=19G,mem_free=9G -o /dev/null -e /dev/null " + CLUSTER_QUEUE,
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
                                        qsub_args="-S /bin/bash -cwd -pe smp 1-12 -l h_vmem=19G,mem_free=9G -o /dev/null -e /dev/null " + CLUSTER_QUEUE,
                                        qstatProgramPath=QSTAT_IMMEDIATE_EXE,
                                        qstatCachedProgramPath=QSTAT_CACHED_EXE))
        elif input_arguments.wfrun == 'ipl_OSX':
            try:
                baw200.write_graph()
            except:
                pass
            print "Running On ipl_OSX"
            baw200.run(plugin=SGEFlavor,
                       plugin_args=dict(template=JOB_SCRIPT,
                                        qsub_args="-S /bin/bash -cwd -pe smp 1-12 -l h_vmem=19G,mem_free=9G -o /dev/null -e /dev/null " + CLUSTER_QUEUE,
                                        qstatProgramPath=QSTAT_IMMEDIATE_EXE,
                                        qstatCachedProgramPath=QSTAT_CACHED_EXE))
        elif input_arguments.wfrun == 'local_4':
            try:
                baw200.write_graph()
            except:
                pass
            print "Running with 4 parallel processes on local machine"
            baw200.run(plugin='MultiProc', plugin_args={'n_procs': 4})
        elif input_arguments.wfrun == 'local_12':
            try:
                baw200.write_graph()
            except:
                pass
            print "Running with 12 parallel processes on local machine"
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
            print "Running sequentially on local machine"
            # baw200.run(updatehash=True)
            baw200.run()
        else:
            print "You must specify the run environment type. [helium_all.q,helium_all.q_graph,ipl_OSX,local_4,local_12,local]"
            print input_arguments.wfrun
            sys.exit(-1)
    except:
        print("ERROR: EXCEPTION CAUGHT IN RUNNING SUBJECT {0}".format(subjectid))
        traceback.print_exc(file=sys.stdout)
        return False
    return True


def MasterProcessingController(argv=None):
    import argparse
    import ConfigParser
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
    parser.add_argument('-doshort', action='store', dest='doshort', default=False, help='If not present, do long')
    parser.add_argument('-rewrite_datasinks', action='store_true', default=False,
                        help='Use if the datasinks should be forced rerun.\nDefault: value in configuration file')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    input_arguments = parser.parse_args()

    expConfig = ConfigParser.ConfigParser()
    expConfig.read(input_arguments.ExperimentConfig)

    # Pipeline-specific information
    GLOBAL_DATA_SINK_REWRITE_FROM_CONFIG = expConfig.getboolean('PIPELINE', 'GLOBAL_DATA_SINK_REWRITE')
    GLOBAL_DATA_SINK_REWRITE = setDataSinkRewriteValue(input_arguments.rewrite_datasinks, GLOBAL_DATA_SINK_REWRITE_FROM_CONFIG)

    # Experiment specific information
    subject_data_file = expConfig.get('EXPERIMENT_DATA', 'SESSION_DB')
    ExperimentName = expConfig.get('EXPERIMENT_DATA', 'EXPERIMENTNAME')
    WORKFLOW_COMPONENTS_STRING = expConfig.get('EXPERIMENT_DATA', 'WORKFLOW_COMPONENTS')
    WORKFLOW_COMPONENTS = eval(WORKFLOW_COMPONENTS_STRING)

    # Platform specific information
    #     Prepend the python search paths
    PYTHON_AUX_PATHS = expConfig.get(input_arguments.processingEnvironment, 'PYTHON_AUX_PATHS')
    PYTHON_AUX_PATHS = PYTHON_AUX_PATHS.split(':')
    PYTHON_AUX_PATHS.extend(sys.path)
    sys.path = PYTHON_AUX_PATHS
    ######################################################################################
    ###### Now ensure that all the required packages can be read in from this custom path
    #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # print sys.path
    from nipype import config  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    config.enable_debug_mode()  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
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

    ## Check to ensure that SimpleITK can be found
    import SimpleITK as sitk
    #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    #####################################################################################
    #  FreeSurfer is extraordinarly finicky and is easily confused and incorrect.
    #  Force that all the FREESURFER env vars are set in subsequent scripts by
    #  ensuring that rough versions of these environmental variables are not
    #  set internal to this script.
    prohibited_env_var_exists = False
    for ENVVAR_TO_CHECK in ['FREESURFER_HOME', 'FSFAST_HOME', 'FSF_OUTPUT_FORMAT', 'SUBJECTS_DIR', 'MNI_DIR', 'FSL_DIR']:
        if ENVVAR_TO_CHECK in os.environ:
            prohibited_env_var_exists = True
            print("ERROR: Environmental Variable {0}={1} exists.  Please unset before continuing.".format(ENVVAR_TO_CHECK, os.environ[ENVVAR_TO_CHECK]))
    if prohibited_env_var_exists:
        sys.exit(-1)

    #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    #####################################################################################
    #     Prepend the shell environment search paths
    PROGRAM_PATHS = expConfig.get(input_arguments.processingEnvironment, 'PROGRAM_PATHS')
    PROGRAM_PATHS = PROGRAM_PATHS.split(':')
    PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))
    os.environ['PATH'] = ':'.join(PROGRAM_PATHS)
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
        file_replace(os.path.join(ATLASPATH, 'ExtendedAtlasDefinition.xml.in'), os.path.join(CACHE_ATLASPATH, 'ExtendedAtlasDefinition.xml'), "@ATLAS_DIRECTORY@", CACHE_ATLASPATH)
    else:
        print("Atlas already exists in experiment cache directory: {0}".format(CACHE_ATLASPATH))
    #  Just to be safe, copy the model file as well
    BCDMODELPATH = expConfig.get(input_arguments.processingEnvironment, 'BCDMODELPATH')
    CACHE_BCDMODELPATH = os.path.join(ExperimentBaseDirectoryCache, os.path.basename(BCDMODELPATH))
    from distutils.file_util import copy_file
    for BCDModelFile in ['LLSModel-2ndVersion.h5', 'T1-2ndVersion.mdl']:
        if BCDModelFile[-2:] == 'h5':
            BCDModelFile = os.path.join('Transforms_h5', BCDModelFile)
        orig = os.path.join(BCDMODELPATH, BCDModelFile)
        new = os.path.join(CACHE_BCDMODELPATH, BCDModelFile)
        new = new.replace('Transforms_h5/', '')  # Flatten back out, even if you needed to get files from subdirectory.
        if not os.path.exists(CACHE_BCDMODELPATH):
            os.mkdir(CACHE_BCDMODELPATH)
        if not os.path.exists(new):
            print("Copying BCD Model file to cache directory: {0}".format(new))
            copy_file(orig, new, preserve_mode=1, preserve_times=1)
        else:
            print("BCD Model exists in cache directory: {0}".format(new))

    CUSTOM_ENVIRONMENT = expConfig.get(input_arguments.processingEnvironment, 'CUSTOM_ENVIRONMENT')
    CUSTOM_ENVIRONMENT = eval(CUSTOM_ENVIRONMENT)
    ## Set custom environmental variables so that subproceses work properly (i.e. for FreeSurfer)
    # print CUSTOM_ENVIRONMENT
    for key, value in CUSTOM_ENVIRONMENT.items():
        # print "SETTING: ", key, value
        os.putenv(key, value)
        os.environ[key] = value
    # print os.environ
    # sys.exit(-1)

    ## If freesurfer is requested, then ensure that a sane environment is available
    if 'FREESURFER' in WORKFLOW_COMPONENTS:
        print "FREESURFER NEEDS TO CHECK FOR SANE ENVIRONMENT HERE."

    CLUSTER_QUEUE = expConfig.get(input_arguments.processingEnvironment, 'CLUSTER_QUEUE')
    CLUSTER_QUEUE_LONG = expConfig.get(input_arguments.processingEnvironment, 'CLUSTER_QUEUE_LONG')
    QSTAT_IMMEDIATE_EXE = expConfig.get(input_arguments.processingEnvironment, 'QSTAT_IMMEDIATE_EXE')
    QSTAT_CACHED_EXE = expConfig.get(input_arguments.processingEnvironment, 'QSTAT_CACHED_EXE')

    ## Setup environment for CPU load balancing of ITK based programs.
    total_CPUS = multiprocessing.cpu_count()
    if input_arguments.wfrun == 'helium_all.q':
        pass
    elif input_arguments.wfrun == 'helium_all.q_graph':
        pass
    elif input_arguments.wfrun == 'ipl_OSX':
        pass
    elif input_arguments.wfrun == 'local_4':
        os.environ['NSLOTS'] = "{0}".format(total_CPUS / 4)
    elif input_arguments.wfrun == 'local_12':
        os.environ['NSLOTS'] = "{0}".format(total_CPUS / 12)
    elif input_arguments.wfrun == 'local':
        os.environ['NSLOTS'] = "{0}".format(total_CPUS / 1)
    elif input_arguments.wfrun == 'ds_runner':
        os.environ['NSLOTS'] = "{0}".format(total_CPUS / 1)
    else:
        print "FAILED RUN: You must specify the run environment type. [helium_all.q,helium_all.q_graph,ipl_OSX,local_4,local_12,local,ds_runner]"
        print input_arguments.wfrun
        sys.exit(-1)

    print "Configuring Pipeline"
    ## Ensure that entire db is built and cached before parallel section starts.
    _ignoreme = OpenSubjectDatabase(ExperimentBaseDirectoryCache, [ "all" ], mountPrefix, subject_data_file)
    to_do_subjects = input_arguments.subject.split(',')
    if to_do_subjects[0] == "all":
        to_do_subjects=_ignoreme.getAllSubjects()
    _ignoreme = None

    ## Create the shell wrapper script for ensuring that all jobs running on remote hosts from SGE
    #  have the same environment as the job submission host.
    JOB_SCRIPT = get_global_sge_script(sys.path, PROGRAM_PATHS, CUSTOM_ENVIRONMENT)
    print JOB_SCRIPT

    # Randomly shuffle to_do_subjects to get max
    import random
    random.shuffle(to_do_subjects)

    ## Make a list of all the arguments to be processed
    sp_args_list = list()
    start_time=time.time()
    subj_index = 1
    for subjectid in to_do_subjects:
        delay = 5*subj_index
        subj_index += 1
        print("START DELAY: {0}".format(delay))
        sp_args=(CACHE_ATLASPATH, CACHE_BCDMODELPATH, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG,QSTAT_IMMEDIATE_EXE,QSTAT_CACHED_EXE,
                                  ExperimentBaseDirectoryCache, ExperimentBaseDirectoryResults, subject_data_file,
                                  GLOBAL_DATA_SINK_REWRITE, JOB_SCRIPT, WORKFLOW_COMPONENTS, input_arguments,
                                  mountPrefix, start_time+delay, subjectid)
        sp_args_list.append(sp_args)

    if 'local' in input_arguments.wfrun:
        print("RUNNING WITHOUT POOL BUILDING")
        for sp_args in sp_args_list:
            DoSingleSubjectProcessing(sp_args)
    else:
        ## Make a pool of workers to submit simultaneously
        from multiprocessing import Pool
        myPool = Pool(processes=64,maxtasksperchild=2)
        all_results=myPool.map_async(DoSingleSubjectProcessing,sp_args_list).get(1e100)

        for indx in range(0,len(sp_args_list)):
            if all_results[indx] == False:
                    print "FAILED for {0}".format(sp_args_list[indx][-1])

    print("THIS RUN OF BAW FOR SUBJS {0} HAS COMPLETED".format(to_do_subjects))
    return 0

if __name__ == "__main__":
    main_status = MasterProcessingController()
    sys.exit(main_status)
