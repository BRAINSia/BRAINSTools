#!/usr/bin/python
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
        return  # input and output files can not match

    # first, see if the pattern is even in the file.
    with open(fname) as f:
        if not any(re.search(pat, line) for line in f):
            print "ERROR:  substitution pattern not found in reference file"
            sys.exit(-1)
            return  # pattern does not occur in file so we are done.

    # pattern is in the file, so perform replace operation.
    with open(fname) as f:
        out = open(out_fname, "w")
        for line in f:
            out.write(re.sub(pat, s_after, line))
        out.close()


def setGlobalDatasinkRewrite(cli, cfg):
    """
    Set the boolean GLOBAL_DATA_SINK_REWRITE

    If the flag '--ignore_datasinks' is set on the command line, pipeline will not rerun the
    pipeline datasinks.  If not, then the configuration file entry 'REWRITE_DATASINKS' under the
    'PIPELINE' heading will control this behavior

    >>> import baw_exp as baw
    >>> baw.setGlobalDatasinkRewrite(True, True); baw.GLOBAL_DATA_SINK_REWRITE
    *** Ignoring datasinks for pipeline rewriting ***
    False
    >>> baw.setGlobalDatasinkRewrite(False, True); baw.GLOBAL_DATA_SINK_REWRITE
    *** Ignoring datasinks for pipeline rewriting ***
    False
    >>> baw.setGlobalDatasinkRewrite(True, False); baw.GLOBAL_DATA_SINK_REWRITE
    *** Ignoring datasinks for pipeline rewriting ***
    False
    >>> baw.setGlobalDatasinkRewrite(False, False); baw.GLOBAL_DATA_SINK_REWRITE
    True

    :param cli: command line value
    :type cli: bool
    :param cfg: configuration file value
    :type cfg: bool

    Sets the variable `GLOBAL_DATA_SINK_REWRITE` constant flag used in :mod:`WorkupT1T2()`

    """
    assert isinstance(cli, bool) and isinstance(cfg, bool), \
      "Inputs are not boolean: {0}, {1}".format(cli, cfg)
    GLOBAL_DATA_SINK_REWRITE=False
    if cli or cfg:
        print "*** Ignoring datasinks for pipeline rewriting ***" # TODO: Use logging
        GLOBAL_DATA_SINK_REWRITE = False
    else:
        GLOBAL_DATA_SINK_REWRITE = True
    return GLOBAL_DATA_SINK_REWRITE


def main(argv=None):
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
    parser.add_argument('--doshort', action='store', dest='doshort', default=False, help='If not present, do long')
    parser.add_argument('--ignore_datasinks', action='store_true',
                        help='Use if the datasinks should not be rerun.\n \
    Default: value in configuration file')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    input_arguments = parser.parse_args()

    expConfig = ConfigParser.ConfigParser()
    expConfig.read(input_arguments.ExperimentConfig)

    # Pipeline-specific information
    ignore_datasinks = expConfig.getboolean('PIPELINE', 'GLOBAL_DATA_SINK_REWRITE')
    GLOBAL_DATA_SINK_REWRITE=setGlobalDatasinkRewrite(input_arguments.ignore_datasinks, ignore_datasinks)

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

    ## Setup environment for CPU load balancing of ITK based programs.
    import multiprocessing
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
    import SessionDB
    subjectDatabaseFile = os.path.join(ExperimentBaseDirectoryCache, 'InternalWorkflowSubjectDB.db')
    subject_list = input_arguments.subject.split(',')
    ## TODO:  Only make DB if db is older than subject_data_file.
    if (not os.path.exists(subjectDatabaseFile)) or (os.path.getmtime(subjectDatabaseFile) < os.path.getmtime(subject_data_file)):
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, subject_list)
        ExperimentDatabase.MakeNewDB(subject_data_file, mountPrefix)
        ExperimentDatabase = None
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, subject_list)
    else:
        print("Using cached database, {0}".format(subjectDatabaseFile))
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, subject_list)
    print "ENTIRE DB for {_subjid}: ".format(_subjid=ExperimentDatabase.getSubjectFilter())
    print "^^^^^^^^^^^^^"
    for row in ExperimentDatabase.getEverything():
        print row
    print "^^^^^^^^^^^^^"

    ## Create the shell wrapper script for ensuring that all jobs running on remote hosts from SGE
    #  have the same environment as the job submission host.
    JOB_SCRIPT = get_global_sge_script(sys.path, PROGRAM_PATHS, CUSTOM_ENVIRONMENT)
    print JOB_SCRIPT

    import WorkupT1T2  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    print "TESTER"
    import ShortWorkupT1T2
    for subjectid in ExperimentDatabase.getAllSubjects():
        if input_arguments.doshort:
            baw200 = ShortWorkupT1T2.ShortWorkupT1T2(subjectid, mountPrefix,
                                                     os.path.join(ExperimentBaseDirectoryCache, str(subjectid)),
                                                     ExperimentBaseDirectoryResults,
                                                     ExperimentDatabase,
                                                     CACHE_ATLASPATH,
                                                     CACHE_BCDMODELPATH,
                                                     GLOBAL_DATA_SINK_REWRITE,
                                                     WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS, CLUSTER_QUEUE=CLUSTER_QUEUE, CLUSTER_QUEUE_LONG=CLUSTER_QUEUE_LONG)
        else:
            baw200 = WorkupT1T2.WorkupT1T2(subjectid, mountPrefix,
                                           os.path.join(ExperimentBaseDirectoryCache, str(subjectid)),
                                           ExperimentBaseDirectoryResults,
                                           ExperimentDatabase,
                                           CACHE_ATLASPATH,
                                           CACHE_BCDMODELPATH,
                                           GLOBAL_DATA_SINK_REWRITE,
                                           WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS, CLUSTER_QUEUE=CLUSTER_QUEUE, CLUSTER_QUEUE_LONG=CLUSTER_QUEUE_LONG, SGE_JOB_SCRIPT=JOB_SCRIPT)
        print "Start Processing"

        SGEFlavor = 'SGE'
        try:
            if input_arguments.wfrun == 'helium_all.q':
                try:
                    baw200.write_graph()
                except:
                    pass
                baw200.run(plugin=SGEFlavor,
                           plugin_args=dict(template=JOB_SCRIPT, qsub_args="-S /bin/bash -cwd -pe smp1 1-12 -l h_vmem=19G,mem_free=9G -o /dev/null -e /dev/null " + CLUSTER_QUEUE))
            elif input_arguments.wfrun == 'helium_all.q_graph':
                try:
                    baw200.write_graph()
                except:
                    pass
                SGEFlavor = 'SGEGraph'  # Use the SGEGraph processing
                baw200.run(plugin=SGEFlavor,
                           plugin_args=dict(template=JOB_SCRIPT, qsub_args="-S /bin/bash -cwd -pe smp1 1-12 -l h_vmem=19G,mem_free=9G -o /dev/null -e /dev/null " + CLUSTER_QUEUE))
            elif input_arguments.wfrun == 'ipl_OSX':
                try:
                    baw200.write_graph()
                except:
                    pass
                print "Running On ipl_OSX"
                baw200.run(plugin=SGEFlavor,
                           plugin_args=dict(template=JOB_SCRIPT, qsub_args="-S /bin/bash -cwd -pe smp1 1-12 -l h_vmem=19G,mem_free=9G -o /dev/null -e /dev/null " + CLUSTER_QUEUE))
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
        except Exception, err:
            print("ERROR: EXCEPTION CAUGHT IN RUNNING SUBJECT {0}".format(subjectid))
            raise err

    print("THIS RUN OF BAW FOR SUBJS {0} HAS COMPLETED".format(ExperimentDatabase.getAllSubjects()))
    return 0

if __name__ == "__main__":
    sys.exit(main())
