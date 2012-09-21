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
def get_global_sge_script(pythonPathsList,binPathsList,customEnvironment={}):
    """This is a wrapper script for running commands on an SGE cluster
so that all the python modules and commands are pathed properly"""

    custEnvString=""
    for key,value in customEnvironment.items():
        custEnvString+="export "+key+"="+value+"\n"

    PYTHONPATH=":".join(pythonPathsList)
    BASE_BUILDS=":".join(binPathsList)
    GLOBAL_SGE_SCRIPT="""#!/bin/bash
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
""".format(PYTHONPATH=PYTHONPATH,BINPATH=BASE_BUILDS,CUSTENV=custEnvString)
    return GLOBAL_SGE_SCRIPT

# From http://stackoverflow.com/questions/1597649/replace-strings-in-files-by-python
def file_replace(fname, out_fname, pat, s_after):
    if fname == out_fname:
        print "ERROR: input and output file names can not match"
        sys.exit(-1)
        return #input and output files can not match
    # first, see if the pattern is even in the file.
    with open(fname) as f:
        if not any(re.search(pat, line) for line in f):
            print "ERROR:  substitution pattern not found in reference file"
            sys.exit(-1)
            return # pattern does not occur in file so we are done.

    # pattern is in the file, so perform replace operation.
    with open(fname) as f:
        out = open(out_fname, "w")
        for line in f:
            out.write(re.sub(pat, s_after, line))
        out.close()

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
    parser.add_argument('--doshort', action='store', dest='doshort', default=True, help='If not present, do long')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    #parser.add_argument('-v', action='store_false', dest='verbose', default=True,
    #                    help='If not present, prints the locations')
    input_arguments = parser.parse_args()


    expConfig = ConfigParser.ConfigParser()
    expConfig.read(input_arguments.ExperimentConfig)

    # Experiment specific information
    subject_data_file=expConfig.get('EXPERIMENT_DATA','SESSION_DB')
    ExperimentName=expConfig.get('EXPERIMENT_DATA','EXPERIMENTNAME')
    WORKFLOW_COMPONENTS_STRING=expConfig.get('EXPERIMENT_DATA','WORKFLOW_COMPONENTS')
    WORKFLOW_COMPONENTS=eval(WORKFLOW_COMPONENTS_STRING)

    # Platform specific information
    #     Prepend the python search paths
    PYTHON_AUX_PATHS=expConfig.get(input_arguments.processingEnvironment,'PYTHON_AUX_PATHS')
    PYTHON_AUX_PATHS=PYTHON_AUX_PATHS.split(':')
    PYTHON_AUX_PATHS.extend(sys.path)
    sys.path=PYTHON_AUX_PATHS
    #print sys.path
    #import SimpleITK
    #     Prepend the shell environment search paths
    PROGRAM_PATHS=expConfig.get(input_arguments.processingEnvironment,'PROGRAM_PATHS')
    PROGRAM_PATHS=PROGRAM_PATHS.split(':')
    PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))
    os.environ['PATH']=':'.join(PROGRAM_PATHS)
    #    Define platform specific output write paths
    mountPrefix=expConfig.get(input_arguments.processingEnvironment,'MOUNTPREFIX')
    BASEOUTPUTDIR=expConfig.get(input_arguments.processingEnvironment,'BASEOUTPUTDIR')
    ExperimentBaseDirectoryPrefix=os.path.realpath(os.path.join(BASEOUTPUTDIR,ExperimentName))
    ExperimentBaseDirectoryCache=ExperimentBaseDirectoryPrefix+"_CACHE"
    ExperimentBaseDirectoryResults=ExperimentBaseDirectoryPrefix +"_Results"
    if not os.path.exists(ExperimentBaseDirectoryCache):
        os.makedirs(ExperimentBaseDirectoryCache)
    if not os.path.exists(ExperimentBaseDirectoryResults):
        os.makedirs(ExperimentBaseDirectoryResults)
    #    Define workup common reference data sets
    #    The ATLAS needs to be copied to the ExperimentBaseDirectoryPrefix
    #    The ATLAS pathing must stay constant
    ATLASPATH=expConfig.get(input_arguments.processingEnvironment,'ATLASPATH')
    if not os.path.exists(ATLASPATH):
        print("ERROR:  Invalid Path for Atlas: {0}".format(ATLASPATH))
        sys.exit(-1)
    CACHE_ATLASPATH=os.path.realpath(os.path.join(ExperimentBaseDirectoryCache,'Atlas'))
    from distutils.dir_util import copy_tree
    if not os.path.exists(CACHE_ATLASPATH):
        print("Copying a reference of the atlas to the experiment cache directory:\n    from: {0}\n    to: {1}".format(ATLASPATH,CACHE_ATLASPATH))
        copy_tree(ATLASPATH,CACHE_ATLASPATH,preserve_mode=1,preserve_times=1)
        ## Now generate the xml file with the correct pathing
        file_replace(os.path.join(ATLASPATH,'ExtendedAtlasDefinition.xml.in'),os.path.join(CACHE_ATLASPATH,'ExtendedAtlasDefinition.xml'),"@ATLAS_DIRECTORY@",CACHE_ATLASPATH)
    else:
        print("Atlas already exists in experiment cache directory: {0}".format(CACHE_ATLASPATH))
    #  Just to be safe, copy the model file as well
    BCDMODELPATH=expConfig.get(input_arguments.processingEnvironment,'BCDMODELPATH')
    CACHE_BCDMODELPATH=os.path.join(ExperimentBaseDirectoryCache,os.path.basename(BCDMODELPATH))
    from distutils.file_util import copy_file
    for BCDModelFiles in ['LLSModel-2ndVersion.h5','T1-2ndVersion.mdl']:
        orig=os.path.join(BCDMODELPATH,BCDModelFiles)
        new=os.path.join(CACHE_BCDMODELPATH,BCDModelFiles)
        if not os.path.exists(CACHE_BCDMODELPATH):
            os.mkdir(CACHE_BCDMODELPATH)
        if not os.path.exists(new):
            print("Copying BCD Model file to cache directory: {0}".format(new))
            copy_file(  orig, new,preserve_mode=1, preserve_times=1)
        else:
            print("BCD Model exists in cache directory: {0}".format(new))


    CUSTOM_ENVIRONMENT=expConfig.get(input_arguments.processingEnvironment,'CUSTOM_ENVIRONMENT')
    CUSTOM_ENVIRONMENT=eval(CUSTOM_ENVIRONMENT)
    ## Set custom environmental variables so that subproceses work properly (i.e. for Freesurfer)
    #print CUSTOM_ENVIRONMENT
    for key,value in CUSTOM_ENVIRONMENT.items():
        #print "SETTING: ", key, value
        os.putenv(key,value)
        os.environ[key]=value
    #print os.environ
    #sys.exit(-1)

    ## If freesurfer is requested, then ensure that a sane environment is available
    if 'FREESURFER' in WORKFLOW_COMPONENTS:
        print "FREESURFER NEEDS TO CHECK FOR SANE ENVIRONMENT HERE."

    CLUSTER_QUEUE=expConfig.get(input_arguments.processingEnvironment,'CLUSTER_QUEUE')

    ## Setup environment for CPU load balancing of ITK based programs.
    import multiprocessing
    total_CPUS=multiprocessing.cpu_count()
    if input_arguments.wfrun == 'helium_all.q':
        pass
    elif input_arguments.wfrun == 'helium_all.q_graph':
        pass
    elif input_arguments.wfrun == 'ipl_OSX':
        pass
    elif input_arguments.wfrun == 'local_4':
        os.environ['NSLOTS']="{0}".format(total_CPUS/4)
    elif input_arguments.wfrun == 'local_12':
        os.environ['NSLOTS']="{0}".format(total_CPUS/12)
    elif input_arguments.wfrun == 'local':
        os.environ['NSLOTS']="{0}".format(total_CPUS/1)
    else:
        print "FAILED RUN: You must specify the run environment type. [helium_all.q,helium_all.q_graph,ipl_OSX,local_4,local_12,local]"
        print input_arguments.wfrun
        sys.exit(-1)

    print "Configuring Pipeline"
    from nipype import config ## NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    config.enable_debug_mode() ## NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified

    import SessionDB
    subjectDatabaseFile=os.path.join( ExperimentBaseDirectoryCache,'InternalWorkflowSubjectDB.db')
    subject_list=input_arguments.subject.split(',')
    ExperimentDatabase=SessionDB.SessionDB(subjectDatabaseFile,subject_list)
    ## TODO:  Only make DB if db is older than subject_data_file.
    if ( not os.path.exists(subjectDatabaseFile) ) or ( os.path.getmtime(subjectDatabaseFile) < os.path.getmtime(subject_data_file) ):
        ExperimentDatabase.MakeNewDB(subject_data_file,mountPrefix)
    else:
        print("Using cached database, {0}".format(subjectDatabaseFile))
    print "ENTIRE DB for {_subjid}: ".format(_subjid=ExperimentDatabase.getSubjectFilter())
    print "^^^^^^^^^^^^^"
    for row in ExperimentDatabase.getEverything():
        print row
    print "^^^^^^^^^^^^^"

    import WorkupT1T2 ## NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    import ShortWorkupT1T2
    for subjectid in ExperimentDatabase.getAllSubjects():
        if input_arguments.doshort:
            baw200=ShortWorkupT1T2.ShortWorkupT1T2(subjectid,mountPrefix,
             os.path.join(ExperimentBaseDirectoryCache,str(subjectid)),
             ExperimentBaseDirectoryResults,
             ExperimentDatabase,
             CACHE_ATLASPATH,
             CACHE_BCDMODELPATH,WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS,CLUSTER_QUEUE=CLUSTER_QUEUE)
        else:
            baw200=WorkupT1T2.WorkupT1T2(subjectid,mountPrefix,
              os.path.join(ExperimentBaseDirectoryCache,str(subjectid)),
              ExperimentBaseDirectoryResults,
              ExperimentDatabase,
              CACHE_ATLASPATH,
              CACHE_BCDMODELPATH,WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS,CLUSTER_QUEUE=CLUSTER_QUEUE)
        print "Start Processing"

        ## Create the shell wrapper script for ensuring that all jobs running on remote hosts from SGE
        #  have the same environment as the job submission host.
        JOB_SCRIPT=get_global_sge_script(sys.path,PROGRAM_PATHS,CUSTOM_ENVIRONMENT)
        #print JOB_SCRIPT

        SGEFlavor='SGE'
        try:
            if input_arguments.wfrun == 'helium_all.q':
                baw200.run(plugin=SGEFlavor,
                    plugin_args=dict(template=JOB_SCRIPT,qsub_args="-S /bin/bash -pe smp1 2-12 -l mem_free=4000M -o /dev/null -e /dev/null "+CLUSTER_QUEUE))
            elif input_arguments.wfrun == 'helium_all.q_graph':
                SGEFlavor='SGEGraph' #Use the SGEGraph processing
                baw200.run(plugin=SGEFlavor,
                    plugin_args=dict(template=JOB_SCRIPT,qsub_args="-S /bin/bash -pe smp1 2-12 -l mem_free=4000M -o /dev/null -e /dev/null "+CLUSTER_QUEUE))
            elif input_arguments.wfrun == 'ipl_OSX':
                baw200.write_graph()
                print "Running On ipl_OSX"
                baw200.run(plugin=SGEFlavor,
                    plugin_args=dict(template=JOB_SCRIPT,qsub_args="-S /bin/bash -pe smp1 2-12 -l mem_free=4000M -o /dev/null -e /dev/null "+CLUSTER_QUEUE))
            elif input_arguments.wfrun == 'local_4':
                baw200.write_graph()
                print "Running with 4 parallel processes on local machine"
                baw200.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
            elif input_arguments.wfrun == 'local_12':
                baw200.write_graph()
                print "Running with 12 parallel processes on local machine"
                baw200.run(plugin='MultiProc', plugin_args={'n_procs' : 12})
            elif input_arguments.wfrun == 'local':
                try:
                    baw200.write_graph()
                except:
                    pass
                print "Running sequentially on local machine"
                baw200.run()
            else:
                print "You must specify the run environment type. [helium_all.q,helium_all.q_graph,ipl_OSX,local_4,local_12,local]"
                print input_arguments.wfrun
                sys.exit(-1)
        except:
            pass

if __name__ == "__main__":
    sys.exit(main())
