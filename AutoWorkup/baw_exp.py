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
echo "With PYTHONPATH={PYTHONPATH}"
echo "With PATH={BINPATH}"
echo "With custom environment:"
echo {CUSTENV}
{CUSTENV}
## NOTE:  nipype inserts the actaul commands that need running below this section.
""".format(PYTHONPATH=PYTHONPATH,BINPATH=BASE_BUILDS,CUSTENV=custEnvString)
    return GLOBAL_SGE_SCRIPT

def main(argv=None):
    import argparse
    import ConfigParser
    import os
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
    group.add_argument('-ExperimentConfig', action="store", dest='ExperimentConfig', required=True,
                       help='The path to the file that describes the entire experiment')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    #parser.add_argument('-v', action='store_false', dest='verbose', default=True,
    #                    help='If not present, prints the locations')
    input_arguments = parser.parse_args()


    expConfig = ConfigParser.ConfigParser()
    expConfig.read(input_arguments.ExperimentConfig)

    # Experiment specific information
    session_db=expConfig.get('EXPERIMENT_DATA','SESSION_DB')
    ExperimentName=expConfig.get('EXPERIMENT_DATA','EXPERIMENTNAME')
    WORKFLOW_COMPONENTS_STRING=expConfig.get('EXPERIMENT_DATA','WORKFLOW_COMPONENTS')
    WORKFLOW_COMPONENTS=eval(WORKFLOW_COMPONENTS_STRING)

    # Platform specific information
    #     Prepend the python search paths
    PYTHON_AUX_PATHS=expConfig.get(input_arguments.processingEnvironment,'PYTHON_AUX_PATHS')
    PYTHON_AUX_PATHS=PYTHON_AUX_PATHS.split(':')
    PYTHON_AUX_PATHS.extend(sys.path)
    sys.path=PYTHON_AUX_PATHS
    #     Prepend the shell environment search paths
    PROGRAM_PATHS=expConfig.get(input_arguments.processingEnvironment,'PROGRAM_PATHS')
    PROGRAM_PATHS=PROGRAM_PATHS.split(':')
    PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))
    os.environ['PATH']=':'.join(PROGRAM_PATHS)
    #    Define platform specific output write paths
    mountPrefix=expConfig.get(input_arguments.processingEnvironment,'MOUNTPREFIX')
    BASEOUTPUTDIR=expConfig.get(input_arguments.processingEnvironment,'BASEOUTPUTDIR')
    FULL_EXPERIMENT_OUTPUTDIR=os.path.join(BASEOUTPUTDIR,ExperimentName)
    if not os.path.exists(FULL_EXPERIMENT_OUTPUTDIR):
        os.makedirs(FULL_EXPERIMENT_OUTPUTDIR)
    ExperimentBaseDirectory=os.path.realpath(os.path.join(BASEOUTPUTDIR,ExperimentName))
    #    Define workup common reference data sets
    ATLASPATH=expConfig.get(input_arguments.processingEnvironment,'ATLASPATH')
    BCDMODELPATH=expConfig.get(input_arguments.processingEnvironment,'BCDMODELPATH')
    CUSTOM_ENVIRONMENT=expConfig.get(input_arguments.processingEnvironment,'CUSTOM_ENVIRONMENT')
    CUSTOM_ENVIRONMENT=eval(CUSTOM_ENVIRONMENT)
    ## Set custom environmental variables so that subproceses work properly (i.e. for Freesurfer)
    #print CUSTOM_ENVIRONMENT
    for key,value in CUSTOM_ENVIRONMENT.items():
        #print "SETTING: ", key, value
        os.putenv(key,value)
        os.environ[key]=value
    print os.environ
    #sys.exit(-1)

    ## If freesurfer is requested, then ensure that a sane environment is available
    if 'FREESURFER' in WORKFLOW_COMPONENTS:
        print "FREESURFER NEEDS TO CHECK FOR SANE ENVIRONMENT HERE."

    CLUSTER_QUEUE=expConfig.get(input_arguments.processingEnvironment,'CLUSTER_QUEUE')

    print "Configuring Pipeline"
    import WorkupT1T2 ## NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    baw200=WorkupT1T2.WorkupT1T2(mountPrefix,
      ExperimentBaseDirectory,
      session_db,
      ATLASPATH,
      BCDMODELPATH,WORKFLOW_COMPONENTS=WORKFLOW_COMPONENTS,CLUSTER_QUEUE=CLUSTER_QUEUE)
    print "Start Processing"

    ## Create the shell wrapper script for ensuring that all jobs running on remote hosts from SGE
    #  have the same environment as the job submission host.
    JOB_SCRIPT=get_global_sge_script(sys.path,PROGRAM_PATHS,CUSTOM_ENVIRONMENT)
    print JOB_SCRIPT
    baw200.write_graph()
    if input_arguments.wfrun == 'helium_all.q':
        baw200.run(plugin='SGE',
            plugin_args=dict(template=JOB_SCRIPT,qsub_args="-S /bin/bash -pe smp1 1-4 -l mem_free=4000M -o /dev/null -e /dev/null "+CLUSTER_QUEUE))
    elif input_arguments.wfrun == 'ipl_OSX':
        print "Running On ipl_OSX"
        baw200.run(plugin='SGE',
            plugin_args=dict(template=JOB_SCRIPT,qsub_args="-S /bin/bash -pe smp1 1-4 -l mem_free=4000M -o /dev/null -e /dev/null "+CLUSTER_QUEUE))
    elif input_arguments.wfrun == 'local_4':
        print "Running with 4 parallel processes on local machine"
        baw200.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
    elif input_arguments.wfrun == 'local_12':
        print "Running with 12 parallel processes on local machine"
        baw200.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
    elif input_arguments.wfrun == 'local':
        print "Running sequentially on local machine"
        baw200.run()
    else:
        print "You must specify the run environment type."


if __name__ == "__main__":
    sys.exit(main())
