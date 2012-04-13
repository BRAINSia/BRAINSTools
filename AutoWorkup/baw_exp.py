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
def get_global_sge_script(pythonPathsList,binPathsList):
    """This is a wrapper script for running commands on an SGE cluster
so that all the python modules and commands are pathed properly"""
    PYTHONPATH=":".join(pythonPathsList)
    BASE_BUILDS=":".join(binPathsList)
    GLOBAL_SGE_SCRIPT="""#!/bin/bash
echo "STARTED at: $(date +'%F-%T')"
echo "Ran on: $(hostname)"
export PATH={BINPATH}
export PYTHONPATH={PYTHONPATH}
echo "With PYTHONPATH={PYTHONPATH}"
echo "With PATH={BINPATH}"
## NOTE:  nipype inserts the actaul commands that need running below this section.
""".format(PYTHONPATH=PYTHONPATH,BINPATH=BASE_BUILDS)
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
    group.add_argument('-processingLevel', action="store", dest='processingLevel', required=False,
                       help='How much processing should be done, 1=basic, 2= include ABC, 3= include ANTS, 4= include Freesurfer',
                       default=2)
    #parser.add_argument('-v', action='store_false', dest='verbose', default=True,
    #                    help='If not present, prints the locations')
    input_arguments = parser.parse_args()


    expConfig = ConfigParser.ConfigParser()
    expConfig.read(input_arguments.ExperimentConfig)
    
    # Experiment specific information
    session_db=expConfig.get('EXPERIMENT_DATA','SESSION_DB')
    ExperimentName=expConfig.get('EXPERIMENT_DATA','EXPERIMENTNAME')
    
    # Platform specific information
    #     Prepend the python search paths
    PYTHON_AUX_PATHS=expConfig.get(input_arguments.processingEnvironment,'PYTHON_AUX_PATHS')
    PYTHON_AUX_PATHS=PYTHON_AUX_PATHS.split(';')
    PYTHON_AUX_PATHS.extend(sys.path)
    sys.path=PYTHON_AUX_PATHS
    #     Prepend the shell environment search paths
    PROGRAM_PATHS=expConfig.get(input_arguments.processingEnvironment,'PROGRAM_PATHS')
    PROGRAM_PATHS=PROGRAM_PATHS.split(';')
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

    print "Configuring Pipeline"
    import WorkupT1T2 ## NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
    baw200=WorkupT1T2.WorkupT1T2(input_arguments.processingLevel, mountPrefix,
      ExperimentBaseDirectory,
      session_db,
      ATLASPATH,
      BCDMODELPATH)
    
    print "Start Processing"

    ## Create the shell wrapper script for ensuring that all jobs running on remote hosts from SGE
    #  have the same environment as the job submission host.
    JOB_SCRIPT=get_global_sge_script(sys.path,PROGRAM_PATHS)
    if input_arguments.wfrun == 'helium_all.q':
        baw200.run(plugin='SGE',
            plugin_args=dict(template=JOB_SCRIPT,qsub_args="-S /bin/bash -q UI -pe smp1 2-4 -o /dev/null -e /dev/null "))
    elif input_arguments.wfrun == 'ipl_OSX':
        baw200.run(plugin='SGE',
            plugin_args=dict(template=JOB_SCRIPT,qsub_args="-S /bin/bash -q OSX -pe smp1 2-4 -o /dev/null -e /dev/null "))
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
	sys.exit(-1)
	
    #baw200.write_graph()

if __name__ == "__main__":
    sys.exit(main())
