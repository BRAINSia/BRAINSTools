import os
import argparse
import textwrap
import csv
import sys

class runOneAW():

    def main(self):
        self.sessionPath = os.path.join(input_arguments.experimentOutputDir, 'session.csv')
        self.configPath = os.path.join(input_arguments.experimentOutputDir, 'localAW.config')
        self.generateSessionCSV()
        self.generateConfigFile()
        self.executeAW()

    def executeAW(self):
        bawCommand = """time python /raid0/homes/jforbes/git/BRAINSStandAlone/AutoWorkup/baw_exp.py \
 -ExperimentConfig {configFile} \
 -pe LOCAL_ENVIRONMENT \
 -wfrun local \
 -subject {subject} \n""".format(configFile=self.configPath, subject=input_arguments.subject)
        print '-'*80
        print '\nExecuting command: \n{bawCommand}'.format(bawCommand=bawCommand)
        os.system(bawCommand)

    def generateSessionCSV(self):
        sessionDict = dict()
        if input_arguments.t1 != []:
            sessionDict['T1-30'] = input_arguments.t1
        if input_arguments.t2 != []:
            sessionDict['T2-30'] = input_arguments.t2
        if sessionDict == dict():
            print 'ERROR: No T1 or T2 images were given as input arguments.'
            sys.exit()
        col_name_list = ["project", "subject", "session", "imagefiles"]
        newFile = csv.writer(open(self.sessionPath, 'wb'), quoting=csv.QUOTE_ALL)
        newFile.writerow(col_name_list)
        line = (input_arguments.project, input_arguments.subject, input_arguments.session, sessionDict)
        newFile.writerow(line)
        print '\nThe session csv file has been generated: {0}\n'.format(self.sessionPath)
        print line

    def generateConfigFile(self):
        configString = """
###  INTENT  ###
# The intent of this configuration file is to define all the information needed to
# completely run an experiment.
#
# The configuration file consists of sections,
# led by a [section] header and followed by name: value entries,
# with continuations in the style of RFC 822 (see section 3.1.1)
#
### How to run an experiment:
#  python ${CHANGEME_PATH_TO_BRAINSSTANDALONE_SOURCE}/AutoWorkup/baw_exp.py \
#               -ExperimentConfig ${1} -pe ${2} -wfrun ${3} -subject ${4}
#          ${1} -- This file
#          ${2} -- One of the defined environments from this file
#          ${3} -- One of the following processing engines:
#                      helium_all.q = only works on Helium cluster with standard SGE
#                      helium_all.q_graph = this is BETA, only on Helium cluster with job-dependent SGE
#                      ipl_OSX = PINC lab with standard SGE
#                      local = local machine with serial node completion
#                      local_4 = local machine with up to 4 simultaneous Nipype processing nodes
#                      local_12 = local machine with up to 12 simultaneous Nipype processing nodes
#          ${4} -- A comma separated list of subjects to process from the SESSION_DB (all others are ignored)
#
#  Conventions of this file:
#         1) variables with leading '_' are optional convenience
#            used to simplyfy other variables and are only
#            internal to this file
#         2) Variables beginging with a letter are mandatory
#            and are read and used by the python baw_exp.py
#            function

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!         Describe the experiment to run                !!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

[EXPERIMENT_DATA]
#   The SESSION_DB is the file path to a formatted csv file
#   describing the data sets one per line, where the different
#   data types are listed in a dictionary.
#     -- Empty lines in the csv file are ignored
#     -- Lines in the file that begin with '#' are ignored
#  'PROJECT','SUBJ','SESSION',"{'T1-30:['{fullpath}/T1_1.nii.gz','{fullpath}/T1_2.nii.gz'],'T2-30':['{fullpath}/T2_1.nii.gz']}"
SESSION_DB=[replaceme_sessionDB]
# The desired output directory for this experiment
EXPERIMENTNAME=TutorialExperimentOutputs
# Components of pipeline to run.  There are some branches of the workflow that are mostly for validation and experimentation.
#WORKFLOW_COMPONENTS=['BASIC','TISSUE_CLASSIFY','SEGMENTATION','FREESURFER','ANTS','AUXLMK']
#WORKFLOW_COMPONENTS=['BASIC','TISSUE_CLASSIFY']
WORKFLOW_COMPONENTS=['BASIC']

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!   Describe the processing environments that will be   !!!!
# !!!!   used to complete this experiment run                !!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

[LOCAL_ENVIRONMENT]
# Where to find Freesurfer base directory
_FREESURFER_HOME=/ipldev/sharedopt/20110601/MacOSX_10.6/freesurfer
_FREESURFER_PATH_EXT=%(_FREESURFER_HOME)s/bin:%(_FREESURFER_HOME)s/fsfast/bin:%(_FREESURFER_HOME)s/tktools:/opt/fsl/bin:%(_FREESURFER_HOME)s/mni/bin
_GRAPHVIZ_BIN_DIR=/usr/local
_BRAINSTOOLS_BUILD_PATH=/raid0/homes/jforbes/git/BRAINSStandAlone-build
# The prefix to add to all image files in the $(SESSION_DB) to account
# for different file system mount points
MOUNTPREFIX=
# The base directory where all experiments of this type go
BASEOUTPUTDIR=[replaceme_outputDir]
#CUSTOM_ENVIRONMENT  This is mostly needed for freesurfer, but may be used for fsl or other programs
CUSTOM_ENVIRONMENT={'FREESURFER_HOME':'%(_FREESURFER_HOME)s', 'FIX_VERTEX_AREA':'', 'FMRI_ANALYSIS_DIR':'%(_FREESURFER_HOME)s/fsfast', 'FSFAST_HOME':'%(_FREESURFER_HOME)s/fsfast', 'FSF_OUTPUT_FORMAT':'nii.gz', 'FSL_BIN':'%(_FREESURFER_HOME)s/bin', 'FSL_DIR':'%(_FREESURFER_HOME)s', 'FS_OVERRIDE':'0', 'FUNCTIONALS_DIR':'%(_FREESURFER_HOME)s/sessions', 'LOCAL_DIR':'%(_FREESURFER_HOME)s/local', 'MINC_BIN_DIR':'%(_FREESURFER_HOME)s/mni/bin', 'MINC_LIB_DIR':'%(_FREESURFER_HOME)s/mni/lib', 'MNI_DATAPATH':'%(_FREESURFER_HOME)s/mni/data', 'MNI_DIR':'%(_FREESURFER_HOME)s/mni', 'MNI_PERL5LIB':'%(_FREESURFER_HOME)s/mni/lib/../System/Library/Perl/5.8.6', 'PERL5LIB':'%(_FREESURFER_HOME)s/mni/lib/../System/Library/Perl/5.8.6', 'SUBJECTS_DIR':'%(_FREESURFER_HOME)s/subjects'}
## The cluster queue to use for submitting jobs.
CLUSTER_QUEUE=-q OSX

############## -- You should not need to modify below here. ###########
# Where to find the Autoworkup scripts.
_BRAINSTOOLS_SCRIPTS=[replaceme_brainToolsDir]/AutoWorkup
# Where SimpleITK, required nipype, and other non-standard python packages reside
## _SIMPLE_ITK_PYTHON_PATH=/${CHANGEME_extra_apps_bin_dir}/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages
_SIMPLE_ITK_PYTHON_PATH=%(_BRAINSTOOLS_BUILD_PATH)s/SimpleITK-build/lib
_NIPYPE_CUSTOM_PATH=%(_BRAINSTOOLS_BUILD_PATH)s/NIPYPE
PYTHON_AUX_PATHS=%(_NIPYPE_CUSTOM_PATH)s:%(_BRAINSTOOLS_SCRIPTS)s:%(_SIMPLE_ITK_PYTHON_PATH)s
# Where to find BRAINSTools
# Paths that need to be configured to find tools needed to run this workflow
PROGRAM_PATHS=%(_BRAINSTOOLS_BUILD_PATH)s/lib:%(_BRAINSTOOLS_BUILD_PATH)s/bin:%(_BRAINSTOOLS_SCRIPTS)s/BRAINSTools:%(_FREESURFER_PATH_EXT)s:%(_GRAPHVIZ_BIN_DIR)s
# The path to the reference atlas spact to be used in this analysis by all BRAINSTools
ATLASPATH=%(_BRAINSTOOLS_BUILD_PATH)s/ReferenceAtlas-build/Atlas/Atlas_20120830
# The path to the model files to be used by BCD.
BCDMODELPATH=%(_BRAINSTOOLS_BUILD_PATH)s/BRAINSTools-build/TestData"""
        firstReplace = configString.replace('[replaceme_sessionDB]',self.sessionPath)
        secondReplace = firstReplace.replace('[replaceme_outputDir]',input_arguments.experimentOutputDir)
        newConfigString = secondReplace.replace('[replaceme_brainToolsDir]',input_arguments.brainsToolsDir)
        handle = open(self.configPath, 'w')
        handle.write(newConfigString)
        handle.close()
        print '\nThe configuration file has been generated: {0}'.format(self.configPath)
        print newConfigString

if __name__ == "__main__":
    # Create and parse input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""
This program is used to run a mini version of AutoWorkup for a specific session and specified T1 and T2 files.

Example Usage:
$ python runOneAW.py -project NAMICHD -subject 10066 -session 100261 \
-t1 /hjohnson/HDNI/AutoWorkupTutorial/10066/100261/ANONRAW/10066_100261_T1-30_2.nii.gz \
/hjohnson/HDNI/AutoWorkupTutorial/10066/100261/ANONRAW/10066_100261_T1-30_22.nii.gz \
-t2 /hjohnson/HDNI/AutoWorkupTutorial/10066/100261/ANONRAW/10066_100261_T2-30_3.nii.gz \
-experimentOutputDir /scratch/autoworkup/runOneAW/test \
-brainsToolsDir /IPLlinux/raid0/homes/jforbes/git/BRAINSStandAlone"""))
    group = parser.add_argument_group('Required')
    group.add_argument('-project', action="store", dest='project', required=True,
                       help='The name of the project to process')
    group.add_argument('-subject', action="store", dest='subject', required=True,
                       help='The name of the subject to process')
    group.add_argument('-session', action="store", dest='session', required=True,
                       help='The name of the session to process')
    group.add_argument('-t1', action="store", dest='t1', nargs='*', default=[],
                       help='The file name(s) of the T1 image(s) to process')
    group.add_argument('-t2', action="store", dest='t2', nargs='*', default=[],
                       help='The file name(s) of the T2 image(s) to process')
    group.add_argument('-experimentOutputDir', action="store", dest='experimentOutputDir', required=True,
                       help='The directory for the experiment output')
    group.add_argument('-brainsToolsDir', action="store", dest='brainsToolsDir', required=True,
                       help='The directory for BRAINSSTANDALONE')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    #parser.add_argument('-v', action='store_false', dest='verbose', default=True,
    #                    help='If not present, prints the locations')
    input_arguments = parser.parse_args()

    runObject = runOneAW()
    runObject.main()
