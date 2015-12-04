BRAINS AutoWorkup Processing pipeline, built on Nipype.

> We recommend to follow NAMICExternalProjects to build BRAINSTools with AutoWorkup.
> The very last part 'Additional Package Installation for BRAINS Auto Workup' of the NAMICExternalProject instruction is required.
> Please refer Readme.md file at the https://github.com/BRAINSia/NAMICExternalProjects

The file ```example/bawConfigurationTemplate.ini``` should be copied and modified to define the
experiment that should be run.

Our goal is to run following command successfully:

```sh
python baw_exp.py -processingLevel 2 -ExperimentConfig Template.config -pe IPL_ENVIRONMENT
```

# BRAINSAutoWorkUp Manual 
***(Under Construction by Regina EY Kim)***
BRAINSAutoWorkUp is developed for batch processing of brain MRI (T1-weighted or T1- and T2-weighted MRI set).
The pipeline is constructed based on NIPYPE.

## Requirements
### Requires NAMICExternalProjects build (https://github.com/BRAINSia/NAMICExternalProjects) which includes
* BRAINSTools build
* SimpleITK build
* NIPYPE

Let's say your NAMICExternalProjects build and located at:
```Shell
% NamicBuild=myNAMIC_build_location
```
### Python Environment
## Preparation
BRAINSAutoWorkup pipeline, located under `BRAINSTools/AutoWorkup/singleSession.py`, requires a configuration 
file with a list file for MRI inputs. 
### A input MRI list file
A list file specifies where T1- (and T2-)weighted images are in the form of python *dictionary*:
```
"projectName","subjectName","MRISessionName","{'T1-30':[t1 repeated scans list],'T2-30':[t2 repeated scans list]}'
```
A rough example of *bawInputList.csv*:
```Python
"projectA","111","12345611","{'T1-30':['somewhere/T1_filename.nii.gz'],'T2-30':['somewhere/T2_COR.nii.gz']}"
"projectB","112","12345614","{'T1-30':['somewhere/T1_filename.nii.gz', 'somewhere/T1_COR_REP1.nii.gz'],'T2-30':['somewhere/T2_COR.nii.gz']}"
"projectA","111","12345698","{'T1-30':['somewhere/T1_COR.nii.gz', 'somewhere/T1_COR_REP1.nii.gz'],'T2-30':['somewhere/T2_COR.nii.gz']}"
```
### Configuration file
A configuration file allows for HPC power utilization as well as local machine running. 
Syntax follows a standard INI configuration file format (Please see https://en.wikipedia.org/wiki/INI_file)
* EXPERIMENT section
```INI
[EXPERIMENT]
SESSION_DB_BASE=path_to_the_list/bawInputList.csv
SESSION_DB_TEMP=%(SESSION_DB_BASE)s
SESSION_DB_LONG=%(SESSION_DB_BASE)s

# A base output directory specification: 
BASE_OUTPUT_DIR=myExperimentDir

# A experiment directory name will be created following the below specifications
# e.g., myExperimentDir/20150715_BAW_Tutorial_base
# e.g., myExperimentDir/20150715_BAW_Tutorial_temp
# e.g., myExperimentDir/20150715_BAW_Tutorial_long
EXPERIMENT_BASE=20150715_BAW_Tutorial_base
EXPERIMENT_TEMP=20150715_BAW_Tutorial_temp
EXPERIMENT_LONG=20150715_BAW_Tutorial_long

# Our current template builder and longitudinal pipeline depend on the results of 
# baseline and template building pipeline. 
EXPERIMENT_TEMP_INPUT=%(EXPERIMENT_BASE)s
EXPERIMENT_LONG_INPUT=%(EXPERIMENT_TEMP)s

# Some, not all, components in the pipeline can be optional and
# they can be specified for each baseline/template builder/longitudinal pipeline:
WORKFLOW_COMPONENTS_BASE=['denoise','landmark','auxlmk','tissue_classify','warp_atlas_to_subject']
WORKFLOW_COMPONENTS_TEMP=[]
WORKFLOW_COMPONENTS_LONG=['denoise','landmark','auxlmk','tissue_classify','segmentation']

# Make sure you have Atlas_20131115 directory under your NAMICExternalProjects build directory.
#    BRAINSAutoWorkUp heavily depend on our atlas, which are highly adjusted
#    for reliable data processing
ATLAS_PATH=$NamicBuild/bin/Atlas/Atlas_20131115
```
* NIPYPE section
```INI
[NIPYPE]
GLOBAL_DATA_SINK_REWRITE=False
```
* DEFAULT section
```INI
[DEFAULT]
MOUNT_PREFIX=
MODULES=
```
* OS specific section:

ALL the OS specifications can be given one configuration file. This allows
developers to test/run/utilize as many computing resources as possible
in one experiment setup. 
      * OSX Example
      ```INI
      [OSX]
      ## The cluster queue to use for submitting "normal running" jobs.
      QUEUE=-q all
      ## The cluster queue to use for submitting "long running" jobs.
      QUEUE_LONG=-q all
      # The QSTAT command for immediate update of values [ use 'qstat' if in doubt ]
      QSTAT_IMMEDIATE=qstat
      QSTAT_CACHED=qstat
      ## Necessary modules to load for jobs
      MODULES=[]

      _GRAPHVIZ_BIN=/usr/bin/graphviz
      VIRTUALENV_DIR=your_virtualEnv_directory/BAWPythonEnv/MacEnv/
      # NAMICExternalProjects build tree
      _BUILD_DIR=your_path_to_NAMIC_built
      ############## -- You may not need to modify below here  ###########
      ```

      and then, you may choose one of the followings (A or B):

      A.
      ```INI
      ############## A -- if you built NAMICExternal without anaconda env (old) ###########
      _BRAINSTOOLS_BIN_DIR=%(_BUILD_DIR)s/bin
      _SIMPLEITK_PYTHON_LIB=%(_BUILD_DIR)s/lib
      _SIMPLEITK_PACKAGE_DIR=%(_BUILD_DIR)s/SimpleITK-build/Wrapping
      _NIPYPE_PACKAGE_DIR=%(_BUILD_DIR)s/NIPYPE
      ############## -- You should not need to modify below here. ###########
      APPEND_PYTHONPATH=%(_NIPYPE_PACKAGE_DIR)s:%(_SIMPLEITK_PYTHON_LIB)s:%(_SIMPLEITK_PACKAGE_DIR)s
      APPEND_PATH=%(_BRAINSTOOLS_BIN_DIR)s:%(_SIMPLEITK_PYTHON_LIB)s:%(_GRAPHVIZ_BIN)s
      ############## end of OSX specification ############## 
      ```

      B.
      ```INI
      ############## -- -- if you built NAMICExternal with anaconda env (newer)  ###########
      _BRAINSTOOLS_BIN_DIR=%(_BUILD_DIR)s/bin
      ############## -- You should not need to modify below here. ###########
      APPEND_PYTHONPATH=
      APPEND_PATH=%(_BRAINSTOOLS_BIN_DIR)s:%(_GRAPHVIZ_BIN)s
      ############## end of OSX specification ############## 
      ```

## Running on Mac OSX
```
env="OSX"

NamicBuild="/Shared/sinapse/scratch/eunyokim/src/NamicExternal/build_OSX_20150619/"
export PATH=$NamicBuild/bin:$PATH
export PYTHONPATH=$NamicBuild/SimpleITK-build/Wrapping/:$NamicBuild/BRAINSTools/AutoWorkup/:$NamicBuild/NIPYPE:$PYTHONPATH

wfrun="local"

python $NamicBuild/BRAINSTools/AutoWorkup/singleSession.py \
      --wfrun $wfrun \
      --workphase atlas-based-reference \
      --pe ${env} \
      --ExperimentConfig ./KIDsHD.config \
      --use-sentinal \
      all

```

## If BCD Fails, how to start the processing anyway!? 

In any failure case of BCD, BAW can be proceeded by 1) fix BCD manually and then 2) specify that manaull fixed file as an input.

- To fix manually, I ran following scripts in order:
cd /raid0/homes/aghayoor/BCD-test/Jira_cases/kidsFailure_ticket_PREDICTIMG-4209
* Step 1
```
#bash 1_reportFailure.sh
#!/bin/bash

llsModel='Atlas/20141004_BCD/LLSModel_50Lmks.h5'
atlasLandmarkWeights='Atlas/20141004_BCD/template_weights_50Lmks.wts'
atlasLandmarks='Atlas/20141004_BCD/template_landmarks_50Lmks.fcsv'
atlasVolume='Atlas/template_t1.nii.gz'
inputTemplateModel='Atlas/20141004_BCD/T1_50Lmks.mdl'
inputVolume='input_t1_to_be_fixed.nii.gz'

BRAINSConstellationDetector  \
  --LLSModel ${llsModel} \
  --acLowerBound 80.000000 \
  --atlasLandmarkWeights ${atlasLandmarkWeights} \
  --atlasLandmarks ${atlasLandmarks} \
  --atlasVolume ${atlasVolume} \
  --houghEyeDetectorMode 1 \
  --inputTemplateModel ${inputTemplateModel} \
  --inputVolume ${inputVolume} \
  --interpolationMode Linear \
  --forceHoughEyeDetectorReportFailure


  # README:
  # Remove "EMSP.fcsv" after this run!
```
* Step 2
```
#bash 2_runBCD_usingEMSP.sh
#!/bin/bash

llsModel='Atlas/20141004_BCD/LLSModel_50Lmks.h5'
atlasLandmarkWeights='Atlas/20141004_BCD/template_weights_50Lmks.wts'
atlasLandmarks='Atlas/20141004_BCD/template_landmarks_50Lmks.fcsv'
atlasVolume='Atlas/template_t1.nii.gz'
inputTemplateModel='Atlas/20141004_BCD/T1_50Lmks.mdl'
inputVolume='input_t1_to_be_fixed.nii.gz'

outputDir=$PWD/HDAdultAtlas_t1_20avg_manually_fixed
mkdir $outputDir
cd $outputDir

BRAINSConstellationDetector  \
  --LLSModel ${llsModel} \
  --acLowerBound 80.000000 \
  --atlasLandmarkWeights ${atlasLandmarkWeights} \
  --atlasLandmarks ${atlasLandmarks} \
  --atlasVolume ${atlasVolume} \
  --houghEyeDetectorMode 1 \
  --inputTemplateModel ${inputTemplateModel} \
  --inputVolume ../EMSP.nrrd \
  --interpolationMode Linear \
  --outputLandmarksInACPCAlignedSpace BCD_ACPC_Landmarks.fcsv \
  --outputLandmarksInInputSpace BCD_Original.fcsv \
  --outputResampledVolume BCD_ACPC.nii.gz \
  --outputTransform BCD_Original2ACPC_transform.h5 \
  --writeBranded2DImage BCD_Branded2DQCimage.png \
  --writedebuggingImagesLevel 3

```
* Step 3:
Copy "EMSP.nrrd"  file to the adjacent of original t1
Change BAW input csv file t1 to above EMSP.nrrd (Name can be anything)
``` 
touch ${inputVolume}_noDenoise
```

BCD should generate correct results if it is run on copied "EMSP.nrrd" file.
