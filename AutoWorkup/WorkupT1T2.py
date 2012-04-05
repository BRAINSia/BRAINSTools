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
import csv
import sys
import string
import argparse
#"""Import necessary modules from nipype."""
#from nipype.utils.config import config
#config.set('logging', 'log_to_file', 'false')
#config.set_log_dir(os.getcwd())
#--config.set('logging', 'workflow_level', 'DEBUG')
#--config.set('logging', 'interface_level', 'DEBUG')
#--config.set('execution','remove_unnecessary_outputs','false')

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import ReconAll

from nipype.utils.misc import package_check
#package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
package_check('numpy', '1.3', 'tutorial1')
package_check('scipy', '0.7', 'tutorial1')
package_check('networkx', '1.0', 'tutorial1')
package_check('IPython', '0.10', 'tutorial1')

from BRAINSTools.BRAINSConstellationDetector import *
from BRAINSTools.BRAINSABC import *
from BRAINSTools.BRAINSDemonWarp import *
from BRAINSTools.BRAINSFit import *
from BRAINSTools.BRAINSMush import *
from BRAINSTools.BRAINSResample import *
from BRAINSTools.BRAINSROIAuto import *
from BRAINSTools.BRAINSLandmarkInitializer import *
from BRAINSTools.BRAINSCut import *
from BRAINSTools.GradientAnisotropicDiffusionImageFilter import *
from BRAINSTools.GenerateSummedGradientImage import *
from BRAINSTools.ANTSWrapper import *
from BRAINSTools.WarpImageMultiTransform import *
from BRAINSTools.WarpAllAtlas import *

#######################  HACK:  Needed to make some global variables for quick
#######################         processing needs
#Generate by running a file system list "ls -1 $AtlasDir *.nii.gz *.xml *.fcsv *.wgts"
atlas_file_list="AtlasPVDefinition.xml ALLPVAIR.nii.gz ALLPVBASALTISSUE.nii.gz ALLPVCRBLGM.nii.gz ALLPVCRBLWM.nii.gz ALLPVCSF.nii.gz ALLPVNOTCSF.nii.gz ALLPVNOTGM.nii.gz ALLPVNOTVB.nii.gz ALLPVNOTWM.nii.gz ALLPVSURFGM.nii.gz ALLPVVB.nii.gz ALLPVWM.nii.gz avg_t1.nii.gz avg_t2.nii.gz tempNOTVBBOX.nii.gz template_ABC_lables.nii.gz template_WMPM2_labels.nii.gz template_WMPM2_labels.txt template_brain.nii.gz template_cerebellum.nii.gz template_class.nii.gz template_headregion.nii.gz template_leftHemisphere.nii.gz template_nac_lables.nii.gz template_nac_lables.txt template_rightHemisphere.nii.gz template_t1.nii.gz template_t1_clipped.nii.gz template_t2.nii.gz template_t2_clipped.nii.gz template_ventricles.nii.gz probabilityMaps/l_caudate_ProbabilityMap.nii.gz probabilityMaps/r_caudate_ProbabilityMap.nii.gz probabilityMaps/l_hippocampus_ProbabilityMap.nii.gz probabilityMaps/r_hippocampus_ProbabilityMap.nii.gz probabilityMaps/l_putamen_ProbabilityMap.nii.gz probabilityMaps/r_putamen_ProbabilityMap.nii.gz probabilityMaps/l_thalamus_ProbabilityMap.nii.gz probabilityMaps/r_thalamus_ProbabilityMap.nii.gz spatialImages/phi.nii.gz spatialImages/rho.nii.gz spatialImages/theta.nii.gz"
atlas_file_names=atlas_file_list.split(' ')
atlas_file_names=["AtlasPVDefinition.xml","ALLPVAIR.nii.gz",
                      "ALLPVBASALTISSUE.nii.gz","ALLPVCRBLGM.nii.gz",
                      "ALLPVCRBLWM.nii.gz","ALLPVCSF.nii.gz","ALLPVNOTCSF.nii.gz",
                      "ALLPVNOTGM.nii.gz","ALLPVNOTVB.nii.gz","ALLPVNOTWM.nii.gz",
                      "ALLPVSURFGM.nii.gz","ALLPVVB.nii.gz","ALLPVWM.nii.gz",
                      "avg_t1.nii.gz","avg_t2.nii.gz","tempNOTVBBOX.nii.gz",
                      "template_ABC_lables.nii.gz","template_WMPM2_labels.nii.gz",
                      "template_WMPM2_labels.txt","template_brain.nii.gz",
                      "template_cerebellum.nii.gz","template_class.nii.gz",
                      "template_headregion.nii.gz","template_leftHemisphere.nii.gz",
                      "template_nac_lables.nii.gz","template_nac_lables.txt",
                      "template_rightHemisphere.nii.gz","template_t1.nii.gz",
                      "template_t1_clipped.nii.gz","template_t2.nii.gz",
                      "template_t2_clipped.nii.gz","template_ventricles.nii.gz",
                      "template_landmarks.fcsv","template_landmark_weights.csv",
"probabilityMaps/l_caudate_ProbabilityMap.nii.gz",
"probabilityMaps/r_caudate_ProbabilityMap.nii.gz",
"probabilityMaps/l_hippocampus_ProbabilityMap.nii.gz",
"probabilityMaps/r_hippocampus_ProbabilityMap.nii.gz",
"probabilityMaps/l_putamen_ProbabilityMap.nii.gz",
"probabilityMaps/r_putamen_ProbabilityMap.nii.gz",
"probabilityMaps/l_thalamus_ProbabilityMap.nii.gz",
"probabilityMaps/r_thalamus_ProbabilityMap.nii.gz",
"spatialImages/phi.nii.gz",
"spatialImages/rho.nii.gz",
"spatialImages/theta.nii.gz"
                      ]


## Remove filename extensions for images, but replace . with _ for other file types
atlas_file_keys=[os.path.basename(fn).replace('.nii.gz','').replace('.','_') for fn in atlas_file_names]
atlas_outputs_filename_match = dict(zip(atlas_file_keys,atlas_file_names))


#############################################################################
#############################################################################
## Utility functions for the pipeline
#############################################################################
#############################################################################
def get_first_T1_and_T2(in_files,T1_count):
    '''
    Returns the first T1 and T2 file in in_files, based on offset in T1_count.
    '''
    return in_files[0],in_files[T1_count]

def GetExtensionlessBaseName(filename):
    '''
    Get the filename without the extension.  Works for .ext and .ext.gz
    '''
    import os
    currBaseName = os.path.basename(filename)
    currExt = os.path.splitext(currBaseName)[1]
    currBaseName = os.path.splitext(currBaseName)[0]
    if currExt == ".gz":
        currBaseName = os.path.splitext(currBaseName)[0]
        currExt = os.path.splitext(currBaseName)[1]
    return currBaseName


def MakeAtlasNode(atlasDirectory):
    BAtlas = pe.Node(interface=nio.DataGrabber(outfields=atlas_file_keys),
                                               name='BAtlas')
    BAtlas.inputs.base_directory = atlasDirectory
    BAtlas.inputs.template = '*'
    ## Prefix every filename with atlasDirectory
    atlas_search_paths=['{0}'.format(fn) for fn in atlas_file_names]
    BAtlas.inputs.field_template = dict(zip(atlas_file_keys,atlas_search_paths))
    ## Give 'atlasDirectory' as the substitution argument
    atlas_template_args_match=[ [[]] for i in atlas_file_keys ] ##build a list of proper lenght with repeated entries
    BAtlas.inputs.template_args = dict(zip(atlas_file_keys,atlas_template_args_match))
    return BAtlas

def create_BRAINSCut_XML(rho,phi,theta,model,
                         r_probabilityMap,l_probabilityMap,
                         atlasT1,atlasBrain,
                         subjT1,subjT2,
                         subjT1GAD,subjT2GAD,subjSGGAD,subjBrain,
                         atlasToSubj,
                         output_dir):
    import re
    import os
    print "*"*80
    print rho
    print phi
    print theta
    print model
    print r_probabilityMap
    print l_probabilityMap
    structure = re.search("r_(\w+)_ProbabilityMap",os.path.basename(r_probabilityMap)).group(1)

    ## The model file name is auto-generated, and needs to be split apart here
    basemodel      =re.search("(.*{structure}Model.*\.txt)(00[0-9]*)".format(structure=structure),model).group(1)
    EpochIterations=re.search("(.*{structure}Model.*\.txt)(00[0-9]*)".format(structure=structure),model).group(2)
    EpochIterations.lstrip('0')

    ## HACK:  Needed to make neural networks easier to use.  This information should be embeded in the model file.
    HiddenNodeDict={'caudate':"14",'putamen':"20",'hippocampus':"20",'thalamus':"14"}

    print "^^"*80
    print "^^"*80
    print "^^"*80
    print "^^"*80
    print basemodel
    print EpochIterations
    NumberOfHiddenNodes=HiddenNodeDict[structure]

    EXTRA_FLAGS=""
    if structure in [ 'putamen','hippocampus']:
        EXTRA_FLAGS="""
     <Image Type="T1GAD" Filename="na" />
     <Image Type="T2GAD" Filename="na" />"""

    XMLSTRING="""<AutoSegProcessDescription>
  <RegistrationConfiguration
         ImageTypeToUse="T1"
         ID="BSpline_ROI"
         BRAINSROIAutoDilateSize="1"
  />
  <ANNParams
         Iterations="{EpochIterations}"
         MaximumVectorsPerEpoch="700000"
         EpochIterations="{EpochIterations}"
         ErrorInterval="1"
         DesiredError="0.000001"
         NumberOfHiddenNodes="{NumberOfHiddenNodes}"
         ActivationSlope="1.0"
         ActivationMinMax="1.0"
  />
  <ApplyModel
         CutOutThresh="0.05"
         MaskThresh="0.4"
         LevelSetImageType="NA"
         GaussianSmoothingSigma="0.0"
  />
  <NeuralNetParams
         MaskSmoothingValue="0.0"
         GradientProfileSize="1"
         TrainingVectorFilename="na"
         TrainingModelFilename="{basemodel}"
         TestVectorFilename="na"
         Normalization="true"
  />
  <ProbabilityMap StructureID="l_{structure}" Gaussian="0.5" GenerateVector="true" Filename="{l_probabilityMap}"/>
  <ProbabilityMap StructureID="r_{structure}" Gaussian="0.5" GenerateVector="true" Filename="{r_probabilityMap}"/>
  <DataSet Type="Atlas" Name="template">
    <Image Type="T1" Filename="{atlasT1}"/>
    <Image Type="T2" Filename="na"/>
    {EXTRA_FLAGS}
    <Image Type="SGGAD" Filename="na"/>

    <SpatialLocation Type="rho"   Filename="{rho}"/>
    <SpatialLocation Type="phi"   Filename="{phi}"/>
    <SpatialLocation Type="theta" Filename="{theta}"/>
  </DataSet>
  <DataSet Name="sessionID" Type="Apply" OutputDir="./">
      <Image Type="T1" Filename="{subjT1}"/>
      <Image Type="T2" Filename="{subjT2}"/>
      <Image Type="T1GAD" Filename="{subjT1GAD}"/>
      <Image Type="T2GAD" Filename="{subjT2GAD}"/>
      <Image Type="SGGAD" Filename="{subjSGGAD}"/>
      <Mask Type="l_{structure}" Filename="{output_dir}/l_{structure}_seg.nii.gz"/>
      <Mask Type="r_{structure}" Filename="{output_dir}/r_{structure}_seg.nii.gz"/>
      <Registration SubjToAtlasRegistrationFilename="na"
                    AtlasToSubjRegistrationFilename="{atlasToSubj}"
                    SubjectBinaryFilename="{subjBrain}"
                    AtlasBinaryFilename="{atlasBrain}"
                    ID="BSpline_ROI"/>
  </DataSet>
</AutoSegProcessDescription>

""".format(structure=structure,rho=rho,phi=phi,theta=theta,
                         basemodel=basemodel,EpochIterations=EpochIterations,NumberOfHiddenNodes=NumberOfHiddenNodes,
                         r_probabilityMap=r_probabilityMap,l_probabilityMap=l_probabilityMap,
                         atlasT1=atlasT1,atlasBrain=atlasBrain,
                         subjT1=subjT1,subjT2=subjT2,
                         subjT1GAD=subjT1GAD,subjT2GAD=subjT2GAD,subjSGGAD=subjSGGAD,subjBrain=subjBrain,
                         EXTRA_FLAGS=EXTRA_FLAGS,
                         atlasToSubj=atlasToSubj,
                         output_dir=output_dir)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    #xml_filename = os.path.join(output_dir,'%s.xml' % structure)
    xml_filename = '%s.xml' % structure
    xml_file = open(xml_filename,'w')
    #xml_file.write(etree.tostring(xml_output, pretty_print=True))
    xml_file.write(XMLSTRING)
    xml_file.close()

    r_struct_fname="{output_dir}/r_{structure}_seg.nii.gz".format(output_dir=output_dir,structure=structure)
    l_struct_fname="{output_dir}/l_{structure}_seg.nii.gz".format(output_dir=output_dir,structure=structure)
    return os.path.realpath(xml_filename), [ r_struct_fname, l_struct_fname ]


def get_list_element( nestedList, index ):
    return nestedList[index]

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def getFirstT1(uid, dbfile):
    from cPickle import load
    with open(dbfile) as fp:
        db = load(fp)
    print("uid:= {0}, dbfile: {1}".format(uid,dbfile))
    print("result:= {0}".format(db[uid]["T1s"]))
    return db[uid]["T1s"][0]

def getT1s(uid, dbfile):
    from cPickle import load
    with open(dbfile) as fp:
        db = load(fp)
    #print("uid:= {0}, dbfile: {1}".format(uid,dbfile))
    #print("result:= {0}".format(db[uid]["T1s"]))
    return db[uid]["T1s"]

def getT1sLength(uid, dbfile):
    from cPickle import load
    with open(dbfile) as fp:
        db = load(fp)
    #print("uid:= {0}, dbfile: {1}".format(uid,dbfile))
    #print("result:= {0}".format(db[uid]["T1s"]))
    return len(db[uid]["T1s"])

def getT2s(uid, dbfile):
    from cPickle import load
    with open(dbfile) as fp:
        db = load(fp)
    #print("uid:= {0}, dbfile: {1}".format(uid,dbfile))
    #print("result:= {0}".format(db[uid]["T1s"]))
    return db[uid]["T2s"]

def getT1sT2s(uid, dbfile,altT1):
    from cPickle import load
    with open(dbfile) as fp:
        db = load(fp)
    #print("uid:= {0}, dbfile: {1}".format(uid,dbfile))
    #print("result:= {0}".format(db[uid]["T1s"]))
    temp=db[uid]["T1s"]
    temp.append(db[uid]["T2s"])
    temp[0]=altT1
    return temp




###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## WorkupT1T2 is the main workflow to be run
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
def WorkupT1T2(processingLevel,mountPrefix,ExperimentBaseDirectory, subject_data_file, atlas_fname_wpath, BCD_model_path,
               InterpolationMode="Linear", Mode=10,DwiList=[] ):
    """
    Run autoworkup on all subjects data defined in the subject_data_file

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectory is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """

    subjectDatabaseFile=os.path.join( ExperimentBaseDirectory,'InternalWorkflowSubjectDB.pickle')
    processingLevel=int(processingLevel)
    subjData=csv.reader(open(subject_data_file,'rb'), delimiter=',', quotechar='"')
    myDB=dict()
    multiLevel=AutoVivification()  #This should be replaced by a more nested dictionary
    nestedDictionary=AutoVivification()
    for row in subjData:
        currDict=dict()
        validEntry=True
        if len(row) == 5:
            site=row[0]
            subj=row[1]
            session=row[2]
            T1s=eval(row[3])
            T2s=eval(row[4])
            fullT1s=[ mountPrefix+i for i in T1s]
            fullT2s=[ mountPrefix+i for i in T2s]
            currDict['T1s']=fullT1s
            currDict['T2s']=fullT2s
            currDict['site']=site
            currDict['subj']=subj
            if len(fullT1s) < 1:
                print("Invalid Entry!  {0}".format(currDict))
                validEntry=False
            if len(fullT2s) < 1:
                print("Invalid Entry!  {0}".format(currDict))
                validEntry=False
            for i in fullT1s:
                if not os.path.exists(i):
                    print("Missing File: {0}".format(i))
                    validEntry=False
            for i in fullT2s:
                if not os.path.exists(i):
                    print("Missing File: {0}".format(i))
                    validEntry=False

            if validEntry == True:
                myDB[session]=currDict
                UNIQUE_ID=site+"_"+subj+"_"+session
                nestedDictionary[site][subj][session]=currDict
                multiLevel[UNIQUE_ID]=currDict
    from cPickle import dump
    dump(multiLevel, open(subjectDatabaseFile,'w'))

    ########### PIPELINE INITIALIZATION #############
    baw200 = pe.Workflow(name="BAW_20120104_workflow")
    baw200.config['execution'] = {
                                     'plugin':'Linear',
                                     #'stop_on_first_crash':'true',
                                     'stop_on_first_crash':'false',
                                     'stop_on_first_rerun': 'false',
                                     'hash_method': 'timestamp',
                                     'single_thread_matlab':'true',
                                     'remove_unnecessary_outputs':'false',
                                     'use_relative_paths':'false',
                                     'remove_node_directories':'false',
                                     'local_hash_check':'true'
                                     }
                                     #'job_finished_timeout':3700
    baw200.config['logging'] = {
      'workflow_level':'DEBUG',
      'filemanip_level':'DEBUG',
      'interface_level':'DEBUG',
      'log_directory': ExperimentBaseDirectory
    }
    baw200.base_dir = ExperimentBaseDirectory

    """TODO: Determine if we want to pass subjectID and scanID, always require full
    paths, get them from the output path, or something else.
    """
    uidSource = pe.Node(interface=IdentityInterface(fields=['uid']),name='99_siteSource')
    uidSource.iterables = ('uid', multiLevel.keys() )

    projSource = pe.Node(interface=IdentityInterface(fields=['proj']),name='99_projSource')
    projSource.iterables = ('proj', multiLevel.keys() )

#    subjSource = pe.Node( Function(function=IdentityInterface, input_names = ['T1List','T2List','altT1'], output_names = ['imagePathList']), run_without_submitting=True, name="99_makeImagePathList")
#    baw200.connect( [ (projSource, makeImagePathList, [(('uid', getT1s, subjectDatabaseFile ), 'T1List')] ), ])
#    baw200.connect( [ (uidSource, makeImagePathList, [(('uid', getT2s, subjectDatabaseFile ), 'T2List')] ), ])

#    subjSource = pe.Node(interface=IdentityInterface(fields=['subj']),name='99_subjSource')
#    subjSource.iterables = ('subj', multiLevel.keys() )

#    sessSource = pe.Node(interface=IdentityInterface(fields=['sess']),name='99_sessSource')
#    sessSource.iterables = ('sess', multiLevel.keys() )


    BAtlas = MakeAtlasNode(atlas_fname_wpath) ## Call function to create node

    ########################################################
    # Run ACPC Detect on first T1 Image - Base Image
    ########################################################
    BCD = pe.Node(interface=BRAINSConstellationDetector(), name="01_BCD")
    ##  Use program default BCD.inputs.inputTemplateModel = T1ACPCModelFile
    ##BCD.inputs.outputVolume =   "BCD_OUT" + "_ACPC_InPlace.nii.gz"                #$# T1AcpcImageList
    BCD.inputs.outputTransform =  "BCD" + "_Original2ACPC_transform.mat"
    BCD.inputs.outputResampledVolume = "BCD" + "_ACPC.nii.gz"
    BCD.inputs.outputLandmarksInInputSpace = "BCD" + "_Original.fcsv"
    BCD.inputs.outputLandmarksInACPCAlignedSpace = "BCD" + "_ACPC_Landmarks.fcsv"
    BCD.inputs.outputMRML = "BCD" + "_ACPC_Scene.mrml"
    BCD.inputs.interpolationMode = InterpolationMode
    BCD.inputs.houghEyeDetectorMode = 1  # Look for dark eyes like on a T1 image, 0=Look for bright eyes like in a T2 image
    BCD.inputs.acLowerBound = 80.0 # Chop the data set 80mm below the AC PC point.
    BCD.inputs.llsModel = os.path.join(BCD_model_path,'LLSModel-2ndVersion.hdf5')
    BCD.inputs.inputTemplateModel = os.path.join(BCD_model_path,'T1-2ndVersion.mdl')

    # Entries below are of the form:
    baw200.connect( [ (uidSource, BCD, [(('uid', getFirstT1, subjectDatabaseFile) , 'inputVolume')] ), ])

    if processingLevel > 0:
        ########################################################
        # Run BLI atlas_to_subject
        ########################################################
        BLI = pe.Node(interface=BRAINSLandmarkInitializer(), name="05_BLI")
        BLI.inputs.outputTransformFilename = "landmarkInitializer_atlas_to_subject_transform.mat"

        baw200.connect([
            (BCD,BLI,[('outputLandmarksInACPCAlignedSpace','inputFixedLandmarkFilename')]),
        ])
        baw200.connect([
            (BAtlas,BLI,[('template_landmarks_fcsv','inputMovingLandmarkFilename')]),
            (BAtlas,BLI,[('template_landmark_weights_csv','inputWeightFilename')])
        ])

    if processingLevel > 0:
        ########################################################
        # Run BLI subject_to_atlas
        ########################################################
        BLI2Atlas = pe.Node(interface=BRAINSLandmarkInitializer(), name="05_BLI2Atlas")
        BLI2Atlas.inputs.outputTransformFilename = "landmarkInitializer_subject_to_atlas_transform.mat"

        baw200.connect([
            (BCD,BLI2Atlas,[('outputLandmarksInInputSpace','inputMovingLandmarkFilename')]),
        ])
        baw200.connect([
            (BAtlas,BLI2Atlas,[('template_landmarks_fcsv','inputFixedLandmarkFilename')]),
            (BAtlas,BLI2Atlas,[('template_landmark_weights_csv','inputWeightFilename')])
        ])
        Resample2Atlas=pe.Node(interface=BRAINSResample(),name="05_Resample2Atlas")
        Resample2Atlas.inputs.interpolationMode = "Linear"
        Resample2Atlas.inputs.outputVolume = "subject2atlas.nii.gz"

        baw200.connect( [ (uidSource, Resample2Atlas, [(('uid', getFirstT1, subjectDatabaseFile ), 'inputVolume')] ), ])
        baw200.connect(BLI2Atlas,'outputTransformFilename',Resample2Atlas,'warpTransform')
        baw200.connect(BAtlas,'template_t1',Resample2Atlas,'referenceVolume')

    if processingLevel > 1:
        ########################################################
        # Run BABC on Multi-modal images
        ########################################################
        def MakeOneFileList(T1List,T2List,altT1):
            """ This funciton uses altT1 for the first T1, and the append the rest of the T1's and T2's """
            imagePathList=list()
            imagePathList.append(altT1)
            for i in T1List[1:]:
                imagePathList.append(i)
            for i in T2List[0:]:
                imagePathList.append(i)
            return imagePathList
        makeImagePathList = pe.Node( Function(function=MakeOneFileList, input_names = ['T1List','T2List','altT1'], output_names = ['imagePathList']), run_without_submitting=True, name="99_makeImagePathList")
        baw200.connect( [ (uidSource, makeImagePathList, [(('uid', getT1s, subjectDatabaseFile ), 'T1List')] ), ])
        baw200.connect( [ (uidSource, makeImagePathList, [(('uid', getT2s, subjectDatabaseFile ), 'T2List')] ), ])
        # -- Standard mode to make 256^3 images
        baw200.connect( BCD,    'outputResampledVolume', makeImagePathList, 'altT1' )

        def MakeOneFileTypeList(T1List,T2List):
            input_types =       ["T1"]*len(T1List)
            input_types.extend( ["T2"]*len(T2List) )
            return ",".join(input_types)
        makeImageTypeList = pe.Node( Function(function=MakeOneFileTypeList, input_names = ['T1List','T2List'], output_names = ['imageTypeList']), run_without_submitting=True, name="99_makeImageTypeList")

        baw200.connect( [ (uidSource, makeImageTypeList, [(('uid', getT1s, subjectDatabaseFile ), 'T1List')] ), ])
        baw200.connect( [ (uidSource, makeImageTypeList, [(('uid', getT2s, subjectDatabaseFile ), 'T2List')] ), ])

        def MakeOutFileList(T1List,T2List):
            def GetExtBaseName(filename):
                '''
                Get the filename without the extension.  Works for .ext and .ext.gz
                '''
                import os
                currBaseName = os.path.basename(filename)
                currExt = os.path.splitext(currBaseName)[1]
                currBaseName = os.path.splitext(currBaseName)[0]
                if currExt == ".gz":
                    currBaseName = os.path.splitext(currBaseName)[0]
                    currExt = os.path.splitext(currBaseName)[1]
                return currBaseName
            all_files=T1List
            all_files.extend(T2List)
            out_corrected_names=[]
            for i in all_files:
                out_name=GetExtBaseName(i)+"_corrected.nii.gz"
                out_corrected_names.append(out_name)
            return out_corrected_names
        makeOutImageList = pe.Node( Function(function=MakeOutFileList, input_names = ['T1List','T2List'], output_names = ['outImageList']), run_without_submitting=True, name="99_makeOutImageList")
        baw200.connect( [ (uidSource, makeOutImageList, [(('uid', getT1s, subjectDatabaseFile ), 'T1List')] ), ])
        baw200.connect( [ (uidSource, makeOutImageList, [(('uid', getT2s, subjectDatabaseFile ), 'T2List')] ), ])

        BABC= pe.Node(interface=BRAINSABC(), name="11_BABC")
        baw200.connect(makeImagePathList,'imagePathList',BABC,'inputVolumes')
        baw200.connect(makeImageTypeList,'imageTypeList',BABC,'inputVolumeTypes')
        baw200.connect(makeOutImageList,'outImageList',BABC,'outputVolumes')
        BABC.inputs.debuglevel = 0
        BABC.inputs.maxIterations = 3
        BABC.inputs.maxBiasDegree = 4
        BABC.inputs.filterIteration = 3
        BABC.inputs.filterMethod = 'GradientAnisotropicDiffusion'
        BABC.inputs.gridSize = [28,20,24]
        BABC.inputs.outputFormat = "NIFTI"
        BABC.inputs.outputLabels = "brain_label_seg.nii.gz"
        BABC.inputs.outputDirtyLabels = "volume_label_seg.nii.gz"
        BABC.inputs.posteriorTemplate = "POSTERIOR_%s.nii.gz"
        BABC.inputs.atlasToSubjectTransform = "atlas_to_subject.mat"
        BABC.inputs.implicitOutputs = ['t1_average_BRAINSABC.nii.gz', 't2_average_BRAINSABC.nii.gz']

        BABC.inputs.resamplerInterpolatorType = InterpolationMode
        ##
        BABC.inputs.outputDir = './'

        baw200.connect(BAtlas,'AtlasPVDefinition_xml',BABC,'atlasDefinition')

        baw200.connect(BLI,'outputTransformFilename',BABC,'atlasToSubjectInitialTransform')
        """
        Get the first T1 and T2 corrected images from BABC
        """
        bfc_files = pe.Node(Function(input_names=['in_files','T1_count'],
                                   output_names=['t1_corrected','t2_corrected'],
                                   function=get_first_T1_and_T2), name='99_bfc_files')

        baw200.connect( [ (uidSource, bfc_files, [(('uid', getT1sLength, subjectDatabaseFile ), 'T1_count')] ), ])
        baw200.connect(BABC,'outputVolumes',bfc_files,'in_files')

        """
        ResampleNACLabels
        """
        ResampleAtlasNACLabels=pe.Node(interface=BRAINSResample(),name="13_ResampleAtlasNACLabels")
        ResampleAtlasNACLabels.inputs.interpolationMode = "NearestNeighbor"
        ResampleAtlasNACLabels.inputs.outputVolume = "atlasToSubjectNACLabels.nii.gz"

        baw200.connect(BABC,'atlasToSubjectTransform',ResampleAtlasNACLabels,'warpTransform')
        baw200.connect(bfc_files,'t1_corrected',ResampleAtlasNACLabels,'referenceVolume')
        baw200.connect(BAtlas,'template_nac_lables',ResampleAtlasNACLabels,'inputVolume')

        """
        BRAINSMush
        """
        BMUSH=pe.Node(interface=BRAINSMush(),name="15_BMUSH")
        BMUSH.inputs.outputVolume = "MushImage.nii.gz"
        BMUSH.inputs.outputMask = "MushMask.nii.gz"
        BMUSH.inputs.lowerThresholdFactor = 1.2
        BMUSH.inputs.upperThresholdFactor = 0.55

        baw200.connect(bfc_files,'t1_corrected',BMUSH,'inputFirstVolume')
        baw200.connect(bfc_files,'t2_corrected',BMUSH,'inputSecondVolume')
        baw200.connect(BABC,'outputLabels',BMUSH,'inputMaskVolume')

        """
        BRAINSROIAuto
        """
        BROI = pe.Node(interface=BRAINSROIAuto(), name="17_BRAINSROIAuto")
        BROI.inputs.closingSize=12
        BROI.inputs.otsuPercentileThreshold=0.01
        BROI.inputs.thresholdCorrectionFactor=1.0
        BROI.inputs.outputROIMaskVolume = "temproiAuto_t1_ACPC_corrected_BRAINSABC.nii.gz"
        baw200.connect(bfc_files,'t1_corrected',BROI,'inputVolume')

        """
        Split the implicit outputs of BABC
        """
        SplitAvgBABC = pe.Node(Function(input_names=['in_files','T1_count'], output_names=['avgBABCT1','avgBABCT2'],
                                 function = get_first_T1_and_T2), run_without_submitting=True, name="99_SplitAvgBABC")
        SplitAvgBABC.inputs.T1_count = 1 ## There is only 1 average T1 image.

        baw200.connect(BABC,'implicitOutputs',SplitAvgBABC,'in_files')


        """
        Gradient Anistropic Diffusion images for BRAINSCut
        """
        GADT1=pe.Node(interface=GradientAnisotropicDiffusionImageFilter(),name="27_GADT1")
        GADT1.inputs.timeStep = 0.025
        GADT1.inputs.conductance = 1
        GADT1.inputs.numberOfIterations = 5
        GADT1.inputs.outputVolume = "GADT1.nii.gz"

        baw200.connect(SplitAvgBABC,'avgBABCT1',GADT1,'inputVolume')

        GADT2=pe.Node(interface=GradientAnisotropicDiffusionImageFilter(),name="27_GADT2")
        GADT2.inputs.timeStep = 0.025
        GADT2.inputs.conductance = 1
        GADT2.inputs.numberOfIterations = 5
        GADT2.inputs.outputVolume = "GADT2.nii.gz"

        def printFullPath(outFileFullPath):
            print("="*80)
            print("="*80)
            print("="*80)
            print("="*80)
            print("{0}".format(outFileFullPath))
            return outFileFullPath
        printOutImage = pe.Node( Function(function=printFullPath, input_names = ['outFileFullPath'], output_names = ['genoutFileFullPath']), run_without_submitting=True, name="99_printOutImage")
        baw200.connect( GADT2, 'outputVolume', printOutImage, 'outFileFullPath' )

        baw200.connect(SplitAvgBABC,'avgBABCT2',GADT2,'inputVolume')

        """
        Sum the gradient images for BRAINSCut
        """
        SGI=pe.Node(interface=GenerateSummedGradientImage(),name="27_SGI")
        SGI.inputs.outputFileName = "SummedGradImage.nii.gz"

        baw200.connect(GADT1,'outputVolume',SGI,'inputVolume1')
        baw200.connect(GADT2,'outputVolume',SGI,'inputVolume2')

        if processingLevel > 1:
            """
            Load the BRAINSCut models & probabiity maps.
            """
            BCM_outputs = ['phi','rho','theta',
                           'r_probabilityMaps','l_probabilityMaps',
                           'models']
            BCM_Models = pe.Node(interface=nio.DataGrabber(input_names=['structures'],
                                                           outfields=BCM_outputs),
                                 name='10_BCM_Models')
            BCM_Models.inputs.base_directory = atlas_fname_wpath
            BCM_Models.inputs.template_args['phi'] = [['spatialImages','phi','nii.gz']]
            BCM_Models.inputs.template_args['rho'] = [['spatialImages','rho','nii.gz']]
            BCM_Models.inputs.template_args['theta'] = [['spatialImages','theta','nii.gz']]
            BCM_Models.inputs.template_args['r_probabilityMaps'] = [['structures']]
            BCM_Models.inputs.template_args['l_probabilityMaps'] = [['structures']]
            BCM_Models.inputs.template_args['models'] = [['structures']]

            BRAINSCut_structures = ['caudate','thalamus','putamen','hippocampus']
            #BRAINSCut_structures = ['caudate','thalamus']
            BCM_Models.iterables = ( 'structures',  BRAINSCut_structures )
            BCM_Models.inputs.template = '%s/%s.%s'
            BCM_Models.inputs.field_template = dict(
                r_probabilityMaps='probabilityMaps/r_%s_ProbabilityMap.nii.gz',
                l_probabilityMaps='probabilityMaps/l_%s_ProbabilityMap.nii.gz',
                models='modelFiles/%sModel*',
                )

            """
            The xml creation and BRAINSCut need to be their own mini-pipeline that gets
            executed once for each of the structures in BRAINSCut_structures.  This can be
            accomplished with a map node and a new pipeline.
            """
            """
            Create xml file for BRAINSCut
            """


            BFitAtlasToSubject = pe.Node(interface=BRAINSFit(),name="30_BFitAtlasToSubject")
            BFitAtlasToSubject.inputs.costMetric="MMI"
            BFitAtlasToSubject.inputs.maskProcessingMode="ROI"
            BFitAtlasToSubject.inputs.numberOfSamples=100000
            BFitAtlasToSubject.inputs.numberOfIterations=[1500,1500]
            BFitAtlasToSubject.inputs.numberOfHistogramBins=50
            BFitAtlasToSubject.inputs.maximumStepLength=0.2
            BFitAtlasToSubject.inputs.minimumStepLength=[0.005,0.005]
            BFitAtlasToSubject.inputs.transformType= ["Affine","BSpline"]
            BFitAtlasToSubject.inputs.maxBSplineDisplacement= 7
            BFitAtlasToSubject.inputs.maskInferiorCutOffFromCenter=65
            BFitAtlasToSubject.inputs.splineGridSize=[28,20,24]
            BFitAtlasToSubject.inputs.outputVolume="Trial_Initializer_Output.nii.gz"
            BFitAtlasToSubject.inputs.outputTransform="Trial_Initializer_Output.mat"
            baw200.connect(SplitAvgBABC,'avgBABCT1',BFitAtlasToSubject,'fixedVolume')
            baw200.connect(BABC,'outputLabels',BFitAtlasToSubject,'fixedBinaryVolume')
            baw200.connect(BAtlas,'template_t1',BFitAtlasToSubject,'movingVolume')
            baw200.connect(BAtlas,'template_brain',BFitAtlasToSubject,'movingBinaryVolume')
            baw200.connect(BLI,'outputTransformFilename',BFitAtlasToSubject,'initialTransform')

            CreateBRAINSCutXML = pe.Node(Function(input_names=['rho','phi','theta',
                                                                  'model',
                                                                  'r_probabilityMap',
                                                                  'l_probabilityMap',
                                                                  'atlasT1','atlasBrain',
                                                                  'subjT1','subjT2',
                                                                  'subjT1GAD','subjT2GAD',
                                                                  'subjSGGAD','subjBrain',
                                                                  'atlasToSubj','output_dir'],
                                                     output_names=['xml_filename','rl_structure_filename_list'],
                                                     function = create_BRAINSCut_XML),
                                            overwrite = True,
                                            name="CreateBRAINSCutXML")

            ## HACK  Makde better directory
            CreateBRAINSCutXML.inputs.output_dir = "." #os.path.join(baw200.base_dir, "BRAINSCut_output")
            baw200.connect(BCM_Models,'models',CreateBRAINSCutXML,'model')
            baw200.connect(BCM_Models,'rho',CreateBRAINSCutXML,'rho')
            baw200.connect(BCM_Models,'phi',CreateBRAINSCutXML,'phi')
            baw200.connect(BCM_Models,'theta',CreateBRAINSCutXML,'theta')
            baw200.connect(BCM_Models,'r_probabilityMaps',CreateBRAINSCutXML,'r_probabilityMap')
            baw200.connect(BCM_Models,'l_probabilityMaps',CreateBRAINSCutXML,'l_probabilityMap')
            baw200.connect(BAtlas,'template_t1',CreateBRAINSCutXML,'atlasT1')
            baw200.connect(BAtlas,'template_brain',CreateBRAINSCutXML,'atlasBrain')
            baw200.connect(SplitAvgBABC,'avgBABCT1',CreateBRAINSCutXML,'subjT1')
            baw200.connect(SplitAvgBABC,'avgBABCT2',CreateBRAINSCutXML,'subjT2')
            baw200.connect(GADT1,'outputVolume',CreateBRAINSCutXML,'subjT1GAD')
            baw200.connect(GADT2,'outputVolume',CreateBRAINSCutXML,'subjT2GAD')
            baw200.connect(SGI,'outputFileName',CreateBRAINSCutXML,'subjSGGAD')
            baw200.connect(BABC,'outputLabels',CreateBRAINSCutXML,'subjBrain')
            baw200.connect(BFitAtlasToSubject,'outputTransform',CreateBRAINSCutXML,'atlasToSubj')
            #CreateBRAINSCutXML.inputs.atlasToSubj = "INTERNAL_REGISTER.mat"
            #baw200.connect(BABC,'atlasToSubjectTransform',CreateBRAINSCutXML,'atlasToSubj')

        if 1 == 1:
            """
            BRAINSCut
            """
            BRAINSCUT = pe.Node(interface=BRAINSCut(),name="BRAINSCUT",
                                   input_names=['netConfiguration'],
                                   output_names=['implicitOutputs'])
            BRAINSCUT.inputs.applyModel = True
            baw200.connect(CreateBRAINSCutXML,'xml_filename',BRAINSCUT,'netConfiguration')
            baw200.connect(CreateBRAINSCutXML,'rl_structure_filename_list',BRAINSCUT,'implicitOutputs')

            """
            BRAINSTalairach
            Not implemented yet.
            """

        ## Make deformed Atlas image space
        if processingLevel > 2:
            print("""
            Run ANTS Registration at processingLevel={0}
            """.format(processingLevel) )
            ComputeAtlasToSubjectTransform = pe.Node(interface=ANTSWrapper(), name="19_ComputeAtlasToSubjectTransform")
            ComputeAtlasToSubjectTransform.inputs.output_prefix = "ANTS_"

            baw200.connect( SplitAvgBABC,'avgBABCT1',ComputeAtlasToSubjectTransform,"fixed_T1_image")
            baw200.connect( SplitAvgBABC,'avgBABCT2',ComputeAtlasToSubjectTransform,"fixed_T2_image")
            baw200.connect( BAtlas,'template_t1',    ComputeAtlasToSubjectTransform,"moving_T1_image")
            baw200.connect( BAtlas,'template_t1',    ComputeAtlasToSubjectTransform,"moving_T2_image")

            WarpAtlas = pe.Node(interface=WarpAllAtlas(), name = "19_WarpAtlas")
            WarpAtlas.inputs.moving_atlas = atlas_fname_wpath
            WarpAtlas.inputs.deformed_atlas = "./"
            baw200.connect( ComputeAtlasToSubjectTransform,'output_affine', WarpAtlas,"affine_transform")
            baw200.connect( ComputeAtlasToSubjectTransform,'output_warp', WarpAtlas,"deformation_field")
            baw200.connect( SplitAvgBABC,'avgBABCT1', WarpAtlas, 'reference_image')

        if processingLevel > 3:
            print("""
            Run Freesurfer ReconAll at processingLevel={0}
            """.format(processingLevel) )
            subj_id = os.path.basename(os.path.dirname(os.path.dirname(baw200.base_dir)))
            scan_id = os.path.basename(os.path.dirname(baw200.base_dir))
            reconall = pe.Node(interface=ReconAll(),name="41_FS510")
            reconall.inputs.subject_id = subj_id+'_'+scan_id
            reconall.inputs.directive = 'all'
            reconall.inputs.subjects_dir = '.'
            baw200.connect(SplitAvgBABC,'avgBABCT1',reconall,'T1_files')
        else:
            print "Skipping freesurfer"

    return baw200

