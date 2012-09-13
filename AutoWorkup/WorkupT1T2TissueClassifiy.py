#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from BRAINSTools.BRAINSABCext import *
"""
    from WorkupT1T2TissueClassify import CreateTissueClassifyWorkflow
    myLocalTCWF= CreateTissueClassifyWorkflow("TissueClassify")
    tissueClassifyWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT1s, subjectDatabaseFile ), 'T1List')] ), ])
    tissueClassifyWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT2s, subjectDatabaseFile ), 'T2List')] ), ])
    tissueClassifyWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT1sLength, subjectDatabaseFile ), 'T1_count')] ), ])
    tissueClassifyWF.connect( BCD,    'outputResampledVolume', myLocalTCWF, 'PrimaryT1' )
    tissueClassifyWF.connect(BAtlas,'ExtendedAtlasDefinition.xml',myLocalTCWF,'atlasDefinition')
    tissueClassifyWF.connect(BLI,'outputTransformFilename',myLocalTCWF,'atlasToSubjectInitialTransform')
"""

def MakeOneFileList(T1List,T2List,PDList,FLList,OtherList,PrimaryT1):
    """ This funciton uses PrimaryT1 for the first T1, and the append the rest of the T1's and T2's """
    imagePathList=list()
    imagePathList.append(PrimaryT1) # Force replacement of the first element
    for i in T1List[1:]:
        imagePathList.append(i) # The reset of the elements
    for i in T2List[0:]:
        imagePathList.append(i)
    for i in PDList[0:]:
        imagePathList.append(i)
    for i in FLList[0:]:
        imagePathList.append(i)
    for i in OtherList[0:]:
        imagePathList.append(i)
    return imagePathList
def MakeOneFileTypeList(T1List,T2List,PDList,FLList,OtherList):
    input_types =       ["T1"]*len(T1List)
    input_types.extend( ["T2"]*len(T2List) )
    input_types.extend( ["PD"]*len(PDList) )
    input_types.extend( ["FL"]*len(FLList) )
    input_types.extend( ["OTHER"]*len(OtherList) )
    return input_types

def MakeOutFileList(T1List,T2List,PDList,FLList,OtherList):
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
    all_files.extend(PDList)
    all_files.extend(FLList)
    all_files.extend(OtherList)
    out_corrected_names=[]
    for i in all_files:
        out_name=GetExtBaseName(i)+"_corrected.nii.gz"
        out_corrected_names.append(out_name)
    return out_corrected_names

def getListIndexOrNoneIfOutOfRange( imageList, index):
    if index < len(imageList):
      return imageList[index]
    else:
      return None
def MakePosteriorDictionaryFunc(posteriorImages):
    posteriorNames=['WM', 'SURFGM', 'ACCUMBEN', 'CAUDATE', 'PUTAMEN', 'GLOBUS', 'THALAMUS', 'HIPPOCAMPUS', 'CRBLGM', 'CRBLWM', 'CSF', 'VB', 'NOTCSF', 'NOTGM', 'NOTWM', 'NOTVB', 'AIR']
    if len(posteriorNames) != len(posteriorImages):
        print "ERROR: ", posteriorNames
        print "ERROR: ", posteriorImages
        return -1
    temp_dictionary=dict(zip(posteriorNames,posteriorImages))
    return temp_dictionary

def CreateTissueClassifyWorkflow(WFname,CLUSTER_QUEUE,InterpolationMode):
    tissueClassifyWF= pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T1List','T2List','PDList','FLList','OtherList',
                                                             'T1_count',
                                                             'PrimaryT1',
        'atlasDefinition','atlasToSubjectInitialTransform']),
        run_without_submitting=True,
        name='InputSpec' )
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['atlasToSubjectTransform','outputLabels','outputHeadLabels',
            #'t1_corrected','t2_corrected',
            't1_average','t2_average','pd_average','fl_average',
            'TissueClassifyOutputDir',
            'posteriorImages'
            ]),
        run_without_submitting=True,
        name='OutputSpec' )


    ########################################################
    # Run BABCext on Multi-modal images
    ########################################################
    makeImagePathList = pe.Node( Function(function=MakeOneFileList,
                                          input_names = ['T1List','T2List','PDList','FLList','OtherList','PrimaryT1'],
                                          output_names = ['imagePathList']), run_without_submitting=True, name="99_makeImagePathList")
    tissueClassifyWF.connect( inputsSpec, 'T1List', makeImagePathList, 'T1List' )
    tissueClassifyWF.connect( inputsSpec, 'T2List', makeImagePathList, 'T2List' )
    tissueClassifyWF.connect( inputsSpec, 'PDList', makeImagePathList, 'PDList' )
    tissueClassifyWF.connect( inputsSpec, 'FLList', makeImagePathList, 'FLList' )
    tissueClassifyWF.connect( inputsSpec, 'OtherList', makeImagePathList, 'OtherList' )
    # -- Standard mode to make 256^3 images
    tissueClassifyWF.connect( inputsSpec, 'PrimaryT1', makeImagePathList, 'PrimaryT1' )

    makeImageTypeList = pe.Node( Function(function=MakeOneFileTypeList,
                                          input_names = ['T1List','T2List','PDList','FLList','OtherList'],
                                          output_names = ['imageTypeList']), run_without_submitting=True, name="99_makeImageTypeList")
    tissueClassifyWF.connect( inputsSpec, 'T1List', makeImageTypeList, 'T1List' )
    tissueClassifyWF.connect( inputsSpec, 'T2List', makeImageTypeList, 'T2List' )
    tissueClassifyWF.connect( inputsSpec, 'PDList', makeImageTypeList, 'PDList' )
    tissueClassifyWF.connect( inputsSpec, 'FLList', makeImageTypeList, 'FLList' )
    tissueClassifyWF.connect( inputsSpec, 'OtherList', makeImageTypeList, 'OtherList' )

    makeOutImageList = pe.Node( Function(function=MakeOutFileList,
                                         input_names = ['T1List','T2List','PDList','FLList','OtherList'],
                                         output_names = ['outImageList']), run_without_submitting=True, name="99_makeOutImageList")
    tissueClassifyWF.connect( inputsSpec, 'T1List', makeOutImageList, 'T1List' )
    tissueClassifyWF.connect( inputsSpec, 'T2List', makeOutImageList, 'T2List' )
    tissueClassifyWF.connect( inputsSpec, 'PDList', makeOutImageList, 'PDList' )
    makeOutImageList.inputs.FLList=[] ## an emptyList HACK
    #HACK tissueClassifyWF.connect( inputsSpec, 'FLList', makeOutImageList, 'FLList' )
    tissueClassifyWF.connect( inputsSpec, 'OtherList', makeOutImageList, 'OtherList' )

    BABCext= pe.Node(interface=BRAINSABCext(), name="BABC")
    #many_cpu_BABC_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 4-12                   -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
    many_cpu_BABC_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 4-12 -l h_vmem=8G,mem_free=8G -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
    BABCext.plugin_args=many_cpu_BABC_options_dictionary
    tissueClassifyWF.connect(makeImagePathList,'imagePathList',BABCext,'inputVolumes')
    tissueClassifyWF.connect(makeImageTypeList,'imageTypeList',BABCext,'inputVolumeTypes')
    tissueClassifyWF.connect(makeOutImageList,'outImageList',BABCext,'outputVolumes')
    BABCext.inputs.debuglevel = 0
    BABCext.inputs.maxIterations = 3
    BABCext.inputs.maxBiasDegree = 4
    BABCext.inputs.filterIteration = 3
    BABCext.inputs.filterMethod = 'GradientAnisotropicDiffusion'
    BABCext.inputs.gridSize = [28,20,24]
    BABCext.inputs.outputFormat = "NIFTI"
    BABCext.inputs.outputLabels = "brain_label_seg.nii.gz"
    BABCext.inputs.outputDirtyLabels = "volume_label_seg.nii.gz"
    BABCext.inputs.posteriorTemplate = "POSTERIOR_%s.nii.gz"
    BABCext.inputs.atlasToSubjectTransform = "atlas_to_subject.mat"
    #BABCext.inputs.implicitOutputs = ['t1_average_BRAINSABC.nii.gz', 't2_average_BRAINSABC.nii.gz']
    BABCext.inputs.interpolationMode = InterpolationMode
    BABCext.inputs.outputDir = './'

    tissueClassifyWF.connect(inputsSpec,'atlasDefinition',BABCext,'atlasDefinition')
    tissueClassifyWF.connect(inputsSpec,'atlasToSubjectInitialTransform',BABCext,'atlasToSubjectInitialTransform')
    """
    Get the first T1 and T2 corrected images from BABCext
    """

    """ HACK:  THIS IS NOT NEEDED!  We should use the averged t1 and averaged t2 images instead!
    def get_first_T1_and_T2(in_files,T1_count):
        '''
        Returns the first T1 and T2 file in in_files, based on offset in T1_count.
        '''
        return in_files[0],in_files[T1_count]
    bfc_files = pe.Node(Function(input_names=['in_files','T1_count'],
                               output_names=['t1_corrected','t2_corrected'],
                               function=get_first_T1_and_T2), run_without_submitting=True, name='99_bfc_files' )
    tissueClassifyWF.connect( inputsSpec, 'T1_count', bfc_files, 'T1_count')
    tissueClassifyWF.connect(BABCext,'outputVolumes',bfc_files, 'in_files')


    tissueClassifyWF.connect(bfc_files,'t1_corrected',outputsSpec,'t1_corrected')
    tissueClassifyWF.connect(bfc_files,'t2_corrected',outputsSpec,'t2_corrected')
    #tissueClassifyWF.connect(bfc_files,'pd_corrected',outputsSpec,'pd_corrected')
    #tissueClassifyWF.connect(bfc_files,'fl_corrected',outputsSpec,'fl_corrected')

    """

    #############
    tissueClassifyWF.connect(BABCext,'atlasToSubjectTransform',outputsSpec,'atlasToSubjectTransform')
    tissueClassifyWF.connect(BABCext,'outputLabels',outputsSpec,'outputLabels')
    tissueClassifyWF.connect(BABCext,'outputDirtyLabels',outputsSpec,'outputHeadLabels')

    tissueClassifyWF.connect( BABCext , 'outputT1AverageImage', outputsSpec, 't1_average')
    tissueClassifyWF.connect( BABCext , 'outputT2AverageImage', outputsSpec, 't2_average')
    tissueClassifyWF.connect( BABCext , 'outputPDAverageImage', outputsSpec, 'pd_average')
    tissueClassifyWF.connect( BABCext , 'outputFLAverageImage', outputsSpec, 'fl_average')
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 0 ), "t1_average")] ), ] )
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 1 ), "t2_average")] ), ] )
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 2 ), "pd_average")] ), ] )

    tissueClassifyWF.connect(BABCext,'outputDir',outputsSpec,'TissueClassifyOutputDir')

    MakePosteriorDictionaryNode = pe.Node( Function(function=MakePosteriorDictionaryFunc,
                                      input_names = ['posteriorImages'],
                                      output_names = ['posteriorDictionary']), run_without_submitting=True, name="99_makePosteriorDictionary")
    tissueClassifyWF.connect(BABCext,'posteriorImages',MakePosteriorDictionaryNode,'posteriorImages')

    tissueClassifyWF.connect(MakePosteriorDictionaryNode,'posteriorDictionary',outputsSpec,'posteriorImages')

    return tissueClassifyWF
