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
    tissueClassifyWF.connect(BAtlas,'AtlasPVDefinition_xml',myLocalTCWF,'atlasDefinition')
    tissueClassifyWF.connect(BLI,'outputTransformFilename',myLocalTCWF,'atlasToSubjectInitialTransform')
"""

def get_first_T1_and_T2(in_files,T1_count):
    '''
    Returns the first T1 and T2 file in in_files, based on offset in T1_count.
    '''
    return in_files[0],in_files[T1_count]


def CreateTissueClassifyWorkflow(WFname,CLUSTER_QUEUE,InterpolationMode):
    tissueClassifyWF= pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T1List','T1_count','T2List','PrimaryT1',
        'atlasDefinition','atlasToSubjectInitialTransform']), name='InputSpec' )



    ########################################################
    # Run BABCext on Multi-modal images
    ########################################################
    def MakeOneFileList(T1List,T2List,PrimaryT1):
        """ This funciton uses PrimaryT1 for the first T1, and the append the rest of the T1's and T2's """
        imagePathList=list()
        imagePathList.append(PrimaryT1)
        for i in T1List[1:]:
            imagePathList.append(i)
        for i in T2List[0:]:
            imagePathList.append(i)
        return imagePathList
    makeImagePathList = pe.Node( Function(function=MakeOneFileList, input_names = ['T1List','T2List','PrimaryT1'], output_names = ['imagePathList']), run_without_submitting=True, name="99_makeImagePathList")
    tissueClassifyWF.connect( inputsSpec, 'T1List', makeImagePathList, 'T1List' )
    tissueClassifyWF.connect( inputsSpec, 'T2List', makeImagePathList, 'T2List' )
    # -- Standard mode to make 256^3 images
    tissueClassifyWF.connect( inputsSpec, 'PrimaryT1', makeImagePathList, 'PrimaryT1' )

    def MakeOneFileTypeList(T1List,T2List):
        input_types =       ["T1"]*len(T1List)
        input_types.extend( ["T2"]*len(T2List) )
        return input_types
    makeImageTypeList = pe.Node( Function(function=MakeOneFileTypeList, input_names = ['T1List','T2List'], output_names = ['imageTypeList']), run_without_submitting=True, name="99_makeImageTypeList")

    tissueClassifyWF.connect( inputsSpec, 'T1List', makeImageTypeList, 'T1List' )
    tissueClassifyWF.connect( inputsSpec, 'T2List', makeImageTypeList, 'T2List' )

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
    tissueClassifyWF.connect( inputsSpec, 'T1List', makeOutImageList, 'T1List' )
    tissueClassifyWF.connect( inputsSpec, 'T2List', makeOutImageList, 'T2List' )

    BABCext= pe.Node(interface=BRAINSABCext(), name="BABC")
    #many_cpu_BABC_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 4-12 -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
    many_cpu_BABC_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 4-12 -l mem_free=8000M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
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
    bfc_files = pe.Node(Function(input_names=['in_files','T1_count'],
                               output_names=['t1_corrected','t2_corrected'],
                               function=get_first_T1_and_T2), name='99_bfc_files')
    tissueClassifyWF.connect( inputsSpec, 'T1_count', bfc_files, 'T1_count')

    #############
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['atlasToSubjectTransform','outputLabels','outputHeadLabels',
            't1_corrected','t2_corrected','outputAverageImages','TissueClassifyOutputDir']), name='OutputSpec' )

    tissueClassifyWF.connect(bfc_files,'t1_corrected',outputsSpec,'t1_corrected')
    tissueClassifyWF.connect(bfc_files,'t2_corrected',outputsSpec,'t2_corrected')

    tissueClassifyWF.connect(BABCext,'outputVolumes',bfc_files, 'in_files')

    tissueClassifyWF.connect(BABCext,'atlasToSubjectTransform',outputsSpec,'atlasToSubjectTransform')
    tissueClassifyWF.connect(BABCext,'outputLabels',outputsSpec,'outputLabels')
    tissueClassifyWF.connect(BABCext,'outputDirtyLabels',outputsSpec,'outputHeadLabels')
    tissueClassifyWF.connect(BABCext,'outputAverageImages',outputsSpec,'outputAverageImages')

    tissueClassifyWF.connect(BABCext,'outputDir',outputsSpec,'TissueClassifyOutputDir')

    return tissueClassifyWF
