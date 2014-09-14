#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from BRAINSABCext import *
from utilities.misc import *
from utilities.distributed import modify_qsub_args

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



def getListIndexOrNoneIfOutOfRange(imageList, index):
    if index < len(imageList):
        return imageList[index]
    else:
        return None


def MakePosteriorDictionaryFunc(posteriorImages):
    from PipeLineFunctionHelpers import POSTERIORS
    if len(POSTERIORS) != len(posteriorImages):
        print "ERROR: ", posteriorNames
        print "ERROR: ", POSTERIORS
        return -1
    temp_dictionary = dict(zip(POSTERIORS, posteriorImages))
    return temp_dictionary


def CreateTissueClassifyWorkflow(WFname, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG, InterpolationMode):
    tissueClassifyWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T1List', 'T2List', 'PDList', 'FLList',
                                                             'OtherList', 'T1_count', 'PrimaryT1',
                                                             'atlasDefinition',
                                                             'atlasToSubjectInitialTransform']),
                         run_without_submitting=True,
                         name='inputspec')
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['atlasToSubjectTransform',
                                                              'atlasToSubjectInverseTransform',
                                                              'outputLabels',
                                                              'outputHeadLabels',  # ???
                                                              #'t1_corrected', 't2_corrected',
                                                              't1_average',
                                                              't2_average',
                                                              'pd_average',
                                                              'fl_average',
                                                              'posteriorImages']),
                          run_without_submitting=True,
                          name='outputspec')

    ########################################################
    # Run BABCext on Multi-modal images
    ########################################################
    makeOutImageList = pe.Node(Function(function=MakeOutFileList,
                                        input_names=['T1List', 'T2List', 'PDList', 'FLList',
                                                     'OtherList','postfix','PrimaryT1'],
                                        output_names=['inImageList','outImageList','imageTypeList']),
                                        run_without_submitting=True, name="99_makeOutImageList")
    tissueClassifyWF.connect(inputsSpec, 'T1List', makeOutImageList, 'T1List')
    tissueClassifyWF.connect(inputsSpec, 'T2List', makeOutImageList, 'T2List')
    tissueClassifyWF.connect(inputsSpec, 'PDList', makeOutImageList, 'PDList')
    tissueClassifyWF.connect(inputsSpec, 'PrimaryT1', makeOutImageList, 'PrimaryT1')
    makeOutImageList.inputs.FLList = []  # an emptyList HACK
    makeOutImageList.inputs.postfix = "_corrected.nii.gz"
    # HACK tissueClassifyWF.connect( inputsSpec, 'FLList', makeOutImageList, 'FLList' )
    tissueClassifyWF.connect(inputsSpec, 'OtherList', makeOutImageList, 'OtherList')


    BABCext = pe.Node(interface=BRAINSABCext(), name="BABC")
    many_cpu_BABC_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,8,8,24), 'overwrite': True}
    BABCext.plugin_args = many_cpu_BABC_options_dictionary
    tissueClassifyWF.connect(makeOutImageList, 'inImageList', BABCext, 'inputVolumes')
    tissueClassifyWF.connect(makeOutImageList, 'imageTypeList', BABCext, 'inputVolumeTypes')
    tissueClassifyWF.connect(makeOutImageList, 'outImageList', BABCext, 'outputVolumes')
    BABCext.inputs.debuglevel = 0
    BABCext.inputs.maxIterations = 3
    BABCext.inputs.maxBiasDegree = 4
    BABCext.inputs.filterIteration = 3
    BABCext.inputs.filterMethod = 'GradientAnisotropicDiffusion'
    BABCext.inputs.atlasToSubjectTransformType = 'SyN'
    # BABCext.inputs.atlasToSubjectTransformType = 'BSpline'
    # BABCext.inputs.gridSize = [28,20,24]
    BABCext.inputs.gridSize = [10, 10, 10]
    BABCext.inputs.outputFormat = "NIFTI"
    BABCext.inputs.outputLabels = "brain_label_seg.nii.gz"
    BABCext.inputs.outputDirtyLabels = "volume_label_seg.nii.gz"
    BABCext.inputs.posteriorTemplate = "POSTERIOR_%s.nii.gz"
    BABCext.inputs.atlasToSubjectTransform = "atlas_to_subject.h5"
    # BABCext.inputs.implicitOutputs = ['t1_average_BRAINSABC.nii.gz', 't2_average_BRAINSABC.nii.gz']
    BABCext.inputs.interpolationMode = InterpolationMode
    BABCext.inputs.outputDir = './'

    tissueClassifyWF.connect(inputsSpec, 'atlasDefinition', BABCext, 'atlasDefinition')
    tissueClassifyWF.connect(inputsSpec, 'atlasToSubjectInitialTransform', BABCext, 'atlasToSubjectInitialTransform')
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
    tissueClassifyWF.connect(BABCext, 'atlasToSubjectTransform', outputsSpec, 'atlasToSubjectTransform')

    def MakeInverseTransformFileName(TransformFileName):
        """### HACK:  This function is to work around a deficiency in BRAINSABCext where the inverse transform name is not being computed properly
          in the list outputs"""
        fixed_inverse_name = TransformFileName.replace(".h5", "_Inverse.h5")
        return [fixed_inverse_name]

    tissueClassifyWF.connect([(BABCext, outputsSpec, [(('atlasToSubjectTransform', MakeInverseTransformFileName), "atlasToSubjectInverseTransform")]), ])
    tissueClassifyWF.connect(BABCext, 'outputLabels', outputsSpec, 'outputLabels')
    tissueClassifyWF.connect(BABCext, 'outputDirtyLabels', outputsSpec, 'outputHeadLabels')

    tissueClassifyWF.connect(BABCext, 'outputT1AverageImage', outputsSpec, 't1_average')
    tissueClassifyWF.connect(BABCext, 'outputT2AverageImage', outputsSpec, 't2_average')
    tissueClassifyWF.connect(BABCext, 'outputPDAverageImage', outputsSpec, 'pd_average')
    tissueClassifyWF.connect(BABCext, 'outputFLAverageImage', outputsSpec, 'fl_average')
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 0 ), "t1_average")] ), ] )
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 1 ), "t2_average")] ), ] )
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 2 ), "pd_average")] ), ] )

    MakePosteriorDictionaryNode = pe.Node(Function(function=MakePosteriorDictionaryFunc,
                                                   input_names=['posteriorImages'],
                                                   output_names=['posteriorDictionary']), run_without_submitting=True, name="99_makePosteriorDictionary")
    tissueClassifyWF.connect(BABCext, 'posteriorImages', MakePosteriorDictionaryNode, 'posteriorImages')

    tissueClassifyWF.connect(MakePosteriorDictionaryNode, 'posteriorDictionary', outputsSpec, 'posteriorImages')

    return tissueClassifyWF
