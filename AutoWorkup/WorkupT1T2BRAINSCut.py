#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from BRAINSTools import *
from BRAINSTools.RF12BRAINSCutWrapper import RF12BRAINSCutWrapper

"""
    from WorkupT1T2BRAINSCutSegmentation import CreateBRAINSCutSegmentationWorkflow
    myLocalcutWF= CreateBRAINSCutSegmentationWorkflow("999999_PersistanceCheckingWorkflow")
    cutWF.connect(SplitAvgBABC,'avgBABCT1',myLocalcutWF,'fixedVolume')
    cutWF.connect(BABC,'outputLabels',myLocalcutWF,'fixedBinaryVolume')
    cutWF.connect(BAtlas,'template_t1',myLocalcutWF,'movingVolume')
    cutWF.connect(BAtlas,'template_brain',myLocalcutWF,'movingBinaryVolume')
    cutWF.connect(BLI,'outputTransformFilename',myLocalcutWF,'initialTransform')
"""
def CreateBRAINSCutWorkflow(WFname,CLUSTER_QUEUE,atlasObject):
    cutWF= pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T1Volume','T2Volume',
        'atlasToSubjectTransform']), name='InputSpec' )

    """
    Gradient Anistropic Diffusion images for BRAINSCut
    """
    GADT1=pe.Node(interface=GradientAnisotropicDiffusionImageFilter(),name="GADT1")
    GADT1.inputs.timeStep = 0.025
    GADT1.inputs.conductance = 1
    GADT1.inputs.numberOfIterations = 5
    GADT1.inputs.outputVolume = "GADT1.nii.gz"

    cutWF.connect(inputsSpec,'T1Volume',GADT1,'inputVolume')

    GADT2=pe.Node(interface=GradientAnisotropicDiffusionImageFilter(),name="GADT2")
    GADT2.inputs.timeStep = 0.025
    GADT2.inputs.conductance = 1
    GADT2.inputs.numberOfIterations = 5
    GADT2.inputs.outputVolume = "GADT2.nii.gz"
    cutWF.connect(inputsSpec,'T2Volume',GADT2,'inputVolume')

    """
    Sum the gradient images for BRAINSCut
    """
    SGI=pe.Node(interface=GenerateSummedGradientImage(),name="SGI")
    SGI.inputs.outputFileName = "SummedGradImage.nii.gz"

    cutWF.connect(GADT1,'outputVolume',SGI,'inputVolume1')
    cutWF.connect(GADT2,'outputVolume',SGI,'inputVolume2')

    """
    BRAINSCut
    """
    RF12BC = pe.Node(interface=RF12BRAINSCutWrapper(),name="RF12_BRAINSCut")
    RF12BC.inputs.trainingVectorFilename = "trainingVectorFilename.txt"
    RF12BC.inputs.xmlFilename = "BRAINSCutSegmentationDefinition.xml"

    RF12BC.inputs.outputBinaryLeftAccumben=     'l_Accumben_seg.nii.gz'
    RF12BC.inputs.outputBinaryRightAccumben=    'r_Accumben_seg.nii.gz'
    RF12BC.inputs.outputBinaryLeftCaudate=      'l_Caudate_seg.nii.gz'
    RF12BC.inputs.outputBinaryRightCaudate=     'r_Caudate_seg.nii.gz'
    RF12BC.inputs.outputBinaryLeftGlobus=       'l_Globus_seg.nii.gz'
    RF12BC.inputs.outputBinaryRightGlobus=      'r_Globus_seg.nii.gz'
    RF12BC.inputs.outputBinaryLeftHippocampus=  'l_Hippocampus_seg.nii.gz'
    RF12BC.inputs.outputBinaryRightHippocampus= 'r_Hippocampus_seg.nii.gz'
    RF12BC.inputs.outputBinaryLeftPutamen=      'l_Putamen_seg.nii.gz'
    RF12BC.inputs.outputBinaryRightPutamen=     'r_Putamen_seg.nii.gz'
    RF12BC.inputs.outputBinaryLeftThalamus=     'l_Thalamus_seg.nii.gz'
    RF12BC.inputs.outputBinaryRightThalamus=    'r_Thalamus_seg.nii.gz'

    cutWF.connect(inputsSpec,'T1Volume',RF12BC,'inputSubjectT1Filename')
    cutWF.connect(inputsSpec,'T2Volume',RF12BC,'inputSubjectT2Filename')
    #cutWF.connect(SGI,'outputVolume',RF12BC,'inputSubjectSGFilename')
    cutWF.connect(SGI,'outputFileName',RF12BC,'inputSubjectSGFilename')
    cutWF.connect(atlasObject,'template_t1',RF12BC,'inputTemplateT1')
    cutWF.connect(atlasObject,'rho',RF12BC,'inputTemplateRhoFilename')
    cutWF.connect(atlasObject,'phi',RF12BC,'inputTemplatePhiFilename')
    cutWF.connect(atlasObject,'theta',RF12BC,'inputTemplateThetaFilename')


    cutWF.connect(atlasObject,'l_accumben_ProbabilityMap',RF12BC,'probabilityMapsLeftAccumben')
    cutWF.connect(atlasObject,'r_accumben_ProbabilityMap',RF12BC,'probabilityMapsRightAccumben')
    cutWF.connect(atlasObject,'l_caudate_ProbabilityMap',RF12BC,'probabilityMapsLeftCaudate')
    cutWF.connect(atlasObject,'r_caudate_ProbabilityMap',RF12BC,'probabilityMapsRightCaudate')
    cutWF.connect(atlasObject,'l_globus_ProbabilityMap',RF12BC,'probabilityMapsLeftGlobus')
    cutWF.connect(atlasObject,'r_globus_ProbabilityMap',RF12BC,'probabilityMapsRightGlobus')
    cutWF.connect(atlasObject,'l_hippocampus_ProbabilityMap',RF12BC,'probabilityMapsLeftHippocampus')
    cutWF.connect(atlasObject,'r_hippocampus_ProbabilityMap',RF12BC,'probabilityMapsRightHippocampus')
    cutWF.connect(atlasObject,'l_putamen_ProbabilityMap',RF12BC,'probabilityMapsLeftPutamen')
    cutWF.connect(atlasObject,'r_putamen_ProbabilityMap',RF12BC,'probabilityMapsRightPutamen')
    cutWF.connect(atlasObject,'l_thalamus_ProbabilityMap',RF12BC,'probabilityMapsLeftThalamus')
    cutWF.connect(atlasObject,'r_thalamus_ProbabilityMap',RF12BC,'probabilityMapsRightThalamus')
    cutWF.connect(atlasObject,'RandomForestAllSubcorticalsBalancedModel_txtD0060NT0060',RF12BC,'modelFilename')

    cutWF.connect(inputsSpec,'atlasToSubjectTransform',RF12BC,'deformationFromTemplateToSubject')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=[
        'outputBinaryLeftAccumben','outputBinaryRightAccumben',
        'outputBinaryLeftCaudate','outputBinaryRightCaudate',
        'outputBinaryLeftGlobus','outputBinaryRightGlobus',
        'outputBinaryLeftHippocampus','outputBinaryRightHippocampus',
        'outputBinaryLeftPutamen','outputBinaryRightPutamen',
        'outputBinaryLeftThalamus','outputBinaryRightThalamus',
        'xmlFilename']), name='OutputSpec' )

    cutWF.connect(RF12BC,'outputBinaryLeftAccumben',outputsSpec,'outputBinaryLeftAccumben')
    cutWF.connect(RF12BC,'outputBinaryRightAccumben',outputsSpec,'outputBinaryRightAccumben')
    cutWF.connect(RF12BC,'outputBinaryLeftCaudate',outputsSpec,'outputBinaryLeftCaudate')
    cutWF.connect(RF12BC,'outputBinaryRightCaudate',outputsSpec,'outputBinaryRightCaudate')
    cutWF.connect(RF12BC,'outputBinaryLeftGlobus',outputsSpec,'outputBinaryLeftGlobus')
    cutWF.connect(RF12BC,'outputBinaryRightGlobus',outputsSpec,'outputBinaryRightGlobus')
    cutWF.connect(RF12BC,'outputBinaryLeftHippocampus',outputsSpec,'outputBinaryLeftHippocampus')
    cutWF.connect(RF12BC,'outputBinaryRightHippocampus',outputsSpec,'outputBinaryRightHippocampus')
    cutWF.connect(RF12BC,'outputBinaryLeftPutamen',outputsSpec,'outputBinaryLeftPutamen')
    cutWF.connect(RF12BC,'outputBinaryRightPutamen',outputsSpec,'outputBinaryRightPutamen')
    cutWF.connect(RF12BC,'outputBinaryLeftThalamus',outputsSpec,'outputBinaryLeftThalamus')
    cutWF.connect(RF12BC,'outputBinaryRightThalamus',outputsSpec,'outputBinaryRightThalamus')
    cutWF.connect(RF12BC,'xmlFilename',outputsSpec,'xmlFilename')

    return cutWF
