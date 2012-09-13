#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from BRAINSTools import *
from BRAINSTools.RF8BRAINSCutWrapper import RF8BRAINSCutWrapper

def GenerateWFName(projectid, subjectid, sessionid,WFName):
    return WFName+'_'+str(subjectid)+"_"+str(sessionid)+"_"+str(projectid)

def CreateLabelMap(listOfImages,LabelImageName,CSVFileName,projectid, subjectid, sessionid):
    """
    A function to create a consolidated label map and a
    csv file of volume measurements.
    """

    """
    subjectANNLabel_l_caudate.nii.gz
    subjectANNLabel_l_hippocampus.nii.gz
    subjectANNLabel_l_putamen.nii.gz
    subjectANNLabel_l_thalamus.nii.gz
    subjectANNLabel_r_caudate.nii.gz
    subjectANNLabel_r_hippocampus.nii.gz
    subjectANNLabel_r_putamen.nii.gz
    subjectANNLabel_r_thalamus.nii.gz
    """
    import SimpleITK as sitk
    import os
    import csv
    orderOfPriority = [
      "l_caudate"     ,
      "r_caudate"     ,
      "l_putamen"     ,
      "r_putamen"     ,
      "l_hippocampus" ,
      "r_hippocampus" ,
      "l_thalamus"    ,
      "r_thalamus"
    ]

    valueDict={
        "l_caudate"     : 1,
        "r_caudate"     : 2,
        "l_putamen"     : 3,
        "r_putamen"     : 4,
        "l_hippocampus" : 5,
        "r_hippocampus" : 6,
        "l_thalamus"    : 7,
        "r_thalamus"    : 8
    }

    labelImage = None
    for segFN in listOfImages:
        im = sitk.ReadImage(segFN)
        im.GetSize()
        remove_pre_postfix=os.path.basename(segFN.replace(".nii.gz","").replace("subjectANNLabel_","").replace("_seg",""))
        structName=remove_pre_postfix.lower()
        if labelImage is None:
            labelImage = im*valueDict[structName]
        else:
            mask=sitk.Not(im)
            ## Clear out an empty space for the next mask to be inserted
            labelImage *= mask
            ## Add in the mask image with it's proper label
            labelImage = labelImage + im*valueDict[structName]
    sitk.WriteImage(labelImage,LabelImageName)

    ls = sitk.LabelStatisticsImageFilter()
    ls.Execute(labelImage,labelImage)
    ImageSpacing=labelImage.GetSpacing()
    csvFile=open(CSVFileName,'w')
    dWriter=csv.DictWriter(csvFile,['projectid', 'subjectid', 'sessionid','Structure','LabelCode','Volume_mm3'],restval='', extrasaction='raise', dialect='excel')
    dWriter.writeheader()
    writeDictionary=dict()
    for name in orderOfPriority:
        value = valueDict[name]
        if ls.HasLabel(value):
            #print "Displaying: ", name, value
            myMeasurementMap = ls.GetMeasurementMap(value)
            dictKeys=myMeasurementMap.GetVectorOfMeasurementNames()
            dictValues=myMeasurementMap.GetVectorOfMeasurementValues()
            measurementDict=dict(zip(dictKeys, dictValues))
            structVolume=ImageSpacing[0]*ImageSpacing[1]*ImageSpacing[2]*measurementDict['Count']
            writeDictionary['Volume_mm3']=structVolume
            writeDictionary['Structure']=name
            writeDictionary['LabelCode']=value
            #writeDictionary['FileName']=os.path.abspath(LabelImageName)
            writeDictionary['projectid']=projectid
            writeDictionary['subjectid']=subjectid
            writeDictionary['sessionid']=sessionid
            dWriter.writerow(writeDictionary)
    return os.path.abspath(LabelImageName),os.path.abspath(CSVFileName)

#==============================================
#==============================================
#==============================================
#==============================================
#==============================================
#==============================================

"""
    from WorkupT1T2BRAINSCutSegmentation import CreateBRAINSCutSegmentationWorkflow
    myLocalcutWF= CreateBRAINSCutSegmentationWorkflow("999999_PersistanceCheckingWorkflow")
    cutWF.connect(SplitAvgBABC,'avgBABCT1',myLocalcutWF,'fixedVolume')
    cutWF.connect(BABC,'outputLabels',myLocalcutWF,'fixedBinaryVolume')
    cutWF.connect(BAtlas,'template_t1',myLocalcutWF,'movingVolume')
    cutWF.connect(BAtlas,'template_brain',myLocalcutWF,'movingBinaryVolume')
    cutWF.connect(BLI,'outputTransformFilename',myLocalcutWF,'initialTransform')
"""
def CreateBRAINSCutWorkflow(projectid, subjectid, sessionid,WFName,CLUSTER_QUEUE,atlasObject):
    cutWF= pe.Workflow(name=GenerateWFName(projectid, subjectid, sessionid,WFName))

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T1Volume','T2Volume',
        'TotalGM','RegistrationROI',
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
    RF8BC = pe.Node(interface=RF8BRAINSCutWrapper(),name="RF8_BRAINSCut")
    many_cpu_RF8BC_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 4-12 -l mem_free=8000M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
    RF8BC.plugin_args=many_cpu_RF8BC_options_dictionary
    RF8BC.inputs.trainingVectorFilename = "trainingVectorFilename.txt"
    RF8BC.inputs.xmlFilename = "BRAINSCutSegmentationDefinition.xml"

    """ HACK These should be l_Caudate_seg.nii.gz
    subjectANNLabel_l_caudate.nii.gz
    subjectANNLabel_l_hippocampus.nii.gz
    subjectANNLabel_l_putamen.nii.gz
    subjectANNLabel_l_thalamus.nii.gz
    subjectANNLabel_r_caudate.nii.gz
    subjectANNLabel_r_hippocampus.nii.gz
    subjectANNLabel_r_putamen.nii.gz
    subjectANNLabel_r_thalamus.nii.gz
    """

    RF8BC.inputs.outputBinaryLeftCaudate=      'subjectANNLabel_l_caudate.nii.gz'
    RF8BC.inputs.outputBinaryRightCaudate=     'subjectANNLabel_r_caudate.nii.gz'
    RF8BC.inputs.outputBinaryLeftHippocampus=  'subjectANNLabel_l_hippocampus.nii.gz'
    RF8BC.inputs.outputBinaryRightHippocampus= 'subjectANNLabel_r_hippocampus.nii.gz'
    RF8BC.inputs.outputBinaryLeftPutamen=      'subjectANNLabel_l_putamen.nii.gz'
    RF8BC.inputs.outputBinaryRightPutamen=     'subjectANNLabel_r_putamen.nii.gz'
    RF8BC.inputs.outputBinaryLeftThalamus=     'subjectANNLabel_l_thalamus.nii.gz'
    RF8BC.inputs.outputBinaryRightThalamus=    'subjectANNLabel_r_thalamus.nii.gz'

    cutWF.connect(inputsSpec,'T1Volume',RF8BC,'inputSubjectT1Filename')
    cutWF.connect(inputsSpec,'T2Volume',RF8BC,'inputSubjectT2Filename')
    cutWF.connect(inputsSpec,'TotalGM',RF8BC,'inputSubjectTotalGMFilename')
    cutWF.connect(inputsSpec,'RegistrationROI',RF8BC,'inputSubjectRegistrationROIFilename')
    # Error cutWF.connect(SGI,'outputVolume',RF8BC,'inputSubjectGadSGFilename')
    cutWF.connect(SGI,'outputFileName',RF8BC,'inputSubjectGadSGFilename')
    cutWF.connect(atlasObject,'template_t1',RF8BC,'inputTemplateT1')
    cutWF.connect(atlasObject,'template_brain',RF8BC,'inputTemplateRegistrationROIFilename')

    cutWF.connect(atlasObject,'rho',RF8BC,'inputTemplateRhoFilename')
    cutWF.connect(atlasObject,'phi',RF8BC,'inputTemplatePhiFilename')
    cutWF.connect(atlasObject,'theta',RF8BC,'inputTemplateThetaFilename')

    cutWF.connect(atlasObject,'l_caudate_ProbabilityMap',RF8BC,'probabilityMapsLeftCaudate')
    cutWF.connect(atlasObject,'r_caudate_ProbabilityMap',RF8BC,'probabilityMapsRightCaudate')
    cutWF.connect(atlasObject,'l_hippocampus_ProbabilityMap',RF8BC,'probabilityMapsLeftHippocampus')
    cutWF.connect(atlasObject,'r_hippocampus_ProbabilityMap',RF8BC,'probabilityMapsRightHippocampus')
    cutWF.connect(atlasObject,'l_putamen_ProbabilityMap',RF8BC,'probabilityMapsLeftPutamen')
    cutWF.connect(atlasObject,'r_putamen_ProbabilityMap',RF8BC,'probabilityMapsRightPutamen')
    cutWF.connect(atlasObject,'l_thalamus_ProbabilityMap',RF8BC,'probabilityMapsLeftThalamus')
    cutWF.connect(atlasObject,'r_thalamus_ProbabilityMap',RF8BC,'probabilityMapsRightThalamus')
    ##TODO:
    cutWF.connect(atlasObject,'RandomForestAllSubcorticalsBalancedModel_txtD0060NT0060_gz',RF8BC,'modelFilename')
    ##HACK: Needs to be fixed
    #RF8BC.inputs.modelFilename='/nfsscratch/PREDICT/TEST_BRAINSCut/20120828ANNModel_Model_RF100.txt'

    cutWF.connect(inputsSpec,'atlasToSubjectTransform',RF8BC,'deformationFromTemplateToSubject')

    mergeAllLabels=pe.Node(interface=Merge(8),name="labelMergeNode")
    # NOTE: Ordering is important
    cutWF.connect(RF8BC,'outputBinaryLeftCaudate',mergeAllLabels,'in1')
    cutWF.connect(RF8BC,'outputBinaryRightCaudate',mergeAllLabels,'in2')
    cutWF.connect(RF8BC,'outputBinaryLeftPutamen',mergeAllLabels,'in3')
    cutWF.connect(RF8BC,'outputBinaryRightPutamen',mergeAllLabels,'in4')
    cutWF.connect(RF8BC,'outputBinaryLeftHippocampus',mergeAllLabels,'in5')
    cutWF.connect(RF8BC,'outputBinaryRightHippocampus',mergeAllLabels,'in6')
    cutWF.connect(RF8BC,'outputBinaryLeftThalamus',mergeAllLabels,'in7')
    cutWF.connect(RF8BC,'outputBinaryRightThalamus',mergeAllLabels,'in8')

    computeOneLabelMap = pe.Node(interface=Function(['listOfImages','LabelImageName','CSVFileName',
         'projectid', 'subjectid', 'sessionid' ],
        ['outputLabelImageName','outputCSVFileName'],
        function=CreateLabelMap),name="ComputeOneLabelMap")
    computeOneLabelMap.inputs.projectid=projectid
    computeOneLabelMap.inputs.subjectid=subjectid
    computeOneLabelMap.inputs.sessionid=sessionid
    computeOneLabelMap.inputs.LabelImageName="allLabels_seg.nii.gz"
    computeOneLabelMap.inputs.CSVFileName = "allLabels_seg.csv"
    cutWF.connect(mergeAllLabels,'out',computeOneLabelMap,'listOfImages')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=[
        'outputBinaryLeftCaudate','outputBinaryRightCaudate',
        'outputBinaryLeftHippocampus','outputBinaryRightHippocampus',
        'outputBinaryLeftPutamen','outputBinaryRightPutamen',
        'outputBinaryLeftThalamus','outputBinaryRightThalamus',
        'outputLabelImageName','outputCSVFileName',
        'xmlFilename']), name='OutputSpec' )

    cutWF.connect(computeOneLabelMap,'outputLabelImageName',outputsSpec,'outputLabelImageName')
    cutWF.connect(computeOneLabelMap,'outputCSVFileName',outputsSpec,'outputCSVFileName')
    cutWF.connect(RF8BC,'outputBinaryLeftCaudate',outputsSpec,'outputBinaryLeftCaudate')
    cutWF.connect(RF8BC,'outputBinaryRightCaudate',outputsSpec,'outputBinaryRightCaudate')
    cutWF.connect(RF8BC,'outputBinaryLeftHippocampus',outputsSpec,'outputBinaryLeftHippocampus')
    cutWF.connect(RF8BC,'outputBinaryRightHippocampus',outputsSpec,'outputBinaryRightHippocampus')
    cutWF.connect(RF8BC,'outputBinaryLeftPutamen',outputsSpec,'outputBinaryLeftPutamen')
    cutWF.connect(RF8BC,'outputBinaryRightPutamen',outputsSpec,'outputBinaryRightPutamen')
    cutWF.connect(RF8BC,'outputBinaryLeftThalamus',outputsSpec,'outputBinaryLeftThalamus')
    cutWF.connect(RF8BC,'outputBinaryRightThalamus',outputsSpec,'outputBinaryRightThalamus')
    cutWF.connect(RF8BC,'xmlFilename',outputsSpec,'xmlFilename')

    return cutWF
