#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from SEMTools import *
from RF12BRAINSCutWrapper import RF12BRAINSCutWrapper

from PipeLineFunctionHelpers import getListIndex


def GenerateWFName(projectid, subjectid, sessionid, WFName):
    return WFName + '_' + str(subjectid) + "_" + str(sessionid) + "_" + str(projectid)


def CreateLabelMap(listOfImages, LabelImageName, CSVFileName, projectid, subjectid, sessionid):
    """
    A function to create a consolidated label map and a
    csv file of volume measurements.
    """

    import SimpleITK as sitk
    import os
    import csv
    orderOfPriority = [
        "l_caudate",
        "r_caudate",
        "l_putamen",
        "r_putamen",
        "l_hippocampus",
        "r_hippocampus",
        "l_thalamus",
        "r_thalamus",
        "l_accumben",
        "r_accumben",
        "l_globus",
        "r_globus"
    ]

    valueDict = {
        "l_caudate": 1,
        "r_caudate": 2,
        "l_putamen": 3,
        "r_putamen": 4,
        "l_hippocampus": 5,
        "r_hippocampus": 6,
        "l_thalamus": 7,
        "r_thalamus": 8,
        "l_accumben": 9,
        "r_accumben": 10,
        "l_globus": 11,
        "r_globus": 12
    }

    labelImage = None
    for segFN in listOfImages:
        im = sitk.ReadImage(segFN)
        im.GetSize()
        remove_pre_postfix = os.path.basename(segFN.replace(".nii.gz", "").replace("subjectANNLabel_", "").replace("_seg", ""))
        structName = remove_pre_postfix.lower()
        if labelImage is None:
            labelImage = im * valueDict[structName]
        else:
            mask = sitk.Not(im)
            ## Clear out an empty space for the next mask to be inserted
            labelImage *= mask
            ## Add in the mask image with it's proper label
            labelImage = labelImage + im * valueDict[structName]
    sitk.WriteImage(labelImage, LabelImageName)

    ls = sitk.LabelStatisticsImageFilter()
    ls.Execute(labelImage, labelImage)
    ImageSpacing = labelImage.GetSpacing()
    csvFile = open(CSVFileName, 'w')
    dWriter = csv.DictWriter(csvFile, ['projectid', 'subjectid', 'sessionid', 'Structure', 'LabelCode', 'Volume_mm3'], restval='', extrasaction='raise', dialect='excel')
    dWriter.writeheader()
    writeDictionary = dict()
    for name in orderOfPriority:
        value = valueDict[name]
        if ls.HasLabel(value):
            # print "Displaying: ", name, value
            myMeasurementMap = ls.GetMeasurementMap(value)
            dictKeys = myMeasurementMap.GetVectorOfMeasurementNames()
            dictValues = myMeasurementMap.GetVectorOfMeasurementValues()
            measurementDict = dict(zip(dictKeys, dictValues))
            structVolume = ImageSpacing[0] * ImageSpacing[1] * ImageSpacing[2] * measurementDict['Count']
            writeDictionary['Volume_mm3'] = structVolume
            writeDictionary['Structure'] = name
            writeDictionary['LabelCode'] = value
            # writeDictionary['FileName']=os.path.abspath(LabelImageName)
            writeDictionary['projectid'] = projectid
            writeDictionary['subjectid'] = subjectid
            writeDictionary['sessionid'] = sessionid
            dWriter.writerow(writeDictionary)
    return os.path.abspath(LabelImageName), os.path.abspath(CSVFileName)

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


def CreateBRAINSCutWorkflow(projectid,
                            subjectid,
                            sessionid,
                            WFName,
                            CLUSTER_QUEUE,
                            CLUSTER_QUEUE_LONG,
                            atlasObject,
                            t1Only=False):
    cutWF = pe.Workflow(name=GenerateWFName(projectid, subjectid, sessionid, WFName))

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T1Volume', 'T2Volume',
                                                             'posteriorDictionary', 'RegistrationROI',
                                                             'atlasToSubjectTransform']), name='inputspec')


    if not t1Only:
        """
        Denoised input for BRAINSCut
        """
        denosingTimeStep=0.0625
        denosingConductance=0.4
        denosingIteration=5

        DenoisedT1 = pe.Node(interface=GradientAnisotropicDiffusionImageFilter(), name="DenoisedT1")
        DenoisedT1.inputs.timeStep = denosingTimeStep
        DenoisedT1.inputs.conductance = denosingConductance
        DenoisedT1.inputs.numberOfIterations = denosingIteration
        DenoisedT1.inputs.outputVolume = "DenoisedT1.nii.gz"

        cutWF.connect(inputsSpec, 'T1Volume', DenoisedT1, 'inputVolume')

        DenoisedT2 = pe.Node(interface=GradientAnisotropicDiffusionImageFilter(), name="DenoisedT2")
        DenoisedT2.inputs.timeStep = denosingTimeStep
        DenoisedT2.inputs.conductance = denosingConductance
        DenoisedT2.inputs.numberOfIterations = denosingIteration
        DenoisedT2.inputs.outputVolume = "DenoisedT2.nii.gz"

        cutWF.connect(inputsSpec, 'T2Volume', DenoisedT2, 'inputVolume')

        """
        Gradient Anistropic Diffusion images for BRAINSCut
        """
        GADT1 = pe.Node(interface=GradientAnisotropicDiffusionImageFilter(), name="GADT1")
        GADT1.inputs.timeStep = 0.025
        GADT1.inputs.conductance = 1
        GADT1.inputs.numberOfIterations = 5
        GADT1.inputs.outputVolume = "GADT1.nii.gz"

        cutWF.connect(inputsSpec, 'T1Volume', GADT1, 'inputVolume')

        GADT2 = pe.Node(interface=GradientAnisotropicDiffusionImageFilter(), name="GADT2")
        GADT2.inputs.timeStep = 0.025
        GADT2.inputs.conductance = 1
        GADT2.inputs.numberOfIterations = 5
        GADT2.inputs.outputVolume = "GADT2.nii.gz"
        cutWF.connect(inputsSpec, 'T2Volume', GADT2, 'inputVolume')

        """
        Sum the gradient images for BRAINSCut
        """
        SGI = pe.Node(interface=GenerateSummedGradientImage(), name="SGI")
        SGI.inputs.outputFileName = "SummedGradImage.nii.gz"

        cutWF.connect(GADT1, 'outputVolume', SGI, 'inputVolume1')
        cutWF.connect(GADT2, 'outputVolume', SGI, 'inputVolume2')


    """
    BRAINSCut
    """
    RF12BC = pe.Node(interface=RF12BRAINSCutWrapper(), name="IQR_NORM_SEP_RF12_BRAINSCut")
    many_cpu_RF12BC_options_dictionary = {'qsub_args': '-S /bin/bash -pe smp1 2-2 -l h_vmem=8G,mem_free=4G -o /dev/null -e /dev/null ' + CLUSTER_QUEUE, 'overwrite': True}
# many_cpu_RF12BC_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 2-8 -l big_mem=true,h_vmem=60G,mem_free=30G -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
# many_cpu_RF12BC_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 4-6 -l big_mem=true,h_vmem=22G,mem_free=22G -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
    RF12BC.plugin_args = many_cpu_RF12BC_options_dictionary
    RF12BC.inputs.trainingVectorFilename = "trainingVectorFilename.txt"
    RF12BC.inputs.xmlFilename = "BRAINSCutSegmentationDefinition.xml"
    RF12BC.inputs.vectorNormalization = "IQR"

    """ HACK These should be l_caudate_seg.nii.gz
    subjectANNLabel_l_caudate.nii.gz
    subjectANNLabel_l_hippocampus.nii.gz
    subjectANNLabel_l_putamen.nii.gz
    subjectANNLabel_l_thalamus.nii.gz
    subjectANNLabel_r_caudate.nii.gz
    subjectANNLabel_r_hippocampus.nii.gz
    subjectANNLabel_r_putamen.nii.gz
    subjectANNLabel_r_thalamus.nii.gz
    """

    """
    RF12BC.inputs.outputBinaryLeftCaudate = 'subjectANNLabel_l_caudate.nii.gz'
    RF12BC.inputs.outputBinaryRightCaudate = 'subjectANNLabel_r_caudate.nii.gz'
    RF12BC.inputs.outputBinaryLeftHippocampus = 'subjectANNLabel_l_hippocampus.nii.gz'
    RF12BC.inputs.outputBinaryRightHippocampus = 'subjectANNLabel_r_hippocampus.nii.gz'
    RF12BC.inputs.outputBinaryLeftPutamen = 'subjectANNLabel_l_putamen.nii.gz'
    RF12BC.inputs.outputBinaryRightPutamen = 'subjectANNLabel_r_putamen.nii.gz'
    RF12BC.inputs.outputBinaryLeftThalamus = 'subjectANNLabel_l_thalamus.nii.gz'
    RF12BC.inputs.outputBinaryRightThalamus = 'subjectANNLabel_r_thalamus.nii.gz'
    RF12BC.inputs.outputBinaryLeftAccumben = 'subjectANNLabel_l_accumben.nii.gz'
    RF12BC.inputs.outputBinaryRightAccumben = 'subjectANNLabel_r_accumben.nii.gz'
    RF12BC.inputs.outputBinaryLeftGlobus = 'subjectANNLabel_l_globus.nii.gz'
    RF12BC.inputs.outputBinaryRightGlobus = 'subjectANNLabel_r_globus.nii.gz'
    """
    "ANNContinuousPredictionl_accumbensubject.nii.gz"
    RF12BC.inputs.outputBinaryLeftCaudate = 'ANNContinuousPredictionl_caudatesubject.nii.gz'
    RF12BC.inputs.outputBinaryRightCaudate = 'ANNContinuousPredictionr_caudatesubject.nii.gz'
    RF12BC.inputs.outputBinaryLeftHippocampus = 'ANNContinuousPredictionl_hippocampussubject.nii.gz'
    RF12BC.inputs.outputBinaryRightHippocampus = 'ANNContinuousPredictionr_hippocampussubject.nii.gz'
    RF12BC.inputs.outputBinaryLeftPutamen = 'ANNContinuousPredictionl_putamensubject.nii.gz'
    RF12BC.inputs.outputBinaryRightPutamen = 'ANNContinuousPredictionr_putamensubject.nii.gz'
    RF12BC.inputs.outputBinaryLeftThalamus = 'ANNContinuousPredictionl_thalamussubject.nii.gz'
    RF12BC.inputs.outputBinaryRightThalamus = 'ANNContinuousPredictionr_thalamussubject.nii.gz'
    RF12BC.inputs.outputBinaryLeftAccumben = 'ANNContinuousPredictionl_accumbensubject.nii.gz'
    RF12BC.inputs.outputBinaryRightAccumben = 'ANNContinuousPredictionr_accumbensubject.nii.gz'
    RF12BC.inputs.outputBinaryLeftGlobus = 'ANNContinuousPredictionl_globussubject.nii.gz'
    RF12BC.inputs.outputBinaryRightGlobus = 'ANNContinuousPredictionr_globussubject.nii.gz'

    cutWF.connect( DenoisedT1, 'outputVolume', RF12BC, 'inputSubjectT1Filename')

    from PipeLineFunctionHelpers import MakeInclusionMaskForGMStructures;
    makeCandidateRegionNode = pe.Node(interface=Function(['posteriorDictionary','candidateRegionFileName'],
                                                    ['outputCandidateRegionFileName'],
                                                    function=MakeInclusionMaskForGMStructures), name="MakeCandidateRegion")
    makeCandidateRegionNode.inputs.candidateRegionFileName = "RF12_CandidateRegionMask.nii.gz"
    cutWF.connect(inputsSpec,'posteriorDictionary',makeCandidateRegionNode,'posteriorDictionary')
    cutWF.connect(makeCandidateRegionNode,'outputCandidateRegionFileName', RF12BC,'candidateRegion')

    if not t1Only:
        cutWF.connect( DenoisedT2, 'outputVolume',  RF12BC, 'inputSubjectT2Filename')
        # cutWF.connect(inputsSpec,'TotalGM',RF12BC,'inputSubjectTotalGMFilename')
        # cutWF.connect(inputsSpec,'RegistrationROI',RF12BC,'inputSubjectRegistrationROIFilename')
        # Error cutWF.connect(SGI,'outputVolume',RF12BC,'inputSubjectGadSGFilename')
        cutWF.connect(SGI, 'outputFileName', RF12BC, 'inputSubjectGadSGFilename')
    cutWF.connect(atlasObject, 'template_t1', RF12BC, 'inputTemplateT1')
    # cutWF.connect(atlasObject,'template_brain',RF12BC,'inputTemplateRegistrationROIFilename')

    cutWF.connect(atlasObject, 'rho', RF12BC, 'inputTemplateRhoFilename')
    cutWF.connect(atlasObject, 'phi', RF12BC, 'inputTemplatePhiFilename')
    cutWF.connect(atlasObject, 'theta', RF12BC, 'inputTemplateThetaFilename')

    cutWF.connect(atlasObject, 'l_caudate_ProbabilityMap', RF12BC, 'probabilityMapsLeftCaudate')
    cutWF.connect(atlasObject, 'r_caudate_ProbabilityMap', RF12BC, 'probabilityMapsRightCaudate')
    cutWF.connect(atlasObject, 'l_hippocampus_ProbabilityMap', RF12BC, 'probabilityMapsLeftHippocampus')
    cutWF.connect(atlasObject, 'r_hippocampus_ProbabilityMap', RF12BC, 'probabilityMapsRightHippocampus')
    cutWF.connect(atlasObject, 'l_putamen_ProbabilityMap', RF12BC, 'probabilityMapsLeftPutamen')
    cutWF.connect(atlasObject, 'r_putamen_ProbabilityMap', RF12BC, 'probabilityMapsRightPutamen')
    cutWF.connect(atlasObject, 'l_thalamus_ProbabilityMap', RF12BC, 'probabilityMapsLeftThalamus')
    cutWF.connect(atlasObject, 'r_thalamus_ProbabilityMap', RF12BC, 'probabilityMapsRightThalamus')
    cutWF.connect(atlasObject, 'l_accumben_ProbabilityMap', RF12BC, 'probabilityMapsLeftAccumben')
    cutWF.connect(atlasObject, 'r_accumben_ProbabilityMap', RF12BC, 'probabilityMapsRightAccumben')
    cutWF.connect(atlasObject, 'l_globus_ProbabilityMap', RF12BC, 'probabilityMapsLeftGlobus')
    cutWF.connect(atlasObject, 'r_globus_ProbabilityMap', RF12BC, 'probabilityMapsRightGlobus')
    # TODO:
    if not t1Only:
        cutWF.connect(atlasObject, 'trainModelFile_txtD0060NT0060_gz', RF12BC, 'modelFilename')
    else:
        ### TODO:  Replace with proper atlasObject name in the future!!! This is a HACK
        ### to avoid changing the hash keys of the input files from the atlas.
        def ChangeModelPathDirectory(multiModalFileName):
                return multiModalFileName.replace('modelFiles', 'T1OnlyModels')
        cutWF.connect([(atlasObject, RF12BC,
                        [(('trainModelFile_txtD0060NT0060_gz', ChangeModelPathDirectory), 'modelFilename')])])

    ## Need to index from next line cutWF.connect(inputsSpec,'atlasToSubjectTransform',RF12BC,'deformationFromTemplateToSubject')
    cutWF.connect([(inputsSpec, RF12BC, [(('atlasToSubjectTransform', getListIndex, 0), 'deformationFromTemplateToSubject')]), ])

    mergeAllLabels = pe.Node(interface=Merge(12), name="labelMergeNode")
    # NOTE: Ordering is important
    cutWF.connect(RF12BC, 'outputBinaryLeftCaudate', mergeAllLabels, 'in1')
    cutWF.connect(RF12BC, 'outputBinaryRightCaudate', mergeAllLabels, 'in2')
    cutWF.connect(RF12BC, 'outputBinaryLeftPutamen', mergeAllLabels, 'in3')
    cutWF.connect(RF12BC, 'outputBinaryRightPutamen', mergeAllLabels, 'in4')
    cutWF.connect(RF12BC, 'outputBinaryLeftHippocampus', mergeAllLabels, 'in5')
    cutWF.connect(RF12BC, 'outputBinaryRightHippocampus', mergeAllLabels, 'in6')
    cutWF.connect(RF12BC, 'outputBinaryLeftThalamus', mergeAllLabels, 'in7')
    cutWF.connect(RF12BC, 'outputBinaryRightThalamus', mergeAllLabels, 'in8')  # HACK:  CHECK ORDERING
    cutWF.connect(RF12BC, 'outputBinaryLeftAccumben', mergeAllLabels, 'in9')
    cutWF.connect(RF12BC, 'outputBinaryRightAccumben', mergeAllLabels, 'in10')
    cutWF.connect(RF12BC, 'outputBinaryLeftGlobus', mergeAllLabels, 'in11')
    cutWF.connect(RF12BC, 'outputBinaryRightGlobus', mergeAllLabels, 'in12')

    computeOneLabelMap = pe.Node(interface=Function(['listOfImages', 'LabelImageName', 'CSVFileName',
                                                     'projectid', 'subjectid', 'sessionid'],
                                                    ['outputLabelImageName', 'outputCSVFileName'],
                                                    function=CreateLabelMap), name="ComputeOneLabelMap")
    computeOneLabelMap.inputs.projectid = projectid
    computeOneLabelMap.inputs.subjectid = subjectid
    computeOneLabelMap.inputs.sessionid = sessionid
    computeOneLabelMap.inputs.LabelImageName = "allLabels.nii.gz"
    computeOneLabelMap.inputs.CSVFileName = "allLabels_seg.csv"
    cutWF.connect(mergeAllLabels, 'out', computeOneLabelMap, 'listOfImages')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=[
        'outputBinaryLeftCaudate', 'outputBinaryRightCaudate',
        'outputBinaryLeftHippocampus', 'outputBinaryRightHippocampus',
        'outputBinaryLeftPutamen', 'outputBinaryRightPutamen',
        'outputBinaryLeftThalamus', 'outputBinaryRightThalamus',
        'outputBinaryLeftAccumben', 'outputBinaryRightAccumben',
        'outputBinaryLeftGlobus', 'outputBinaryRightGlobus',
        'outputLabelImageName', 'outputCSVFileName',
        'xmlFilename']), name='outputspec')

    cutWF.connect(computeOneLabelMap, 'outputLabelImageName', outputsSpec, 'outputLabelImageName')
    cutWF.connect(computeOneLabelMap, 'outputCSVFileName', outputsSpec, 'outputCSVFileName')
    cutWF.connect(RF12BC, 'outputBinaryLeftCaudate', outputsSpec, 'outputBinaryLeftCaudate')
    cutWF.connect(RF12BC, 'outputBinaryRightCaudate', outputsSpec, 'outputBinaryRightCaudate')
    cutWF.connect(RF12BC, 'outputBinaryLeftHippocampus', outputsSpec, 'outputBinaryLeftHippocampus')
    cutWF.connect(RF12BC, 'outputBinaryRightHippocampus', outputsSpec, 'outputBinaryRightHippocampus')
    cutWF.connect(RF12BC, 'outputBinaryLeftPutamen', outputsSpec, 'outputBinaryLeftPutamen')
    cutWF.connect(RF12BC, 'outputBinaryRightPutamen', outputsSpec, 'outputBinaryRightPutamen')
    cutWF.connect(RF12BC, 'outputBinaryLeftThalamus', outputsSpec, 'outputBinaryLeftThalamus')
    cutWF.connect(RF12BC, 'outputBinaryRightThalamus', outputsSpec, 'outputBinaryRightThalamus')
    cutWF.connect(RF12BC, 'outputBinaryLeftAccumben', outputsSpec, 'outputBinaryLeftAccumben')
    cutWF.connect(RF12BC, 'outputBinaryRightAccumben', outputsSpec, 'outputBinaryRightAccumben')
    cutWF.connect(RF12BC, 'outputBinaryLeftGlobus', outputsSpec, 'outputBinaryLeftGlobus')
    cutWF.connect(RF12BC, 'outputBinaryRightGlobus', outputsSpec, 'outputBinaryRightGlobus')
    cutWF.connect(RF12BC, 'xmlFilename', outputsSpec, 'xmlFilename')

    return cutWF
