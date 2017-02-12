#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from nipype.interfaces.utility import Merge, Function, IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.semtools import *
from PipeLineFunctionHelpers import getListIndex

from .RF12BRAINSCutWrapper import RF12BRAINSCutWrapper


def GenerateWFName(projectid, subjectid, sessionid, WFName):
    return WFName + '_' + str(subjectid) + "_" + str(sessionid) + "_" + str(projectid)


def CreateLabelMap(listOfImages, LabelImageName, CSVFileName, posteriorDictionary, projectid, subjectid, sessionid):
    """
    A function to create a consolidated label map and a
    csv file of volume measurements.
    """

    import SimpleITK as sitk
    import os
    import csv

    def CleanUpSegmentationsWithExclusionProbabilityMaps(initial_seg, probMapOfExclusion, percentageThreshold=0.85):
        """This function is used to clean up grey matter sub-cortical segmentations
    by removing tissue that is more than 85% chance of being either WM or CSF
    The inputs are the initial segmentation, the WM Probability, and the CSF Probability
    """
        seg = sitk.Cast(initial_seg, sitk.sitkUInt8)
        print("AA", initial_seg)
        # print "BB", dict(sitk.Statistics(seg))
        exclude_Mask = sitk.Cast(sitk.BinaryThreshold(probMapOfExclusion, percentageThreshold, 1.0, 0, 1), sitk.sitkUInt8)
        # print "CC", dict(sitk.Statistics(exclude_Mask))
        cleanedUpSeg = seg * exclude_Mask
        # print "DD", dict(sitk.Statistics(cleanedUpSeg))
        return cleanedUpSeg

    def CleanUpGMSegmentationWithWMCSF(initial_seg_fn, posteriorDictionary, WMThreshold, CSFThreshold):
        initial_seg = sitk.Cast(sitk.ReadImage(initial_seg_fn.encode('ascii','replace')), sitk.sitkUInt8)

        #WM_FN = posteriorDictionary['WM']
        #WM_PROB = sitk.ReadImage(WM_FN.encode('ascii','replace'))
        #WM_removed = CleanUpSegmentationsWithExclusionProbabilityMaps(initial_seg, WM_PROB, WMThreshold)

        CSF_FN = posteriorDictionary['CSF']
        CSF_PROB = sitk.ReadImage(CSF_FN.encode('ascii','replace'))
        CSF_removed = CleanUpSegmentationsWithExclusionProbabilityMaps(initial_seg, CSF_PROB, CSFThreshold)
        return CSF_removed

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

    cleaned_labels_map = dict()
    labelImage = None
    print("ZZZ")
    x = 0
    for segFN in listOfImages:
        x = x + 1
        print(x, segFN)
        ## Clean up the segmentations
        curr_segROI = CleanUpGMSegmentationWithWMCSF(segFN, posteriorDictionary, 0.85, 0.85)
        print("Y")
        curr_segROI.GetSize()
        remove_pre_postfix = segFN.replace(".nii.gz", "")
        remove_pre_postfix = os.path.basename(remove_pre_postfix.replace("subjectANNLabel_", "").replace("_seg", ""))
        remove_pre_postfix = os.path.basename(remove_pre_postfix.replace("ANNContinuousPrediction", "").replace("subject", ""))
        structName = remove_pre_postfix.lower()
        cleaned_fileName = os.path.join(os.path.dirname(segFN), "cleaned_" + structName + "_seg.nii.gz")
        print("=" * 20, structName, " ", cleaned_fileName)
        cleaned_labels_map[structName] = cleaned_fileName
        sitk.WriteImage(curr_segROI, cleaned_fileName.encode('ascii','replace'))
        if labelImage is None:
            labelImage = curr_segROI * valueDict[structName]
        else:
            not_mask = sitk.Not(curr_segROI)
            ## Clear out an empty space for the next mask to be inserted
            labelImage *= not_mask
            ## Add in the mask image with it's proper label
            labelImage = labelImage + curr_segROI * valueDict[structName]
    sitk.WriteImage(labelImage, LabelImageName.encode('ascii','replace'))

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
            # myMeasurementMap = ls.GetMeasurementMap(value)
            # dictKeys = myMeasurementMap.GetVectorOfMeasurementNames()
            # dictValues = myMeasurementMap.GetVectorOfMeasurementValues()
            # measurementDict = dict(zip(dictKeys, dictValues))
            structVolume = ImageSpacing[0] * ImageSpacing[1] * ImageSpacing[2] * ls.GetCount(value) # measurementDict['Count']
            writeDictionary['Volume_mm3'] = structVolume
            writeDictionary['Structure'] = name
            writeDictionary['LabelCode'] = value
            # writeDictionary['FileName']=os.path.abspath(LabelImageName)
            writeDictionary['projectid'] = projectid
            writeDictionary['subjectid'] = subjectid
            writeDictionary['sessionid'] = sessionid
            dWriter.writerow(writeDictionary)

    CleanedLeftCaudate = cleaned_labels_map['l_caudate']
    CleanedRightCaudate = cleaned_labels_map['r_caudate']
    CleanedLeftHippocampus = cleaned_labels_map['l_hippocampus']
    CleanedRightHippocampus = cleaned_labels_map['r_hippocampus']
    CleanedLeftPutamen = cleaned_labels_map['l_putamen']
    CleanedRightPutamen = cleaned_labels_map['r_putamen']
    CleanedLeftThalamus = cleaned_labels_map['l_thalamus']
    CleanedRightThalamus = cleaned_labels_map['r_thalamus']
    CleanedLeftAccumben = cleaned_labels_map['l_accumben']
    CleanedRightAccumben = cleaned_labels_map['r_accumben']
    CleanedLeftGlobus = cleaned_labels_map['l_globus']
    CleanedRightGlobus = cleaned_labels_map['r_globus']
    return os.path.abspath(LabelImageName), os.path.abspath(CSVFileName), CleanedLeftCaudate, CleanedRightCaudate, CleanedLeftHippocampus, CleanedRightHippocampus, CleanedLeftPutamen, CleanedRightPutamen, CleanedLeftThalamus, CleanedRightThalamus, CleanedLeftAccumben, CleanedRightAccumben, CleanedLeftGlobus, CleanedRightGlobus

#==============================================
#==============================================
#==============================================
#==============================================
#==============================================
#==============================================

def CreateBRAINSCutWorkflow(projectid,
                            subjectid,
                            sessionid,
                            CLUSTER_QUEUE,
                            CLUSTER_QUEUE_LONG,
                            WFName,
                            t1Only):
    cutWF = pe.Workflow(name=GenerateWFName(projectid, subjectid, sessionid, WFName))

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T1Volume', 'T2Volume',
                                                             'posteriorDictionary', 'RegistrationROI',
                                                             'atlasToSubjectTransform', 'template_t1_denoised_gaussian',
                                                             'rho', 'phi', 'theta',
                                                             'l_caudate_ProbabilityMap', 'r_caudate_ProbabilityMap',
                                                             'l_hippocampus_ProbabilityMap', 'r_hippocampus_ProbabilityMap',
                                                             'l_putamen_ProbabilityMap', 'r_putamen_ProbabilityMap',
                                                             'l_thalamus_ProbabilityMap', 'r_thalamus_ProbabilityMap',
                                                             'l_accumben_ProbabilityMap', 'r_accumben_ProbabilityMap',
                                                             'l_globus_ProbabilityMap', 'r_globus_ProbabilityMap',
                                                             'trainModelFile_txtD0060NT0060_gz']),
                                                              name='inputspec')

    #Denoised T1 input for BRAINSCut
    denosingTimeStep = 0.0625
    denosingConductance = 0.4
    denosingIteration = 5

    DenoisedT1 = pe.Node(interface=GradientAnisotropicDiffusionImageFilter(), name="DenoisedT1")
    DenoisedT1.inputs.timeStep = denosingTimeStep
    DenoisedT1.inputs.conductance = denosingConductance
    DenoisedT1.inputs.numberOfIterations = denosingIteration
    DenoisedT1.inputs.outputVolume = "DenoisedT1.nii.gz"

    cutWF.connect(inputsSpec, 'T1Volume', DenoisedT1, 'inputVolume')

    #Gradient Anistropic Diffusion T1 images for BRAINSCut
    GADT1 = pe.Node(interface=GradientAnisotropicDiffusionImageFilter(), name="GADT1")
    GADT1.inputs.timeStep = 0.025
    GADT1.inputs.conductance = 1
    GADT1.inputs.numberOfIterations = 5
    GADT1.inputs.outputVolume = "GADT1.nii.gz"

    cutWF.connect(inputsSpec, 'T1Volume', GADT1, 'inputVolume')

    if not t1Only:
        #Denoised T1 input for BRAINSCut
        DenoisedT2 = pe.Node(interface=GradientAnisotropicDiffusionImageFilter(), name="DenoisedT2")
        DenoisedT2.inputs.timeStep = denosingTimeStep
        DenoisedT2.inputs.conductance = denosingConductance
        DenoisedT2.inputs.numberOfIterations = denosingIteration
        DenoisedT2.inputs.outputVolume = "DenoisedT2.nii.gz"

        cutWF.connect(inputsSpec, 'T2Volume', DenoisedT2, 'inputVolume')

        #Gradient Anistropic Diffusion T1 images for BRAINSCut
        GADT2 = pe.Node(interface=GradientAnisotropicDiffusionImageFilter(), name="GADT2")
        GADT2.inputs.timeStep = 0.025
        GADT2.inputs.conductance = 1
        GADT2.inputs.numberOfIterations = 5
        GADT2.inputs.outputVolume = "GADT2.nii.gz"
        cutWF.connect(inputsSpec, 'T2Volume', GADT2, 'inputVolume')

        #Sum the gradient images for BRAINSCut
        SGI = pe.Node(interface=GenerateSummedGradientImage(), name="SGI")
        SGI.inputs.outputFileName = "SummedGradImage.nii.gz"

        cutWF.connect(GADT1, 'outputVolume', SGI, 'inputVolume1')
        cutWF.connect(GADT2, 'outputVolume', SGI, 'inputVolume2')

    #BRAINSCut
    RF12BC = pe.Node(interface=RF12BRAINSCutWrapper(), name="IQR_NORM_SEP_RF12_BRAINSCut")
    # HACK
    # import os
    # RF12BC.inputs.environ = dict(os.environ)
    # many_cpu_RF12BC_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,4,2,2), 'overwrite': True}
    # RF12BC.plugin_args = many_cpu_RF12BC_options_dictionary
    # END HACK
    RF12BC.inputs.trainingVectorFilename = "trainingVectorFilename.txt"
    RF12BC.inputs.xmlFilename = "BRAINSCutSegmentationDefinition.xml"
    RF12BC.inputs.vectorNormalization = "IQR"

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

    cutWF.connect(DenoisedT1, 'outputVolume', RF12BC, 'inputSubjectT1Filename')

    from PipeLineFunctionHelpers import MakeInclusionMaskForGMStructures
    makeCandidateRegionNode = pe.Node(interface=Function(['posteriorDictionary', 'candidateRegionFileName'],
                                                         ['outputCandidateRegionFileName'],
                                                         function=MakeInclusionMaskForGMStructures), name="MakeCandidateRegion")
    makeCandidateRegionNode.inputs.candidateRegionFileName = "RF12_CandidateRegionMask.nii.gz"
    cutWF.connect(inputsSpec, 'posteriorDictionary', makeCandidateRegionNode, 'posteriorDictionary')
    cutWF.connect(makeCandidateRegionNode, 'outputCandidateRegionFileName', RF12BC, 'candidateRegion')

    cutWF.connect([(inputsSpec, RF12BC, [('template_t1_denoised_gaussian', 'inputTemplateT1'),
                                         # ('template_brain', 'inputTemplateRegistrationROIFilename'),
                                         ('rho', 'inputTemplateRhoFilename'),
                                         ('phi', 'inputTemplatePhiFilename'),
                                         ('theta', 'inputTemplateThetaFilename'),
                                         ('l_caudate_ProbabilityMap', 'probabilityMapsLeftCaudate'),
                                         ('r_caudate_ProbabilityMap', 'probabilityMapsRightCaudate'),
                                         ('l_hippocampus_ProbabilityMap', 'probabilityMapsLeftHippocampus'),
                                         ('r_hippocampus_ProbabilityMap', 'probabilityMapsRightHippocampus'),
                                         ('l_putamen_ProbabilityMap', 'probabilityMapsLeftPutamen'),
                                         ('r_putamen_ProbabilityMap', 'probabilityMapsRightPutamen'),
                                         ('l_thalamus_ProbabilityMap', 'probabilityMapsLeftThalamus'),
                                         ('r_thalamus_ProbabilityMap', 'probabilityMapsRightThalamus'),
                                         ('l_accumben_ProbabilityMap', 'probabilityMapsLeftAccumben'),
                                         ('r_accumben_ProbabilityMap', 'probabilityMapsRightAccumben'),
                                         ('l_globus_ProbabilityMap', 'probabilityMapsLeftGlobus'),
                                         ('r_globus_ProbabilityMap', 'probabilityMapsRightGlobus'),
                 ])])

    # TODO:
    if not t1Only:
        cutWF.connect(DenoisedT2, 'outputVolume', RF12BC, 'inputSubjectT2Filename')
        # cutWF.connect(inputsSpec,'TotalGM',RF12BC,'inputSubjectTotalGMFilename')
        # cutWF.connect(inputsSpec,'RegistrationROI',RF12BC,'inputSubjectRegistrationROIFilename')
        # Error cutWF.connect(SGI,'outputVolume',RF12BC,'inputSubjectGadSGFilename')
        cutWF.connect(SGI, 'outputFileName', RF12BC, 'inputSubjectGadSGFilename')
        cutWF.connect(inputsSpec, 'trainModelFile_txtD0060NT0060_gz', RF12BC, 'modelFilename')
    else:
        ### TODO:  Replace with proper atlas file name in the future!!! This is a HACK
        ### to avoid changing the hash keys of the input files from the atlas.
        def ChangeModelPathDirectory(multiModalFileName):
            return multiModalFileName.replace('modelFiles', 'T1OnlyModels')
        cutWF.connect([(inputsSpec, RF12BC,
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
                                                     'posteriorDictionary',
                                                     'projectid', 'subjectid', 'sessionid'],
                                                    ['outputLabelImageName', 'outputCSVFileName',
                                                     'CleanedLeftCaudate',
                                                     'CleanedRightCaudate',
                                                     'CleanedLeftHippocampus',
                                                     'CleanedRightHippocampus',
                                                     'CleanedLeftPutamen',
                                                     'CleanedRightPutamen',
                                                     'CleanedLeftThalamus',
                                                     'CleanedRightThalamus',
                                                     'CleanedLeftAccumben',
                                                     'CleanedRightAccumben',
                                                     'CleanedLeftGlobus',
                                                     'CleanedRightGlobus'
                                                     ],
                                                    function=CreateLabelMap), name="ComputeOneLabelMap")
    computeOneLabelMap.inputs.projectid = projectid
    computeOneLabelMap.inputs.subjectid = subjectid
    computeOneLabelMap.inputs.sessionid = sessionid
    computeOneLabelMap.inputs.LabelImageName = "allLabels.nii.gz"
    computeOneLabelMap.inputs.CSVFileName = "allLabels_seg.csv"
    cutWF.connect(inputsSpec, 'posteriorDictionary', computeOneLabelMap, 'posteriorDictionary')
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
    cutWF.connect(computeOneLabelMap, 'CleanedLeftCaudate', outputsSpec, 'outputBinaryLeftCaudate')
    cutWF.connect(computeOneLabelMap, 'CleanedRightCaudate', outputsSpec, 'outputBinaryRightCaudate')
    cutWF.connect(computeOneLabelMap, 'CleanedLeftHippocampus', outputsSpec, 'outputBinaryLeftHippocampus')
    cutWF.connect(computeOneLabelMap, 'CleanedRightHippocampus', outputsSpec, 'outputBinaryRightHippocampus')
    cutWF.connect(computeOneLabelMap, 'CleanedLeftPutamen', outputsSpec, 'outputBinaryLeftPutamen')
    cutWF.connect(computeOneLabelMap, 'CleanedRightPutamen', outputsSpec, 'outputBinaryRightPutamen')
    cutWF.connect(computeOneLabelMap, 'CleanedLeftThalamus', outputsSpec, 'outputBinaryLeftThalamus')
    cutWF.connect(computeOneLabelMap, 'CleanedRightThalamus', outputsSpec, 'outputBinaryRightThalamus')
    cutWF.connect(computeOneLabelMap, 'CleanedLeftAccumben', outputsSpec, 'outputBinaryLeftAccumben')
    cutWF.connect(computeOneLabelMap, 'CleanedRightAccumben', outputsSpec, 'outputBinaryRightAccumben')
    cutWF.connect(computeOneLabelMap, 'CleanedLeftGlobus', outputsSpec, 'outputBinaryLeftGlobus')
    cutWF.connect(computeOneLabelMap, 'CleanedRightGlobus', outputsSpec, 'outputBinaryRightGlobus')

    cutWF.connect(RF12BC, 'xmlFilename', outputsSpec, 'xmlFilename')

    return cutWF
