## \author Ali Ghayoor
##
## This workflow computes the statistics of RISs over regions of interest from the input label map.
##

import os
import SimpleITK as sitk
import nipype
from nipype.interfaces import ants
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/oS
import nipype.pipeline.engine as pe  # pypeline engine
from SEMTools import *

def CreateMeasurementWorkflow(WFname, LABELS_CONFIG_FILE):
    #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    ###### UTILITY FUNCTIONS #######
    # This function returns a label map that only covers the FOV of the input DWI scan
    def CreateDWILabelMap(T2LabelMapVolume,DWIBrainMask):
        import os
        import SimpleITK as sitk
        T2LabelMapVolume = sitk.ReadImage(T2LabelMapVolume,sitk.sitkUInt8)
        DWIBrainMask = sitk.ReadImage(DWIBrainMask)
        # 1- Dilate input DWI mask
        dilateFilter = sitk.BinaryDilateImageFilter()
        dilateFilter.SetKernelRadius(1)
        dilated_mask = dilateFilter.Execute( DWIBrainMask )
        # 2- Resample dilated mask to the space of T2LabelMap (1x1x1)
        # Use Linear interpolation + thresholding
        resFilt = sitk.ResampleImageFilter()
        resFilt.SetReferenceImage(T2LabelMapVolume)
        resFilt.SetOutputPixelType(sitk.sitkFloat32)
        resFilt.SetInterpolator(sitk.sitkLinear)
        resampled_dilated_mask = resFilt.Execute(dilated_mask)
        # Thresholding by 0
        threshFilt = sitk.BinaryThresholdImageFilter()
        thresh_resampled_dilated_mask = threshFilt.Execute(resampled_dilated_mask,0.0001,1.0,1,0)
        # 3- Multiply this binary mask to the T2 labelmap volume
        DWILabelMapVolume = thresh_resampled_dilated_mask * T2LabelMapVolume
        outputVolume = os.path.realpath('DWILabelMapVolume.nrrd')
        sitk.WriteImage(DWILabelMapVolume, outputVolume)
        return outputVolume

    def MakeResamplerInFileList(FAImage,MDImage,RDImage,FrobeniusNormImage,Lambda1Image,Lambda2Image,Lambda3Image):
        RISsList = [FAImage,MDImage,RDImage,FrobeniusNormImage,Lambda1Image,Lambda2Image,Lambda3Image]
        return RISsList

    # This functions computes statistics of each input RIS volume over all input labels
    # and writes the results as a CSV file
    def ComputeStatistics(inputVolume,T2LabelMapVolume,DWILabelMapVolume,labelCodesFile):
        import os
        import SimpleITK as sitk
        #### Util Funcs ####
        def createLabelsDictionary(labelCodesFile):
            import csv
            labelsDictionary={}
            with open(labelCodesFile) as lf:
                reader = csv.reader(lf, delimiter=',')
                for line in reader:
                  if line[0][0] == "#":
                     continue
                  else:
                     labelsDictionary[line[0]] = line[1]
            return labelsDictionary

        def computeVoxelVolume(inputVolume):
            import operator
            return reduce(operator.mul, inputVolume.GetSpacing())

        def ReturnStatisticsList(labelID,voxelVolume,resampledRISVolume,DWILabelMap,T2LabelMap):
            statFilter = sitk.LabelStatisticsImageFilter()
            # RIS stats over input label ID
            statFilter.Execute(resampledRISVolume, DWILabelMap)
            mean = statFilter.GetMean(labelID)
            std = statFilter.GetSigma(labelID)
            maximum = statFilter.GetMaximum(labelID)
            minimum = statFilter.GetMinimum(labelID)
            median = statFilter.GetMedian(labelID)
            effectiveVolume = statFilter.GetCount(labelID)*voxelVolume
            # compute total volume of input label ID in the non-cropped labelmap (T2LabelMap)
            statFilter.Execute(resampledRISVolume, T2LabelMap)
            totalVolume = statFilter.GetCount(labelID)*voxelVolume
            # if effectiveVolume is 0, the label is missed in dwi scan, or it doesn't exists in current labelmaps
            # in both cases we need zero confidence coefficient for that.
            if effectiveVolume == 0:
                confidence_coeficient = 0
                maximum = 0
                minimum = 0
            else:
                confidence_coeficient=effectiveVolume/totalVolume
            # Now create statistics list
            statsList = [format(mean,'.4f'),
                         format(std,'.4f'),
                         format(maximum,'.4f'),
                         format(minimum,'.4f'),
                         format(median,'.4f'),
                         effectiveVolume,
                         totalVolume,
                         format(confidence_coeficient,'.3f')]
            return statsList

        def writeLabelStatistics(filename,statisticsDictionary):
            import csv
            with open(filename, 'wb') as lf:
                headerdata = [['#Label', 'mean', 'std', 'max', 'min', 'median', 'effective_volume', 'total_volume', 'confidence_coeficient']]
                wr = csv.writer(lf, delimiter=',')
                wr.writerows(headerdata)
                for key, value in sorted(statisticsDictionary.items()):
                    wr.writerows([[key] + value])
        #### #### #### ####
        resampledRISVolume = sitk.ReadImage(inputVolume)
        T2LabelMap = sitk.ReadImage(T2LabelMapVolume)
        DWILabelMap = sitk.ReadImage(DWILabelMapVolume)
        labelsDictionary = createLabelsDictionary(labelCodesFile)
        statisticsDictionary={}
        voxelVolume=computeVoxelVolume(resampledRISVolume)
        for key in labelsDictionary:
            labelID = int(key)
            statisticsDictionary[labelsDictionary[key]] = ReturnStatisticsList(labelID,voxelVolume,resampledRISVolume,DWILabelMap,T2LabelMap)
        inputBaseName = os.path.basename(inputVolume)
        inputName = os.path.splitext(inputBaseName)[0]
        RISName = inputName.split('_',1)[0]
        CSVStatisticsFile = os.path.realpath(RISName + '_statistics.csv')
        writeLabelStatistics(CSVStatisticsFile,statisticsDictionary)
        return CSVStatisticsFile

    # This function helps to pick desirable output from the output list
    def pickFromList(inputlist,item):
        return inputlist[item]
    #################################
    #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    MeasurementWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T2LabelMapVolume','DWIBrainMask','LabelsConfigFile',
                                                             'FAImage','MDImage','RDImage','FrobeniusNormImage',
                                                             'Lambda1Image','Lambda2Image','Lambda3Image']),
                         name='inputsSpec')
    inputsSpec.inputs.LabelsConfigFile = LABELS_CONFIG_FILE

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['FA_stats','MD_stats','RD_stats','FrobeniusNorm_stats',
                                                              'Lambda1_stats','Lambda2_stats','Lambda3_stats']),
                          name='outputsSpec')

    # Step1: Create the labelmap volume for DWI scan
    CreateDWILabelMapNode = pe.Node(interface=Function(function = CreateDWILabelMap,
                                                       input_names=['T2LabelMapVolume','DWIBrainMask'],
                                                       output_names=['DWILabelMapVolume']),
                                    name="CreateDWILabelMap")
    MeasurementWF.connect(inputsSpec, 'T2LabelMapVolume', CreateDWILabelMapNode, 'T2LabelMapVolume')
    MeasurementWF.connect(inputsSpec, 'DWIBrainMask', CreateDWILabelMapNode, 'DWIBrainMask')

    # Now we have two labelmap volumes (both have 1x1x1 voxel lattice):
    # (1) T2LabelMap: Used to compute total_volume for each label
    # (2) DWILabelMap: It is probably cropped and missed some part of labels,
    #                  and is used to compute all stats like [mean,std,max,min,median,effective_volume].

    # Step2: Resample each RIS to T2LabelmapVolume voxel lattice
    MakeResamplerInFilesListNode = pe.Node(Function(function=MakeResamplerInFileList,
                                                    input_names=['FAImage','MDImage','RDImage','FrobeniusNormImage',
                                                                 'Lambda1Image','Lambda2Image','Lambda3Image'],
                                                    output_names=['RISsList']),
                                           name="MakeResamplerInFilesListNode")
    MeasurementWF.connect([(inputsSpec, MakeResamplerInFilesListNode, [('FAImage','FAImage'),
                                                                       ('MDImage','MDImage'),
                                                                       ('RDImage','RDImage'),
                                                                       ('FrobeniusNormImage','FrobeniusNormImage'),
                                                                       ('Lambda1Image','Lambda1Image'),
                                                                       ('Lambda2Image','Lambda2Image'),
                                                                       ('Lambda3Image','Lambda3Image')])])


    ResampleRISsNode = pe.MapNode(interface=BRAINSResample(), name="ResampleRISs",
                                  iterfield=['inputVolume', 'outputVolume'])
    ResampleRISsNode.inputs.interpolationMode = 'Linear'
    ResampleRISsNode.inputs.pixelType = 'float'
    ResampleRISsNode.inputs.outputVolume = ['FA_res.nrrd','MD_res.nrrd','RD_res.nrrd','frobenius_norm_res.nrrd',
                                            'lambda1_res.nrrd','lambda2_res.nrrd','lambda3_res.nrrd']
    MeasurementWF.connect(inputsSpec,'T2LabelMapVolume',ResampleRISsNode,'referenceVolume')
    MeasurementWF.connect(MakeResamplerInFilesListNode,'RISsList',ResampleRISsNode,'inputVolume')

    # Step3: Computes statistics of each resampled RIS over all input labels
    # and writes the results as a CSV file (a csv file for each RIS)
    ComputeStatisticsNode = pe.MapNode(interface=Function(function = ComputeStatistics,
                                                          input_names=['inputVolume','T2LabelMapVolume','DWILabelMapVolume','labelCodesFile'],
                                                          output_names=['CSVStatisticsFile']),
                                       name="ComputeStatistics",
                                       iterfield=['inputVolume'])
    MeasurementWF.connect(ResampleRISsNode, 'outputVolume', ComputeStatisticsNode, 'inputVolume')
    MeasurementWF.connect(inputsSpec, 'T2LabelMapVolume', ComputeStatisticsNode, 'T2LabelMapVolume')
    MeasurementWF.connect(CreateDWILabelMapNode, 'DWILabelMapVolume', ComputeStatisticsNode, 'DWILabelMapVolume')
    MeasurementWF.connect(inputsSpec, 'LabelsConfigFile', ComputeStatisticsNode, 'labelCodesFile')
    MeasurementWF.connect(ComputeStatisticsNode, ('CSVStatisticsFile', pickFromList, 0), outputsSpec, 'FA_stats')
    MeasurementWF.connect(ComputeStatisticsNode, ('CSVStatisticsFile', pickFromList, 1), outputsSpec, 'MD_stats')
    MeasurementWF.connect(ComputeStatisticsNode, ('CSVStatisticsFile', pickFromList, 2), outputsSpec, 'RD_stats')
    MeasurementWF.connect(ComputeStatisticsNode, ('CSVStatisticsFile', pickFromList, 3), outputsSpec, 'FrobeniusNorm_stats')
    MeasurementWF.connect(ComputeStatisticsNode, ('CSVStatisticsFile', pickFromList, 4), outputsSpec, 'Lambda1_stats')
    MeasurementWF.connect(ComputeStatisticsNode, ('CSVStatisticsFile', pickFromList, 5), outputsSpec, 'Lambda2_stats')
    MeasurementWF.connect(ComputeStatisticsNode, ('CSVStatisticsFile', pickFromList, 6), outputsSpec, 'Lambda3_stats')

    return MeasurementWF
