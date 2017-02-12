#! /usr/bin/env python
"""
DWIWorkFlow.py
==============

The purpose of this pipeline is to complete all the pre-processing steps needed to turn diffusion-weighted images into FA images that will be used to build a template diffusion tensor atlas for fiber tracking.

Usage:
  DWIWorkFlow.py --inputDWIScan DWISCAN --inputT2Scan T2SCAN --inputBrainLabelsMapImage BLMImage --program_paths PROGRAM_PATHS --python_aux_paths PYTHON_AUX_PATHS [--workflowCacheDir CACHEDIR] [--resultDir RESULTDIR]
  DWIWorkFlow.py -v | --version
  DWIWorkFlow.py -h | --help

Options:
  -h --help                                 Show this help and exit
  -v --version                              Print the version and exit
  --inputDWIScan DWISCAN                    Path to the input DWI scan for further processing
  --inputT2Scan T2SCAN                      Path to the input T2 scan
  --inputBrainLabelsMapImage BLMImage       Path to the input brain labels map image
  --program_paths PROGRAM_PATHS             Path to the directory where binary files are places
  --python_aux_paths PYTHON_AUX_PATHS       Path to the AutoWorkup directory
  --workflowCacheDir CACHEDIR               Base directory that cache outputs of workflow will be written to (default: ./)
  --resultDir RESULTDIR                     Outputs of dataSink will be written to a sub directory under the resultDir named by input scan sessionID (default: CACHEDIR)
"""
from __future__ import print_function

#############################  UTILITY FUNCTIONS  #####################################
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# remove the skull from the T2 volume
def ExtractBRAINFromHead(RawScan, BrainLabels):
    import os
    import SimpleITK as sitk
    # Remove skull from the head scan
    assert os.path.exists(RawScan), "File not found: %s" % RawScan
    assert os.path.exists(BrainLabels), "File not found: %s" % BrainLabels
    headImage = sitk.ReadImage(RawScan.encode('ascii','replace'))
    labelsMap = sitk.ReadImage(BrainLabels.encode('ascii','replace'))
    rs = sitk.ResampleImageFilter()
    rs.SetInterpolator(sitk.sitkLinear)
    rs.SetTransform(sitk.Transform(3, sitk.sitkIdentity))
    rs.SetReferenceImage(labelsMap)
    resampledHead = rs.Execute(headImage)

    label_mask = labelsMap>0
    brainImage = sitk.Cast(resampledHead,sitk.sitkInt16) * sitk.Cast(label_mask,sitk.sitkInt16)
    outputVolume = os.path.realpath('T2Stripped.nrrd')
    sitk.WriteImage(brainImage, outputVolume.encode('ascii','replace'))
    return outputVolume

def MakeResamplerInFileList(inputT2, inputLabelMap):
    imagesList = [inputT2, inputLabelMap]
    return imagesList

# This function helps to pick desirable output from the outputVolume list
def pickFromList(inlist,item):
    return inlist[item]

# Create registration mask for ANTs from resampled label map image
def CreateAntsRegistrationMask(brainMask):
    import os
    import SimpleITK as sitk
    assert os.path.exists(brainMask), "File not found: %s" % brainMask
    labelsMap = sitk.ReadImage(brainMask.encode('ascii','replace'))
    label_mask = labelsMap>0
    # dilate the label mask
    dilateFilter = sitk.BinaryDilateImageFilter()
    dilateFilter.SetKernelRadius(12)
    dilated_mask = dilateFilter.Execute( label_mask )
    regMask = dilated_mask
    registrationMask = os.path.realpath('registrationMask.nrrd')
    sitk.WriteImage(regMask, registrationMask.encode('ascii','replace'))
    return registrationMask

# Save direction cosine for the input volume
def SaveDirectionCosineToMatrix(inputVolume):
    import os
    import SimpleITK as sitk
    assert os.path.exists(inputVolume), "File not found: %s" % inputVolume
    t2 = sitk.ReadImage(inputVolume.encode('ascii','replace'))
    directionCosine = t2.GetDirection()
    return directionCosine

def MakeForceDCFilesList(inputB0, inputT2, inputLabelMap):
    import os
    assert os.path.exists(inputB0), "File not found: %s" % inputB0
    assert os.path.exists(inputT2), "File not found: %s" % inputT2
    assert os.path.exists(inputLabelMap), "File not found: %s" % inputLabelMap
    imagesList = [inputB0, inputT2, inputLabelMap]
    return imagesList

# Force DC to ID
def ForceDCtoID(inputVolume):
    import os
    import SimpleITK as sitk
    inImage = sitk.ReadImage(inputVolume.encode('ascii','replace'))
    inImage.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    outputVolume = os.path.realpath('IDDC_'+ os.path.basename(inputVolume))
    sitk.WriteImage(inImage, outputVolume.encode('ascii','replace'))
    return outputVolume

def pickCompositeTransfromFromList(composite_transform_as_list):
    if isinstance(composite_transform_as_list, basestring):
        return composite_transform_as_list;
    return composite_transform_as_list[0]

def RestoreDCFromSavedMatrix(inputVolume, inputDirectionCosine):
    import os
    import SimpleITK as sitk
    inImage = sitk.ReadImage(inputVolume.encode('ascii','replace'))
    inImage.SetDirection(inputDirectionCosine)
    outputVolume = os.path.realpath('CorrectedDWI.nrrd')
    sitk.WriteImage(inImage, outputVolume.encode('ascii','replace'))
    return outputVolume

def GetRigidTransformInverse(inputTransform):
    import os
    import SimpleITK as sitk
    inputTx = sitk.ReadTransform(inputTransform)
    versorRigidTx = sitk.VersorRigid3DTransform()
    versorRigidTx.SetFixedParameters(inputTx.GetFixedParameters())
    versorRigidTx.SetParameters(inputTx.GetParameters())
    invTx = versorRigidTx.GetInverse()
    inverseTransform = os.path.realpath('Inverse_'+ os.path.basename(inputTransform))
    sitk.WriteTransform(invTx, inverseTransform)
    return inverseTransform
#######################################################################################
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

def runMainWorkflow(DWI_scan, T2_scan, labelMap_image, BASE_DIR, dataSink_DIR):
    print("Running the workflow ...")

    sessionID = os.path.basename(os.path.dirname(DWI_scan))
    subjectID = os.path.basename(os.path.dirname(os.path.dirname(DWI_scan)))
    siteID = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(DWI_scan))))

    #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
    ####### Workflow ###################
    WFname = 'DWIWorkflow_CACHE_' + sessionID
    DWIWorkflow = pe.Workflow(name=WFname)
    DWIWorkflow.base_dir = BASE_DIR

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T2Volume', 'DWIVolume','LabelMapVolume']),
                         name='inputspec')

    inputsSpec.inputs.DWIVolume = DWI_scan
    inputsSpec.inputs.T2Volume = T2_scan
    inputsSpec.inputs.LabelMapVolume = labelMap_image

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['CorrectedDWI','CorrectedDWI_in_T2Space','tensor_image','DWIBrainMask',
                                                              'FAImage','MDImage','RDImage','FrobeniusNormImage',
                                                              'Lambda1Image','Lambda2Image','Lambda3Image',
                                                              'ukfTracks','ukf2ndTracks' ]),
                          name='outputsSpec')

    # Step0: remove the skull from the T2 volume
    ExtractBRAINFromHeadNode = pe.Node(interface=Function(function = ExtractBRAINFromHead,
                                                          input_names=['RawScan','BrainLabels'],
                                                          output_names=['outputVolume']),
                                       name="ExtractBRAINFromHead")

    DWIWorkflow.connect(inputsSpec, 'T2Volume', ExtractBRAINFromHeadNode, 'RawScan')
    DWIWorkflow.connect(inputsSpec, 'LabelMapVolume', ExtractBRAINFromHeadNode, 'BrainLabels')

    # Step1: extract B0 from DWI volume
    EXTRACT_B0 = pe.Node(interface=extractNrrdVectorIndex(),name="EXTRACT_B0")
    EXTRACT_B0.inputs.vectorIndex = 0
    EXTRACT_B0.inputs.outputVolume = 'B0_Image.nrrd'
    DWIWorkflow.connect(inputsSpec,'DWIVolume',EXTRACT_B0,'inputVolume')

    # Step2: Register T2 to B0 space using BRAINSFit
    BFit_T2toB0 = pe.Node(interface=BRAINSFit(), name="BFit_T2toB0")
    BFit_T2toB0.inputs.costMetric = "MMI"
    BFit_T2toB0.inputs.numberOfSamples = 100000
    BFit_T2toB0.inputs.numberOfIterations = [1500]
    BFit_T2toB0.inputs.numberOfHistogramBins = 50
    BFit_T2toB0.inputs.maximumStepLength = 0.2
    BFit_T2toB0.inputs.minimumStepLength = [0.00005]
    BFit_T2toB0.inputs.useRigid = True
    BFit_T2toB0.inputs.useAffine = True
    BFit_T2toB0.inputs.maskInferiorCutOffFromCenter = 65
    BFit_T2toB0.inputs.maskProcessingMode = "ROIAUTO"
    BFit_T2toB0.inputs.ROIAutoDilateSize = 13
    BFit_T2toB0.inputs.backgroundFillValue = 0.0
    BFit_T2toB0.inputs.initializeTransformMode = 'useCenterOfHeadAlign'
    BFit_T2toB0.inputs.strippedOutputTransform = "T2ToB0_RigidTransform.h5"
    BFit_T2toB0.inputs.writeOutputTransformInFloat = True
    DWIWorkflow.connect(EXTRACT_B0, 'outputVolume', BFit_T2toB0, 'fixedVolume')
    DWIWorkflow.connect(ExtractBRAINFromHeadNode, 'outputVolume', BFit_T2toB0, 'movingVolume')

    # Step3: Use T_rigid to "resample" T2 and label map images to B0 image space
    MakeResamplerInFilesListNode = pe.Node(Function(function=MakeResamplerInFileList,
                                                    input_names=['inputT2','inputLabelMap'],
                                                    output_names=['imagesList']),
                                           name="MakeResamplerInFilesListNode")

    DWIWorkflow.connect([(ExtractBRAINFromHeadNode,MakeResamplerInFilesListNode,[('outputVolume','inputT2')]),
                         (inputsSpec,MakeResamplerInFilesListNode,[('LabelMapVolume','inputLabelMap')])])

    ResampleToB0Space = pe.MapNode(interface=BRAINSResample(), name="ResampleToB0Space",
                                   iterfield=['inputVolume', 'pixelType', 'outputVolume'])
    ResampleToB0Space.inputs.interpolationMode = 'Linear'
    ResampleToB0Space.inputs.outputVolume = ['T2toB0.nrrd','BRAINMaskToB0.nrrd']
    ResampleToB0Space.inputs.pixelType = ['ushort','binary']
    DWIWorkflow.connect(BFit_T2toB0,'strippedOutputTransform',ResampleToB0Space,'warpTransform')
    DWIWorkflow.connect(EXTRACT_B0,'outputVolume',ResampleToB0Space,'referenceVolume')
    DWIWorkflow.connect(MakeResamplerInFilesListNode,'imagesList',ResampleToB0Space,'inputVolume')

    # Step4: Create registration mask from resampled label map image
    CreateRegistrationMask = pe.Node(interface=Function(function = CreateAntsRegistrationMask,
                                                        input_names=['brainMask'],
                                                        output_names=['registrationMask']),
                                     name="CreateAntsRegistrationMask")

    DWIWorkflow.connect(ResampleToB0Space, ('outputVolume', pickFromList, 1),
                        CreateRegistrationMask, 'brainMask')

    # Step5: Save direction cosine for the resampled T2 image
    SaveDirectionCosineToMatrixNode = pe.Node(interface=Function(function = SaveDirectionCosineToMatrix,
                                                                 input_names=['inputVolume'],
                                                                 output_names=['directionCosine']),
                                              name="SaveDirectionCosineToMatrix")

    DWIWorkflow.connect(ResampleToB0Space, ('outputVolume', pickFromList, 0),
                        SaveDirectionCosineToMatrixNode, 'inputVolume')

    # Step6: Force DC to ID
    MakeForceDCFilesListNode = pe.Node(Function(function=MakeForceDCFilesList,
                                                input_names=['inputB0','inputT2','inputLabelMap'],
                                                output_names=['imagesList']),
                                       name="MakeForceDCFilesListNode")

    DWIWorkflow.connect([(EXTRACT_B0,MakeForceDCFilesListNode,[('outputVolume','inputB0')]),
                         (ResampleToB0Space,MakeForceDCFilesListNode,[(('outputVolume', pickFromList, 0),'inputT2')]),
                         (CreateRegistrationMask,MakeForceDCFilesListNode,[('registrationMask','inputLabelMap')])])

    ForceDCtoIDNode = pe.MapNode(interface=Function(function = ForceDCtoID,
                                                    input_names=['inputVolume'],
                                                    output_names=['outputVolume']),
                                 name="ForceDCtoID",
                                 iterfield=['inputVolume'])

    DWIWorkflow.connect(MakeForceDCFilesListNode, 'imagesList', ForceDCtoIDNode, 'inputVolume')

    # Step7: Run antsRegistration in one direction
    antsReg_B0ToTransformedT2 = pe.Node(interface=ants.Registration(), name="antsReg_B0ToTransformedT2")
    antsReg_B0ToTransformedT2.inputs.interpolation = "Linear"
    antsReg_B0ToTransformedT2.inputs.dimension = 3
    antsReg_B0ToTransformedT2.inputs.transforms = ["SyN"]
    antsReg_B0ToTransformedT2.inputs.transform_parameters = [(0.25, 3.0, 0.0)]
    antsReg_B0ToTransformedT2.inputs.metric = ['MI']
    antsReg_B0ToTransformedT2.inputs.sampling_strategy = [None]
    antsReg_B0ToTransformedT2.inputs.sampling_percentage = [1.0]
    antsReg_B0ToTransformedT2.inputs.metric_weight = [1.0]
    antsReg_B0ToTransformedT2.inputs.radius_or_number_of_bins = [32]
    antsReg_B0ToTransformedT2.inputs.number_of_iterations = [[70, 50, 40]]
    antsReg_B0ToTransformedT2.inputs.convergence_threshold = [1e-6]
    antsReg_B0ToTransformedT2.inputs.convergence_window_size = [10]
    antsReg_B0ToTransformedT2.inputs.use_histogram_matching = [True]
    antsReg_B0ToTransformedT2.inputs.shrink_factors = [[3, 2, 1]]
    antsReg_B0ToTransformedT2.inputs.smoothing_sigmas = [[2, 1, 0]]
    antsReg_B0ToTransformedT2.inputs.sigma_units = ["vox"]
    antsReg_B0ToTransformedT2.inputs.use_estimate_learning_rate_once = [False]
    antsReg_B0ToTransformedT2.inputs.write_composite_transform = True
    antsReg_B0ToTransformedT2.inputs.collapse_output_transforms = False
    antsReg_B0ToTransformedT2.inputs.initialize_transforms_per_stage = False
    antsReg_B0ToTransformedT2.inputs.output_transform_prefix = 'Tsyn'
    antsReg_B0ToTransformedT2.inputs.winsorize_lower_quantile = 0.01
    antsReg_B0ToTransformedT2.inputs.winsorize_upper_quantile = 0.99
    antsReg_B0ToTransformedT2.inputs.float = True
    antsReg_B0ToTransformedT2.inputs.num_threads = -1
    antsReg_B0ToTransformedT2.inputs.args = '--restrict-deformation 0x1x0'

    DWIWorkflow.connect(ForceDCtoIDNode, ('outputVolume', pickFromList, 1), antsReg_B0ToTransformedT2, 'fixed_image')
    DWIWorkflow.connect(ForceDCtoIDNode, ('outputVolume', pickFromList, 2), antsReg_B0ToTransformedT2, 'fixed_image_mask')
    DWIWorkflow.connect(ForceDCtoIDNode, ('outputVolume', pickFromList, 0), antsReg_B0ToTransformedT2, 'moving_image')

    # Step8: Now, all necessary transforms are acquired. It's a time to
    #        transform input DWI image into T2 image space
    # {DWI} --> ForceDCtoID --> gtractResampleDWIInPlace(using SyN transfrom)
    # --> Restore DirectionCosine From Saved Matrix --> gtractResampleDWIInPlace(inverse of T_rigid from BFit)
    # --> {CorrectedDW_in_T2Space}

    DWI_ForceDCtoIDNode = pe.Node(interface=Function(function = ForceDCtoID,
                                                     input_names=['inputVolume'],
                                                     output_names=['outputVolume']),
                                  name='DWI_ForceDCtoIDNode')

    DWIWorkflow.connect(inputsSpec,'DWIVolume',DWI_ForceDCtoIDNode,'inputVolume')

    gtractResampleDWI_SyN = pe.Node(interface=gtractResampleDWIInPlace(),
                                    name="gtractResampleDWI_SyN")

    DWIWorkflow.connect(DWI_ForceDCtoIDNode,'outputVolume',gtractResampleDWI_SyN,'inputVolume')
    DWIWorkflow.connect(antsReg_B0ToTransformedT2,('composite_transform',pickCompositeTransfromFromList),gtractResampleDWI_SyN,'warpDWITransform')
    DWIWorkflow.connect(ForceDCtoIDNode,('outputVolume', pickFromList, 1),gtractResampleDWI_SyN,'referenceVolume') # fixed image of antsRegistration
    gtractResampleDWI_SyN.inputs.outputVolume = 'IDDC_correctedDWI.nrrd'

    RestoreDCFromSavedMatrixNode = pe.Node(interface=Function(function = RestoreDCFromSavedMatrix,
                                                              input_names=['inputVolume','inputDirectionCosine'],
                                                              output_names=['outputVolume']),
                                           name='RestoreDCFromSavedMatrix')

    DWIWorkflow.connect(gtractResampleDWI_SyN,'outputVolume',RestoreDCFromSavedMatrixNode,'inputVolume')
    DWIWorkflow.connect(SaveDirectionCosineToMatrixNode,'directionCosine',RestoreDCFromSavedMatrixNode,'inputDirectionCosine')
    DWIWorkflow.connect(RestoreDCFromSavedMatrixNode,'outputVolume', outputsSpec, 'CorrectedDWI')

    GetRigidTransformInverseNode = pe.Node(interface=Function(function = GetRigidTransformInverse,
                                                              input_names=['inputTransform'],
                                                              output_names=['inverseTransform']),
                                           name='GetRigidTransformInverse')

    DWIWorkflow.connect(BFit_T2toB0,'strippedOutputTransform',GetRigidTransformInverseNode,'inputTransform')

    gtractResampleDWIInPlace_Trigid = pe.Node(interface=gtractResampleDWIInPlace(),
                                              name="gtractResampleDWIInPlace_Trigid")

    DWIWorkflow.connect(RestoreDCFromSavedMatrixNode,'outputVolume',gtractResampleDWIInPlace_Trigid,'inputVolume')
    DWIWorkflow.connect(GetRigidTransformInverseNode,'inverseTransform',gtractResampleDWIInPlace_Trigid,'inputTransform') #Inverse of rigid transform from BFit
    gtractResampleDWIInPlace_Trigid.inputs.outputVolume = 'CorrectedDWI_in_T2Space_estimate.nrrd'
    gtractResampleDWIInPlace_Trigid.inputs.outputResampledB0 = 'CorrectedDWI_in_T2Space_estimate_B0.nrrd'

    # Setp9: An extra registration step to tune the alignment between the CorrecetedDWI_in_T2Space image and T2 image.
    BFit_TuneRegistration = pe.Node(interface=BRAINSFit(), name="BFit_TuneRegistration")
    BFit_TuneRegistration.inputs.costMetric = "MMI"
    BFit_TuneRegistration.inputs.numberOfSamples = 100000
    BFit_TuneRegistration.inputs.numberOfIterations = [1500]
    BFit_TuneRegistration.inputs.numberOfHistogramBins = 50
    BFit_TuneRegistration.inputs.maximumStepLength = 0.2
    BFit_TuneRegistration.inputs.minimumStepLength = [0.00005]
    BFit_TuneRegistration.inputs.useRigid = True
    BFit_TuneRegistration.inputs.useAffine = True
    BFit_TuneRegistration.inputs.maskInferiorCutOffFromCenter = 65
    BFit_TuneRegistration.inputs.maskProcessingMode = "ROIAUTO"
    BFit_TuneRegistration.inputs.ROIAutoDilateSize = 13
    BFit_TuneRegistration.inputs.backgroundFillValue = 0.0
    BFit_TuneRegistration.inputs.initializeTransformMode = 'useCenterOfHeadAlign'
    BFit_TuneRegistration.inputs.strippedOutputTransform = "CorrectedB0inT2Space_to_T2_RigidTransform.h5"
    BFit_TuneRegistration.inputs.writeOutputTransformInFloat = True
    DWIWorkflow.connect(ExtractBRAINFromHeadNode, 'outputVolume', BFit_TuneRegistration, 'fixedVolume') #T2 brain volume
    DWIWorkflow.connect(gtractResampleDWIInPlace_Trigid, 'outputResampledB0', BFit_TuneRegistration, 'movingVolume') # CorrectedB0_in_T2Space

    gtractResampleDWIInPlace_TuneRigidTx = pe.Node(interface=gtractResampleDWIInPlace(),
                                                   name="gtractResampleDWIInPlace_TuneRigidTx")
    DWIWorkflow.connect(gtractResampleDWIInPlace_Trigid,'outputVolume',gtractResampleDWIInPlace_TuneRigidTx,'inputVolume')
    DWIWorkflow.connect(BFit_TuneRegistration,'strippedOutputTransform',gtractResampleDWIInPlace_TuneRigidTx,'inputTransform')
    gtractResampleDWIInPlace_TuneRigidTx.inputs.outputVolume = 'CorrectedDWI_in_T2Space.nrrd'
    gtractResampleDWIInPlace_TuneRigidTx.inputs.outputResampledB0 = 'CorrectedDWI_in_T2Space_B0.nrrd'

    # Finally we pass the outputs of the gtractResampleDWIInPlace_TuneRigidTx to the outputsSpec
    DWIWorkflow.connect(gtractResampleDWIInPlace_TuneRigidTx, 'outputVolume', outputsSpec, 'CorrectedDWI_in_T2Space')

    # Step10: Create brain mask from the input labelmap
    DWIBRAINMASK = pe.Node(interface=BRAINSResample(), name='DWIBRAINMASK')
    DWIBRAINMASK.inputs.interpolationMode = 'Linear'
    DWIBRAINMASK.inputs.outputVolume = 'BrainMaskForDWI.nrrd'
    DWIBRAINMASK.inputs.pixelType = 'binary'
    DWIWorkflow.connect(gtractResampleDWIInPlace_TuneRigidTx,'outputResampledB0',DWIBRAINMASK,'referenceVolume')
    DWIWorkflow.connect(inputsSpec,'LabelMapVolume',DWIBRAINMASK,'inputVolume')
    DWIWorkflow.connect(DWIBRAINMASK, 'outputVolume', outputsSpec, 'DWIBrainMask')

    # Step11: DTI estimation
    DTIEstim = pe.Node(interface=dtiestim(), name="DTIEstim")
    DTIEstim.inputs.method = 'wls'
    DTIEstim.inputs.tensor_output = 'DTI_Output.nrrd'
    DTIEstim.inputs.threshold = 0
    DWIWorkflow.connect(gtractResampleDWIInPlace_TuneRigidTx, 'outputVolume', DTIEstim, 'dwi_image')
    DWIWorkflow.connect(DWIBRAINMASK, 'outputVolume', DTIEstim, 'brain_mask')
    DWIWorkflow.connect(DTIEstim, 'tensor_output', outputsSpec, 'tensor_image')

    # Step12: DTI process
    DTIProcess = pe.Node(interface=dtiprocess(), name='DTIProcess')
    DTIProcess.inputs.fa_output = 'FA.nrrd'
    DTIProcess.inputs.md_output = 'MD.nrrd'
    DTIProcess.inputs.RD_output = 'RD.nrrd'
    DTIProcess.inputs.frobenius_norm_output = 'frobenius_norm_output.nrrd'
    DTIProcess.inputs.lambda1_output = 'lambda1_output.nrrd'
    DTIProcess.inputs.lambda2_output = 'lambda2_output.nrrd'
    DTIProcess.inputs.lambda3_output = 'lambda3_output.nrrd'
    DTIProcess.inputs.scalar_float = True

    DWIWorkflow.connect(DTIEstim, 'tensor_output', DTIProcess, 'dti_image')
    DWIWorkflow.connect(DTIProcess, 'fa_output', outputsSpec, 'FAImage')
    DWIWorkflow.connect(DTIProcess, 'md_output', outputsSpec, 'MDImage')
    DWIWorkflow.connect(DTIProcess, 'RD_output', outputsSpec, 'RDImage')
    DWIWorkflow.connect(DTIProcess, 'frobenius_norm_output', outputsSpec, 'FrobeniusNormImage')
    DWIWorkflow.connect(DTIProcess, 'lambda1_output', outputsSpec, 'Lambda1Image')
    DWIWorkflow.connect(DTIProcess, 'lambda2_output', outputsSpec, 'Lambda2Image')
    DWIWorkflow.connect(DTIProcess, 'lambda3_output', outputsSpec, 'Lambda3Image')
    """
    # Step13: UKF Processing
    UKFNode = pe.Node(interface=UKFTractography(), name= "UKFRunRecordStates")
    UKFNode.inputs.tracts = "ukfTracts.vtk"
    #UKFNode.inputs.tractsWithSecondTensor = "ukfSecondTensorTracks.vtk"
    UKFNode.inputs.numTensor = '2'
    UKFNode.inputs.freeWater = True ## default False
    UKFNode.inputs.recordFA = True ## default False
    UKFNode.inputs.recordTensors = True ## default False
    #UKFNode.inputs.recordCovariance = True ## default False
    #UKFNode.inputs.recordState = True ## default False
    #UKFNode.inputs.recordFreeWater = True ## default False
    #UKFNode.inputs.recordTrace = True ## default False
    #UKFNode.inputs.recordNMSE = True ## default False

    DWIWorkflow.connect(gtractResampleDWIInPlace_TuneRigidTx, 'outputVolume', UKFNode, 'dwiFile')
    DWIWorkflow.connect(DWIBRAINMASK, 'outputVolume', UKFNode, 'maskFile')
    DWIWorkflow.connect(UKFNode,'tracts',outputsSpec,'ukfTracks')
    #DWIWorkflow.connect(UKFNode,'tractsWithSecondTensor',outputsSpec,'ukf2ndTracks')
    """

    ## Write all outputs with DataSink
    DWIDataSink = pe.Node(interface=nio.DataSink(), name='DWIDataSink')
    DWIDataSink.inputs.base_directory = dataSink_DIR
    DWIDataSink.inputs.container = sessionID

    DWIWorkflow.connect(outputsSpec, 'ukfTracks', DWIDataSink, 'Outputs.@ukfTracks')
    #DWIWorkflow.connect(outputsSpec, 'ukf2ndTracks', DWIDataSink, 'Outputs.@ukf2ndTracks')
    DWIWorkflow.connect(outputsSpec, 'CorrectedDWI', DWIDataSink, 'Outputs.@CorrectedDWI')
    DWIWorkflow.connect(outputsSpec, 'CorrectedDWI_in_T2Space', DWIDataSink, 'Outputs.@CorrectedDWI_in_T2Space')
    DWIWorkflow.connect(outputsSpec, 'tensor_image', DWIDataSink, 'Outputs.@tensor_image')
    DWIWorkflow.connect(outputsSpec, 'DWIBrainMask', DWIDataSink, 'Outputs.@DWIBrainMask')
    DWIWorkflow.connect(outputsSpec, 'FAImage', DWIDataSink, 'Outputs.@FAImage')
    DWIWorkflow.connect(outputsSpec, 'MDImage', DWIDataSink, 'Outputs.@MDImage')
    DWIWorkflow.connect(outputsSpec, 'RDImage', DWIDataSink, 'Outputs.@RDImage')
    DWIWorkflow.connect(outputsSpec, 'FrobeniusNormImage', DWIDataSink, 'Outputs.@FrobeniusNormImage')
    DWIWorkflow.connect(outputsSpec, 'Lambda1Image', DWIDataSink, 'Outputs.@Lambda1Image')
    DWIWorkflow.connect(outputsSpec, 'Lambda2Image', DWIDataSink, 'Outputs.@Lambda2Image')
    DWIWorkflow.connect(outputsSpec, 'Lambda3Image', DWIDataSink, 'Outputs.@Lambda3Image')

    DWIWorkflow.write_graph()
    DWIWorkflow.run()


if __name__ == '__main__':
  import os
  import glob
  import sys

  from docopt import docopt
  argv = docopt(__doc__, version='1.0')
  print(argv)

  DWISCAN = argv['--inputDWIScan']
  assert os.path.exists(DWISCAN), "Input DWI scan is not found: %s" % DWISCAN

  T2SCAN = argv['--inputT2Scan']
  assert os.path.exists(T2SCAN), "Input T2 scan is not found: %s" % T2SCAN

  LabelMapImage = argv['--inputBrainLabelsMapImage']
  assert os.path.exists(LabelMapImage), "Input Brain labels map image is not found: %s" % LabelMapImage

  PROGRAM_PATHS = argv['--program_paths']

  PYTHON_AUX_PATHS = argv['--python_aux_paths']

  if argv['--workflowCacheDir'] == None:
      print("*** workflow cache directory is set to current working directory.")
      CACHEDIR = os.getcwd()
  else:
      CACHEDIR = argv['--workflowCacheDir']
      assert os.path.exists(CACHEDIR), "Cache directory is not found: %s" % CACHEDIR

  if argv['--resultDir'] == None:
      print("*** data sink result directory is set to the cache directory.")
      RESULTDIR = CACHEDIR
  else:
      RESULTDIR = argv['--resultDir']
      assert os.path.exists(RESULTDIR), "Results directory is not found: %s" % RESULTDIR

  print('=' * 100)

  #\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
  #####################################################################################
  #     Prepend the shell environment search paths
  PROGRAM_PATHS = PROGRAM_PATHS.split(':')
  PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))
  os.environ['PATH'] = ':'.join(PROGRAM_PATHS)

  CUSTOM_ENVIRONMENT=dict()

  # Platform specific information
  #     Prepend the python search paths
  PYTHON_AUX_PATHS = PYTHON_AUX_PATHS.split(':')
  PYTHON_AUX_PATHS.extend(sys.path)
  sys.path = PYTHON_AUX_PATHS

  import SimpleITK as sitk
  import nipype
  from nipype.interfaces import ants
  from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
  from nipype.interfaces.base import traits, isdefined, BaseInterface
  from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
  import nipype.interfaces.io as nio   # Data i/oS
  import nipype.pipeline.engine as pe  # pypeline engine
  from nipype.interfaces.freesurfer import ReconAll
  from nipype.interfaces.semtools import *
  #####################################################################################
  exit = runMainWorkflow(DWISCAN, T2SCAN, LabelMapImage, CACHEDIR, RESULTDIR)

  sys.exit(exit)
