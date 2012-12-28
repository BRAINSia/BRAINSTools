

## start {
def WFPerSubjectDef ( inputListOfSubjectVolumes, 
                      inputTemplateDir,
                      outputCacheDir,
                      outputResultsDir):

    #
    # workflow definition
    #
    import nipype.pipeline.engine as pe
    
    WFPerSubject = pe.Workflow(name="subject")
    WFPerSubject.base_dir = outputCacheDir
    
    import nipype.interfaces.io as nio
    datasink = pe.Node( nio.DataSink(), name = 'sinker' )
    datasink.inputs.base_directory = ( outputResultsDir ) 
    
    datasink.inputs.regexp_substitutions = [ ('_inputVolume.*BABC..', '') ,
                                             ('_Cache',''),
                                             ('.nii.gz', '') ]
    
    # --------------------------------------------------------------------------------------- #
    # 1. Deform atlas to subject with really gross transform
    #
    from BRAINSFit import BRAINSFit
    
    inputTemplateT1 = inputTemplateDir + "/template_t1.nii.gz"
    
    BFitAtlasToSubject = pe.Node( interface=BRAINSFit(),
                                  name="01_AtlasToSubjectRegistration")
    
    BFitAtlasToSubject.inputs.fixedVolume           = inputListOfSubjectVolumes['t1']
    BFitAtlasToSubject.inputs.movingVolume          = inputTemplateT1
    BFitAtlasToSubject.inputs.outputVolume          = "atlas_to_subject_warped.nii.gz"  
    BFitAtlasToSubject.inputs.initialTransform      = inputListOfSubjectVolumes['transform']
    BFitAtlasToSubject.inputs.outputTransform       = "atlas_to_subject.h5"  
    
    ## fixed parameter specification
    BFitAtlasToSubject.inputs.costMetric = "MMI"
    BFitAtlasToSubject.inputs.maskProcessingMode = "ROIAUTO"
    BFitAtlasToSubject.inputs.numberOfSamples = 100000
    BFitAtlasToSubject.inputs.numberOfIterations = [1500,1500,1500,1500]
    BFitAtlasToSubject.inputs.numberOfHistogramBins = 50
    BFitAtlasToSubject.inputs.maximumStepLength = 0.2
    BFitAtlasToSubject.inputs.minimumStepLength = [0.005,0.005,0.005,0.005]
    BFitAtlasToSubject.inputs.transformType = "ScaleVersor3D,ScaleSkewVersor3D,Affine,BSpline" 
    BFitAtlasToSubject.inputs.relaxationFactor = 0.5  
    BFitAtlasToSubject.inputs.translationScale = 1000
    BFitAtlasToSubject.inputs.reproportionScale = 1  
    BFitAtlasToSubject.inputs.skewScale = 1
    BFitAtlasToSubject.inputs.useExplicitPDFDerivativesMode = "AUTO"  
    BFitAtlasToSubject.inputs.useCachingOfBSplineWeightsMode = "ON"  
    BFitAtlasToSubject.inputs.maxBSplineDisplacement = 7 
    BFitAtlasToSubject.inputs.projectedGradientTolerance = 1e-05 
    BFitAtlasToSubject.inputs.costFunctionConvergenceFactor = 1e+09 
    BFitAtlasToSubject.inputs.backgroundFillValue = 0 
    BFitAtlasToSubject.inputs.maskInferiorCutOffFromCenter = 65  
    BFitAtlasToSubject.inputs.splineGridSize = [56,40,48] 
    
    ## add to workflow
    WFPerSubject.add_nodes( [BFitAtlasToSubject])
    
    WFPerSubject.connect( BFitAtlasToSubject, 'outputTransform',
                          datasink, '01_BFitAtlasToSubject')
    
    
    # --------------------------------------------------------------------------------------- #
    # 2. Warp Probability Maps
    #
    
    rois     = [ "l_accumben" , "l_caudate" ,"l_putamen" ,"l_globus" ,"l_thalamus" ,"l_hippocampus",
                 "r_accumben" , "r_caudate" ,"r_putamen" ,"r_globus" ,"r_thalamus" ,"r_hippocampus" ]
    
    import MyUtilities  
    from nipype.interfaces.utility import Function
    GetProbabilityMapFilename = pe.Node( name = "utilGetProbabilityMapFilename", 
                           interface = Function( input_names = ["inputROI", "inputDir"],
                                                 output_names = ["outputFilename"],
                                                 function = MyUtilities.GetProbabilityFilename )
                          )
    GetProbabilityMapFilename.iterables = ( "inputROI", rois )                        
    GetProbabilityMapFilename.inputs.inputDir  = inputTemplateDir
    
    from BRAINSResample import BRAINSResample 
    WarpProbabilityMap = pe.Node( interface=BRAINSResample(), 
                                  name = "02_WarpProabilityMap",
                                )
    
    WarpProbabilityMap.inputs.outputVolume    = "outputTemplateWarpedToSubject.nii.gz"
    WarpProbabilityMap.inputs.referenceVolume = inputListOfSubjectVolumes['t1']
    
    ## add to the workflow and connect
    WFPerSubject.add_nodes([ WarpProbabilityMap])
    WFPerSubject.connect( GetProbabilityMapFilename , 'outputFilename',
                          WarpProbabilityMap,'inputVolume' )
    WFPerSubject.connect( BFitAtlasToSubject, 'outputTransform',
                          WarpProbabilityMap,'warpTransform' )
    WFPerSubject.connect( WarpProbabilityMap, 'outputVolume',
                          datasink, '02_WarpedProbability' )

    # --------------------------------------------------------------------------------------- #
    # Smoothing 
    #

    # dummy node for sigmal iteration
    from nipype.interfaces.utility import IdentityInterface
    sigmaInputNode = pe.Node( interface = IdentityInterface( fields = ['sigma'] ),
                              name="sigmaInput" )
    sigmaInputNode.iterables = ( 'sigma' ,[0,0.5,1,1.5,2])
    WFPerSubject.add_nodes([ sigmaInputNode] )


    SmoothingProbMap = pe.Node( name = "SmoothingProbMap",
                                interface = Function( input_names = ["inputVolume",
                                                                     "inputSigma",
                                                                     "outputVolume"],
                                                      output_names = ["outputSmoothProbabilityMap"],
                                                      function = MyUtilities.SmoothProbabilityMap )
                              )
    WFPerSubject.add_nodes([ SmoothingProbMap ] )
    WFPerSubject.connect( sigmaInputNode, 'sigma' ,
                          SmoothingProbMap, 'inputSigma' )
    WFPerSubject.connect( WarpProbabilityMap, 'outputVolume' ,
                          SmoothingProbMap, 'inputVolume' )

    SmoothingProbMap.inputs.outputVolume= "outputSmoothProbabilityMap.nii.gz"

    # --------------------------------------------------------------------------------------- #
    # 3. Create Mask by Threshold
    #
    CreateMask = pe.Node( name = "03_CreateMask",
                          interface = Function( input_names = ["inputVolume",
                                                               "lowerThreshold",
                                                               "upperThreshold",
                                                               "outputFilename"],
                                                output_names = ["outputMask"],
                                                function =  MyUtilities.ThresholdProbabilityMap )
                          )
    
    CreateMask.inputs.lowerThreshold = 0.05
    CreateMask.inputs.upperThreshold = 0.95
    CreateMask.inputs.outputFilename = "roiMask.nii.gz"
    
    WFPerSubject.add_nodes( [CreateMask ] )
    WFPerSubject.connect( SmoothingProbMap, 'outputSmoothProbabilityMap',
                          CreateMask, 'inputVolume')
    
    WFPerSubject.connect( CreateMask, 'outputMask',
                          datasink, '03_ROI' )
#    
    # --------------------------------------------------------------------------------------- #
    # 4. Compute Statistice 
    #
    LabelStatistics = pe.Node( name = "04_LabelStatistics",
                               interface = Function( input_names = ["inputLabel",
                                                                    "inputVolume",
                                                                    "outputCSVFilename"],
                                                     output_names = ["outputCSVFilename",
                                                                     "outputDictionarySet"],
                                                     function = MyUtilities.LabelStatistics ),
                              )
    LabelStatistics.iterables  = ("inputVolume", [ inputListOfSubjectVolumes['t1'], 
                                                   inputListOfSubjectVolumes['t2']]
                                 )
    LabelStatistics.inputs.outputCSVFilename = "labelStatistics.csv"
    
    WFPerSubject.add_nodes( [LabelStatistics] )
    WFPerSubject.connect( CreateMask, "outputMask",
                          LabelStatistics, "inputLabel" )
    WFPerSubject.connect( LabelStatistics, 'outputCSVFilename',
                          datasink, 'labelStatistics')
    # --------------------------------------------------------------------------------------- #
    # 5. Linear Transform based on local statistics  and get statistics 
    #
    # Statistics has to be computed only for the ROI that applies 

    # dummy node for normalization iteration
    from nipype.interfaces.utility import IdentityInterface
    normalizationListNode= pe.Node( interface = IdentityInterface( fields = ['normalization'] ),
                              name="normalization" )
    NormalizationMethods = ['zScore', 
                            'MAD',
                            'DoubleSigmoid',
                            'Sigmoid',
                            'QEstimator',
                            'Linear']
    normalizationListNode.iterables = ( 'normalization' , NormalizationMethods)
    WFPerSubject.add_nodes([ normalizationListNode ] )


    NormalizeAndGetStatofROI_Node  = pe.Node( name = "05_NormalizeStats",
                                              interface = Function( input_names = ['inputSet_LabelStat',
                                                                                   'inputMethod',
                                                                                   'outputVolume',
                                                                                   'outputCSVFilename' ],
                                                                    output_names = ['outputCSV',
                                                                                    'outputVolume'],
                                                          function = MyUtilities.NormalizeAndComputeStatOfROI )
                                            )
    # inputs

    NormalizeAndGetStatofROI_Node.inputs.outputVolume = "NormalizedVolume.nii.gz"
    NormalizeAndGetStatofROI_Node.inputs.outputCSVFilename = "NormalizedVolumeStats.csv"
    # connect to the work flow
    WFPerSubject.add_nodes( [NormalizeAndGetStatofROI_Node] )
    WFPerSubject.connect( LabelStatistics, 'outputDictionarySet',
                          NormalizeAndGetStatofROI_Node, 'inputSet_LabelStat')
    WFPerSubject.connect( normalizationListNode, 'normalization',
                          NormalizeAndGetStatofROI_Node, 'inputMethod')
    
    # datasink
    WFPerSubject.connect( NormalizeAndGetStatofROI_Node, 'outputCSV',
                          datasink, 'labelStatistics.@Stat')

    WFPerSubject.connect( NormalizeAndGetStatofROI_Node, 'outputVolume',
                          datasink, 'labelStatistics.@Volume')
    WFPerSubject.run()
    WFPerSubject.write_graph(graph2use='orig')

    return LabelStatistics.outputs.outputCSVFilename
    ## end }



## main
#
# command line argument specification
#
#import argparse
#
#robustStatArgParser = argparse.ArgumentParser( description = "robust stat. argument")
#
#robustStatArgParser.add_argument('--inputSubjectT1', help='input subject average t1', 
#                                 required=True)
#robustStatArgParser.add_argument('--inputSubjectT2', help='input subject average t2', 
#                                 required=True)
#robustStatArgParser.add_argument('--inputInitialTransform', help='input initial trasnform from atlas to subject', 
#                                 required=True)
#robustStatArgParser.add_argument('--inputTemplateDir', help='input template directory', 
#                                 required=True)
#robustStatArgParser.add_argument('--outputDir',      help="output directory for the workflow", 
#                                 required=True)
#
#cmdArgs = robustStatArgParser.parse_args()
#
#WFPerSubnect( cmdArgs.inputSubjectT1, 
#              cmdArgs.inputSubjectT2,
#              cmdArgs.inputInitialTransform,
#              cmdArgs.inputTemplateDir,
#              cmdArgs.outputDir )

