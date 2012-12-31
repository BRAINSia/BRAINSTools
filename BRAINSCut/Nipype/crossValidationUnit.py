def getProbabilityMapFilename( roiList):
    probabilityMapFilename = {}
    for roi in roiList:
        probabilityMapFilename[roi] = roi + "_probabilityMap.nii.gz" 
    return probabilityMapFilename


def unitWorkUp ( configurationFilename, 
                 doApply = False,
                 baseDir = "."):
    import os
    import sys
    import nipype.pipeline.engine as pe
    from nipype.interfaces.utility import Function
    import ConfigurationParser
    import crossValidationUnit as this
    
    from nipype import config
    config.enable_debug_mode()
    
    workflow = pe.Workflow( name = 'balancedTraning' )
    workflow.base_dir = baseDir
    
    configurationMap = ConfigurationParser.ConfigurationSectionMap( configurationFilename) 
    Options          = configurationMap[ 'Options' ]
    roiDict          = Options[ 'roiBooleanCreator'.lower() ]

    #
    #-------------------------------- filenameGeneratorND is dummy node
    # to create proper probability file location for nipype
    #

    filenameGeneratorND = pe.Node( name      = "filenameGeneratorND",
                                   interface = Function( 
                                      input_names  = ['roiList',
                                                      'gaussianSigma'],
                                      output_names = ['probabilityMapFilename'],
                                      function     = this.getProbabilityMapFilename )
                                 )
    filenameGeneratorND.inputs.roiList = roiDict.keys()

    #
    #--------------------------------  start from generate probability
    #
    probabilityMapGeneratorND = pe.Node( name = "probabilityMapGeneratorND",
                                         interface = Function( 
                                             input_names = ['configurationFilename',
                                                            'probabilityMapDict',
                                                            'gaussianSigma',
                                                            'outputXmlFilename'],
                                             output_names = [ 'probabilityMapDict',
                                                              'outputXmlFilename',
                                                              'outputConfigurationFilename'],
                                             function     = ConfigurationParser.BRAINSCutGenerateProbabilityMap )
                                       )
    
    probabilityMapGeneratorND.inputs.outputXmlFilename = 'netConfiguration.xml'
    probabilityMapGeneratorND.inputs.configurationFilename = configurationFilename 
    probabilityMapGeneratorND.inputs.gaussianSigma = Options[ 'gaussianSigma'.lower() ]
    
    workflow.connect( filenameGeneratorND, 'probabilityMapFilename',
                      probabilityMapGeneratorND, 'probabilityMapDict' )
    
    #
    #--------------------------------  create vectors for each ROI
    #
    configFileND = pe.Node( name = "configFileND",
                            interface = Function(
                                input_names = ['originalFilename',
                                               'editedFilenamePrefix' ],
                                output_names = ['editiedFilenames'],
                                function     = ConfigurationParser.ConfigurationFileEditor ) 
                          )
    
    configFileND.inputs.originalFilename = configurationFilename  
    configFileND.inputs.editedFilenamePrefix = 'ROI'
    workflow.add_nodes( [ configFileND ] )
    
    vectorCreatorND = pe.MapNode( name = "vectorCreatorND", 
                                  interface = Function(
                                      input_names = ['configurationFilename',
                                                     'probabilityMapDict',
                                                     'normalization',
                                                     'outputXmlFilename',
                                                     'outputVectorFilename'],
                                      output_names = ['outputVectorFilename',
                                                      'outputVectorHdrFilename',
                                                      'outputNormalization',
                                                      'outputXmlFilename'],
                                      function     = ConfigurationParser.BRAINSCutCreateVector ),
                                  iterfield = [ 'configurationFilename']
                                )
    vectorCreatorND.inputs.outputVectorFilename = 'oneROIVectorFile.txt'
    vectorCreatorND.inputs.outputXmlFilename = 'oneROICreateVectorNetConfiguration.xml'
    import ast
    normalizationOption = Options[ 'normalization'.lower()]  
    #normalizationOption = ast.literal_eval( Options[ 'normalization'.lower()]  )
    print( """Normalization Option: {str}
           """.format( str=normalizationOption ) )
    vectorCreatorND.iterables = ( 'normalization', normalizationOption )
    #
    #--------------------------------  workflow connections
    #
    workflow.connect( configFileND, 'editiedFilenames',
                      vectorCreatorND, 'configurationFilename' )
    workflow.connect( probabilityMapGeneratorND, 'probabilityMapDict',
                      vectorCreatorND, 'probabilityMapDict' )
    
    #
    #--------------------------------  balance and combine each ROI vectors
    #
    balaceND = pe.Node( name = "balanceND",
                        interface = Function(
                            input_names = ['inputVectorFilenames'],
                            output_names = ['outputVectorFilenames',
                                            'outputVectorHdrFilenames'],
                            function = ConfigurationParser.BalanceInputVectors )
                      )
    workflow.connect( vectorCreatorND, 'outputVectorFilename',
                      balaceND, 'inputVectorFilenames' )
    
    combineND = pe.Node( name = "combineND",
                         interface = Function(
                            input_names = ['inputVectorFilenames',
                                           'outputVectorFilename'],
                            output_names = ['outputVectorFilename',
                                            'outputVectorHdrFilename'],
                            function = ConfigurationParser.CombineInputVectors )
                       )
    workflow.connect( balaceND, 'outputVectorFilenames',
                      combineND, 'inputVectorFilenames')
    
    combineND.inputs.outputVectorFilename = 'allCombinedVector.txtANN'
    
    #
    #--------------------------------  train
    #
    trainND = pe.Node( name = "trainND", 
                       interface = Function( 
                           input_names = ['configurationFilename',
                                          'inputVectorFilename',
                                          'outputModelFilenamePrefix',
                                          'outputXmlFilename',
                                          'methodParameter'],
                           output_names = ['outputTrainedModelFilename',
                                           'outputMethodParameter'],
                           function = ConfigurationParser.BRAINSCutTrainModel ),
                     )
    #methodParameter = { '--method': 'RandomForest',
    #                    '--numberOfTrees': 60,
    #                    '--randomTreeDepth ': 60 }
    import ast
    methodFromConfiguFile =  Options['modelParameter'.lower()] 
    trainND.iterables= ( 'methodParameter', methodFromConfiguFile ) 
    trainND.inputs.outputXmlFilename = 'trianNetConfiguration.xml'
    trainND.inputs.outputModelFilenamePrefix = 'trainModelFile.txt'
    
    workflow.connect( probabilityMapGeneratorND, 'outputConfigurationFilename',
                      trainND, 'configurationFilename')
    workflow.connect( combineND, 'outputVectorFilename',
                      trainND, 'inputVectorFilename')
    
    #
    #--------------------------------  apply
    #
    # make output dir for each subject as a 
    if doApply:
        applyND = pe.Node( name = "applyND", 
                           interface = Function( 
                               input_names = ['configurationFilename',
                                              'probabilityMapDict',
                                              'normalization',
                                              'inputModelFilename',
                                              'methodParameter',
                                              'outputXmlFilename'
                                              ],
                               output_names = ['outputLabelDict'],
                               function = ConfigurationParser.BRAINSCutApplyModel )
                         )
        #methodParameter = { '--method': 'RandomForest',
        #                    '--numberOfTrees': 60,
        #                    '--randomTreeDepth ': 60 }
        applyND.inputs.outputXmlFilename = 'applyConfiguration.xml'
        workflow.connect( probabilityMapGeneratorND, 'outputConfigurationFilename',
                          applyND, 'configurationFilename')
        workflow.connect( vectorCreatorND, 'outputNormalization',
                          applyND, 'normalization' )
        workflow.connect( probabilityMapGeneratorND, 'probabilityMapDict',
                          applyND, 'probabilityMapDict' )
        workflow.connect( trainND, 'outputTrainedModelFilename',
                          applyND, 'inputModelFilename' )
        workflow.connect( trainND, 'outputMethodParameter',
                          applyND, 'methodParameter' )
        #
        # analysis
        #
        #analysisND = pe.Node( name = "analysisND",
        #                      interface = Function(
        #                          input_names['inputImageDict',
        #                                      'inputManualDict',
        #                                      'outputCSVFilename'],
        #                          output_names['outputCSVFilename'],
        #                          function = analysis.similarityFromApplyOutput )
        #                    )

    #
    #
    ##workflow.run(updatehash=True)
    workflow.run()
    
def main(argv=None):
    #-------------------------------- argument parser
    import argparse
    argParser = argparse.ArgumentParser( description ='10-cross validation command line argument parser')
    
    argParser.add_argument( '--configurationFilename',    help="""configurationFilename
                                                               Configuration file name with FULL PATH""", 
                            dest='configurationFilename', required=True )
    argParser.add_argument( '--doApply',    help="""doApply
                                                    """,
                            action = 'store_true',
                            dest='doApply' )
    argParser.add_argument( '--baseDir',    help="""baseDir
                                                    """, 
                            dest='baseDir', required=False, default="." )
    args = argParser.parse_args()
    
    unitWorkUp( args.configurationFilename,
                args.doApply,
                args.baseDir)

if __name__ == "__main__":
    import sys
    sys.exit(main())
