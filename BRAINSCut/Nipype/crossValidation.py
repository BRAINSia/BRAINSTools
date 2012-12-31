##############################################################################
def getDictionaryValues( inputDictionary ):
    print( inputDictionary )
    returnList = []
    for session in inputDictionary.iterkeys():
        sessionDict = inputDictionary[ session ]
        for label in sessionDict.iterkeys():
            for img in sessionDict[label]:
                returnList.append( img )
    return returnList


##############################################################################
def getProbabilityMapFilename( roiList):
    print("""getProbabilityMapFilename""")
    probabilityMapFilename = {}
    for roi in roiList:
        probabilityMapFilename[roi] = roi + "_probabilityMap.nii.gz" 
    return probabilityMapFilename
##############################################################################
def get_global_sge_script(pythonPathsList,binPathsList,customEnvironment={}):
    print("""get_global_sge_script""")
    """This is a wrapper script for running commands on an SGE cluster
so that all the python modules and commands are pathed properly"""

    custEnvString=""
    for key,value in customEnvironment.items():
        custEnvString+="export "+key+"="+value+"\n"

    PYTHONPATH=":".join(pythonPathsList)
    BASE_BUILDS=":".join(binPathsList)
    GLOBAL_SGE_SCRIPT="""#!/bin/bash
echo "STARTED at: $(date +'%F-%T')"
echo "Ran on: $(hostname)"
export PATH={BINPATH}
export PYTHONPATH={PYTHONPATH}

echo "========= CUSTOM ENVIORNMENT SETTINGS =========="
echo "export PYTHONPATH={PYTHONPATH}"
echo "export PATH={BINPATH}"
echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

echo "With custom environment:"
echo {CUSTENV}
{CUSTENV}
## NOTE:  nipype inserts the actual commands that need running below this section.
""".format(PYTHONPATH=PYTHONPATH,BINPATH=BASE_BUILDS,CUSTENV=custEnvString)
    return GLOBAL_SGE_SCRIPT
##############################################################################
def writeConfigFile( originalFilename,
                     outputConfigFilename, 
                     outputAdditionaListFiles ):
    
    print("""****************************
          writeConfigFile
          """)
    import ConfigParser

    inConfigParser = ConfigParser.ConfigParser()
    try:
        print( """read
               {fn}""".format( fn = originalFilename ))
        inConfigParser.read( originalFilename );
    except:
        print( """ERROR
               fail to read file {fn}
               """.format( fn = originalFilename ))


    outConfigParser = ConfigParser.RawConfigParser()

    for section in inConfigParser.sections():
        outConfigParser.add_section( section )
        for option in inConfigParser.options( section ):
             outConfigParser.set( section, option, inConfigParser.get( section ,option) )

    print( outputAdditionaListFiles )
    for option in outputAdditionaListFiles.iterkeys():
        print( outputAdditionaListFiles[ option ] )
        print( outputAdditionaListFiles[ option ] )
        print( outputAdditionaListFiles[ option ] )
        print( outputAdditionaListFiles[ option ] )
        outConfigParser.set( 'ListFiles', 
                             option, 
                             outputAdditionaListFiles[ option ] )
    
    with open( outputConfigFilename, 'wb' ) as outConfigFile:
        outConfigParser.write( outConfigFile )
    import os
    return os.path.abspath( outputConfigFilename)


#############################################################################
def writeListFile( sessionDict, 
                   outFilenameDict,
                   tagsToWrite):
    print("""****************************
          writeListFile
          """)
    import csv
    for outTag in outFilenameDict.iterkeys():
        outFile = open( outFilenameDict[ outTag ], "wb")

        writer = csv.DictWriter( outFile, 
                                 sessionDict[ sessionDict.keys()[0] ].keys() )
        writer.writeheader()

        print( tagsToWrite )
        for session in sessionDict.iterkeys():
            sessionRow = sessionDict[ session ]
            if tagsToWrite[ session ] == outTag:
                print( "Add {s} ".format( s = sessionRow['sessionID'] ))
                writer.writerow( sessionRow)
            else:
                print( "Drop {s} ".format( s = sessionRow['sessionID'] ))
        outFile.close()

    import os.path
    returnOutFilenameDict = {}
    for fn in outFilenameDict:
        returnOutFilenameDict = os.path.abspath( outFilenameDict[ fn] )
    return returnOutFilenameDict 

#############################################################################
def getStartAndEndIndex( p_iTh,
                         p_numberOfElementsPerSubset):
    print("""****************************
          getStartAndEndIndex
          """)
    if p_iTh == 0:
        StartIndex = 0
    else:
        StartIndex = sum( p_numberOfElementsPerSubset[0:p_iTh]) 
    EndIndx = StartIndex + p_numberOfElementsPerSubset[ p_iTh ] -1 
    print( p_numberOfElementsPerSubset )
    print( "({s},{e})".format( s=StartIndex, e=EndIndx ))
    return StartIndex, EndIndx

#############################################################################
def readListFileBySessionID( inputFilename,
                             totalNumber = -1 ):
    import csv
    listDict = {}

    try: 
        print("""Read
              {fn}""".format( fn = inputFilename ))
        with open(  inputFilename, "r") as inFile:
            reader=csv.reader( inFile, delimiter=",", skipinitialspace=True)
            header = reader.next()
            print( header )
            for row in reader:
                rowWithHeader = zip( header, row)
                rowDict = {}
                for ( name, value ) in rowWithHeader:
                    rowDict[ name ] = value.strip()
                listDict[ rowDict[ 'sessionID' ] ]=  rowDict 
    except:
        print( """ERROR
               fail to read file {fn}
               """.format( fn = inputFilename))
    import sys
    if totalNumber > 0 and len( listDict ) != totalNumber:
        print("""ERROR
              Total number of feature images are not equal to the main list.
              n( inputFilename ) = {n} != {t}
              """.format( n=len( listDict ) ,t=totalNumber ))
        sys.exit()

    return listDict
#############################################################################
def getTags( sessionList,
             nTest, 
             numberOfElementInSubset,
             randomize=False):
    print("""****************************
          getTags 
          """)
    if randomize:
        sessionOrder=getRandomizedSessionOrder(sessionList)
    else:
        orderIndex=0
        sessionOrder = {}
        for session in sessionList:
            sessionOrder[ session ] = orderIndex
            orderIndex = orderIndex +1
    applyStart, applyEnd = getStartAndEndIndex( nTest, numberOfElementInSubset)
    tags = {}
    for session in sessionList:
        if applyStart <= sessionOrder[ session ] and sessionOrder[ session ] <= applyEnd:
            tags[session]='Apply'
        else:
            tags[session]='Train'
    print("""Generate tags:::
             {t}""".format( t=tags))
    return tags 
#############################################################################
def getRandomizedSessionOrder( sessionList):
    print("""****************************
          getRandomizedSessionOrder 
          """)
    ## randomize the order 
    import random
    print ("""input list:::
           {s}""".format( s=sessionList ))
    random.shuffle( sessionList, random.random )
    print( """randomized list:::
           {s}""".format( s=sessionList ))
    orderIndex = 0
    sessionOrder = {}
    for s in sessionList:
        print( """ assign sessionOrder[{s}] = {orderIndex}""".format( s=s, orderIndex=orderIndex))
        sessionOrder[ s ] = orderIndex 
        orderIndex = orderIndex + 1
    return sessionOrder

#############################################################################
def generateNewFilenames( nTest, 
                          featureList,
                          outputPrefix):
    print("""****************************
          generateNewFilenames 
          """)
    returnConfigFilename = outputPrefix + "_Test" + str(nTest) + "_configuration.config"
    returnMainListFilename = { 'Train':outputPrefix + "_Test" + str(nTest) + "_mainTrainList.csv",
                               'Apply':outputPrefix + "_Test" + str(nTest) + "_mainApplyList.csv"}
    returnFeatureListFilenameDict = {}
    for ft in featureList:
        currentPefix=outputPrefix + "_" + "Test" +str(nTest) + "_" + str(ft)
        returnFeatureListFilenameDict[ft] = {'Train': currentPefix + "_featureTrainList.csv",
                                             'Apply': currentPefix + "_featureApplyList.csv" }
    print("""
          returnMainListFilename: {fn1}
          returnFeatureListFilenameDict: {fn2}
          """.format( fn1=returnMainListFilename, fn2=returnFeatureListFilenameDict))
    return returnConfigFilename, returnMainListFilename, returnFeatureListFilenameDict
    
#############################################################################
def createConfigurationFileForCrossValidationUnitTest( inputConfigurationFilename,
                                                       outputConfigurationFilenamePrefix):
    print("""****************************
          createConfigurationFileForCrossValidationUnitTest
          """)
    import os.path
    outputConfigurationFilenamePrefix = os.path.abspath( outputConfigurationFilenamePrefix )

    import ConfigurationParser
    m_configurationMap =  ConfigurationParser.ConfigurationSectionMap( inputConfigurationFilename )

    # get list filenames
    import crossValidation as this
    listFilenames = m_configurationMap[ 'ListFiles' ]
    mainListFilename = listFilenames[ 'subjectListFilename'.lower() ] 
    featureListFilenamesDict = listFilenames[ 'featureListFileDictionary'.lower() ]
    numberOfElementsInSubset = listFilenames[ 'numberOfElementInSubset'.lower() ]
    numberOfTotalSession = sum(numberOfElementsInSubset)

    # read files into sessionID -> data
    mainSessionDict = this.readListFileBySessionID( mainListFilename, 
                                                    numberOfTotalSession)
    featureSessionDict = {}
    if  len(featureListFilenamesDict) > 0:
        for ft in featureListFilenamesDict.iterkeys():
            featureSessionDict[ft] = this.readListFileBySessionID( featureListFilenamesDict[ft],
                                                                   numberOfTotalSession)

    #{ iterate throug subsets
    outputConfigFilenameDict= {}
    for nTest in range ( 0, len( numberOfElementsInSubset ) ): 
        trainApplyTagList = this.getTags( mainSessionDict.keys(),
                                          nTest, 
                                          numberOfElementsInSubset);
        newConfigFilename, newMainFilename, newFeatureFilenameDict = this.generateNewFilenames( nTest, 
                                   featureListFilenamesDict.keys(),
                                   outputConfigurationFilenamePrefix)
        this.writeListFile( mainSessionDict, 
                            newMainFilename,
                            trainApplyTagList)
        trainFeatureStr = {}
        applyFeatureStr = {}
        print( "++++++++++++++++++++++++++++++++newFeatureFilenameDict++++++++++++++++++++++++++++++++")
        print( newFeatureFilenameDict )
        if len( featureSessionDict ) > 0 :
            for ft in featureSessionDict.iterkeys():
                this.writeListFile( featureSessionDict[ft],
                                    newFeatureFilenameDict[ft],
                                    trainApplyTagList)
                trainFeatureStr[ft]=newFeatureFilenameDict[ft]['Train']
                applyFeatureStr[ft]=newFeatureFilenameDict[ft]['Apply']
            
        print( newMainFilename['Train'] )
        print( newMainFilename['Apply'] )
        print( trainFeatureStr )
        print( applyFeatureStr )
        this.writeConfigFile( inputConfigurationFilename,
                              newConfigFilename, 
                              {'subjectListFilename':newMainFilename['Train'],
                               'applySubjectListFilename':newMainFilename['Apply'],
                               'featureListFileDictionary':str( trainFeatureStr ),
                               'applyFeatureListFileDictionary':str( applyFeatureStr )})
        outputConfigFilenameDict[ "Test"+ str(nTest) ] = os.path.abspath( newConfigFilename )
    return outputConfigFilenameDict
#############################################################################

def extractConfigFile ( configurationFiledict ):
    print("""****************************
          extractConfigFile
          """)
    print( configurationFiledict )
    return configurationFiledict.values()

def crossValidationWorkUp( crossValidationConfigurationFilename,
                           baseDir,
                           runOption,
                           PythonBinDir,
                           BRAINSStandAloneSrcDir,
                           BRAINSStandAloneBuildDir):
    print("""****************************
          crossValidationWorkUp
          """)
    from nipype import config
    config.enable_debug_mode()

    import crossValidation as this
    import ConfigurationParser
    myConfigurationMap = ConfigurationParser.ConfigurationSectionMap( 
                              crossValidationConfigurationFilename )
    
    import nipype.pipeline.engine as pe
    from nipype.interfaces.utility import Function
    import ast
    print( """ before
           createeachvalidationunitnd
           """)
    createConfigurationFiles = pe.Node( name = "createConfigurationFiles",
                                        interface = Function(  
                                           input_names = ['inputConfigurationFilename',
                                                          'outputConfigurationFilenamePrefix'],
                                           output_names = ['outputConfigFilenameDict'],
                                           function = this.createConfigurationFileForCrossValidationUnitTest )
                                      )

    preprocessing = pe.Workflow( name = 'Preprocessing' )
    preprocessing.base_dir = baseDir +"/PreprocessingDir"

    createConfigurationFiles.inputs.inputConfigurationFilename = crossValidationConfigurationFilename
    createConfigurationFiles.inputs.outputConfigurationFilenamePrefix = 'createConfigurationFiles'
    
    extractConfigurationFileListND = pe.Node( name = "extractConfigurationFileListND",
                                              interface = Function(
                                                  input_names = ['configurationFiledict'],
                                                  output_names = ['configurationFileList'],
                                                  function = this.extractConfigFile )
                                            )
    preprocessing.connect( createConfigurationFiles, 'outputConfigFilenameDict',
                           extractConfigurationFileListND, 'configurationFiledict')

    preprocessing.run()

    #------------------------------------------------------------------------------------
    # Data graber for outputs 
    #
    import nipype.interfaces.io as nio
    dg = nio.DataGrabber()
    dg.inputs.base_directory = baseDir + "/PreprocessingDir/Preprocessing/createConfigurationFiles/"
    dg.inputs.template = "*config"
    mainConfigFiles = dg.run()

    print( mainConfigFiles.outputs.outfiles  )
    print( mainConfigFiles.outputs.outfiles  )
    print( mainConfigFiles.outputs.outfiles  )
    print( mainConfigFiles.outputs.outfiles  )
    
    #------------------------------------------------------------------------------------
    workflow = pe.Workflow( name = 'crossValidationWF' )
    workflow.base_dir = baseDir

    #------------------------------------------------------------------------------------
    # Generate Probability Map 
    #
    Options = myConfigurationMap[ 'Options' ]
    roiDict = Options[ 'roiBooleanCreator'.lower() ]


    #-------------------------------- probMapFilenameGenerator is dummy node
    # to create proper probability file location for nipype
    #
    print("""************************
          probMapFilenameGenerator
          """)

    probMapFilenameGenerator = pe.Node( name      = "probMapFilenameGenerator",
                                        interface = Function( 
                                           input_names  = ['roiList'],
                                           output_names = ['probabilityMapFilename'],
                                           function     = this.getProbabilityMapFilename )
                                      )
    print( roiDict)
    probMapFilenameGenerator.inputs.roiList = roiDict.keys()
    print("""************************
          probabilityMapGeneratorND
          """)

    
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

    gaussianSigmaParam = ast.literal_eval( Options[ 'gaussianSigma'.lower() ] ) 
    print ( gaussianSigmaParam )
    probabilityMapGeneratorND.iterables =  ('configurationFilename', mainConfigFiles.outputs.outfiles)
    probabilityMapGeneratorND.inputs.gaussianSigma = gaussianSigmaParam

    workflow.connect( probMapFilenameGenerator, 'probabilityMapFilename',
                      probabilityMapGeneratorND, 'probabilityMapDict' )
    
    #
    #--------------------------------  create vectors for each ROI
    #
    print("""************************
          configFileND 
          """)
    configFileND = pe.Node( name = "configFileND",
                            interface = Function(
                                input_names = ['originalFilename',
                                               'editedFilenamePrefix' ],
                                output_names = ['editedFilenames'],
                                function     = ConfigurationParser.ConfigurationFileEditor ) 
                          )
    
    configFileND.inputs.editedFilenamePrefix = 'ROI'
    workflow.connect( probabilityMapGeneratorND, 'outputConfigurationFilename',
                      configFileND, 'originalFilename' )
    
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
                                  iterfield = ['configurationFilename']
                                )
    vectorCreatorND.inputs.outputVectorFilename = 'oneROIVectorFile.txt'
    vectorCreatorND.inputs.outputXmlFilename = 'oneROICreateVectorNetConfiguration.xml'
    normalizationOption = Options[ 'normalization'.lower()]  
    print( """Normalization Option: {str}
           """.format( str=normalizationOption ) )
    vectorCreatorND.iterables = ( 'normalization', normalizationOption )
    #
    #--------------------------------  workflow connections
    #
    workflow.connect( configFileND, 'editedFilenames',
                      vectorCreatorND, 'configurationFilename' )
    workflow.connect( probabilityMapGeneratorND, 'probabilityMapDict',
                      vectorCreatorND, 'probabilityMapDict' )
    
    #
    #--------------------------------  balance and combine each ROI vectors
    #
    print("""************************
          balanceND
          """)
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
    print("""************************
          trainND 
          """)
    trainND = pe.Node( name = "trainND", 
                       interface = Function( 
                           input_names = ['configurationFilename',
                                          'inputVectorFilename',
                                          'outputModelFilenamePrefix',
                                          'outputXmlFilename',
                                          'methodParameter'],
                           output_names = ['outputTrainedModelFilename',
                                           'outputMethodParameter'],
                           function = ConfigurationParser.BRAINSCutTrainModel )
                     )
    #methodParameter = { '--method': 'RandomForest',
    #                    '--numberOfTrees': 60,
    #                    '--randomTreeDepth ': 60 }
    methodFromConfiguFile = Options['modelParameter'.lower()] 
    trainND.iterables= ( 'methodParameter', methodFromConfiguFile)
                         
    trainND.inputs.outputXmlFilename = 'trianNetConfiguration.xml'
    trainND.inputs.outputModelFilenamePrefix = 'trainModelFile.txt'
    
    workflow.connect( probabilityMapGeneratorND, 'outputConfigurationFilename',
                      trainND, 'configurationFilename')
    workflow.connect( combineND, 'outputVectorFilename',
                      trainND, 'inputVectorFilename')
    #
    #--------------------------------  apply
    #
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

    #####################################################################################
    # Data Sink
    #
    import os
    LabelsDS = pe.Node( nio.DataSink(), name='LabelDS')
    LabelsDS.inputs.base_directory = os.path.join( baseDir , "Result" )
    LabelsDS.inputs.regexp_substitutions = [ ('/_', '/'),
                                             ('configurationFilename.*_Test','Test'),
                                             ('_configuration.config/normalization_','/'),
                                             ('methodParameter_--method',''),
                                             ('RandomForest','RF/'),
                                             ('.--randomTreeDepth','TreeDepth'),
                                             ('.--numberOfTrees','_TreeNumber'),
                                             ('ANNContinuousPrediction(?P<roi>.+)(?P<session>\d\d\d\d\d).nii.gz',r'\g<session>_\g<roi>_ANNContinuous.nii.gz')
                                             ]
    #ANNContinuousPredictionl_accumben77478

    workflow.connect( [ ( applyND, LabelsDS,
                          [ ( ( 'outputLabelDict', getDictionaryValues), 'Labels')] ) ] )

    #####################################################################################
    # analysis 
    #

    #####################################################################################
    # Running
    #
    if runOption == "cluster":
        ############################################
        # Platform specific information
        #     Prepend the python search paths
        pythonPath = BRAINSStandAloneSrcDir + "/BRAINSCut/BRAINSFeatureCreators/RobustStatisticComputations:" + BRAINSStandAloneSrcDir + "/AutoWorkup/:" + BRAINSStandAloneSrcDir + "/AutoWorkup/BRAINSTools/:" + BRAINSStandAloneBuildDir + "/SimpleITK-build/bin/" + BRAINSStandAloneBuildDir + "/SimpleITK-build/lib:" + PythonBinDir
        binPath = BRAINSStandAloneBuildDir + "/bin:" + BRAINSStandAloneBuildDir+ "/lib"

        PYTHON_AUX_PATHS= pythonPath
        PYTHON_AUX_PATHS=PYTHON_AUX_PATHS.split(':')                                                                                  
        PYTHON_AUX_PATHS.extend(sys.path)                                                                                             
        sys.path=PYTHON_AUX_PATHS                                                                                                     
        #print sys.path
        import SimpleITK as sitk
        #     Prepend the shell environment search paths
        PROGRAM_PATHS= binPath
        PROGRAM_PATHS=PROGRAM_PATHS.split(':')
        import os
        PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))                                                                           
        os.environ['PATH']=':'.join(PROGRAM_PATHS)

        Cluster_Script = get_global_sge_script( PYTHON_AUX_PATHS, 
                                                PROGRAM_PATHS,
                                                {}
                                              )
        workflow.run( plugin='SGE',
                      plugin_args = dict( template = Cluster_Script, 
                                          qsub_args = "-S /bin/bash -pe smp1 4-8 -o /dev/null "))
    else:
        print("""************************
              run 
              """)
        try:
            workflow.write_graph(graph2use='flat')
        except:
            pass
        workflow.run()


def main(argv=None):
    import os
    import sys
    import nipype.pipeline.engine as pe
    from nipype.interfaces.utility import Function
    import ConfigurationParser
    
    from nipype import config
    config.enable_debug_mode()
    
    workflow = pe.Workflow( name = 'crossValidation' )
    workflow.base_dir = '.'
    
    #-------------------------------- argument parser
    import argparse
    argParser = argparse.ArgumentParser( description ="""****************************
        10-cross validation command line argument parser
        """)
    # workup arguments
    argWfGrp = argParser.add_argument_group( 'argWfGrp', """****************************
        auto workflow arguments for cross validation
        """)
    argWfGrp.add_argument( '--crossValidationConfigurationFilename',    
        help="""configurationFilename
        Configuration file name with FULL PATH""", 
        dest='crossValidationConfigurationFilename', required=True )
    argWfGrp.add_argument( '--baseDir',    help="""baseDir
        """, 
        dest='baseDir', required=False, default="." )
    argWfGrp.add_argument( '--runOption',    help="""runOption [local/cluster]
        """, 
        dest='runOption', required=False, default="local" )
    argWfGrp.add_argument( '--PythonBinDir',    help="""PythonBinDir [local/cluster]
        """, 
        dest='PythonBinDir', required=False, default="NA" )
    argWfGrp.add_argument( '--BRAINSStandAloneSrcDir',    help="""BRAINSStandAloneSrcDir [local/cluster]
        """, 
        dest='BRAINSStandAloneSrcDir', required=False, default="NA" )
    argWfGrp.add_argument( '--BRAINSStandAloneBuildDir',    help="""BRAINSStandAloneBuildDir [local/cluster]
        """, 
        dest='BRAINSStandAloneBuildDir', required=False, default="NA" )

    # test arguments
    argTestGrp = argParser.add_argument_group( 'argTestGrp', """****************************
        arguments for testing
        """)
    argTestGrp.add_argument( '--unitTest', action='store_true',
        dest='unitTest', help="""****************************
        List of test function name
        """)
    args = argParser.parse_args()

    #-------------------------------- 
    if not args.unitTest: 
        crossValidationWorkUp ( args.crossValidationConfigurationFilename,
                                args.baseDir,
                                args.runOption,
                                args.PythonBinDir,
                                args.BRAINSStandAloneSrcDir, 
                                args.BRAINSStandAloneBuildDir)
    
    #-------------------------------- 
    if args.unitTest:
        testElementPerSubject = [ 3, 4, 5 ]
        getStartAndEndIndex ( 0, testElementPerSubject )
        getStartAndEndIndex ( 1, testElementPerSubject )
        getStartAndEndIndex ( 2, testElementPerSubject )
        
        featureDict = {'GadSG':'testGadFeatureList.csv',
                       't2':'t2FeatureList.csv'}


        sessionList = ["s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12"]
        getRandomizedSessionOrder( sessionList )
        myTag = getTags( sessionList, 
                         2,
                         testElementPerSubject)
        featureFilenameDict = {'f1':'f1.csv', 'f2':'f2.csv'}
        configFilename,mainFilenameDict, featureFilenameDict = generateNewFilenames( 3,
                                                   featureFilenameDict.keys(),
                                                   "outputPrefix")
        import ConfigurationParser
        m_configurationMap =  ConfigurationParser.ConfigurationSectionMap( args.crossValidationConfigurationFilename )

        listFiles = m_configurationMap[ 'ListFiles' ]
        mainListFilename = listFiles['subjectListFilename'.lower() ]  
        sessionDict = readListFileBySessionID( mainListFilename )
        myTag = getTags( sessionDict.keys(),
                         2,
                         listFiles[  'numberOfElementInSubset'.lower() ] )
        writeListFile( sessionDict ,
                       mainFilenameDict,
                       myTag )

    
import sys

if __name__ == "__main__":
    sys.exit(main())
