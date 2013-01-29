## input: segmentations from BCut
##        manual traces
## output: a set of similarity measures
##         mean/median of them
##         ICC measures
##############################################################################
##############################################################################
def computeSummary( rObject ):
    import rpy2.robjects as robjects
    rObject.r('''
    require( psy )                                                                              
    require( xtable )

    writeLatexTable <- function(data, outputFilename, doAppend )
    {
      latexText <- xtable( data  )
      print( paste( "writeLatexTable::", outputFilename ) )
      print( latexText, file  = outputFilename,
                        append= doAppend)
    }

    numericCols <- c(                                                                           
        "FP", "SimilarityIndex", "totalSearchVol", "RelativeOverlap", "union", "Sensitivity",   
        "HausdorffAvg", "Precision",   "TP", "refVol",   "TN", "Hausdorff",      "beta", "sessionID",
        "FScore", "Specificity",       "alpha", "intersection",  "FN", "autoVol" )              
        
    dataColumns <- c( numericCols, 'roi' )                                                      
    
    error.bar <- function(x, y, upper, lower=upper, length=0.1,...){                            
      if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))   
      stop("vectors must be same length")
      arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)                           
    } 
    computeSummary <- function( csvFilename, 
                                outputFilePrefix
                                )
    {
      print( outputFilePrefix )                                                                 
      dt<-read.csv( csvFilename, header=T )                                                     
      print( head( dt ) )
      
      for ( cROI in levels( factor( dt$roi )) )                                                 
      {
        subResult <- data.frame( matrix( nrow=0,                                                
                                         ncol= length( numericCols ) ,                          
                                         dimnames = list( NULL, col.names=numericCols ) ))      
        print( paste("processing ", cROI ) )
        roiDT <- subset( dt, dt$roi == cROI , select=numericCols )                              
        print(head(roiDT))
        subResult[ 'min', ]       <- apply( roiDT, 2, min )                                     
        subResult[ 'which.min', ] <- roiDT$sessionID[ apply( roiDT, 2, which.min ) ]            
        subResult[ 'mean', ]      <- apply( roiDT, 2, mean ) 
        subResult[ 'max', ]       <- apply( roiDT, 2, max )                                     
        subResult[ 'which.max', ] <- roiDT$sessionID[ apply( roiDT, 2, which.max ) ]            
        
        print( subResult )                                                                      
        
        print( "write tex file" )                                                               
        currentFilePrefix = paste( outputFilePrefix, '_', cROI, sep = "" )                      
        write.csv( subResult, paste( currentFilePrefix, '_summary.csv', sep=""),                
                   quote = FALSE )
        writeLatexTable( subset( subResult, select=c( "FP", "TP", "TN", "FN", 
                                                  "FScore", "Precision", "Specificity", "Sensitivity", "alpha", "beta")),
                     paste( currentFilePrefix, '_summary.tex', sep=""),
                     FALSE )
        writeLatexTable( subset( subResult, select=c( "SimilarityIndex", "RelativeOverlap", "Hausdorff", "HausdorffAvg", 
                                                  "totalSearchVol", "autoVol", "refVol" ) ),
                         paste( currentFilePrefix, '_summary.tex', sep=""),
                         TRUE)
        
        print( "compute ICCs")
        iccDT <- subset( roiDT, select=c(autoVol, refVol))                                      
        iccResult <- icc( iccDT )
        print( iccResult )
        write.csv( iccResult, paste( currentFilePrefix, '_icc.csv', sep=""),                    
                   quote = FALSE, 
                   row.names = FALSE )                                                          

        print( "write ICC tex file" )
        iccTable <- data.frame( x=c("icc.agreement", "icc.consistency"),
                                v = c ( iccResult$icc.agreement, iccResult$icc.consistency )
                              )
        print( iccTable )

        writeLatexTable( iccTable, 
                         paste( currentFilePrefix, '_icc.tex', sep=""),
                         FALSE )
                   
        print( "ICC graph" )
        iccMin = min(iccDT)
        iccMax = max(iccDT)
        range = iccMax - iccMin 
        iccMin = iccMin - range *0.01
        iccMax = iccMax + range *0.01

        print( paste( "iccMin = ", iccMin ) )
        print( paste( "iccMax = ", iccMax ) )

        pdf( paste( currentFilePrefix, '_icc.pdf', sep="") )
        plot( iccDT , pch=19,
              xlim= c( iccMin, iccMax),
              ylim= c( iccMin, iccMax) )
        abline( a=0, b=1, col="red"  )
        legend( "topleft", 
                c( paste( "ICC(c):", round( iccResult$icc.consistency, 2) ),                    
                   paste( "ICC(a):", round( iccResult$icc.agreement, 2 ) )  ),                  
                bty="n",
                cex=2
              ) 
        dev.off() 
        pdf( paste( currentFilePrefix, '_summary.pdf', sep=""))
        print( paste( "Write ", currentFilePrefix, '_summary.pdf', sep="") )                    
        print( c( subResult[ 'mean', "SimilarityIndex"],
                  subResult[ 'mean', "RelativeOverlap" ],                                       
                  subResult[ 'mean', "Sensitivity" ],
                  subResult[ 'mean', "Specificity" ] ) )                                        
        print( subResult[ 'mean', "SimilarityIndex"] + subResult[ 'mean', "RelativeOverlap" ] ) 
        xbar <- barplot( c( subResult[ 'mean', "SimilarityIndex"],
                            subResult[ 'mean', "RelativeOverlap" ],                             
                            subResult[ 'mean', "Sensitivity" ],
                            subResult[ 'mean', "Specificity" ] ),                               
                          ylim=c(0,1)) 
        abline( h=c(0.2, 0.4, 0.6, 0.8), 
                col="blue",
                lty=2 )
        axis( 1, at = xbar,
                 labels = c( 'Similarity', 'RelativeOverlap', 'Sensitivity', 'Specificity' ) )  
        error.bar( xbar, c( subResult[ 'mean', "SimilarityIndex"],
                         subResult[ 'mean', "RelativeOverlap" ],
                         subResult[ 'mean', "Sensitivity" ],
                         subResult[ 'mean', "Specificity" ] ),
                         c( subResult[ 'min', "SimilarityIndex"],
                         subResult[ 'min', "RelativeOverlap" ],
                         subResult[ 'min', "Sensitivity" ],
                         subResult[ 'min', "Specificity" ] ),
                      c( subResult[ 'max', "SimilarityIndex"],
                         subResult[ 'max', "RelativeOverlap" ],
                         subResult[ 'max', "Sensitivity" ],
                         subResult[ 'max', "Specificity" ] ) )
        text( xbar, c( subResult[ 'min', "SimilarityIndex"],
                       subResult[ 'min', "RelativeOverlap" ],
                       subResult[ 'min', "Sensitivity" ],
                       subResult[ 'min', "Specificity" ] ),
                    c( subResult[ 'which.min', "SimilarityIndex"],
                       subResult[ 'which.min', "RelativeOverlap" ],
                       subResult[ 'which.min', "Sensitivity" ],
                       subResult[ 'which.min', "Specificity" ] ),
                    cex=1.3)
        text( xbar, c( subResult[ 'max', "SimilarityIndex"],
                       subResult[ 'max', "RelativeOverlap" ],
                       subResult[ 'max', "Sensitivity" ],
                       subResult[ 'max', "Specificity" ] ),
                    c( subResult[ 'which.max', "SimilarityIndex"],
                       subResult[ 'which.max', "RelativeOverlap" ],
                       subResult[ 'which.max', "Sensitivity" ],
                       subResult[ 'which.max', "Specificity" ] ),
                    cex=1.3)
        dev.off()
      }

    }
    ''')
    return rObject

#########################################################################################
def computeSummaryFromCSV( inputCSVFilename,
                           outputCSVPrefix
                           ):
    import rpy2.robjects as robjects
    import analysis as this
    robjects = this.computeSummary( robjects )
    rComputeSummary = robjects.globalenv['computeSummary']

    import os
    outputCSVPrefix = os.path.abspath( outputCSVPrefix )
    res = rComputeSummary( inputCSVFilename , 
                           outputCSVPrefix
                           )

    import glob
    outputCSVList = glob.glob( outputCSVPrefix + "*csv" ) 
    outputCSVList.append( glob.glob( outputCSVPrefix + "*pdf" ) )
    outputCSVList.append( glob.glob( outputCSVPrefix + "*tex" ) )
    return outputCSVList 


#########################################################################################
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
#########################################################################################
def getDefMask( img, tolerance ):
    import SimpleITK as sitk
    lowerThreshold = tolerance
    upperThreshold = 1.0 - tolerance
    binary =  sitk.BinaryThreshold( img, lowerThreshold, upperThreshold)
    return binary

#########################################################################################
def getLabelVolume( img, label=1):
    import SimpleITK as sitk
    binary =  sitk.BinaryThreshold( img, label, label)
    stat = sitk.LabelStatisticsImageFilter()
    stat.Execute( binary, binary )
    try:
        count = stat.GetCount( 1 )
    except:
        count = 0
        pass
    volume = count *(img.GetSpacing()[0]*img.GetSpacing()[1]*img.GetSpacing()[2])
    print( """Computed volume is
            {vl} mm^3""".format( vl=volume ))
    return volume

#########################################################################################
def printImageInfo( img ):
    import SimpleITK as sitk
    print("""Image info:::
          spacing: {sp}
          pixelID: {pid}
          dimension: {d}
          """.format( sp= img.GetSpacing(), 
                      pid = img.GetPixelIDValue(), 
                      d = img.GetDimension() )) 
#########################################################################################
def computeSimilarity( autoFilename, defFilename, refFilename, autoLabel, roi, session):
    import SimpleITK as sitk
    import os
    import analysis as this
    floatTolerance=0.01

    print( """ compute similarity of label : 
           {l}""".format( l = autoLabel ))

    autoImg = sitk.BinaryThreshold( sitk.ReadImage( autoFilename ), autoLabel, autoLabel )
    defImg = sitk.ReadImage( defFilename )
    refImg = sitk.BinaryThreshold( sitk.ReadImage( refFilename ), 1 )

    this.printImageInfo( autoImg )
    this.printImageInfo( defImg )
    this.printImageInfo( refImg )
    
    OUT = {}
    OUT[ 'roi' ] = roi
    OUT[ 'sessionID' ] = session 
    OUT['autoVol'] = this.getLabelVolume( autoImg )
    OUT['refVol'] = this.getLabelVolume( refImg )
    
    defMsk = this.getDefMask( defImg, floatTolerance )
    OUT['totalSearchVol'] = this.getLabelVolume( defMsk )
    
    OUT['union'] = this.getLabelVolume( autoImg | refImg )
    OUT['intersection'] = this.getLabelVolume( autoImg & refImg )
    OUT['TP'] = OUT['union']
    OUT['FP'] = this.getLabelVolume( autoImg - refImg, 1 )
    
    autoNeg = sitk.BinaryThreshold( (defMsk - autoImg), 1, 1)
    refNeg  = sitk.BinaryThreshold( (defMsk - refImg ), 1, 1)
    
    OUT['FN'] = this.getLabelVolume( autoNeg & refImg )
    OUT['TN'] = this.getLabelVolume( autoNeg & refNeg )
    
    OUT['alpha'] = OUT['FP'] / ( OUT['FP'] + OUT['TN'] )
    OUT['beta']  = OUT['FN'] / ( OUT['TP'] + OUT['FN'] )
    OUT['Sensitivity'] = OUT['TP'] / ( OUT['TP'] + OUT['FN'] )
    OUT['Specificity'] = OUT['TN'] / ( OUT['FP'] + OUT['TN'] )
    OUT['Precision'] = OUT['TP'] / (OUT['TP'] + OUT['FP'] )
    OUT['FScore'] = 2 * OUT['Precision'] * OUT['Sensitivity'] / ( OUT['Precision'] + OUT['Sensitivity'] )
    OUT['RelativeOverlap'] = OUT['intersection']/ OUT['union']
    OUT['SimilarityIndex'] = 2 * OUT['intersection'] / ( OUT['autoVol'] + OUT['refVol'] )
    
    if OUT['autoVol'] != 0:
        hausdorffFilter = sitk.HausdorffDistanceImageFilter()
        hausdorffFilter.Execute( autoImg, refImg )
        OUT['Hausdorff'] = hausdorffFilter.GetHausdorffDistance ()
        OUT['HausdorffAvg'] = hausdorffFilter.GetAverageHausdorffDistance()
    else:
        OUT['Hausdorff'] = -1
        OUT['HausdorffAvg'] = -1
    


    for ele in OUT.iterkeys():
        print( "{e} = {v}".format( e = ele, v = OUT[ele]))
    return OUT

#########################################################################################
def computeICCs( raterA, raterB):
    return ICCs

#########################################################################################
def writeCSV( dataDictList, 
              outputCSVFilename):
    import csv
    f = open( outputCSVFilename, 'wb')
    w = csv.DictWriter( f, dataDictList[0].keys() )
    w.writeheader()
    for row in dataDictList:
        w.writerow( row )
    f.close()

    return outputCSVFilename
#########################################################################################
def getData( ResultDir, 
             normalization, 
             methodParameter, 
             sessionID, 
             optionalString = ''): # ex = ANNLabel_seg.nii.gz
    import nipype.interfaces.io as nio
    #
    # get label map                                                                                                       
    #
    DG = nio.DataGrabber( infields = ['normalization','method','sessionID'],                                              
                          outfields = ['outputList'] )                                                                    
    DG.inputs.base_directory = ResultDir
    DG.inputs.template = 'Test*/%s/RF/%s/%s*' + optionalString
    DG.inputs.template_args = dict( outputList = [[ 'normalization', 'method', 'sessionID' ]]) 
    DG.inputs.normalization = normalization
    DG.inputs.method = methodParameter
    DG.inputs.sessionID =sessionID
    dt = DG.run()  
    print( """Data grabber with {opt}:::                                                                                  
           {str}
           """.format( opt = optionalString, str=dt.outputs.outputList ))          

    return dt.outputs.outputList
    

#########################################################################################
def experimentAnalysis( resultDir, 
                 outputCSVFilename, 
                 normalization,
                 methodParameter, 
                 manualDict,
                 roiList,
                 sessionIDList,
                 doRawComparison ):
    import nipype.interfaces.io as nio
    import ast
    import analysis as this

    roiLabel = {}
    label = 1
    for roi in sorted( set( roiList) ):
        roiLabel [ roi ] = label
        print( """{r} ===> {l}""".format( r=roi, l=label))
        label = label +1


    autoFileList = []
    defFileList = []
    refFileList = []
    autoLabelList = []
    identityDictList = []
    roiSeqList= []
    sessionList = []
    for sessionID in sessionIDList:
        #
        # get label map
        #

        roiManualDict = ast.literal_eval(  manualDict[ sessionID]['roiList']  )
        identityDict = { } 
        for roi in roiList:
            identityDict[ 'sessionID' ] = sessionID
            identityDict[ 'roi' ] = roi
            identityDictList.append( identityDict )
            roiSeqList.append( roi )
            sessionList.append( sessionID )

            if doRawComparison:
                labelDT = this.getData( resultDir,     
                                        normalization,
                                        methodParameter,
                                        sessionID,
                                        roi+"_ANNContinuous.nii.gz" )
                autoLabelList.append( 1 ) 

            else :
                labelDT = this.getData( resultDir,       
                                        normalization,
                                        methodParameter, 
                                        sessionID,
                                        "ANNLabel_seg.nii.gz" )
                autoLabelList.append( roiLabel[ roi ] ) 

            autoFileList.append( labelDT )

            print( """compute roi of ::
                   {s}
                   """.format( s = roi ))
            #
            # get get deformation 
            #
            defFileList.append( this.getData( resultDir,       normalization,
                                         methodParameter, sessionID,
                                         roi+"*.nii.gzdef.nii.gz" )  
                              )
            # read manual image
            refFileList.append( roiManualDict[ roi ] )


    print( """LENGTH:::
        length( roiSeqList) = {roiL}
        length( sessionList ) = {sessionL}
        length( autoFileList ) = {autuL}
        length( refFilename ) = {refL}
    """.format( roiL = len( roiSeqList), sessionL = len(sessionList), autuL=len(autoFileList), 
                refL = len( refFileList)  ))

    import nipype.pipeline.engine as pe
    import os
    workFlowName = 'experimentAnalysis'

    exp = pe.Workflow( name = workFlowName )
    outputCSVFilename = os.path.abspath( outputCSVFilename ) 
    exp.base_dir = os.path.dirname( outputCSVFilename )

    
    from nipype.interfaces.utility import Function 
    computeSimilarityND = pe.MapNode( name = "computeSimilarityND",
                                      interface = Function( input_names =[ 'autoFilename',
                                                                           'defFilename',
                                                                           'refFilename',
                                                                           'autoLabel',
                                                                           'roi',
                                                                           'session'] ,
                                                            output_names = [ 'outDict' ],
                                                            function = this.computeSimilarity ),
                                      iterfield = ['autoFilename', 
                                                   'defFilename',
                                                   'refFilename',
                                                   'autoLabel',
                                                   'roi',
                                                   'session']
                                    )
    computeSimilarityND.inputs.autoFilename = autoFileList
    computeSimilarityND.inputs.defFilename = defFileList
    computeSimilarityND.inputs.refFilename = refFileList
    computeSimilarityND.inputs.autoLabel = autoLabelList
    computeSimilarityND.inputs.roi = roiSeqList 
    computeSimilarityND.inputs.session = sessionList

    exp.add_nodes( [ computeSimilarityND ] )

    writeCSVFileND = pe.Node( name = 'writeCSVFileND',
                              interface = Function( input_names = ['dataDictList',
                                                                   'outputCSVFilename'],
                                                    output_names = ['outputCSVFilename'],
                                                    function = this.writeCSV )
                            )
    writeCSVFileND.inputs.outputCSVFilename = outputCSVFilename
    exp.connect( computeSimilarityND, 'outDict',
                      writeCSVFileND, 'dataDictList')
   
    exp.run()
    return outputCSVFilename 
     
#########################################################################################
def similarityComputeWorkflow( ResultDir, 
                               OutputDir,
                               ExperimentalConfigurationFile,
                               runOption,
                               PythonBinDir,
                               BRAINSStandAloneSrcDir,
                               BRAINSStandAloneBuildDir):
    
    import sys
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    #
    # get normalization option from experimental configuration file
    #
    import ConfigurationParser as configParser
    import analysis as this
    configMap = configParser.ConfigurationSectionMap( ExperimentalConfigurationFile )
    normalizationOptions = configMap[ 'Options'][ 'normalization' ]
    print(""" Normalization Option:::
          {str}
          """.format( str = normalizationOptions ))
  
    #
    # get methods
    #
    import ast
    methodOptionsDictList = configMap['Options'][ 'modelParameter'.lower() ]
    methodOptions = []
    print( methodOptionsDictList )
    for option in methodOptionsDictList:
        methodStr = 'TreeDepth'+str( option['--randomTreeDepth']) +'_TreeNumber'+str(option['--numberOfTrees'])
        methodOptions.append( methodStr )
    print(""" Method Option:::
          {str}
          """.format( str = methodStr ))

    #
    # get roiList
    #
    roiList = configMap[ 'Options']['roiBooleanCreator'.lower() ].keys()
    print(""" ROIList:::
          {str}
          """.format( str = roiList ))

    #
    # get sessionList and manualDict
    #
    import XMLConfigurationGenerator
    subjectListFilename = configMap['ListFiles'][ 'subjectListFilename'.lower() ]
    manualDict  = XMLConfigurationGenerator.combineCSVs( subjectListFilename, {} )
    sessionList = manualDict.keys()

    #
    # workflow
    #
    workFlowName = 'outputDataCollector'
    workflow = pe.Workflow( name = workFlowName  )
    workflow.base_dir = OutputDir 
    
    from nipype.interfaces.utility import Function
    experimentalND = pe.Node( name = 'experimentalND',
                              interface = Function( input_names=['resultDir',
                                                                 'outputCSVFilename',
                                                                 'normalization',
                                                                 'methodParameter',
                                                                 'manualDict',
                                                                 'roiList',
                                                                 'sessionIDList',
                                                                 'doRawComparison'],
                                                    output_names='outputCSVFilename',
                                                    function = this.experimentAnalysis 
                                                   )
                              )
    experimentalND.inputs.resultDir = ResultDir
    experimentalND.inputs.outputCSVFilename = 'experimentalResult.csv' 
    experimentalND.inputs.roiList = roiList
    experimentalND.inputs.manualDict = manualDict
    experimentalND.inputs.sessionIDList = sessionList
    #experimentalND.inputs.doRawComparison = doRawComparison
    experimentalND.iterables = [ ( 'normalization', normalizationOptions),
                                 ( 'methodParameter', methodOptions ),
                                 ( 'doRawComparison', [ True, False] )
                               ]
    workflow.add_nodes( [ experimentalND ] )

    summaryND = pe.Node( name = 'summaryND',
                         interface = Function( input_names = ['inputCSVFilename',
                                                              'outputCSVPrefix'
                                                              ],
                                               output_names = ['outputCSVList'],
                                               function = this.computeSummaryFromCSV )
                        )

    summaryND.inputs.outputCSVPrefix = 'summaryOutput'
    workflow.connect( experimentalND, 'outputCSVFilename',
                      summaryND, 'inputCSVFilename' )

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
        workflow.run()

def main(argv=None):
    import os
    import sys
    
    from nipype import config
    config.enable_debug_mode()
    
    #-------------------------------- argument parser
    import argparse
    argParser = argparse.ArgumentParser( description ="""****************************
        10-cross validation analysis 
        """)
    # workup arguments
    argWfGrp = argParser.add_argument_group( 'argWfGrp', """****************************
        auto workflow arguments for cross validation
        """)
    argWfGrp.add_argument( '--experimentalConfigurationFile',    
        help="""experimentalConfigurationFile
        Configuration file name with FULL PATH""", 
        dest='experimentalConfigurationFile', required=True )
    argWfGrp.add_argument( '--expDir',    help="""expDir
        """, 
        dest='expDir', required=False, default="." )
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

    args =argParser.parse_args()
    similarityComputeWorkflow( args.expDir, 
                               args.baseDir,
                               args.experimentalConfigurationFile,
                               args.runOption,
                               args.PythonBinDir,
                               args.BRAINSStandAloneSrcDir,
                               args.BRAINSStandAloneBuildDir)
                               
    

import sys

if __name__ == "__main__":
    sys.exit(main())
#########################################################################################
# unit test
#
#ResultDir = "/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec/Experiment_20121222/Result/Labels/"
#outputDir = '/ipldev/scratch/eunyokim/src/BRAINSStandAlone/BRAINSStandAlone/BRAINSCut/Nipype/'
#
#outputCSVFilename = '/ipldev/scratch/eunyokim/src/BRAINSStandAlone/BRAINSStandAlone/BRAINSCut/Nipype/output.csv'
#normalization = 'Linear'
#methodParameter = 'TreeDepth50_TreeNumber50'
#configFilename = '/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec/Dec22/model.config'
#
#similarityComputeWorkflow( ResultDir, outputDir, configFilename )
#
#import ConfigurationParser as configParser
#import XMLConfigurationGenerator
#configMap = configParser.ConfigurationSectionMap( configFilename )  
#subjectListFilename = configMap[ 'ListFiles' ]['subjectListFilename'.lower() ]
#manualDict  = XMLConfigurationGenerator.combineCSVs( subjectListFilename, {} )
#sessionIDList = manualDict.keys()
#print( sessionIDList )
#
#testDictList = [ {'roi':'l_accumben','vol':1000,'NF':3},
#                 {'roi':'r_accumben','vol':100,'NF':3} ]
#writeCSV(  testDictList, 'out.csv' )
#roiList =  [
#            'l_thalamus'
#           ]
#experimentAnalysis( ResultDir, 
#                    outputCSVFilename, 
#                    normalization, 
#                    methodParameter, 
#                    manualDict, 
#                    roiList, 
#                    sessionIDList )
#
