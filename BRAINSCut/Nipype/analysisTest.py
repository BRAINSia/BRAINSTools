#########################################################################################
def computeSummary( rObject ):
    import rpy2.robjects as robjects
    rObject.r('''
    require( psy )
    numericCols <- c(
        "FP", "SimilarityIndex", "totalSearchVol", "RelativeOverlap", "union", "Sensitivity",
        "HausdorffAvg", "Precision",   "TP", "refVol",   "TN", "Hausdorff",      "beta", "sessionID",
        "FScore", "Specificity",       "alpha", "intersection",  "FN", "autoVol" )
    
    dataColumns <- c( numericCols, 'roi' )

    computeSummary <- function( csvFilename, outputFilePrefix )
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
        subResult[ 'which.min', ] <- as.character( roiDT$sessionID[ apply( roiDT, 2, which.min ) ] )
        subResult[ 'mean', ]      <- apply( roiDT, 2, mean ) 
        subResult[ 'max', ]       <- apply( roiDT, 2, max ) 
        subResult[ 'which.max', ] <- as.character( roiDT$sessionID[ apply( roiDT, 2, which.max ) ] )

        print( subResult )

        print( outputFilePrefix )
        write.csv( subResult, paste( outputFilePrefix, '_', cROI, '_summary.csv', sep=""),
                   quote = FALSE )
        
        iccDT <- subset( roiDT, select=c(autoVol, refVol))
        iccResult <- icc( iccDT )
        write.csv( iccResult, paste( outputFilePrefix, '_', cROI, '_icc.csv', sep=""),
                   quote = FALSE, 
                   row.names = FALSE )


      }
    }
    ''')
    return rObject


  
#########################################################################################
def computeSummaryFromCSV( inputCSVFilename,
                           outputCSVPrefix):
    import rpy2.robjects as robjects
    robjects = computeSummary( robjects )
    rComputeSummary = robjects.globalenv['computeSummary']
    res = rComputeSummary( inputCSVFilename , outputCSVPrefix )




#computeSummaryFromCSV( '/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec/Dec29/outputDataCollector/_methodParameter_TreeDepth50_TreeNumber50_normalization_DoubleSigmoid_Q01/experimentalND/experimentalResult.csv' )
 
#computeSummaryFromCSV( '/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec/Experiment_20121222/outputDataCollector/_methodParameter_TreeDepth50_TreeNumber50_normalization_Linear/experimentalND/experimentalResult.csv')
#computeSummaryFromCSV('test.csv')
