##############################################################################
## Utility functions to compute similarities between label maps
## Regina Kim
##############################################################################
##############################################################################
def getLabelVolume(img, label=1):
    import SimpleITK as sitk
    binary = sitk.BinaryThreshold(img, label, label)
    stat = sitk.LabelStatisticsImageFilter()
    stat.Execute(binary, binary)
    try:
        count = stat.GetCount(1)
    except:
        count = 0
        pass
    volume = count * (img.GetSpacing()[0] * img.GetSpacing()[1] * img.GetSpacing()[2])
    print( """Computed volume is
            {vl} mm^3""".format( vl=volume ))
    return volume

#########################################################################################


def printImageInfo(img):
    import SimpleITK as sitk
    print("""Image info:::
          spacing: {sp}
          pixelID: {pid}
          dimension: {d}
          """.format( sp=img.GetSpacing(),
                      pid=img.GetPixelIDValue(),
                      d=img.GetDimension()))
#########################################################################################


def getDefMask(img, tolerance):
    import SimpleITK as sitk
    lowerThreshold = tolerance
    upperThreshold = 1.0 - tolerance
    binary = sitk.BinaryThreshold(img, lowerThreshold, upperThreshold)
    return binary
#########################################################################################

def computeSimilarity(autoFilename, refFilename, autoLabel, roi, session, defFilename = NULL):
    import SimpleITK as sitk
    import os
    import analysis as this
    floatTolerance = 0.01

    print( """ compute similarity of label :
           {l}""".format( l=autoLabel ))

    autoImg = sitk.BinaryThreshold(sitk.ReadImage(autoFilename), autoLabel, autoLabel)

    refImg = sitk.BinaryThreshold(sitk.ReadImage(refFilename), 1)

    this.printImageInfo(autoImg)
    this.printImageInfo(refImg)

    OUT = {}
    OUT['roi'] = roi
    OUT['sessionID'] = session
    OUT['autoVol'] = this.getLabelVolume(autoImg)
    OUT['refVol'] = this.getLabelVolume(refImg)

    if defFilename not NULL:
        defImg = sitk.ReadImage(defFilename)
        this.printImageInfo(defImg)
        defMsk = this.getDefMask(defImg, floatTolerance)

        OUT['totalSearchVol'] = this.getLabelVolume(defMsk)

        autoNeg = sitk.BinaryThreshold((defMsk - autoImg), 1, 1)
        refNeg = sitk.BinaryThreshold((defMsk - refImg), 1, 1)

        OUT['FN'] = this.getLabelVolume(autoNeg & refImg)
        OUT['TN'] = this.getLabelVolume(autoNeg & refNeg)

    OUT['union'] = this.getLabelVolume(autoImg | refImg)
    OUT['intersection'] = this.getLabelVolume(autoImg & refImg)
    OUT['TP'] = OUT['union']
    OUT['FP'] = this.getLabelVolume(autoImg - refImg, 1)

    OUT['alpha'] = OUT['FP'] / (OUT['FP'] + OUT['TN'])
    OUT['beta'] = OUT['FN'] / (OUT['TP'] + OUT['FN'])
    OUT['Sensitivity'] = OUT['TP'] / (OUT['TP'] + OUT['FN'])
    OUT['Specificity'] = OUT['TN'] / (OUT['FP'] + OUT['TN'])
    OUT['Precision'] = OUT['TP'] / (OUT['TP'] + OUT['FP'])
    OUT['FScore'] = 2 * OUT['Precision'] * OUT['Sensitivity'] / (OUT['Precision'] + OUT['Sensitivity'])
    OUT['RelativeOverlap'] = OUT['intersection'] / OUT['union']
    OUT['SimilarityIndex'] = 2 * OUT['intersection'] / (OUT['autoVol'] + OUT['refVol'])

    if OUT['autoVol'] != 0:
        hausdorffFilter = sitk.HausdorffDistanceImageFilter()
        hausdorffFilter.Execute(autoImg, refImg)
        OUT['Hausdorff'] = hausdorffFilter.GetHausdorffDistance()
        OUT['HausdorffAvg'] = hausdorffFilter.GetAverageHausdorffDistance()
    else:
        OUT['Hausdorff'] = -1
        OUT['HausdorffAvg'] = -1

    for ele in OUT.iterkeys():
        print("{e} = {v}".format(e=ele, v=OUT[ele]))
    return OUT

#########################################################################################
def computeSummary(rObject):
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


def computeSummaryFromCSV(inputCSVFilename,
                          outputCSVPrefix
                          ):
    import rpy2.robjects as robjects
    import analysis as this
    robjects = this.computeSummary(robjects)
    rComputeSummary = robjects.globalenv['computeSummary']

    import os
    outputCSVPrefix = os.path.abspath(outputCSVPrefix)
    res = rComputeSummary(inputCSVFilename,
                          outputCSVPrefix
                          )

    import glob
    outputCSVList = glob.glob(outputCSVPrefix + "*csv")
    outputCSVList.append(glob.glob(outputCSVPrefix + "*pdf"))
    outputCSVList.append(glob.glob(outputCSVPrefix + "*tex"))
    return outputCSVList
