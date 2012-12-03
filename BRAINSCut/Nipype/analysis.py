## input: segmentations from BCut
##        manual traces
## output: a set of similarity measures
##         mean/median of them
##         ICC measures


def computeSimilarity( imageA, imageB):
    import SimpleITK as sitk
    imgA = sitk.ReadImage( imageA )
    imgB = sitk.ReadImage( imageB )

    binaryLowerThreshold = 1
    binaryUpperThreshold = 1
    binaryImgA = sitk.BinaryThreshold( imgA, 
                                       binaryLowerThreshold, 
                                       binaryUpperThreshold)
    binaryImgB = sitk.BinaryThreshold( imgB, 
                                       binaryLowerThreshold, 
                                       binaryUpperThreshold)

    outSI = sitk.SimilarityIndexImageFilter( 
    
                                       

    return volumeA, volumeB, RO, SI, HD
def computeICCs( raterA, raterB):
    return ICCs
def writeCSV( dataDict, 
              outputFilename):
    return outputFilename
#########################################################################################
def similarityFromApplyOutput( BCutResultDict, 
                               ManualDict,
                               outputCSVFilename ):
    
    import sys
    if not set( BCutResultDict.keys() ) in set( ManualDict.keys() ):
        print("""ERROR
              BCutResult dictionary has some ROI that does not have 
              reference volumes
              """)
        sys.exit()

    returnSimilarityDict = {}
    for session in  BCutResultDict.iterkeys(): 
        sessionResultDict= {} 
        for roi in session.iterkeys():
            sessionResultDict[ roi ] = computeSimilarity( BCutResultDict[ roi ],
                                                          ManualDict[ roi] )
        returnSimilarityDict[ session ] =sessionResultDict

    returnOutputCSVFilename = writeCSV( returnSimilarityDict, outputCSVFilename )

    return returnOutputCSVFilename

#########################################################################################
# unit test
#
 
