
# --------------------------------------------------------------------------------------- #
def GetProbabilityFilename  ( inputROI,
                              inputDir
                              ):

  outputFilename = inputDir + "/" + "/probabilityMaps/" + inputROI + "_ProbabilityMap.nii.gz"
  return outputFilename

# --------------------------------------------------------------------------------------- #
def ThresholdProbabilityMap ( inputVolume, 
                              lowerThreshold,
                              upperThreshold,
                              outputFilename):
  import os
  import sys
  import SimpleITK as sitk

  inImg = sitk.Cast( sitk.ReadImage( inputVolume ),
                     sitk.sitkFloat32 )
  outImg = sitk.BinaryThreshold( inImg, 
                                 lowerThreshold,
                                 upperThreshold )

  sitk.WriteImage( outImg, outputFilename )
  returnFile = os.path.realpath( outputFilename )

  return returnFile

# --------------------------------------------------------------------------------------- #
def SmoothProbabilityMap ( inputVolume, 
                           inputSigma,
                           outputVolume):
  import os
  import sys
  import SimpleITK as sitk

  inImg = sitk.Cast( sitk.ReadImage( inputVolume ),
                     sitk.sitkFloat32 )

  normalizeAcrossScale = False

  if inputSigma <= 0:
    outImg = inImg
  else :
    smoother = sitk.SmoothingRecursiveGaussianImageFilter()
    outImg = smoother.Execute( inImg,
                               inputSigma,
                               normalizeAcrossScale )

  sitk.WriteImage( outImg, outputVolume )
  returnFile = os.path.realpath( outputVolume )

  return returnFile


# --------------------------------------------------------------------------------------- #
def LabelStatistics ( inputLabel,
                      inputVolume,
                      outputCSVFilename):
    labelValue=1

    import SimpleITK as sitk

    statCalculator = sitk.LabelStatisticsImageFilter()
    inImg = sitk.Cast( sitk.ReadImage( inputVolume ),
                       sitk.sitkFloat32 )
    inMsk = sitk.ReadImage( inputLabel) 


    statCalculator.Execute( inImg, inMsk )

    ## volume
    imageSpacing = inMsk.GetSpacing()
    Volume = imageSpacing[0] * imageSpacing[1] * imageSpacing[2] * statCalculator.GetCount( labelValue )

    Mean = statCalculator.GetMean( labelValue )
    Median = statCalculator.GetMedian( labelValue )
    Minimum = statCalculator.GetMinimum( labelValue )
    Maximum = statCalculator.GetMaximum( labelValue )
    Sigma = statCalculator.GetSigma( labelValue )
    Variance = statCalculator.GetVariance( labelValue )

    outputDictionary = dict()

    outputDictionary[ 'Mean' ] = Mean
    outputDictionary[ 'Median' ] = Median
    outputDictionary[ 'Minimum' ] = Minimum
    outputDictionary[ 'Maximum' ] = Maximum
    outputDictionary[ 'Sigma' ] = Sigma
    outputDictionary[ 'Variance' ] = Variance
    print "################################################## Mean:: "     + str(Mean) 
    print "################################################## Median:: "   + str(Median)
    print "################################################## Minimum:: "  + str(Minimum)
    print "################################################## Maximum:: "  + str(Maximum)
    print "################################################## Sigma:: "    + str(Sigma)
    print "################################################## Variance:: " + str(Variance)

    ## TODO 25/75 quantiles
    LowerHalfMsk = sitk.BinaryThreshold( inImg, Minimum, Median )
    Quantile25Calculator = sitk.LabelStatisticsImageFilter()
    Quantile25Calculator.Execute( inImg, inMsk*LowerHalfMsk  )
    Quantile25 = Quantile25Calculator.GetMedian( labelValue )
    outputDictionary[ 'Quantile25' ] = Quantile25
    print "################################################## Quantile25:: " + str(Quantile25)

    UpperHalfMsk = sitk.BinaryThreshold( inImg, Median, Maximum )
    Quantile75Calculator = sitk.LabelStatisticsImageFilter()
    Quantile75Calculator.Execute( inImg, inMsk*UpperHalfMsk )
    Quantile75 = Quantile75Calculator.GetMedian( labelValue )
    outputDictionary[ 'Quantile75' ] = Quantile75
    print "################################################## Quantile75:: " + str(Quantile75)
    
    ## TODO MAD 
    AbsoluteFilter = sitk.AbsImageFilter()
    AbsImg = AbsoluteFilter.Execute( inImg-Median )

    MADCalculator = sitk.LabelStatisticsImageFilter()
    MADCalculator.Execute( AbsImg, inMsk )
    MAD = MADCalculator.GetMedian( labelValue )
    outputDictionary[ 'MAD' ] = MAD
    print "################################################## MAD:: " + str(MAD)

    import csv
    csvFile = open( outputCSVFilename, 'w')
    dWriter = csv.DictWriter( csvFile, outputDictionary.keys() )
    dWriter.writeheader()
    dWriter.writerow( outputDictionary )

    import os
    import sys
    returnFile = os.path.realpath( outputCSVFilename )

    outputDictionarySet = dict( volume = inputVolume,
                                labelVolume = inputLabel , 
                                statDict = outputDictionary )
    #returnDict = { 'outputStatDictionary': outputDictionary,
    #               'outputCSVFilename': returnFile }
    #return  returnDict
    return returnFile, outputDictionarySet

# --------------------------------------------------------------------------------------- #
def NormalizeInputVolume ( inputVolume,
                           inputMethod, 
                           inputStats,
                           outputVolume ):
  import SimpleITK as sitk
  import os
  from math import e
  
  inImg = sitk.Cast( sitk.ReadImage( inputVolume ),
                     sitk.sitkFloat32 )

  IQR = inputStats['Quantile75']-inputStats['Quantile25']

  # create exp image
  #expImg = sitk.Image ( inImg.GetSize(), sitk.sitkFloat32 );
  #expImg.CopyInformation( inImg )
  expImg = inImg
  expImg = expImg * 0 # make zero image
  expImg = expImg + e # make exp. image

  if inputMethod == 'zScore':
    print "zScore Normalization"
    outImg = (inImg - inputStats['Mean']) / inputStats['Sigma']
  elif inputMethod == 'MAD':
    print "MAD Normalization"
    outImg = (inImg - inputStats['Median']) / inputStats['MAD']
  elif inputMethod == 'Sigmoid':
    print "Sigmoid Normalization"
    outImg = 1 / ( 1 + expImg ** ( -2 * ( inImg - inputStats['Median'] )/IQR ) )
  elif inputMethod == 'QEstimator':
    print "QEstimator Normalization"
    outImg = (inImg - inputStats['Median']) / IQR
  elif inputMethod == 'Linear':
    print "Linear Normalization"
    outImg = ( inImg - inputStats['Minimum'] ) / ( inputStats['Maximum']  - inputStats['Minimum'] )
  elif inputMethod == 'DoubleSigmoid':
    print "Double Sigmoid Normalization"
    outMsk1 = sitk.BinaryThreshold( inImg, inputStats['Minimum'], inputStats['Median'] )
    outImg1 = 1 / ( 1 + expImg  ** ( -2 * ( inImg - ( inputStats['Median'] ) )/ 
                        (  inputStats['Median'] - inputStats['Quantile25'] ) ) )

    outMsk2 = sitk.BinaryThreshold( inImg, ( inputStats['Median']+0.00001) , inputStats['Maximum'] )
    outImg2 = 1 / ( 1 + expImg  ** ( -2 * ( inImg - inputStats['Median']  )/ 
                        ( inputStats['Quantile75']- (inputStats['Median']+0.00001  )) ) )
    outImg = outImg1* sitk.Cast( outMsk1, sitk.sitkFloat32 ) + outImg2* sitk.Cast( outMsk2, sitk.sitkFloat32 )

  sitk.WriteImage( outImg, outputVolume )
  returnFile = os.path.realpath( outputVolume )

  return returnFile

# --------------------------------------------------------------------------------------- #
def NormalizeAndComputeStatOfROI( inputSet_LabelStat,
                                  inputMethod,
                                  outputVolume,
                                  outputCSVFilename ):
  from MyUtilities import NormalizeInputVolume 
  m_outputVolume = NormalizeInputVolume( inputSet_LabelStat['volume'], 
                                         inputMethod, 
                                         inputSet_LabelStat['statDict'],
                                         outputVolume )

  from MyUtilities import LabelStatistics
  m_outputCSV, m_outputDictSet = LabelStatistics( inputSet_LabelStat['labelVolume'], 
                                                  m_outputVolume, 
                                                  outputCSVFilename )
  return m_outputCSV, m_outputVolume

