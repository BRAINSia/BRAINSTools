"""
usage: atlasSmallIslandCleanup.py --inputAtlasPath=<argument> --outputAtlasPath=<argument> --inputT1Path=<argument> [--inputT2Path=<argument>] [--includeLabelsList=<argument> | --excludeLabelsList=<argument>] --maximumIslandVoxelCount=<argument> [--useFullyConnectedInConnectedComponentFilter] [--forceSuspiciousLabelChange] [--noDilation]
atlasSmallIslandCleanup.py -h | --help
"""

import SimpleITK as sitk
import math

class DustCleanup():

  def __init__(self, arguments):
    self.inputAtlasPath = arguments['--inputAtlasPath']
    self.outputAtlasPath = arguments['--outputAtlasPath']
    self.inputT1Path = arguments['--inputT1Path']
    self.inputT2Path = arguments['--inputT2Path']
    self.includeLabelsList = self.evalInputListArg(arguments['--includeLabelsList'])
    self.excludeLabelsList = self.evalInputListArg(arguments['--excludeLabelsList'])
    self.maximumIslandVoxelCount = int(arguments['--maximumIslandVoxelCount'])
    self.useFullyConnectedInConnectedComponentFilter = arguments['--useFullyConnectedInConnectedComponentFilter']
    self.forceSuspiciousLabelChange = arguments['--forceSuspiciousLabelChange']
    self.noDilation = arguments['--noDilation']
    self.islandStatistics = {'Total': {'numberOfIslandsCleaned': 0, 'numberOfIslands': 0}}

  def evalInputListArg(self, inputArg):
    if inputArg:
      return map(int, inputArg.split(','))
    else:
      return None

  def main(self):
    labelImage = sitk.Cast(sitk.ReadImage(self.inputAtlasPath), sitk.sitkInt16)
    inputT1VolumeImage = sitk.ReadImage(self.inputT1Path)
    if self.inputT2Path:
      inputT2VolumeImage = sitk.ReadImage(self.inputT2Path)
    else:
      inputT2VolumeImage = None
    labelsList = self.getLabelsList(inputT1VolumeImage, labelImage)
    for label in labelsList:
      labelImage = self.relabelCurrentLabel(labelImage, inputT1VolumeImage, inputT2VolumeImage, label)
    self.printIslandStatistics()
    sitk.WriteImage(labelImage, self.outputAtlasPath)

  def getLabelsList(self, volumeImage, labelImage):
    labelStatsObject = self.getLabelStatsObject(volumeImage, labelImage)
    labelsList = self.getLabelListFromLabelStatsObject(labelStatsObject)
    if self.excludeLabelsList:
      return self.removeLabelsFromLabelsList(labelsList, self.excludeLabelsList)
    if self.includeLabelsList:
      return self.verifyIncludeLabelsList(labelsList, self.includeLabelsList)
    return labelsList

  def removeLabelsFromLabelsList(self, labelsList, excludeList):
    for val in excludeList:
      try:
        labelsList.remove(val)
      except ValueError:
        #print "WARNING: Can not remove Label value", val, "is NOT a valid label in the input atlas:", self.inputAtlasPath
        pass
    return labelsList

  def verifyIncludeLabelsList(self, labelsList, includeList):
    verifiedList = list()
    for val in includeList:
      if val in labelsList:
        verifiedList.append(val)
      else:
        #print "WARNING: Can not validate Label value", val, "is NOT a valid label in the input atlas:", self.inputAtlasPath
        pass
    return verifiedList

  def printIslandStatistics(self):
    print "-"*50
    print "Label, numberOfIslandsCleaned, numberOfIslands, IslandVoxelCount, numberOfIslandsCleanedForIslandVoxelCount"
    for val in sorted(self.islandStatistics):
      labelStats = [str(val), str(self.islandStatistics[val]['numberOfIslandsCleaned']),
                     str(self.islandStatistics[val]['numberOfIslands'])]
      if val != 'Total':
        for i in range(1, self.maximumIslandVoxelCount + 1):
          labelStats.extend([str(i), str(self.islandStatistics[val][i])])
      print ','.join(labelStats)

  def relabelCurrentLabel(self, labelImage, inputT1VolumeImage, inputT2VolumeImage, label):

    self.islandStatistics[label] = {'numberOfIslandsCleaned': 0}

    for currentIslandSize in range(1, self.maximumIslandVoxelCount + 1):
      maskForCurrentLabel = sitk.BinaryThreshold(labelImage, label, label)
      relabeledConnectedRegion = self.getRelabeldConnectedRegion(maskForCurrentLabel, currentIslandSize)
      labelStatsT1WithRelabeledConnectedRegion = self.getLabelStatsObject(inputT1VolumeImage, relabeledConnectedRegion)
      if inputT2VolumeImage:
        labelStatsT2WithRelabeledConnectedRegion = self.getLabelStatsObject(inputT2VolumeImage, relabeledConnectedRegion)
      labelList = self.getLabelListFromLabelStatsObject(labelStatsT1WithRelabeledConnectedRegion)
      labelList.remove(0)  #remove background label from labelList
      labelList.reverse()

      if currentIslandSize == 1: #use island size 1 to get # of islands since this label map is not dilated
        self.islandStatistics[label]['numberOfIslands'] = len(labelList)
        self.islandStatistics['Total']['numberOfIslands'] += len(labelList)

      numberOfIslandsCleaned = 0

      for currentLabel in labelList:
        islandVoxelCount = labelStatsT1WithRelabeledConnectedRegion.GetCount(currentLabel)
        if islandVoxelCount < currentIslandSize:
          continue
        elif islandVoxelCount == currentIslandSize and currentLabel != 1: #stop if you reach largest island
          meanT1Intensity = labelStatsT1WithRelabeledConnectedRegion.GetMean(currentLabel)
          if inputT2VolumeImage:
            meanT2Intensity = labelStatsT2WithRelabeledConnectedRegion.GetMean(currentLabel)
          else:
            meanT2Intensity = None
          targetLabels = self.getTargetLabels(labelImage, relabeledConnectedRegion, inputT1VolumeImage, currentLabel)
          diffDict = self.calculateLabelIntensityDifferenceValue(meanT1Intensity, meanT2Intensity,
                                                                 targetLabels, inputT1VolumeImage,
                                                                 inputT2VolumeImage, labelImage)
          if self.forceSuspiciousLabelChange:
            diffDict.pop(label)
          sortedLabelList = self.getDictKeysListSortedByValue(diffDict)
          currentLabelBinaryThresholdImage = sitk.BinaryThreshold(relabeledConnectedRegion, currentLabel, currentLabel)
          labelImage = self.relabelImage(labelImage, currentLabelBinaryThresholdImage, sortedLabelList[0])
          numberOfIslandsCleaned += 1
        else:
          break

      self.islandStatistics[label][currentIslandSize] = numberOfIslandsCleaned
      self.islandStatistics[label]['numberOfIslandsCleaned'] += numberOfIslandsCleaned
      self.islandStatistics['Total']['numberOfIslandsCleaned'] += numberOfIslandsCleaned

    return labelImage

  def getRelabeldConnectedRegion(self, maskForCurrentLabel, currentIslandSize):
    if (currentIslandSize > 1) and (not self.noDilation):
      dilationKernelRadius = self.calcDilationKernelRadius(currentIslandSize)
      dilatedMaskForCurrentLabel = self.dilateLabelMap(maskForCurrentLabel, dilationKernelRadius)
      relabeledConnectedLabelMap = self.runConnectedComponentsAndRelabel(dilatedMaskForCurrentLabel)
      return sitk.Mask(relabeledConnectedLabelMap, maskForCurrentLabel, outsideValue=0)
    else:
      return self.runConnectedComponentsAndRelabel(maskForCurrentLabel)

  def calcDilationKernelRadius(self, currentIslandSize):
    # use the equation for the volume of a sphere to calculate the kernel radius value
    return int(math.ceil(math.pow(currentIslandSize/((4./3.)*math.pi), (1./3.))))

  def runConnectedComponentsAndRelabel(self, binaryImage):
    if not self.useFullyConnectedInConnectedComponentFilter:
      connectedRegion = sitk.ConnectedComponent(binaryImage, fullyConnected=False)
    else:
      connectedRegion = sitk.ConnectedComponent(binaryImage, fullyConnected=True)
    relabeledConnectedRegion = sitk.RelabelComponent(connectedRegion)
    return relabeledConnectedRegion

  def getLabelStatsObject(self, volumeImage, labelImage):
    labelStatsObject = sitk.LabelStatisticsImageFilter()
    labelStatsObject.Execute(volumeImage, labelImage)

    return labelStatsObject

  def getLabelListFromLabelStatsObject(self, labelStatsObject):
    if sitk.Version().MajorVersion() > 0 or sitk.Version().MinorVersion() >= 9:
      compontentLabels = labelStatsObject.GetLabels()
    else: #if sitk version < 0.9 then use older function call GetValidLabels
      compontentLabels = labelStatsObject.GetValidLabels()
    return list(compontentLabels)

  def getTargetLabels(self, labelImage, relabeledConnectedRegion, inputVolumeImage, currentLabel):
    currentLabelBinaryThresholdImage = sitk.BinaryThreshold(relabeledConnectedRegion, currentLabel, currentLabel)
    castedCurrentLabelBinaryThresholdImage = sitk.Cast(currentLabelBinaryThresholdImage, sitk.sitkInt16)

    dilatedBinaryLabelMap = self.dilateLabelMap(castedCurrentLabelBinaryThresholdImage, 1)
    outsideValue = -1
    reducedLabelMapImage = sitk.Mask(labelImage, dilatedBinaryLabelMap, outsideValue=outsideValue)

    reducedLabelMapT1LabelStats = self.getLabelStatsObject(inputVolumeImage, reducedLabelMapImage)
    targetLabels = self.getLabelListFromLabelStatsObject(reducedLabelMapT1LabelStats)
    targetLabels = self.removeOutsideValueFromTargetLabels(targetLabels, outsideValue)
    return targetLabels

  def removeOutsideValueFromTargetLabels(self, targetLabels, outsideValue):
    if outsideValue in targetLabels:
      targetLabels.remove(outsideValue)
    return targetLabels

  def dilateLabelMap(self, inputLabelImage, kernelRadius):
    myFilter = sitk.BinaryDilateImageFilter()
    myFilter.SetBackgroundValue(0.0)
    myFilter.SetBoundaryToForeground(False)
    myFilter.SetDebug(False)
    myFilter.SetForegroundValue(1.0)
    myFilter.SetKernelRadius((kernelRadius, kernelRadius, kernelRadius))
    myFilter.SetKernelType(2)  # Kernel Type=Box
    myFilter.SetNumberOfThreads(8)
    output = myFilter.Execute(inputLabelImage)
    castedOutput = sitk.Cast(output, sitk.sitkInt16)

    return castedOutput

  def calculateLabelIntensityDifferenceValue(self, averageT1IntensitySuspiciousLabel,
                                             averageT2IntensitySuspiciousLabel,
                                             targetLabels, inputT1VolumeImage,
                                             inputT2VolumeImage, inputLabelImage):
    """
    Calculates a measurement for each label that is on the border of the suspicious label.
    This value is the square root of the sum of the squared difference in the average T1
    intensity values and the squared difference in the average T2 intensity values of the
    two islands in the comparison. The calculated value for each border label will later be
    sorted in ascending order - meaning that the smallest value has the "closest" average
    intensity to the suspicious label.
    """

    squareRootDiffLabelDict = dict()
    labelStatsT1WithInputLabelImage = self.getLabelStatsObject(inputT1VolumeImage, inputLabelImage)
    if inputT2VolumeImage:
      labelStatsT2WithInputLabelImage = self.getLabelStatsObject(inputT2VolumeImage, inputLabelImage)

    for targetLabel in targetLabels:
      averageT1IntensityTargetLabel = labelStatsT1WithInputLabelImage.GetMean(targetLabel)
      squareDiffAverageT1 = math.pow(averageT1IntensitySuspiciousLabel -
                                     averageT1IntensityTargetLabel, 2)
      if inputT2VolumeImage:
        averageT2IntensityTargetLabel = labelStatsT2WithInputLabelImage.GetMean(targetLabel)
        squareDiffAverageT2 = math.pow(averageT2IntensitySuspiciousLabel -
                                       averageT2IntensityTargetLabel, 2)
      else:
        squareDiffAverageT2 = 0
      squareRootDiff = math.sqrt(squareDiffAverageT1 + squareDiffAverageT2)

      squareRootDiffLabelDict[int(targetLabel)] = squareRootDiff

    return squareRootDiffLabelDict

  def relabelImage(self, labelImage, newRegion, newLabel):
    castedLabelImage = sitk.Cast(labelImage, sitk.sitkInt16)
    castedNewRegion = sitk.Cast(newRegion, sitk.sitkInt16)
    negatedMask = sitk.BinaryNot(castedNewRegion)
    negatedImage = sitk.Mask(castedLabelImage, negatedMask)
    maskTimesNewLabel = sitk.Multiply(castedNewRegion, newLabel)
    relabeledImage = sitk.Add(negatedImage, maskTimesNewLabel)
    return relabeledImage

  def getDictKeysListSortedByValue(self, val):
    return sorted(val, key=val.get)

if __name__ == '__main__':
  from docopt import docopt
  arguments = docopt(__doc__)
  print arguments
  print "-"*50
  Object = DustCleanup(arguments)
  Object.main()
