"""
atlasSmallIslandCleanup.py
============================
Description:

Author:

Usage:

usage: atlasSmallIslandCleanup.py --inputAtlasPath=<argument> --outputAtlasPath=<argument> --inputT1Path=<argument> [--inputT2Path=<argument>] [--includeLabelsList=<argument> | --excludeLabelsList=<argument>] --maximumIslandVoxelCount=<argument> [--useFullyConnectedInConnectedComponentFilter] [--forceSuspiciousLabelChange] [--noDilation]
atlasSmallIslandCleanup.py -h | --help
"""

import SimpleITK as sitk
import math


class DustCleanup:
    """
    This class represents a...
    """

    def __init__(self, arguments):
        from collections import OrderedDict

        self.inputAtlasPath = arguments["--inputAtlasPath"]
        self.outputAtlasPath = arguments["--outputAtlasPath"]
        self.inputT1Path = arguments["--inputT1Path"]
        self.inputT2Path = arguments["--inputT2Path"]
        self.includeLabelsList = self.eval_input_list_arg(arguments["--includeLabelsList"])
        self.excludeLabelsList = self.eval_input_list_arg(arguments["--excludeLabelsList"])
        self.maximumIslandVoxelCount = int(arguments["--maximumIslandVoxelCount"])
        self.useFullyConnectedInConnectedComponentFilter = arguments[
            "--useFullyConnectedInConnectedComponentFilter"
        ]
        self.forceSuspiciousLabelChange = arguments["--forceSuspiciousLabelChange"]
        self.noDilation = arguments["--noDilation"]
        self.islandStatistics = OrderedDict(
            {"Total": {"numberOfIslandsCleaned": 0, "numberOfIslands": 0}}
        )

    def eval_input_list_arg(self, inputArg):
        """
        This function...
        :param inputArg:
        :return: list(map(int, inputArg.split(','))) OR None
        """
        if inputArg:
            return list(map(int, inputArg.split(",")))
        else:
            return None

    def main(self):
        """
        This is the main function...
        """
        labelImage = sitk.Cast(sitk.ReadImage(self.inputAtlasPath), sitk.sitkInt16)
        inputT1VolumeImage = sitk.ReadImage(self.inputT1Path)
        if self.inputT2Path:
            inputT2VolumeImage = sitk.ReadImage(self.inputT2Path)
        else:
            inputT2VolumeImage = None
        labelsList = self.get_labels_list(inputT1VolumeImage, labelImage)
        for label in labelsList:
            labelImage = self.relabel_current_label(
                labelImage, inputT1VolumeImage, inputT2VolumeImage, label
            )
        self.print_island_statistics()
        sitk.WriteImage(labelImage, self.outputAtlasPath)

    def get_labels_list(self, volumeImage, labelImage):
        """
        This function....
        :param volumeImage:
        :param labelImage:
        :return: labelsList
        """
        labelStatsObject = self.get_label_stats_object(volumeImage, labelImage)
        labelsList = self.get_label_list_from_label_stats_object(labelStatsObject)
        if self.excludeLabelsList:
            return self.remove_labels_from_labels_list(labelsList, self.excludeLabelsList)
        if self.includeLabelsList:
            return self.verify_include_labels_list(labelsList, self.includeLabelsList)
        return labelsList

    def remove_labels_from_labels_list(self, labelsList, excludeList):
        """
        This function...
        :param labelsList:
        :param excludeList:
        :return: labelsList
        """
        for val in excludeList:
            try:
                labelsList.remove(val)
            except ValueError:
                # print "WARNING: Can not remove Label value", val, "is NOT a valid label in the input atlas:", self.inputAtlasPath
                pass
        return labelsList

    def verify_include_labels_list(self, labelsList, includeList):
        """
        This function...
        :param labelsList:
        :param includeList:
        :return: verifiedList
        """
        verifiedList = list()
        for val in includeList:
            if val in labelsList:
                verifiedList.append(val)
            else:
                # print "WARNING: Can not validate Label value", val, "is NOT a valid label in the input atlas:", self.inputAtlasPath
                pass
        return verifiedList

    def print_island_statistics(self):
        """
        This function prints...
        """
        print()
        "-" * 50
        print()
        "Label, numberOfIslandsCleaned, numberOfIslands, IslandVoxelCount, numberOfIslandsCleanedForIslandVoxelCount"
        for islandName in sorted(self.islandStatistics.keys()):
            labelStats = [
                str(islandName),
                str(self.islandStatistics[islandName]["numberOfIslandsCleaned"]),
                str(self.islandStatistics[islandName]["numberOfIslands"]),
            ]
            if islandName != "Total":
                for i in range(1, self.maximumIslandVoxelCount + 1):
                    labelStats.extend(
                        [str(i), str(self.islandStatistics[islandName][i])]
                    )
            print()
            ",".join(labelStats)

    def relabel_current_label(
        self, labelImage, inputT1VolumeImage, inputT2VolumeImage, label_key
    ):
        """
        This function...
        :param labelImage:
        :param inputT1VolumeImage:
        :param inputT2VolumeImage:
        :param label_key:
        :return: labelImage
        """
        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        label_key = str(label_key)  # all keys must be strings in order to sort
        self.islandStatistics[label_key] = OrderedDict({"numberOfIslandsCleaned": 0})
        label_value = int(label_key)
        for currentIslandSize in range(1, self.maximumIslandVoxelCount + 1):
            maskForCurrentLabel = sitk.BinaryThreshold(
                labelImage, label_value, label_value
            )
            relabeledConnectedRegion = self.get_relabeld_connected_region(
                maskForCurrentLabel, currentIslandSize
            )
            labelStatsT1WithRelabeledConnectedRegion = self.get_label_stats_object(
                inputT1VolumeImage, relabeledConnectedRegion
            )
            if inputT2VolumeImage:
                labelStatsT2WithRelabeledConnectedRegion = self.get_label_stats_object(
                    inputT2VolumeImage, relabeledConnectedRegion
                )
            labelList = self.get_label_list_from_label_stats_object(
                labelStatsT1WithRelabeledConnectedRegion
            )
            labelList.remove(0)  # remove background label from labelList
            labelList.reverse()

            if (
                currentIslandSize == 1
            ):  # use island size 1 to get # of islands since this label map is not dilated
                self.islandStatistics[label_key]["numberOfIslands"] = len(labelList)
                self.islandStatistics["Total"]["numberOfIslands"] += len(labelList)

            numberOfIslandsCleaned = 0

            for currentLabel in labelList:
                islandVoxelCount = labelStatsT1WithRelabeledConnectedRegion.GetCount(
                    currentLabel
                )
                if islandVoxelCount < currentIslandSize:
                    continue
                elif (
                    islandVoxelCount == currentIslandSize and currentLabel != 1
                ):  # stop if you reach largest island
                    meanT1Intensity = labelStatsT1WithRelabeledConnectedRegion.GetMean(
                        currentLabel
                    )
                    if inputT2VolumeImage:
                        meanT2Intensity = labelStatsT2WithRelabeledConnectedRegion.GetMean(
                            currentLabel
                        )
                    else:
                        meanT2Intensity = None
                    targetLabels = self.get_target_labels(
                        labelImage,
                        relabeledConnectedRegion,
                        inputT1VolumeImage,
                        currentLabel,
                    )
                    diffDict = self.calculate_label_intensity_difference_value(
                        meanT1Intensity,
                        meanT2Intensity,
                        targetLabels,
                        inputT1VolumeImage,
                        inputT2VolumeImage,
                        labelImage,
                    )
                    if self.forceSuspiciousLabelChange:
                        diffDict.pop(label_key)
                    sortedLabelList = [
                        int(x) for x in self.get_dict_keys_list_sorted_by_value(diffDict)
                    ]
                    currentLabelBinaryThresholdImage = sitk.BinaryThreshold(
                        relabeledConnectedRegion, currentLabel, currentLabel
                    )
                    labelImage = self.relabel_image(
                        labelImage, currentLabelBinaryThresholdImage, sortedLabelList[0]
                    )
                    numberOfIslandsCleaned += 1
                else:
                    break

            self.islandStatistics[label_key][currentIslandSize] = numberOfIslandsCleaned
            self.islandStatistics[label_key][
                "numberOfIslandsCleaned"
            ] += numberOfIslandsCleaned
            self.islandStatistics["Total"][
                "numberOfIslandsCleaned"
            ] += numberOfIslandsCleaned

        return labelImage

    def get_relabeld_connected_region(self, maskForCurrentLabel, currentIslandSize):
        """
        This function...
        :param maskForCurrentLabel:
        :param currentIslandSize:
        :return: sitk.Mask(relabeledConnectedLabelMap, maskForCurrentLabel, outsideValue=0) OR self.run_connected_components_and_relabel(maskForCurrentLabel)
        """
        if (currentIslandSize > 1) and (not self.noDilation):
            dilationKernelRadius = self.calc_dilation_kernel_radius(currentIslandSize)
            dilatedMaskForCurrentLabel = self.dilate_label_map(
                maskForCurrentLabel, dilationKernelRadius
            )
            relabeledConnectedLabelMap = self.run_connected_components_and_relabel(
                dilatedMaskForCurrentLabel
            )
            return sitk.Mask(
                relabeledConnectedLabelMap, maskForCurrentLabel, outsideValue=0
            )
        else:
            return self.run_connected_components_and_relabel(maskForCurrentLabel)

    def calc_dilation_kernel_radius(self, currentIslandSize):
        """
        This function...
        :param currentIslandSize:
        :return: int(math.ceil(math.pow(currentIslandSize / ((4. / 3.) * math.pi), (1. / 3.))))
        """
        # use the equation for the volume of a sphere to calculate the kernel radius value
        return int(
            math.ceil(
                math.pow(currentIslandSize / ((4.0 / 3.0) * math.pi), (1.0 / 3.0))
            )
        )

    def run_connected_components_and_relabel(self, binaryImage):
        """
        This function...
        :param binaryImage:
        :return: relabeledConnectedRegion
        """
        if not self.useFullyConnectedInConnectedComponentFilter:
            connectedRegion = sitk.ConnectedComponent(binaryImage, False)
        else:
            connectedRegion = sitk.ConnectedComponent(binaryImage, True)
        relabeledConnectedRegion = sitk.RelabelComponent(connectedRegion)
        return relabeledConnectedRegion

    def get_label_stats_object(self, volumeImage, labelImage):
        """
        This function...
        :param volumeImage:
        :param labelImage:
        :return: labelStatsObject
        """
        labelStatsObject = sitk.LabelStatisticsImageFilter()
        labelStatsObject.Execute(volumeImage, labelImage)

        return labelStatsObject

    def get_label_list_from_label_stats_object(self, labelStatsObject):
        """
        This function...
        :param labelStatsObject:
        :return: list(compontentLabels)
        """
        if sitk.Version().MajorVersion() > 0 or sitk.Version().MinorVersion() >= 9:
            compontentLabels = labelStatsObject.GetLabels()
        else:  # if sitk version < 0.9 then use older function call GetValidLabels
            compontentLabels = labelStatsObject.GetValidLabels()
        return list(compontentLabels)

    def get_target_labels(
        self, labelImage, relabeledConnectedRegion, inputVolumeImage, currentLabel
    ):
        """
        This function...
        :param labelImage:
        :param relabeledConnectedRegion:
        :param inputVolumeImage:
        :param currentLabel:
        :return: targetLabels
        """
        currentLabelBinaryThresholdImage = sitk.BinaryThreshold(
            relabeledConnectedRegion, currentLabel, currentLabel
        )
        castedCurrentLabelBinaryThresholdImage = sitk.Cast(
            currentLabelBinaryThresholdImage, sitk.sitkInt16
        )

        dilatedBinaryLabelMap = self.dilate_label_map(
            castedCurrentLabelBinaryThresholdImage, 1
        )
        outsideValue = -1
        reducedLabelMapImage = sitk.Mask(
            labelImage, dilatedBinaryLabelMap, outsideValue=outsideValue
        )

        reducedLabelMapT1LabelStats = self.get_label_stats_object(
            inputVolumeImage, reducedLabelMapImage
        )
        targetLabels = self.get_label_list_from_label_stats_object(
            reducedLabelMapT1LabelStats
        )
        targetLabels = self.remove_outside_value_from_target_labels(
            targetLabels, outsideValue
        )
        return targetLabels

    def remove_outside_value_from_target_labels(self, targetLabels, outsideValue):
        """
        This function...
        :param targetLabels:
        :param outsideValue:
        :return: targetLabels
        """
        if outsideValue in targetLabels:
            targetLabels.remove(outsideValue)
        return targetLabels

    def dilate_label_map(self, inputLabelImage, kernelRadius):
        """
        This function...
        :param inputLabelImage:
        :param kernelRadius:
        :return: castedOutput
        """
        myFilter = sitk.BinaryDilateImageFilter()
        myFilter.SetBackgroundValue(0.0)
        myFilter.SetBoundaryToForeground(False)
        myFilter.SetDebug(False)
        myFilter.SetForegroundValue(1.0)
        myFilter.SetKernelRadius((kernelRadius, kernelRadius, kernelRadius))
        myFilter.SetKernelType(2)  # Kernel Type=Box
        # myFilter.SetNumberOfWorkUnits(8)
        output = myFilter.Execute(inputLabelImage)
        castedOutput = sitk.Cast(output, sitk.sitkInt16)

        return castedOutput

    def calculate_label_intensity_difference_value(
        self,
        averageT1IntensitySuspiciousLabel,
        averageT2IntensitySuspiciousLabel,
        targetLabels,
        inputT1VolumeImage,
        inputT2VolumeImage,
        inputLabelImage,
    ):
        """
        Calculates a measurement for each label that is on the border of the suspicious label.
        This value is the square root of the sum of the squared difference in the average T1
        intensity values and the squared difference in the average T2 intensity values of the
        two islands in the comparison. The calculated value for each border label will later be
        sorted in ascending order - meaning that the smallest value has the "closest" average
        intensity to the suspicious label.
        :param averageT1IntensitySuspiciousLabel:
        :param averageT2IntensitySuspiciousLabel:
        :param targetLabels:
        :param inputT1VolumeImage:
        :param inputT2VolumeImage:
        :param inputLabelImage:
        :return: squareRootDiffLabelDict
        """
        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        squareRootDiffLabelDict = OrderedDict()
        labelStatsT1WithInputLabelImage = self.get_label_stats_object(
            inputT1VolumeImage, inputLabelImage
        )
        if inputT2VolumeImage:
            labelStatsT2WithInputLabelImage = self.get_label_stats_object(
                inputT2VolumeImage, inputLabelImage
            )

        for targetLabel in targetLabels:
            averageT1IntensityTargetLabel = labelStatsT1WithInputLabelImage.GetMean(
                targetLabel
            )
            squareDiffAverageT1 = math.pow(
                averageT1IntensitySuspiciousLabel - averageT1IntensityTargetLabel, 2
            )
            if inputT2VolumeImage:
                averageT2IntensityTargetLabel = labelStatsT2WithInputLabelImage.GetMean(
                    targetLabel
                )
                squareDiffAverageT2 = math.pow(
                    averageT2IntensitySuspiciousLabel - averageT2IntensityTargetLabel, 2
                )
            else:
                squareDiffAverageT2 = 0
            squareRootDiff = math.sqrt(squareDiffAverageT1 + squareDiffAverageT2)

            squareRootDiffLabelDict[str(targetLabel)] = squareRootDiff

        return squareRootDiffLabelDict

    def relabel_image(self, labelImage, newRegion, newLabel):
        """
        This function...
        :param labelImage:
        :param newRegion:
        :param newLabel:
        :return: relabeledImage
        """
        castedLabelImage = sitk.Cast(labelImage, sitk.sitkInt16)
        castedNewRegion = sitk.Cast(newRegion, sitk.sitkInt16)
        negatedMask = sitk.BinaryNot(castedNewRegion)
        negatedImage = sitk.Mask(castedLabelImage, negatedMask)
        maskTimesNewLabel = sitk.Multiply(castedNewRegion, newLabel)
        relabeledImage = sitk.Add(negatedImage, maskTimesNewLabel)
        return relabeledImage

    def get_dict_keys_list_sorted_by_value(self, val):
        """
        This function...
        :param val:
        :return: sorted(val, key=val.get)
        """
        return sorted(val, key=val.get)


if __name__ == "__main__":
    from docopt import docopt

    arguments = docopt(__doc__)
    print()
    arguments
    print()
    "-" * 50
    Object = DustCleanup(arguments)
    Object.main()
