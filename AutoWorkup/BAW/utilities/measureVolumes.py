#!/usr/bin/env python

"""
measureVolumes.py
============================
Description:
    The purpose of this is to..

Usage:

"""


from builtins import range
from builtins import str

from utilities.distributed import modify_qsub_args
from utilities.misc import *


def make_label_dictionary(inputColorLookUpTableFilename):
    """
    Construct dictionary:
        label No.: label name

    :param inputColorLookUpTableFilename:
    :return:
    """
    # inputColorLookUpTableFilename="/Shared/johnsonhj/HDNI/ReferenceData/20150709_HDAdultAtlas/BAWHDAdultAtlas_FreeSurferConventionColorLUT_20150709.txt"
    # import csv
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    labelDictionary = OrderedDict()
    with open(inputColorLookUpTableFilename) as f:
        contents = f.readlines()
        for line in contents:
            currentline = line.split()
            # print(currentline)
            if len(currentline) > 0 and currentline[0].isdigit():
                # print( currentline )
                labelNo = int(currentline[0])
                labelName = currentline[1]
                labelDictionary[labelNo] = labelName
                # print(labelDictionary)
    return labelDictionary


"""
#Unit test:
inputColorLookUpTableFilename="/Shared/johnsonhj/HDNI/ReferenceData/20150709_HDAdultAtlas/BAWHDAdultAtlas_FreeSurferConventionColorLUT_20150709.txt"

labelDict=make_label_dictionary(inputColorLookUpTableFilename)
print(labelDict)
"""


def get_label_volumes(labelVolume, RefVolume, labelDictionary):
    """
    Get label volumes using
    1. reference volume and
    2. labeldictionary

    :param labelVolume:
    :param RefVolume:
    :param labelDictionary:
    :return:
    """
    import SimpleITK as sitk
    import os
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    labelImg = sitk.ReadImage(labelVolume, sitk.sitkInt64)
    RefImg = sitk.ReadImage(RefVolume, sitk.sitkFloat64)
    labelStatFilter = sitk.LabelStoutputatisticsImageFilter()
    labelStatFilter.Execute(RefImg, labelImg)
    ImageSpacing = RefImg.GetSpacing()

    outputLabelVolumes = list()
    for value in labelStatFilter.GetLabels():
        structVolume = (
            ImageSpacing[0]
            * ImageSpacing[1]
            * ImageSpacing[2]
            * labelStatFilter.GetCount(value)
        )
        labelVolDict = OrderedDict()
        labelVolDict["Volume_mm3"] = structVolume

        if value in list(labelDictionary.keys()):
            print(("{0} --> {1}".format(value, labelDictionary[value])))
            labelVolDict["LabelName"] = labelDictionary[value]
        else:
            print(("** Caution: {0} --> No name exists!".format(value)))
            labelVolDict["LabelName"] = "NA"
        labelVolDict["LabelCode"] = value
        labelVolDict["FileName"] = os.path.abspath(labelVolume)
        outputLabelVolumes.append(labelVolDict)
    return outputLabelVolumes


"""
#Unit test::
labelName="/Shared/sinapse/CACHE/20160405_PREDICTHD_long_Results/PHD_024/0138/49757/TissueClassify/JointFusion_HDAtlas20_2015_label.nii.gz"
t1Name="/Shared/sinapse/CACHE/20160405_PREDICTHD_long_Results/PHD_024/0138/49757/TissueClassify/t1_average_BRAINSABC.nii.gz"

labelVolDict= get_label_volumes(labelName, t1Name, labelDict)
print(labelVolDict)
"""


def write_dictionary_to_csv(inputList, outputFilename):
    """
    This function...

    :param inputList:
    :param outputFilename:
    :return:
    """
    import csv
    import os

    csvFile = open(outputFilename, "w")
    dWriter = csv.DictWriter(
        csvFile,
        ["LabelCode", "LabelName", "Volume_mm3", "FileName"],
        restval="",
        extrasaction="raise",
        dialect="excel",
    )
    dWriter.writeheader()
    for line in inputList:
        dWriter.writerow(line)
    return os.path.abspath(outputFilename)


"""
#Unit test::
write_dictionary_to_csv(labelVolDict, "~/Desktop/test.csv")

"""


def write_dictionary_to_json(inputList, outputFilename):
    """
    This function..

    :param inputList:
    :param outputFilename:
    :return:
    """
    import json

    with open(outputFilename, "w") as fp:
        json.dump(inputList, fp)
    import os

    outputFilename = os.path.abspath(outputFilename)
    return outputFilename


"""
#Unit test::
write_dictionary_to_json(labelVolDict, "~/Desktop/test.json")
"""


def volume_measure(
    inputColorLookUpTableFilename,
    labelFilename,
    inputReferenceFilename,
    outputFileBasename,
):
    """
    This function...

    :param inputColorLookUpTableFilename:
    :param labelFilename:
    :param inputReferenceFilename:
    :param outputFileBasename:
    :return:
    """
    labelDict = make_label_dictionary(inputColorLookUpTableFilename)
    measurementsList = get_label_volumes(labelFilename, inputReferenceFilename, labelDict)
    csvFilename = outputFileBasename + "CSV.csv"
    jsonFilename = outputFileBasename + "JSON.json"
    write_dictionary_to_csv(measurementsList, csvFilename)
    write_dictionary_to_json(measurementsList, jsonFilename)

    import os

    csvFilename = os.path.abspath(csvFilename)
    jsonFilename = os.path.abspath(jsonFilename)
    return (csvFilename, jsonFilename)


import sys
import getopt


def main():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "hc:l:r:o:",
            [
                "help",
                "colorTable=",
                "labelFilename=",
                "referenceFilename=",
                "outputFileBasename=",
            ],
        )
    except getopt.GetoptError as err:
        print((str(err)))
        print(
            "WorkupComputeLabelVolume.py -c <colorTable> -l <labelFilename>  -r <referenceFilename> -o <outputFileBasename>"
        )
        sys.exit(2)
    colorTable = ""
    labelFilename = ""
    outputFileBasename = ""
    referenceFilename = ""
    for opt, arg in opts:
        if opt == "-h":
            print(
                "WorkupComputeLabelVolume.py -c <colorTable> -l <labelFilename>  -r <referenceFilename> -o <outputFileBasename>"
            )
            sys.ext()
        elif opt in ("-c", "--colorTable"):
            colorTable = arg
        elif opt in ("-l", "--labelFilename"):
            labelFilename = arg
        elif opt in ("-r", "--referenceFilename"):
            referenceFilename = arg
        elif opt in ("-o", "--outputFileBasename"):
            outputFileBasename = arg

    if colorTable and labelFilename and outputFileBasename:
        print(
            (
                """ Arguments:
        color table: {0}
        labelFile: {1}
        referenceFilename: {2}
        outputFileBasename: {3}""".format(
                    colorTable, labelFilename, referenceFilename, outputFileBasename
                )
            )
        )

        outputFiles = volume_measure(
            colorTable, labelFilename, referenceFilename, outputFileBasename
        )
        print(outputFiles)
    else:
        print(
            "WorkupComputeLabelVolume.py -c <colorTable> -l <labelFilename> -o <outputFileBasename>"
        )


if __name__ == "__main__":
    main()
