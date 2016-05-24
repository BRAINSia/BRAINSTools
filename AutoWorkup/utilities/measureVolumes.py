#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range

from utilities.misc import *
from utilities.distributed import modify_qsub_args

def MakeLabelDictionary(inputColorLookUpTableFilename):
    """
    Construct dictionary:
        label No.: label name
    """
    #inputColorLookUpTableFilename="/Shared/johnsonhj/HDNI/ReferenceData/20150709_HDAdultAtlas/BAWHDAdultAtlas_FreeSurferConventionColorLUT_20150709.txt"
    import csv
    labelDictionary=dict()
    with open(inputColorLookUpTableFilename) as f:
        contents = f.readlines()
        for line in contents:
            currentline=line.split()
            #print(currentline)
            if (len(currentline)>0 and currentline[0].isdigit()):
                #print currentline
                labelNo=int(currentline[0])
                labelName=currentline[1]
                labelDictionary[labelNo] = labelName
    return labelDictionary
"""
#Unit test:
inputColorLookUpTableFilename="/Shared/johnsonhj/HDNI/ReferenceData/20150709_HDAdultAtlas/BAWHDAdultAtlas_FreeSurferConventionColorLUT_20150709.txt"

labelDict=MakeLabelDictionary(inputColorLookUpTableFilename)
print(labelDict)
"""


def GetLabelVolumes(labelVolume, RefVolume, labelDictionary):
    """
    Get label volumes using
    1. reference volume and
    2. labeldictionary
    """
    import SimpleITK as sitk
    import os
    labelImg = sitk.ReadImage(labelVolume, sitk.sitkInt64)
    RefImg = sitk.ReadImage(RefVolume, sitk.sitkFloat64)
    labelStatFilter = sitk.LabelStatisticsImageFilter()
    labelStatFilter.Execute(RefImg, labelImg)
    ImageSpacing = RefImg.GetSpacing()

    outputLabelVolumes=list()
    for value in labelStatFilter.GetLabels():
        structVolume = ImageSpacing[0] * ImageSpacing[1] * ImageSpacing[2] * labelStatFilter.GetCount(value)
        labelVolDict=dict()
        labelVolDict['Volume_mm3'] = structVolume

        if value in labelDictionary.keys():
            print("{0} --> {1}".format(value,labelDictionary[value]))
            labelVolDict['LabelName']=labelDictionary[value]
        else:
            print("** Caution: {0} --> No name exists!".format(value))
            labelVolDict['LabelName']='NA'
        labelVolDict['LabelCode'] = value
        labelVolDict['FileName']=os.path.abspath(labelVolume)
        outputLabelVolumes.append(labelVolDict)
    return outputLabelVolumes
"""
#Unit test::
labelName="/Shared/sinapse/CACHE/20160405_PREDICTHD_long_Results/PHD_024/0138/49757/TissueClassify/JointFusion_HDAtlas20_2015_label.nii.gz"
t1Name="/Shared/sinapse/CACHE/20160405_PREDICTHD_long_Results/PHD_024/0138/49757/TissueClassify/t1_average_BRAINSABC.nii.gz"

labelVolDict= GetLabelVolumes(labelName, t1Name, labelDict)
print(labelVolDict)
"""

def WriteDictionaryToCSV(inputList, outputFilename):
    import csv
    import os
    csvFile=open(outputFilename, 'w')
    dWriter=csv.DictWriter(csvFile,
                           ['LabelCode','LabelName','Volume_mm3','FileName'],
                          restval='',
                          extrasaction='raise',
                          dialect='excel')
    dWriter.writeheader()
    for line in inputList:
        dWriter.writerow(line)
    return os.path.abspath(outputFilename)

"""
#Unit test::
WriteDictionaryToCSV(labelVolDict, "~/Desktop/test.csv")
"""


def WriteDictionaryToJson(inputList, outputFilename):
    import json
    with open(outputFilename, 'w') as fp:
        json.dump(inputList, fp)
    import os
    outputFilename = os.path.abspath(outputFilename)
    return (outputFilename)

"""
#Unit test::
WriteDictionaryToJson(labelVolDict, "~/Desktop/test.json")
"""

def VolumeMeasure(inputColorLookUpTableFilename,
                  labelFilename,
                  inputReferenceFilename,
                  outputFileBasename):
    labelDict = MakeLabelDictionary(inputColorLookUpTableFilename)
    measurementsList = GetLabelVolumes(labelFilename, inputReferenceFilename, labelDict)
    csvFilename=outputFileBasename+"CSV.csv"
    jsonFilename=outputFileBasename+"JSON.json"
    WriteDictionaryToCSV( measurementsList, csvFilename)
    WriteDictionaryToJson( measurementsList, jsonFilename)

    import os
    csvFilename = os.path.abspath( csvFilename )
    jsonFilename = os.path.abspath( jsonFilename )
    return( csvFilename , jsonFilename )

import sys
import getopt
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:l:r:o:",
                                   ["help",
                                    "colorTable=",
                                    "labelFilename=",
                                    "referenceFilename=",
                                    "outputFileBasename="])
    except getopt.GetoptError as err:
        print (str(err))
        print ("WorkupComputeLabelVolume.py -c <colorTable> -l <labelFilename>  -r <referenceFilename> -o <outputFileBasename>")
        sys.exit(2)
    colorTable = ""
    labelFilename = ""
    outputFileBasename = ""
    referenceFilename = ""
    for opt, arg in opts:
        if opt == '-h':
            print ("WorkupComputeLabelVolume.py -c <colorTable> -l <labelFilename>  -r <referenceFilename> -o <outputFileBasename>")
            sys.ext()
        elif opt in ("-c","--colorTable"):
            colorTable=arg
        elif opt in ("-l","--labelFilename"):
            labelFilename=arg
        elif opt in ("-r","--referenceFilename"):
            referenceFilename=arg
        elif opt in ("-o","--outputFileBasename"):
            outputFileBasename=arg

    if( colorTable and labelFilename and outputFileBasename ):
        print (""" Arguments:
        color table: {0}
        labelFile: {1}
        referenceFilename: {2}
        outputFileBasename: {3}""".format(colorTable, labelFilename, referenceFilename, outputFileBasename))

        outputFiles=VolumeMeasure( colorTable, labelFilename, referenceFilename, outputFileBasename )
        print (outputFiles)
    else:
        print ("WorkupComputeLabelVolume.py -c <colorTable> -l <labelFilename> -o <outputFileBasename>")

if __name__ == "__main__":
    main()
