"""
BRAINSTalairachext.py
======

This script creates a landmarks fcsv file (that includes AC, PC, SLA and IRP points) to be passed to the input of BRAINSTalairach.

This scripts gets an input (T1/T2) volume and a label map image to find coordinates of SLA and IRP points. These two points are used by BRAINSTalairach to find a box around "cerebrum" part of the human brain.
Also, this script uses the output landmarks file of BCD to find the coordinates of AC and PC points.

USAGE:

  BRAINSTalairachext.py \
  --inputVolume <input T1/T2 Volume> \
  --inputLabelsImage <input label map image> \
  --inputLandmarksFile <input landmarks file (output of BCD)> \
  --outputTalairachLandmarksFile <output landmarks file (including AC,PC,SLA,IRP)>
"""
from __future__ import print_function

import sys, getopt
import os.path
import csv

import SimpleITK as sitk
print(sitk.Version())

def csv_file_writer(filename,data):
  with open(filename, 'w') as lf:
    headerdata1 = [['#Fiducial', 'List', 'file', filename],
                   ['#numPoints', '=', len(data)]]
    headerdata2 = [['#symbolScale = 5'],
                  ['#visibility = 1'],
                  ['#textScale = 4.5'],
                  ['#color = 0.4','1','1'],
                  ['#selectedColor = 1','0.5','0.5'],
                  ['#label','x','y','z','sel','vis']]
    wr = csv.writer(lf, delimiter=' ')
    wr.writerows(headerdata1)
    wr = csv.writer(lf, delimiter=',')
    wr.writerows(headerdata2)
    wr.writerows(data)

def csv_file_reader(filename,dataList):
  import csv
  with open(filename) as lf:
    reader = csv.reader(lf, delimiter=',')
    for line in reader:
      lmkName = line[0]
      if lmkName == "AC":
        dataList.append(line)
      elif lmkName == "PC":
        dataList.append(line)
      else:
        continue

def main(argv):
  inputVolume = ''
  inputLabelsImage = ''
  inputLandmarksFile = ''
  outputTalairachLandmarksFile = ''

  try:
    opts, args = getopt.getopt(argv,"hi:m:l:o:",["inputVolume=","inputLabelsImage=","inputLandmarksFile=","outputTalairachLandmarksFile="])
  except getopt.GetoptError:
    print('BRAINSTalairachext.py -i <inputVolume> -m <inputLabelsImage> -l <inputLandmarksFile> -o <outputTalairachLandmarksFile>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print('BRAINSTalairachext.py -i <inputVolume> -m <inputLabelsImage> -l <inputLandmarksFile> -o <outputTalairachLandmarksFile>')
      sys.exit()
    elif opt in ("-i", "--inputVolume"):
      inputVolume = arg
    elif opt in ("-m", "--inputLabelsImage"):
      inputLabelsImage = arg
    elif opt in ("-l", "--inputLandmarksFile"):
      inputLandmarksFile = arg
    elif opt in ("-o", "--outputTalairachLandmarksFile"):
      outputTalairachLandmarksFile = arg

  print('Input Volume is "', inputVolume)
  print('Input Labels Image is "', inputLabelsImage)
  print('Input Landmarks File is "', inputLandmarksFile)
  print('Output Talairach Landmarks File is "', outputTalairachLandmarksFile)

  input_img = sitk.ReadImage(inputVolume.encode('ascii','replace'))
  img_labels = sitk.ReadImage(inputLabelsImage.encode('ascii','replace'))

  exclusionLabels=((img_labels == 11) +
                   (img_labels == 35) +
                   (img_labels == 38) +
                   (img_labels == 39) +
                   (img_labels == 40) +
                   (img_labels == 41) +
                   (img_labels == 51) +
                   (img_labels == 52) +
                   (img_labels == 71) +
                   (img_labels == 72) +
                   (img_labels == 73) +
                   (img_labels == 230)+
                   (img_labels == 255))

  important_labels = img_labels*(1-exclusionLabels)
  unified_important_labels = important_labels>0

  # debug
  #sitk.WriteImage(unified_important_labels,'./unified_important_labels.nrrd'.encode('ascii','replace'))
  ##

  labelStatFilter = sitk.LabelStatisticsImageFilter()
  labelStatFilter.Execute(input_img, unified_important_labels)
  bbox = labelStatFilter.GetBoundingBox(1)

  I=bbox[4]
  R=bbox[0]
  P=bbox[3]
  indx_IRP=[R,P,I]

  S=bbox[5]
  L=bbox[1]
  A=bbox[2]
  indx_SLA=[L,A,S]

  print(("IRP indeces: ",indx_IRP))
  print(("SLA indeces: ",indx_SLA))

  itk_IRP = unified_important_labels.TransformIndexToPhysicalPoint(indx_IRP)
  itk_SLA = unified_important_labels.TransformIndexToPhysicalPoint(indx_SLA)

  IRP=[-itk_IRP[0],-itk_IRP[1],itk_IRP[2]]
  SLA=[-itk_SLA[0],-itk_SLA[1],itk_SLA[2]]

  print(("IRP: ",IRP))
  print(("SLA: ",SLA))

  data=[['SLA', SLA[0], SLA[1], SLA[2],1,1],
        ['IRP', IRP[0], IRP[1], IRP[2],1,1]]

  csv_file_reader(inputLandmarksFile,data)
  csv_file_writer(outputTalairachLandmarksFile,data)

if __name__ == "__main__":
  main(sys.argv[1:])
