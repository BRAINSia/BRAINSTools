#! /usr/bin/env python
"""
computeIPD.py
==============

This program takes a list of input fcsv files and computes inter-pupillary distance (IPD) for each input fcsv file. Then, it places the results into a CSV file with 2 columns:  sessionID,IPD.

Usage:
  computeIPD.py --inputFilesList INPUTLIST --outputIPDsList IPDsLIST
  computeIPD.py -v | --version
  computeIPD.py -h | --help

Options:
  -h --help                     Show this help and exit
  -v --version                  Print the version and exit
  --inputFilesList INPUTLIST    List of input fcsv files
  --outputIPDsList IPDsLIST     List of output IPDs

Example:
  computeIPD.py --inputFilesList TRACKDataFCSVList.txt --outputIPDsList TRACKIPDs.txt
  computeIPD.py --inputFilesList PREDICTDataFCSVList.txt --outputIPDsList PREDICTIPDs.txt
"""


def csv_file_reader(fcsvFile, dataList):
  import csv
  import numpy
  sessionID = os.path.basename(os.path.dirname(os.path.dirname(fcsvFile)))
  with open(fcsvFile) as lf:
       reader = csv.reader(lf, delimiter=',')
       for line in reader:
         lmkName = line[0]
         if lmkName == "RE":
           REpoint=numpy.array((float(line[1]), float(line[2]), float(line[3])))
         elif lmkName == "LE":
           LEpoint=numpy.array((float(line[1]), float(line[2]), float(line[3])))
         else:
           continue
  IPD = numpy.linalg.norm(REpoint - LEpoint)
  dataList.append([sessionID, IPD])

def csv_file_writer(outputCSVFile, data):
  with open(outputCSVFile, 'w') as lf:
      headerdata = [['#sessionID', 'IPD']]
      wr = csv.writer(lf, delimiter=',')
      wr.writerows(headerdata)
      wr.writerows(data)

if __name__ == '__main__':
  import os
  import glob
  import sys
  import csv

  from docopt import docopt
  argv = docopt(__doc__, version='1.0')
  print(argv)

  INPUTLIST = argv['--inputFilesList']
  assert os.path.exists(INPUTLIST), "Input files list is not found: %s" % INPUTLIST

  IPDsLIST = argv['--outputIPDsList']

  print(('=' * 100))

  dataList=[]
  with open(INPUTLIST) as lf:
      reader = csv.reader(lf)
      for line in reader:
        fcsvFile = line[0]
        csv_file_reader(fcsvFile, dataList)

  csv_file_writer(IPDsLIST, dataList)

  sys.exit(0)
