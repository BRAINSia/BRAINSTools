from __future__ import print_function


import sys, getopt
import os
import ast

import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.semtools import BRAINSSnapShotWriter



def print_usage():
    print("python ./BAWScreenShots.py \\\n"
          "-i  ./small_list.csv \\\n"
          "-d /Shared/sinapse/CACHE/20160202_PREDICTHD_base_Results/  \\\n"
          "-c /scratch/eunyokim/20160202_PREDICTHD_base_Snapshot/")
    print("BAWScreenShots.py\n"
          "              -i <inputfile>\n"
          "              -d <inputDirectory>\n"
          "              -c <cacheDirectory>\n"
          "              -s <inputSubDirectory>")

def readInputFile(inputFilename):
    inputList = []
    with open(inputFilename) as infile:
        for line in infile:
            sessionInfo = ast.literal_eval(line)
            # print sessionInfo
            if type(sessionInfo) is list:
                pass
            else:
                sessionDict = dict()
                sessionDict['project'] = sessionInfo[0]
                sessionDict['subject'] = sessionInfo[1]
                sessionDict['session'] = sessionInfo[2]
                sessionDict['imagefiles'] = ast.literal_eval(sessionInfo[3])

            inputList.append(sessionDict)
    """
    printing for debugging
    """
    # for line in inputList:
    #  for key,value in line.items():
    #    # print key,value
    #    if key == 'imagefiles':
    #      for key,value in sessionDict['imagefiles'] .items():
    #         print key,value
    return inputList


def main(argv=None):
    inputfile = ''
    inputDirectory = ''
    inputSubDirectory = 'TissueClassify'
    searchVolumeList = ['t1_average_BRAINSABC.nii.gz',
                        't2_average_BRAINSABC.nii.gz']
    searchBinaryVolumeList = ['JointFusion_HDAtlas20_2015_fs_standard_label.nii.gz']
    try:
        opts, _ = getopt.getopt(argv, "hi:d:c:s:", ["inputFilename=",
                                                       "inputDirectory=",
                                                       "cacheDirectory=",
                                                       "inputSubDirectory="])
    except getopt.GetoptError:
        print_usage()

        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit()
        elif opt in ("-i", "--inputFilename"):
            inputfile = arg
        elif opt in ("-d", "--inputDirectory"):
            inputDirectory = arg
        elif opt in ("-c", "--cacheDirectory"):
            cacheDirectory = arg
        elif opt in ("-s", "--inputSubDirectory"):
            inputSubDirectory = arg

    print ('Input file is {0}'.format(inputfile))

    inputDictList = readInputFile(inputfile)
    snapShotWF = pe.Workflow(name="BAWSnapshots")  # templage generate work flow
    snapShotWF.base_dir = os.path.abspath(cacheDirectory)

    for oneSession in inputDictList:
        sessionid = oneSession['session']
        sessionSubj = oneSession['subject']
        sessionProj = oneSession['project']
        sessionInputDir = os.path.join(inputDirectory,
                                       sessionProj, sessionSubj, sessionid,
                                       inputSubDirectory)
        anatomicalImageList = []
        print ("For session = {0} ".format(sessionid))
        for filename in searchVolumeList:
            fullPathFilename = os.path.join(sessionInputDir, filename)
            if os.path.isfile(fullPathFilename):
                anatomicalImageList.append(fullPathFilename)
            else:
                print
                "    {0} does not exist. Skip.".format(fullPathFilename)

        binaryImageList = []
        for filename in searchBinaryVolumeList:
            fullPathBinaryFilename = os.path.join(sessionInputDir, filename)
            if os.path.isfile(fullPathBinaryFilename):
                binaryImageList.append(fullPathBinaryFilename)
            else:
                print
                "    {0} does not exist. Skip.".format(fullPathBinaryFilename)

        if len(anatomicalImageList) > 0:
            print("  Generating SnapShot Writer for {0} ... ... ".format(sessionid))
            ## SnapShotWriter for Segmented result checking:
            SnapShotWriterNodeName = "SnapShotWriter_" + str(sessionid)
            SnapShotWriter = pe.Node(interface=BRAINSSnapShotWriter(), name=SnapShotWriterNodeName)

            SnapShotWriter.inputs.outputFilename = 'snapShot' + str(sessionid) + '.png'  # output specification
            SnapShotWriter.inputs.inputPlaneDirection = [2, 1, 1, 1, 1, 0, 0]
            SnapShotWriter.inputs.inputSliceToExtractInPhysicalPoint = [-3, -7, -3, 5, 7, 22, -22]
            SnapShotWriter.inputs.inputVolumes = anatomicalImageList
            if binaryImageList > 0:
                SnapShotWriter.inputs.inputBinaryVolumes = binaryImageList

            snapShotWF.add_nodes([SnapShotWriter])
    snapShotWF.run()

if __name__ == "__main__":
    main(sys.argv[1:])
