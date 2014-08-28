#!/usr/bin/env python
usage="""
   Usage::
   python WorkupAddsonBrainStem.py --landmarkFilename --brainStemFilename --tissueLabelFilename
"""

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine



########################################################
# brainstem computation from white matter mask
########################################################
def brainStem( tissueLabelFilename,
               landmarkFilename,
               brainStemFilename):
    import csv
    import os
    import SimpleITK as sitk
    myLandmark=dict()
    with open(landmarkFilename, 'r') as myLandmarkFile:
        landmark = myLandmarkFile.readlines()[8:]
        for row in landmark:
            lmk = row.split(',')
            lmkID = lmk[0]
            lmkX = float(lmk[1])
            lmkY = float(lmk[2])
            lmkZ = float(lmk[3])
            #print (lmkID, lmkX, lmkY, lmkZ)
            #RAS to LPS
            myLandmark[lmkID] = (-lmkX,-lmkY,lmkZ)

    ## read a label file 'TissueClassify/fixed_brainlabels_seg.nii.gz'
    brainLbl = sitk.ReadImage(tissueLabelFilename)

    ## get index locations
    myLmkIndex = dict()
    for lmk in myLandmark.keys():
        myLmkIndex[ lmk ] = brainLbl.TransformPhysicalPointToIndex( myLandmark[lmk] )

    imageSize=brainLbl.GetSize()
    print( imageSize)

    #print ("imageSize : "+str(imageSize))
    #cropLower = (myLmkIndex['lat_right'][0],myLmkIndex['mid_lat'][1],0)
    cropLower = (myLmkIndex['lat_right'][0],
                 myLmkIndex['mid_lat'][1],
                 myLmkIndex['dens_axis'][2])
    #cropUpper   = (myLmkIndex['lat_left'][0],myLmkIndex['mid_prim_sup'][1],myLmkIndex['PC'][2])
    cropUpper   = (imageSize[0]-myLmkIndex['lat_left'][0],
                   imageSize[1]-myLmkIndex['mid_prim_sup'][1],
                   imageSize[2]-myLmkIndex['PC'][2])

    print ("cropLower : "+str(cropLower))
    print ("cropUpper : "+str(cropUpper))
    brainStem_area = sitk.Crop(brainLbl, cropLower, cropUpper)
    #print brainStem_area.GetSize()

    brainStem = sitk.BinaryThreshold( brainStem_area, 1,1)
    #sitk.WriteImage(brainStem,"./brainStem.nii.gz")

    brainStem_connected=sitk.ConnectedComponent(brainStem)
    brainStem_largest_connected=sitk.BinaryThreshold(brainStem_connected,1,1)
    #sitk.WriteImage(brainStem_connected,"./brainStem_connected.nii.gz")
    print brainStemFilename
    sitk.WriteImage(brainStem_largest_connected,brainStemFilename)

    outputBrainstemFilename = os.path.abspath( brainStemFilename )
    return outputBrainstemFilename

def CreateBrainstemWorkflow(WFname, CLUSTER_QUEUE, outputFilename):
    brainstemWF = pe.Workflow( name=WFname )

    inputsSpec = pe.Node( interface=IdentityInterface(fields=['inputTissueLabelFilename',
                                                              'inputLandmarkFilename'] ),
                          run_without_submitting=True,
                          name='inputspec')
    outputSpec = pe.Node( interface=IdentityInterface(fields=['outputBrainstemFilename']),
                          run_without_submitting=True,
                          name='outputspec')

    generateBrainStemNode = pe.Node( Function( function=brainStem,
                                               input_names=['tissueLabelFilename',
                                                            'landmarkFilename',
                                                            'brainStemFilename'],
                                               output_names=['outputBrainstemFilename']),
                                     run_without_submitting=False,
                                     name='brainStem' )

    brainstemWF.connect( inputsSpec, 'inputTissueLabelFilename',
                         generateBrainStemNode, 'tissueLabelFilename' )
    brainstemWF.connect( inputsSpec, 'inputLandmarkFilename',
                         generateBrainStemNode, 'landmarkFilename' )
    generateBrainStemNode.inputs.brainStemFilename = outputFilename

    brainstemWF.connect( generateBrainStemNode, 'outputBrainstemFilename',
                         outputSpec, 'outputBrainstemFilename')

    return brainstemWF

import sys
import getopt
def main(argv=None):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hlbt:v", ["help",
                                                            "landmarkFilename=",
                                                            "brainStemFilename=",
                                                            "tissueLabelFilename="])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            print usage
            sys.exit(0)
        elif o == '--landmarkFilename':
            lmkFilename = a
        elif o == '--brainStemFilename':
            bsFilename = a
        elif o == '--tissueLabelFilename':
            tlFilename = a
    print lmkFilename,bsFilename,tlFilename

    brainStem( tlFilename, lmkFilename, bsFilename )

if __name__ == "__main__":
    sys.exit(main())
