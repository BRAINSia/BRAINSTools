#!/usr/bin/env python
from __future__ import print_function
usage = """
   Usage::
   python WorkupAddsonBrainStem.py --landmarkFilename --brainStemFilename --tissueLabelFilename
"""

from nipype.interfaces.utility import Function, IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine


# #######################################################
# brainstem computation from white matter mask
########################################################
def brainStem(tissueLabelFilename,
              landmarkFilename,
              brainStemFilename,
              ouputTissuelLabelFilename):
    import os
    import SimpleITK as sitk

    def cropAndResampleInPlace(inputBrainLabelFilename,
                               lower, upper, thresholdUpper, thresholdLower,
                               outputImageFilename):
        brainLbl = sitk.ReadImage(inputBrainLabelFilename)
        brainStem_area = sitk.Crop(brainLbl, lower, upper)

        # HACK: WM label value should be given
        brainStem = sitk.BinaryThreshold(brainStem_area, thresholdUpper, thresholdLower)
        brainStem_connected = sitk.ConnectedComponent(brainStem)
        brainStem_largest_connected = sitk.BinaryThreshold(brainStem_connected, 1, 1)

        ## Fill Hole
        radius = 3
        kernel = sitk.sitkBox
        foregroundValue = 1
        brainStem_fillHole = sitk.BinaryFillhole(brainStem_largest_connected > 0, True, foregroundValue)
        ## Resample In Place
        rs = sitk.ResampleImageFilter()
        rs.SetInterpolator(sitk.sitkNearestNeighbor)
        rs.SetTransform(sitk.Transform(3, sitk.sitkIdentity))
        rs.SetReferenceImage(brainLbl)
        brainStemInPlace = rs.Execute(brainStem_fillHole)

        sitk.WriteImage(brainStemInPlace, outputImageFilename)

        import os

        return os.path.abspath(outputImageFilename)


    myLandmark = dict()
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
            myLandmark[lmkID] = (-lmkX, -lmkY, lmkZ)

    ## read a label file 'TissueClassify/fixed_brainlabels_seg.nii.gz'
    brainLbl = sitk.ReadImage(tissueLabelFilename)

    """
    brain stem
    """
    ## get index locations
    myLmkIndex = dict()
    for lmk in myLandmark.keys():
        myLmkIndex[lmk] = brainLbl.TransformPhysicalPointToIndex(myLandmark[lmk])

    imageSize = brainLbl.GetSize()
    cropLower = [myLmkIndex['lat_right'][0],
                 myLmkIndex['mid_lat'][1],
                 myLmkIndex['dens_axis'][2]]
    cropLower = [int(x) for x in cropLower]
    cropUpper = [imageSize[0] - myLmkIndex['lat_left'][0],
                 imageSize[1] - myLmkIndex['mid_prim_sup'][1],
                 imageSize[2] - myLmkIndex['PC'][2]]
    cropUpper = [int(x) for x in cropUpper]
    wmLabelNo = 1
    brainStem = cropAndResampleInPlace(tissueLabelFilename, cropLower, cropUpper, wmLabelNo, wmLabelNo,
                                       brainStemFilename)

    """
    below the brain stem
    """
    cropUpper = [0, 0, 0]
    cropLower = [0, 0, myLmkIndex['dens_axis'][2]]

    brainStemExtraBottom = cropAndResampleInPlace(tissueLabelFilename, cropLower, cropUpper, 0, 255,
                                                  brainStemFilename + "_InValid.nii.gz")

    brainStemBinary = ( sitk.ReadImage(brainStem) > 0 )
    brainStemExtraBottomBinary = 1 - (sitk.ReadImage(brainStemExtraBottom) > 0)

    outputTissueLabel = brainLbl * (1 - (brainStemBinary > 0) ) * (
    1 - (brainStemExtraBottomBinary > 0)) + brainStemBinary * 30

    errod_brain_mask = sitk.ErodeObjectMorphology(( outputTissueLabel > 0 ), 2)
    LargestComponentCode = 1
    one_region_mask = ( sitk.RelabelComponent(sitk.ConnectedComponent(errod_brain_mask)) == LargestComponentCode )
    dilate_one_region = ( sitk.DilateObjectMorphology(one_region_mask, 3) > 0 )
    cleanedOutputTissueLabel = outputTissueLabel * dilate_one_region

    full_output_path = os.path.abspath(ouputTissuelLabelFilename)
    sitk.WriteImage(cleanedOutputTissueLabel, full_output_path)
    return full_output_path


def CreateBrainstemWorkflow(WFname, CLUSTER_QUEUE, outputFilename):
    brainstemWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['inputTissueLabelFilename',
                                                             'inputLandmarkFilename']),
                         run_without_submitting=True,
                         name='inputspec')
    outputSpec = pe.Node(interface=IdentityInterface(fields=['ouputTissuelLabelFilename']),
                         run_without_submitting=True,
                         name='outputspec')

    generateBrainStemNode = pe.Node(Function(function=brainStem,
                                             input_names=['tissueLabelFilename',
                                                          'landmarkFilename',
                                                          'brainStemFilename',
                                                          'ouputTissuelLabelFilename'],
                                             output_names=['ouputTissuelLabelFilename']),
                                    run_without_submitting=False,
                                    name='brainStem')

    brainstemWF.connect(inputsSpec, 'inputTissueLabelFilename',
                        generateBrainStemNode, 'tissueLabelFilename')
    brainstemWF.connect(inputsSpec, 'inputLandmarkFilename',
                        generateBrainStemNode, 'landmarkFilename')
    generateBrainStemNode.inputs.brainStemFilename = outputFilename + "_brainStem.nii.gz"
    generateBrainStemNode.inputs.ouputTissuelLabelFilename = outputFilename

    brainstemWF.connect(generateBrainStemNode, 'ouputTissuelLabelFilename',
                        outputSpec, 'ouputTissuelLabelFilename')

    return brainstemWF


import sys
import getopt


def main(argv=None):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hlbt:v", ["help",
                                                            "landmarkFilename=",
                                                            "brainStemFilename=",
                                                            "tissueLabelFilename="])
    except getopt.error as msg:
        print(msg)
        print("for help use --help")
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            print(usage)
            sys.exit(0)
        elif o == '--landmarkFilename':
            lmkFilename = a
        elif o == '--brainStemFilename':
            bsFilename = a
        elif o == '--tissueLabelFilename':
            tlFilename = a
    print(lmkFilename, bsFilename, tlFilename)

    brainStem(tlFilename, lmkFilename, bsFilename)


if __name__ == "__main__":
    sys.exit(main())
