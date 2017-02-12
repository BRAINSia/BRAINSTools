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
                               physBB1, physBB2, thresholdUpper, thresholdLower,
                               outputImageFilename):

        brainLbl = sitk.ReadImage(inputBrainLabelFilename.encode('ascii','replace'))

        roiBBStart_index = brainLbl.TransformPhysicalPointToIndex(physBB1)
        roiBBStop_index  = brainLbl.TransformPhysicalPointToIndex(physBB2)
        del(physBB1)
        del(physBB2)
        roiBBStart_final_index =[ min(roiBBStart_index[i],roiBBStop_index[i]) for i in range(0,3) ]
        roiBBStop_final_index = [ max(roiBBStart_index[i],roiBBStop_index[i]) for i in range(0,3) ]
        roiBBSize_final_index = [ roiBBStop_final_index[i] - roiBBStart_final_index[i] for i in range(0,3) ]
        del(roiBBStart_index)
        del(roiBBStop_index)

#        print( "XX"*30)
#        print( brainLbl.GetSize() )
#        print( roiBBStart_final_index )
#        print( roiBBSize_final_index )

        brainStem_area = sitk.RegionOfInterest(brainLbl,roiBBSize_final_index,roiBBStart_final_index)

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

        sitk.WriteImage(brainStemInPlace, outputImageFilename.encode('ascii','replace'))

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
            #RAS to LPS
            myLandmark[lmkID] = (-lmkX, -lmkY, lmkZ)

    """
    brain stem
    """
    roiBBStart = [myLandmark['lat_right'][0],
                 myLandmark['mid_lat'][1],
                 myLandmark['dens_axis'][2]]

    roiBBStop = [myLandmark['lat_left'][0],
                 myLandmark['mid_prim_sup'][1],
                 myLandmark['PC'][2]]

    wmLabelNo = 1
    brainStem = cropAndResampleInPlace(tissueLabelFilename, roiBBStart, roiBBStop, wmLabelNo, wmLabelNo,
                                       brainStemFilename)

    """
    below the brain stem
    """
    ## read a label file 'TissueClassify/fixed_brainlabels_seg.nii.gz'
    brainLbl = sitk.ReadImage(tissueLabelFilename.encode('ascii','replace'))
    brainLblMaxIndex = [ x-1 for x in brainLbl.GetSize() ]
    fov1 = brainLbl.TransformIndexToPhysicalPoint([0,0,0])
    fov2 = brainLbl.TransformIndexToPhysicalPoint(brainLblMaxIndex)

    roiBBStart = [min(fov1[0],fov2[0]), min(fov1[1],fov2[1]) , min(fov1[2],fov2[2]) ]
    roiBBStop  =  [max(fov1[0],fov2[0]), max(fov1[1],fov2[1]) , myLandmark['dens_axis'][2] ]

    noExtraBottomBrainStem = cropAndResampleInPlace(tissueLabelFilename, roiBBStart, roiBBStop, 0, 255,
                                                  brainStemFilename + "_InValid.nii.gz")

    brainStemBinary = ( sitk.ReadImage(brainStem.encode('ascii','replace')) > 0 )
    noExtraBottomBrainStemBinary = (sitk.ReadImage(noExtraBottomBrainStem.encode('ascii','replace')) > 0)

    outputTissueLabel = brainLbl * (1 - (brainStemBinary > 0) ) * (
    1 - (noExtraBottomBrainStemBinary > 0)) + brainStemBinary * 30

    errod_brain_mask = sitk.ErodeObjectMorphology(( outputTissueLabel > 0 ), 2)
    LargestComponentCode = 1
    one_region_mask = ( sitk.RelabelComponent(sitk.ConnectedComponent(errod_brain_mask)) == LargestComponentCode )
    dilate_one_region = ( sitk.DilateObjectMorphology(one_region_mask, 3) > 0 )
    cleanedOutputTissueLabel = outputTissueLabel * dilate_one_region

    full_output_path = os.path.abspath(ouputTissuelLabelFilename)
    sitk.WriteImage(cleanedOutputTissueLabel, full_output_path.encode('ascii','replace'))
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
