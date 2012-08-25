## NOTE: THERE CAN NOT BE ANY GLOBAL imports in this file
##       NIPYPE pipeline functions must be self contained
##       and any import needed for a function must be
##       included in the function itself.

## This file contains misc SimpleITK based functions for use in nipype
## nodes.

## AVOID REFORMATTING THIS FILE, it causes the hash to change in
## nipype and that require re-running the function.


def ClipT1ImageWithBrainMask(t1_image,brain_labels,clipped_file_name):
    import os
    import sys
    import SimpleITK as sitk
    ## Now clean up the posteriors based on anatomical knowlege.
    ## sometimes the posteriors are not relevant for priors
    ## due to anomolies around the edges.
    t1=sitk.Cast(sitk.ReadImage(t1_image),sitk.sitkFloat32)
    bl=sitk.Cast(sitk.ReadImage(brain_labels),sitk.sitkFloat32)
    clipped=t1*bl
    sitk.WriteImage(clipped,clipped_file_name)
    clipped_file=os.path.realpath(clipped_file_name)
    return clipped_file
