"""
FixLabelMapsTools.py
========================
Description:
    These are nodes used to fix-up label maps by merging information from various different sources.

Author:
    Hans J. Johnson
    
Usage:

"""


def fix_label_map_from_neuromorphemetrics_2012(
    fusionFN, FixedHeadFN, posteriorListOfTuples, LeftHemisphereFN, outFN, OUT_DICT
):
    """
    This funciton...

    :param fusionFN:
    :param FixedHeadFN:
    :param posteriorListOfTuples:
    :param LeftHemisphereFN:
    :param outFN:
    :param OUT_DICT:
    :return:
    """
    import SimpleITK as sitk
    import os
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    posterior_dict = OrderedDict(posteriorListOfTuples)

    def force_mask_insert(inlabels, newmask, newmaskvalue):
        """
        This function...

        :param inlabels:
        :param newmask:
        :param newmaskvalue:
        :return:
        """
        inlabels = sitk.Cast(inlabels, sitk.sitkUInt32)
        newmask = sitk.Cast((newmask > 0), sitk.sitkUInt32)
        outlabels = inlabels * sitk.Cast((1 - newmask), sitk.sitkUInt32)
        outlabels = outlabels + newmask * newmaskvalue
        return sitk.Cast(outlabels, sitk.sitkUInt32)

    ## TODO: get_largest_label is copied from elsewhere
    def get_largest_label(inputMask, UseErosionCleaning):
        """
        This function...

        :param inputMask:
        :param UseErosionCleaning:
        :return:
        """
        LargestComponentCode = 1
        if UseErosionCleaning:
            erosionMask = sitk.ErodeObjectMorphology(inputMask, 1)
        else:
            erosionMask = inputMask
        CC = sitk.ConnectedComponent(erosionMask)
        Rlabel = sitk.RelabelComponent(CC)
        largestMask = Rlabel == LargestComponentCode
        if UseErosionCleaning:
            dilateMask = sitk.DilateObjectMorphology(largestMask, 1)
        else:
            dilateMask = largestMask

        return largestMask * dilateMask > 0

    def recode_nonlargest(outlabels, keepCode, UNKNOWN_LABEL_CODE):
        """
        This function...

        :param outlabels:
        :param keepCode:
        :param UNKNOWN_LABEL_CODE:
        :return:
        """
        orig_mask = outlabels == keepCode
        connected_mask = get_largest_label(orig_mask, False)
        small_regions = orig_mask - connected_mask
        outlabels = force_mask_insert(outlabels, connected_mask, keepCode)
        outlabels = force_mask_insert(outlabels, small_regions, UNKNOWN_LABEL_CODE)
        return outlabels

    def minimize_size_of_image(outlabels):
        """This function will find the largest integer value in the labelmap, and
        cast the image to the smallest possible integer size so that no loss of data
        results.

        :param outlabels:
        :return:
        """

        measureFilt = sitk.StatisticsImageFilter()
        measureFilt.Execute(outlabels)
        imgMin = measureFilt.GetMinimum()
        imgMax = measureFilt.GetMaximum()
        if imgMax < (2 ** 8) - 1:
            outlabels = sitk.Cast(outlabels, sitk.sitkUInt8)
        elif imgMax < (2 ** 16) - 1:
            outlabels = sitk.Cast(outlabels, sitk.sitkUInt16)
        elif imgMax < (2 ** 32) - 1:
            outlabels = sitk.Cast(outlabels, sitk.sitkUInt32)
        elif imgMax < (2 ** 64) - 1:
            outlabels = sitk.Cast(outlabels, sitk.sitkUInt64)
        return outlabels

    fusionIm = sitk.Cast(sitk.ReadImage(fusionFN), sitk.sitkUInt32)
    FixedHead = sitk.Cast(sitk.ReadImage(FixedHeadFN), sitk.sitkUInt32)

    BRAINSABC_DICT = OrderedDict({"BRAINSTEM": 30, "CSF": 4, "BLOOD": 5})

    ## Intialize by cloning
    outlabels = sitk.Image(fusionIm)
    lbl_orig_mask = outlabels > 0
    lbl_outter_ring = sitk.BinaryDilate(lbl_orig_mask, 2) - lbl_orig_mask
    # DEBUG sitk.WriteImage(lbl_outter_ring,"/tmp/lbl_outter_ring.nii.gz")

    vb_post = sitk.ReadImage(posterior_dict["VB"])
    ring_vb = lbl_outter_ring * sitk.BinaryThreshold(
        vb_post, 0.5, 1.01, 1, 0
    )  # just outside mask
    # DEBUG sitk.WriteImage(ring_vb,"/tmp/ring_vb.nii.gz")
    inner_vb = lbl_orig_mask * sitk.BinaryThreshold(
        vb_post, 0.85, 1.01, 1, 0
    )  # inside mask, but very high probability
    # DEBUG sitk.WriteImage(inner_vb,"/tmp/inner_vb.nii.gz")
    # background = 0 , suspicous = 999
    ## Add blood from BRAINSABC to mask as as value OUT_DICT['BLOOD']
    blood_labels = (
        (FixedHead == BRAINSABC_DICT["BLOOD"]) * (outlabels == 0 | outlabels == 999)
        | ring_vb
        | inner_vb
    )
    # DEBUG sitk.WriteImage(blood_labels,"/tmp/blood_labels.nii.gz")
    outlabels = force_mask_insert(outlabels, blood_labels, OUT_DICT["BLOOD"])

    csf_post = sitk.ReadImage(posterior_dict["CSF"])
    ring_csf = lbl_outter_ring * sitk.BinaryThreshold(
        csf_post, 0.5, 1.01, 1, 0
    )  # just outside mask
    # DEBUG sitk.WriteImage(ring_csf,"/tmp/ring_csf.nii.gz")
    inner_csf = lbl_orig_mask * sitk.BinaryThreshold(
        csf_post, 0.85, 1.01, 1, 0
    )  # inside mask, but very high probability
    # DEBUG sitk.WriteImage(inner_csf,"/tmp/inner_csf.nii.gz")
    ## Add CSF from BRAINSABC to mask as as value OUT_DICT['RH_CSF']
    csf_labels = (
        (FixedHead == BRAINSABC_DICT["CSF"]) * (outlabels == 0 | outlabels == 999)
        | ring_csf
        | inner_csf
    )
    # DEBUG sitk.WriteImage(csf_labels,"/tmp/csf_labels.nii.gz")
    outlabels = force_mask_insert(outlabels, csf_labels, OUT_DICT["RH_CSF"])
    # DEBUG sitk.WriteImage(outlabels,"/tmp/outlabels.nii.gz")

    ## Now split CSF based on LeftHemisphereMask
    if LeftHemisphereFN != None:
        LeftHemisphereIm = sitk.Cast(sitk.ReadImage(LeftHemisphereFN), sitk.sitkUInt32)
        left_hemi_pre = outlabels == OUT_DICT["LH_CSF"]
        outlabels = force_mask_insert(
            outlabels, left_hemi_pre, OUT_DICT["RH_CSF"]
        )  ## Make all CSF Right hemisphere
        left_hemi_post = (
            LeftHemisphereIm
            * sitk.Cast((outlabels == OUT_DICT["RH_CSF"]), sitk.sitkUInt32)
            > 0
        )  # SplitCSF with LeftHemisphereMask
        outlabels = force_mask_insert(
            outlabels, left_hemi_post, OUT_DICT["LH_CSF"]
        )  ## Make all CSF Right hemisphere
    ## Now extend brainstem lower
    ## BrainStem often mislabled to Cerebellum WM (label 7 and 46)
    ## Fix for brainstem for the mislabeld Cerebellum as well.
    misLabelDict = OrderedDict({"none": 0, "leftCrblWM": 7, "rightCrblWM": 46})
    for misLabel in misLabelDict:
        brain_stem = (FixedHead == BRAINSABC_DICT["BRAINSTEM"]) * (
            outlabels == misLabelDict[misLabel]
        )
        outlabels = force_mask_insert(
            outlabels, brain_stem, OUT_DICT["BRAINSTEM"]
        )  ## Make all CSF Right hemisphere

    VALID_REGION = sitk.Cast(
        (FixedHead > 0) | ring_csf | inner_csf | ring_vb | inner_vb, sitk.sitkUInt32
    )
    outlabels = outlabels * VALID_REGION
    ## Caudate = 36 37
    ## Putamen = 57 58
    ## Pallidus = 55,56
    ## Thalamus = 59,60
    ## Hippocampus = 47,48
    ## Accumbens  = 23,30
    UNKNOWN_LABEL_CODE = OUT_DICT["UNKNOWN"]
    labels_to_ensure_connected = OUT_DICT["CONNECTED"]
    for keepCode in labels_to_ensure_connected:
        outlabels = recode_nonlargest(outlabels, keepCode, UNKNOWN_LABEL_CODE)

    ## FILL IN HOLES
    unkown_holes = (VALID_REGION > 0) * (outlabels == 0)
    outlabels = force_mask_insert(
        outlabels, unkown_holes, UNKNOWN_LABEL_CODE
    )  ## Fill unkown regions with unkown code
    outlabels = minimize_size_of_image(outlabels)

    fixedFusionLabelFN = os.path.realpath(outFN)
    sitk.WriteImage(outlabels, fixedFusionLabelFN)
    # print("\n\n\n\n\n\n{0}\n\n\n\nXXXXXXXX".format(fixedFusionLabelFN))
    return fixedFusionLabelFN


def recode_label_map(InputFileName, OutputFileName, RECODE_TABLE):
    """
    This funciton...

    :param InputFileName:
    :param OutputFileName:
    :param RECODE_TABLE:
    :return:
    """
    import SimpleITK as sitk
    import os

    def minimize_size_of_image(outlabels):
        """This function will find the largest integer value in the labelmap, and
        cast the image to the smallest possible integer size so that no loss of data
        results.

        :param outlabels:
        :return:
        """
        measureFilt = sitk.StatisticsImageFilter()
        measureFilt.Execute(outlabels)
        imgMin = measureFilt.GetMinimum()
        imgMax = measureFilt.GetMaximum()
        if imgMax < (2 ** 8) - 1:
            outlabels = sitk.Cast(outlabels, sitk.sitkUInt8)
        elif imgMax < (2 ** 16) - 1:
            outlabels = sitk.Cast(outlabels, sitk.sitkUInt16)
        elif imgMax < (2 ** 32) - 1:
            outlabels = sitk.Cast(outlabels, sitk.sitkUInt32)
        elif imgMax < (2 ** 64) - 1:
            outlabels = sitk.Cast(outlabels, sitk.sitkUInt64)
        return outlabels

    LabelImage = sitk.Cast(sitk.ReadImage(InputFileName), sitk.sitkUInt32)
    for (old, new) in RECODE_TABLE:
        LabelImage = (
            sitk.Cast((LabelImage == old), sitk.sitkUInt32) * (new - old) + LabelImage
        )
    LabelImage = minimize_size_of_image(LabelImage)
    recodedFN = os.path.realpath(OutputFileName)
    sitk.WriteImage(LabelImage, recodedFN)
    return recodedFN
