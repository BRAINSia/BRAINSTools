## \author Hans J. Johnson
## These are nodes used to fix-up label maps by merging information from
## various different sources.
##

def FixLabelMapFromNeuromorphemetrics2012(fusionFN,FixedHeadFN,posterior_dict,LeftHemisphereFN,outFN, OUT_DICT):
    import SimpleITK as sitk
    import os

    def ForceMaskInsert(inlabels,newmask,newmaskvalue):
        inlabels = sitk.Cast(inlabels,sitk.sitkUInt32)
        newmask = sitk.Cast( (newmask>0) , sitk.sitkUInt32)
        outlabels=inlabels*sitk.Cast( (1-newmask), sitk.sitkUInt32)
        outlabels = outlabels + newmask*newmaskvalue
        return sitk.Cast(outlabels,sitk.sitkUInt32)
    ## TODO: GetLargestLabel is copied from elsewhere
    def GetLargestLabel(inputMask, UseErosionCleaning):
        LargestComponentCode = 1
        if UseErosionCleaning:
            erosionMask = sitk.ErodeObjectMorphology(inputMask, 1)
        else:
            erosionMask = inputMask
        CC = sitk.ConnectedComponent(erosionMask)
        Rlabel = sitk.RelabelComponent(CC)
        largestMask = ( Rlabel == LargestComponentCode)
        if UseErosionCleaning:
            dilateMask = sitk.DilateObjectMorphology(largestMask, 1)
        else:
            dilateMask = largestMask

        return (largestMask * dilateMask > 0)

    def RecodeNonLargest(outlabels,keepCode,UNKNOWN_LABEL_CODE):
        orig_mask = (outlabels ==  keepCode)
        connected_mask = GetLargestLabel(orig_mask,False)
        small_regions = ( orig_mask - connected_mask )
        outlabels = ForceMaskInsert(outlabels,connected_mask,keepCode)
        outlabels = ForceMaskInsert(outlabels,small_regions,UNKNOWN_LABEL_CODE)
        return outlabels

    def MinimizeSizeOfImage(outlabels):
        """This function will find the largest integer value in the labelmap, and
        cast the image to the smallest possible integer size so that no loss of data
        results."""
        measureFilt  = sitk.StatisticsImageFilter()
        measureFilt.Execute(outlabels)
        imgMin=measureFilt.GetMinimum()
        imgMax=measureFilt.GetMaximum()
        if imgMax < (2**8)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt8 )
        elif imgMax < (2**16)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt16 )
        elif imgMax < (2**32)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt32 )
        elif imgMax < (2**64)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt64 )
        return outlabels

    fusionIm=sitk.Cast(sitk.ReadImage(fusionFN),sitk.sitkUInt32)
    FixedHead=sitk.Cast(sitk.ReadImage(FixedHeadFN),sitk.sitkUInt32)

    BRAINSABC_DICT = { 'BRAINSTEM': 30, 'CSF': 4 , 'BLOOD': 5 }

    ## Intialize by cloning
    outlabels = sitk.Image(fusionIm)
    lbl_orig_mask = ( outlabels > 0 )
    lbl_outter_ring = sitk.BinaryDilate(lbl_orig_mask,2) - lbl_orig_mask
    # DEBUG sitk.WriteImage(lbl_outter_ring,"/tmp/lbl_outter_ring.nii.gz")


    vb_post = sitk.ReadImage(posterior_dict["VB"])
    ring_vb = lbl_outter_ring * sitk.BinaryThreshold(vb_post, 0.5,1.01,1,0) # just outside mask
    # DEBUG sitk.WriteImage(ring_vb,"/tmp/ring_vb.nii.gz")
    inner_vb = lbl_orig_mask * sitk.BinaryThreshold(vb_post , 0.85, 1.01, 1, 0) # inside mask, but very high probability
    # DEBUG sitk.WriteImage(inner_vb,"/tmp/inner_vb.nii.gz")
    # background = 0 , suspicous = 999
    ## Add blood from BRAINSABC to mask as as value OUT_DICT['BLOOD']
    blood_labels=(FixedHead == BRAINSABC_DICT['BLOOD']) * (outlabels == 0 | outlabels == 999) | ring_vb | inner_vb
    # DEBUG sitk.WriteImage(blood_labels,"/tmp/blood_labels.nii.gz")
    outlabels = ForceMaskInsert(outlabels,blood_labels,OUT_DICT['BLOOD'])

    csf_post = sitk.ReadImage(posterior_dict["CSF"])
    ring_csf = lbl_outter_ring * sitk.BinaryThreshold(csf_post, 0.5, 1.01, 1, 0) # just outside mask
    # DEBUG sitk.WriteImage(ring_csf,"/tmp/ring_csf.nii.gz")
    inner_csf = lbl_orig_mask * sitk.BinaryThreshold(csf_post, 0.85, 1.01, 1, 0) # inside mask, but very high probability
    # DEBUG sitk.WriteImage(inner_csf,"/tmp/inner_csf.nii.gz")
    ## Add CSF from BRAINSABC to mask as as value OUT_DICT['RH_CSF']
    csf_labels=(FixedHead == BRAINSABC_DICT['CSF'] ) * (outlabels == 0 | outlabels == 999) | ring_csf | inner_csf
    # DEBUG sitk.WriteImage(csf_labels,"/tmp/csf_labels.nii.gz")
    outlabels= ForceMaskInsert(outlabels,csf_labels,OUT_DICT['RH_CSF'])
    # DEBUG sitk.WriteImage(outlabels,"/tmp/outlabels.nii.gz")

    ## Now split CSF based on LeftHemisphereMask
    if LeftHemisphereFN != None:
        LeftHemisphereIm=sitk.Cast(sitk.ReadImage(LeftHemisphereFN),sitk.sitkUInt32)
        left_hemi_pre = ( outlabels == OUT_DICT['LH_CSF'] )
        outlabels = ForceMaskInsert(outlabels,left_hemi_pre,OUT_DICT['RH_CSF'])  ## Make all CSF Right hemisphere
        left_hemi_post =  (LeftHemisphereIm * sitk.Cast ( ( outlabels == OUT_DICT['RH_CSF'] ),sitk.sitkUInt32) > 0 ) # SplitCSF with LeftHemisphereMask
        outlabels = ForceMaskInsert(outlabels,left_hemi_post,OUT_DICT['LH_CSF'])  ## Make all CSF Right hemisphere
    ## Now extend brainstem lower
    ## BrainStem often mislabled to Cerebellum WM (label 7 and 46)
    ## Fix for brainstem for the mislabeld Cerebellum as well.
    misLabelDict = { "none":0, "leftCrblWM":7, "rightCrblWM":46}
    for misLabel in misLabelDict:
        brain_stem = (FixedHead == BRAINSABC_DICT['BRAINSTEM']) * ( outlabels == misLabelDict[misLabel])
        outlabels = ForceMaskInsert(outlabels,brain_stem,OUT_DICT['BRAINSTEM'])  ## Make all CSF Right hemisphere

    VALID_REGION = sitk.Cast( (FixedHead > 0) | ring_csf | inner_csf | ring_vb | inner_vb ,sitk.sitkUInt32)
    outlabels = outlabels * VALID_REGION
    ## Caudate = 36 37
    ## Putamen = 57 58
    ## Pallidus = 55,56
    ## Thalamus = 59,60
    ## Hippocampus = 47,48
    ## Accumbens  = 23,30
    UNKNOWN_LABEL_CODE=OUT_DICT['UNKNOWN']
    labels_to_ensure_connected = OUT_DICT['CONNECTED']
    for keepCode in labels_to_ensure_connected:
        outlabels = RecodeNonLargest(outlabels,keepCode,UNKNOWN_LABEL_CODE)

    ## FILL IN HOLES
    unkown_holes = ( VALID_REGION > 0 ) * ( outlabels == 0 )
    outlabels = ForceMaskInsert(outlabels,unkown_holes,UNKNOWN_LABEL_CODE)  ## Fill unkown regions with unkown code
    outlabels = MinimizeSizeOfImage(outlabels)

    fixedFusionLabelFN=os.path.realpath(outFN)
    sitk.WriteImage(outlabels,fixedFusionLabelFN)
    #print("\n\n\n\n\n\n{0}\n\n\n\nXXXXXXXX".format(fixedFusionLabelFN))
    return fixedFusionLabelFN

def RecodeLabelMap(InputFileName,OutputFileName,RECODE_TABLE):
    import SimpleITK as sitk
    import os

    def MinimizeSizeOfImage(outlabels):
        """This function will find the largest integer value in the labelmap, and
        cast the image to the smallest possible integer size so that no loss of data
        results."""
        measureFilt  = sitk.StatisticsImageFilter()
        measureFilt.Execute(outlabels)
        imgMin=measureFilt.GetMinimum()
        imgMax=measureFilt.GetMaximum()
        if imgMax < (2**8)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt8 )
        elif imgMax < (2**16)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt16 )
        elif imgMax < (2**32)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt32 )
        elif imgMax < (2**64)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt64 )
        return outlabels
    LabelImage=sitk.Cast(sitk.ReadImage(InputFileName),sitk.sitkUInt32)
    for (old,new) in RECODE_TABLE:
        LabelImage = sitk.Cast((LabelImage == old),sitk.sitkUInt32)*(new - old)+LabelImage
    LabelImage = MinimizeSizeOfImage(LabelImage)
    recodedFN=os.path.realpath(OutputFileName)
    sitk.WriteImage(LabelImage,recodedFN)
    return recodedFN
