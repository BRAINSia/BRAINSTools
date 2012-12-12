## NOTE: THERE CAN NOT BE ANY GLOBAL imports in this file
##       NIPYPE pipeline functions must be self contained
##       and any import needed for a function must be
##       included in the function itself.

## This file contains misc SimpleITK based functions for use in nipype
## nodes.

## AVOID REFORMATTING THIS FILE, it causes the hash to change in
## nipype and that require re-running the function.

# Globals
POSTERIORS = ['WM', 'SURFGM', 'ACCUMBEN', 'CAUDATE', 'PUTAMEN', 'GLOBUS', 'THALAMUS',
              'HIPPOCAMPUS', 'CRBLGM', 'CRBLWM', 'CSF', 'VB', 'NOTCSF', 'NOTGM', 'NOTWM',
              'NOTVB', 'AIR']

def getListIndex( imageList, index):
    return imageList[index]

def ClipT1ImageWithBrainMask(t1_image,brain_labels,clipped_file_name):
    import os
    import sys
    import SimpleITK as sitk
    ## Now clean up the posteriors based on anatomical knowlege.
    ## sometimes the posteriors are not relevant for priors
    ## due to anomolies around the edges.
    t1=sitk.Cast(sitk.ReadImage(t1_image),sitk.sitkFloat32)
    bl=sitk.Cast(sitk.ReadImage(brain_labels),sitk.sitkFloat32)
    bl_binary=sitk.Cast(sitk.BinaryThreshold(bl,1,1000000),sitk.sitkFloat32)
    clipped=t1*bl_binary
    sitk.WriteImage(clipped,clipped_file_name)
    clipped_file=os.path.realpath(clipped_file_name)
    return clipped_file

def UnwrapPosteriorImagesFromDictionaryFunction(postDict):
    return postDict.values()



def FixWMPartitioning(brainMask,PosteriorsList):
    """"There were some errors in mis-classifications for WM/NON_WM"""
    import SimpleITK as sitk
    import os

    def FillHolePreserveEdge(inputMask,HOLE_FILL_SIZE):
        """This function fills holes and tries to preserve
           the exterior topology.  Holes that are within 3 units
           of the exterior topology may not be completely filled.
           Any voxel in the original mask will be guanteed to be
           in the returned mask."""
        return sitk.BinaryThreshold(
                  inputMask +
                     sitk.ErodeObjectMorphology(
                                      sitk.VotingBinaryHoleFilling(BM,[HOLE_FILL_SIZE,HOLE_FILL_SIZE,HOLE_FILL_SIZE])
                     ,HOLE_FILL_SIZE)
                ,1,10000)

    BM=sitk.BinaryThreshold(sitk.ReadImage(brainMask),1,1000)
    BM_FILLED= FillHolePreserveEdge(BM,3)

    NOTCSV_index=-1
    CSF_index=-1
    NOTGM_index=-1
    GM_index=-1
    NOTWM_index=-1
    WM_index=-1
    NOTVB_index=-1
    VB_index=-1
    AIR_index=-1
    for i in range(0,len(PosteriorsList)):
        if os.path.basename(PosteriorsList[i]) == 'POSTERIOR_NOTCSF.nii.gz':
            NOTCSF_index=i
        elif os.path.basename(PosteriorsList[i]) == 'POSTERIOR_CSF.nii.gz':
            CSF_index=i
        elif os.path.basename(PosteriorsList[i]) == 'POSTERIOR_NOTGM.nii.gz':
            NOTGM_index=i
        elif os.path.basename(PosteriorsList[i]) == 'POSTERIOR_SURFGM.nii.gz':
            GM_index=i
        elif os.path.basename(PosteriorsList[i]) == 'POSTERIOR_NOTWM.nii.gz':
            NOTWM_index=i
        elif os.path.basename(PosteriorsList[i]) == 'POSTERIOR_WM.nii.gz':
            WM_index=i
        elif os.path.basename(PosteriorsList[i]) == 'POSTERIOR_NOTVB.nii.gz':
            NOTVB_index=i
        elif os.path.basename(PosteriorsList[i]) == 'POSTERIOR_VB.nii.gz':
            VB_index=i
        elif os.path.basename(PosteriorsList[i]) == 'POSTERIOR_AIR.nii.gz':
            AIR_index=i
                
    def ShiftValueForHardPartition(BM_FILLED,PosteriorsList,NOTREGION_index,REGION_index,REGION_NAME,NOTREGION_NAME):
        NOTREGION=sitk.ReadImage(PosteriorsList[NOTREGION_index])
        REGION=sitk.ReadImage(PosteriorsList[REGION_index])
        ALL_REGION=NOTREGION+REGION
        NEW_REGION=ALL_REGION*sitk.Cast(BM_FILLED,sitk.sitkFloat32)
        NEW_NOTREGION=ALL_REGION*sitk.Cast((1-BM_FILLED),sitk.sitkFloat32)
        NEW_REGION_FN=os.path.realpath('POSTERIOR_{0}.nii.gz'.format(REGION_NAME))
        NEW_NOTREGION_FN=os.path.realpath('POSTERIOR_{0}.nii.gz'.format(NOTREGION_NAME))
        sitk.WriteImage(NEW_REGION,NEW_REGION_FN)
        sitk.WriteImage(NEW_NOTREGION,NEW_NOTREGION_FN)
        PosteriorsList[NOTREGION_index]=NEW_NOTREGION_FN
        PosteriorsList[REGION_index]=NEW_REGION_FN
        return PosteriorsList

    UpdatedPosteriorsList=PosteriorsList[:]
    UpdatedPosteriorsList = ShiftValueForHardPartition(BM_FILLED,UpdatedPosteriorsList,NOTCSV_index,CSF_index,'CSF','NOTCSF')
    UpdatedPosteriorsList = ShiftValueForHardPartition(BM_FILLED,UpdatedPosteriorsList,NOTGM_index,GM_index,'SURFGM','NOTGM')
    UpdatedPosteriorsList = ShiftValueForHardPartition(BM_FILLED,UpdatedPosteriorsList,NOTWM_index,WM_index,'WM','NOTWM')
    UpdatedPosteriorsList = ShiftValueForHardPartition(BM_FILLED,UpdatedPosteriorsList,NOTVB_index,VB_index,'VB','NOTVB')
    
    AirMask=sitk.BinaryThreshold( sitk.ReadImage(PosteriorsList[AIR_index]),0.50,1000000)
    nonAirMask=sitk.Cast(1-AirMask,sitk.sitkUInt8)
    nonAirRegionMask=os.path.realpath('NonAirMask.nii.gz')
    sitk.WriteImage(nonAirMask,nonAirRegionMask)
    
    POSTERIOR_LABELS=dict()                          #(FG,Label) 
    POSTERIOR_LABELS["POSTERIOR_ACCUMBEN.nii.gz"] =    (1,20)
    POSTERIOR_LABELS["POSTERIOR_AIR.nii.gz"] =         (0,0)
    POSTERIOR_LABELS["POSTERIOR_CAUDATE.nii.gz"] =     (1,21)
    POSTERIOR_LABELS["POSTERIOR_CRBLGM.nii.gz"] =      (1,11)
    POSTERIOR_LABELS["POSTERIOR_CRBLWM.nii.gz"] =      (1,12)
    POSTERIOR_LABELS["POSTERIOR_CSF.nii.gz"] =         (1,4)
    POSTERIOR_LABELS["POSTERIOR_GLOBUS.nii.gz"] =      (1,23)
    POSTERIOR_LABELS["POSTERIOR_HIPPOCAMPUS.nii.gz"] = (1,25)
    POSTERIOR_LABELS["POSTERIOR_NOTCSF.nii.gz"] =      (0,6)
    POSTERIOR_LABELS["POSTERIOR_NOTGM.nii.gz"] =       (0,7)
    POSTERIOR_LABELS["POSTERIOR_NOTVB.nii.gz"] =       (0,9)
    POSTERIOR_LABELS["POSTERIOR_NOTWM.nii.gz"] =       (0,8)
    POSTERIOR_LABELS["POSTERIOR_PUTAMEN.nii.gz"] =     (1,22)
    POSTERIOR_LABELS["POSTERIOR_SURFGM.nii.gz"] =      (1,2)
    POSTERIOR_LABELS["POSTERIOR_THALAMUS.nii.gz"] =    (1,24)
    POSTERIOR_LABELS["POSTERIOR_VB.nii.gz"] =          (1,5)
    POSTERIOR_LABELS["POSTERIOR_WM.nii.gz"] =          (1,1)

    MatchingFGCodeList=list()
    MatchingLabelList=list()
    for full_post_path_fn in UpdatedPosteriorsList:
        post_key=os.path.basename(full_post_path_fn)
        MatchingFGCodeList.append(POSTERIOR_LABELS[post_key][0])
        MatchingLabelList.append(POSTERIOR_LABELS[post_key][1])
        
    return UpdatedPosteriorsList,MatchingFGCodeList,MatchingLabelList,nonAirRegionMask
