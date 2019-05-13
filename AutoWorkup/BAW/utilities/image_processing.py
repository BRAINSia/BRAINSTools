"""
image_processing.py
============================
Description:
    The purpose of this is to image process functions for pipelines
    
Author:

Usage:

"""

# Image processing functions for pipelines
# from builtins import range


def fix_wm_partitioning(brainMask, PosteriorsList):
    """"There were some errors in mis-classifications for WM/NON_WM"""
    """
    This Function takes in...

    :param brainMask:
    :param PosteriorsList:
    :return: UpdatedPosteriorsList, MatchingFGCodeList, MatchingLabelList, nonAirRegionMask
    """
    import SimpleITK as sitk
    import os

    def fill_hole_preserved_edge(inputMask, HOLE_FILL_SIZE):
        """This function fills holes and tries to preserve
           the exterior topology.  Holes that are within 3 units
           of the exterior topology may not be completely filled.
           Any voxel in the original mask will be guanteed to be
           in the returned mask.

           :param inputMask:
           :param HOLE_FILL_SIZE:
           :return:
           """

        return sitk.BinaryThreshold(
            inputMask
            + sitk.ErodeObjectMorphology(
                sitk.VotingBinaryHoleFilling(
                    BM, [HOLE_FILL_SIZE, HOLE_FILL_SIZE, HOLE_FILL_SIZE]
                ),
                HOLE_FILL_SIZE,
            ),
            1,
            10000,
        )

    print(("Reading {0} of type {1}".format(brainMask, type(brainMask))))
    BM = sitk.BinaryThreshold(sitk.ReadImage(brainMask), 1, 1000)
    BM_FILLED = fill_hole_preserved_edge(BM, 3)

    NOTCSF_index = None  # Note: Purposfully using '-1' as it will force an error.
    CSF_index = None
    NOTGM_index = None
    GM_index = None
    NOTWM_index = None
    WM_index = None
    NOTVB_index = None
    VB_index = None
    AIR_index = None
    bnames = [os.path.basename(fname) for fname in PosteriorsList]
    for i in range(0, len(PosteriorsList)):
        if bnames[i] == "POSTERIOR_NOTCSF.nii.gz":
            NOTCSF_index = i
        elif bnames[i] == "POSTERIOR_CSF.nii.gz":
            CSF_index = i
        elif bnames[i] == "POSTERIOR_NOTGM.nii.gz":
            NOTGM_index = i
        elif bnames[i] == "POSTERIOR_SURFGM.nii.gz":
            GM_index = i
        elif bnames[i] == "POSTERIOR_NOTWM.nii.gz":
            NOTWM_index = i
        elif bnames[i] == "POSTERIOR_WM.nii.gz":
            WM_index = i
        elif bnames[i] == "POSTERIOR_NOTVB.nii.gz":
            NOTVB_index = i
        elif bnames[i] == "POSTERIOR_VB.nii.gz":
            VB_index = i
        elif bnames[i] == "POSTERIOR_AIR.nii.gz":
            AIR_index = i

    def shift_value_for_hard_partition(
        BM_FILLED,
        ShiftPosteriorsList,
        NOTREGION_index,
        REGION_index,
        REGION_NAME,
        NOTREGION_NAME,
    ):
        """
        This function takes in...

        :param BM_FILLED:
        :param ShiftPosteriorsList:
        :param NOTREGION_index:
        :param REGION_index:
        :param REGION_NAME:
        :param NOTREGION_NAME:
        :return:
        """
        print(
            (
                "Reading {0} of type {1}".format(
                    ShiftPosteriorsList[NOTREGION_index],
                    type(ShiftPosteriorsList[NOTREGION_index]),
                )
            )
        )
        NOTREGION = sitk.ReadImage(ShiftPosteriorsList[NOTREGION_index])
        print(
            (
                "Reading {0} of type {1}".format(
                    ShiftPosteriorsList[REGION_index],
                    type(ShiftPosteriorsList[REGION_index]),
                )
            )
        )
        REGION = sitk.ReadImage(ShiftPosteriorsList[REGION_index])
        ALL_REGION = NOTREGION + REGION
        NEW_REGION = ALL_REGION * sitk.Cast(BM_FILLED, sitk.sitkFloat32)
        NEW_NOTREGION = ALL_REGION * sitk.Cast((1 - BM_FILLED), sitk.sitkFloat32)
        NEW_REGION_FN = os.path.realpath("POSTERIOR_{0}.nii.gz".format(REGION_NAME))
        NEW_NOTREGION_FN = os.path.realpath(
            "POSTERIOR_{0}.nii.gz".format(NOTREGION_NAME)
        )
        sitk.WriteImage(NEW_REGION, NEW_REGION_FN)
        sitk.WriteImage(NEW_NOTREGION, NEW_NOTREGION_FN)
        ShiftPosteriorsList[NOTREGION_index] = NEW_NOTREGION_FN
        ShiftPosteriorsList[REGION_index] = NEW_REGION_FN
        return ShiftPosteriorsList

    UpdatedPosteriorsList = list(PosteriorsList)
    UpdatedPosteriorsList = shift_value_for_hard_partition(
        BM_FILLED, UpdatedPosteriorsList, NOTCSF_index, CSF_index, "CSF", "NOTCSF"
    )
    UpdatedPosteriorsList = shift_value_for_hard_partition(
        BM_FILLED, UpdatedPosteriorsList, NOTGM_index, GM_index, "SURFGM", "NOTGM"
    )
    UpdatedPosteriorsList = shift_value_for_hard_partition(
        BM_FILLED, UpdatedPosteriorsList, NOTWM_index, WM_index, "WM", "NOTWM"
    )
    UpdatedPosteriorsList = shift_value_for_hard_partition(
        BM_FILLED, UpdatedPosteriorsList, NOTVB_index, VB_index, "VB", "NOTVB"
    )

    print(
        (
            "Reading {0} of type {1}".format(
                PosteriorsList[AIR_index], type(PosteriorsList[AIR_index])
            )
        )
    )
    AirMask = sitk.BinaryThreshold(
        sitk.ReadImage(PosteriorsList[AIR_index]), 0.50, 1000000
    )
    nonAirMask = sitk.Cast(1 - AirMask, sitk.sitkUInt8)
    nonAirRegionMask = os.path.realpath("NonAirMask.nii.gz")
    sitk.WriteImage(nonAirMask, nonAirRegionMask)

    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    # input order is maintained across runs.
    # This prevents "Outdated cache" errors
    # That were causing the rest of
    # the pipeline to be constantly re-run
    POSTERIOR_LABELS = OrderedDict()  # (FG,Label)
    POSTERIOR_LABELS["POSTERIOR_WM.nii.gz"] = (1, 1)
    POSTERIOR_LABELS["POSTERIOR_SURFGM.nii.gz"] = (1, 2)
    POSTERIOR_LABELS["POSTERIOR_BASAL.nii.gz"] = (1, 21)
    POSTERIOR_LABELS["POSTERIOR_GLOBUS.nii.gz"] = (1, 23)
    POSTERIOR_LABELS["POSTERIOR_THALAMUS.nii.gz"] = (1, 24)
    POSTERIOR_LABELS["POSTERIOR_HIPPOCAMPUS.nii.gz"] = (1, 25)
    POSTERIOR_LABELS["POSTERIOR_CRBLGM.nii.gz"] = (1, 11)
    POSTERIOR_LABELS["POSTERIOR_CRBLWM.nii.gz"] = (1, 12)
    POSTERIOR_LABELS["POSTERIOR_CSF.nii.gz"] = (1, 4)
    POSTERIOR_LABELS["POSTERIOR_VB.nii.gz"] = (1, 5)
    POSTERIOR_LABELS["POSTERIOR_NOTCSF.nii.gz"] = (0, 6)
    POSTERIOR_LABELS["POSTERIOR_NOTGM.nii.gz"] = (0, 7)
    POSTERIOR_LABELS["POSTERIOR_NOTWM.nii.gz"] = (0, 8)
    POSTERIOR_LABELS["POSTERIOR_NOTVB.nii.gz"] = (0, 9)
    POSTERIOR_LABELS["POSTERIOR_AIR.nii.gz"] = (0, 0)

    MatchingFGCodeList = list()
    MatchingLabelList = list()
    for full_post_path_fn in UpdatedPosteriorsList:
        post_key = os.path.basename(full_post_path_fn)
        MatchingFGCodeList.append(POSTERIOR_LABELS[post_key][0])
        MatchingLabelList.append(POSTERIOR_LABELS[post_key][1])

    return (
        UpdatedPosteriorsList,
        MatchingFGCodeList,
        MatchingLabelList,
        nonAirRegionMask,
    )
