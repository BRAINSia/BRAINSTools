"""
partials.py
=================
Description:

Author:

Usage:

"""
import os.path
from warnings import warn

import SimpleITK as sitk
import numpy as np

from ..common import check_file
from ..config import _config

partials = [
    "accumben",
    "air",
    "caudate",
    "crblgm",
    "crblwm",
    "csf",
    "globus",
    "hippocampus",
    "notcsf",
    "notgm",
    "notvb",
    "notwm",
    "putamen",
    "surfgm",
    "thalamus",
    "vb",
    "wm",
]
accumulated = [
    "background_total",
    "gm_total",
    "csf_total",
    "vb_total",
    "globus_total",
    "wm_total",
]

_isAccumulated = False
_tolerance = [0.51, 1.01]


def _format_partial_assert_string():
    """
    Returns assertion string for label

    :return: assertString
    """
    assertString = "Partial label is not recognized: %s\nValid labels are:"
    for p in partials + accumulated:
        assertString = "\n".join([assertString, p.upper() + ","])
    assertString = assertString[:-1]
    return assertString


def _check_label(label):
    """
    Verifies string is in the list of valid labels for partial volume resuls

    :param label:
    """
    errorString = _format_partial_assert_string()
    if label.lower() == "icv":
        pass
    else:
        assert label.lower() in partials + accumulated, errorString % label


def _set_if_accumulated(label):
    """
    Sets _isAccumulated to True if label is in accumulated list

    :param label:
    """
    global _isAccumulated
    if label.lower() in accumulated + ["icv"]:
        _isAccumulated = True
    else:
        _isAccumulated = False
        # if label == "ICV":
        #     print "_isAccumulated: ", _isAccumulated


def calculate_binary_volume(dirname, label, _isAccumulated=True, tolerance=_tolerance):
    """
    This function...

    :param dirname:
    :param label:
    :param _isAccumulated:
    :param tolerance:
    :return:
    """
    label = label.upper()

    maskSum = 0.0
    if _isAccumulated:
        fileDir = _config.get("Results", "accumulated")
    else:
        # print "Not accumulated: ", label
        fileDir = _config.get("Results", "partials")

    if label == "ICV":
        for sublabel in accumulated:
            if sublabel == "background_total":
                continue
            else:
                # print "sublabel: ", sublabel, calculate_binary_volume(dirname, sublabel, True)
                maskSum += calculate_binary_volume(dirname, sublabel, True)
        return maskSum

    labelFile = os.path.join(dirname, fileDir, "POSTERIOR_" + label + ".nii.gz")
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    lowerTol, upperTol = tolerance
    binary = sitk.BinaryThreshold(image, lowerTol, upperTol, 1, 0)
    nda = sitk.GetArrayFromImage(binary)
    maskSum += nda.sum()
    ## print maskSum
    size = image.GetSpacing()
    ## print size
    return maskSum * size[0] * size[1] * size[2]


def calculate_partial_volume(dirname, label, _isAccumulated=True):
    """
    This function...

    :param dirname:
    :param label:
    :param _isAccumulated:
    :return:
    """
    label = label.upper()

    maskSum = 0.0
    if _isAccumulated:
        fileDir = _config.get("Results", "accumulated")
    else:
        fileDir = _config.get("Results", "partials")

    if label == "ICV":
        for sublabel in accumulated:
            if sublabel == "background_total":
                continue
            else:
                # print "sublabel: ", sublabel, calculate_partial_volume(dirname, sublabel, True)
                maskSum += calculate_partial_volume(dirname, sublabel, True)
        return maskSum

    labelFile = os.path.join(dirname, fileDir, "POSTERIOR_" + label + ".nii.gz")
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    nda = sitk.GetArrayFromImage(image)
    maskSum = nda.sum()
    ## print maskSum
    size = image.GetSpacing()
    ## print size
    return maskSum * size[0] * size[1] * size[2]


def get_posterior_volume(*args, **kwds):
    """
    This function...

    :param args:
    :param kwds:
    :return:
    """
    dirname = labels = project = subject = session = experimentDir = None
    experimentDir = _config.get("Results", "directory")
    # parse keywords
    for key, value in list(kwds.items()):
        if key == "dirname":
            dirname = check_file(value)
        elif key == "labels":
            labels = value
        elif key == "project":
            project = value
        elif key == "subject":
            subject = value
        elif key == "session":
            session = value
        elif key == "binary":
            binary = value
    # Set keyword-only defaults
    if binary is None:
        binary = True
    # parse ordered arguments
    args = list(args)
    if len(args) > 0:
        if session is None:
            session = args.pop()
        if subject is None:
            subject = args.pop()
        if project is None:
            project = args.pop()
        if labels is None:
            labels = args
        if dirname is None:
            try:
                dirname = check_file(
                    os.path.join(experimentDir, project, subject, session)
                )
            except Exception as err:
                raise err
    assert dirname is not None
    volume = 0.0
    if isinstance(labels, str):
        labels = [labels]

    for label in labels:
        label = label.upper()
        _check_label(label)
        _set_if_accumulated(label)
        if binary:
            volume += calculate_binary_volume(dirname, label, _isAccumulated)
        else:
            volume += calculate_partial_volume(dirname, label, _isAccumulated)
    return volume
