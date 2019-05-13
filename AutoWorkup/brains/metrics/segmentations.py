"""
segmentations.py
=================
Description:

Author:

Usage:

"""
from itertools import starmap, zip_longest

import os.path
from builtins import map
from builtins import range

import SimpleITK as sitk
import numpy as np

from .partials import calcutateBinaryVolume
from ..common import check_file
from ..config import _config
from collections import (
    OrderedDict,
)  # Need OrderedDict internally to ensure consistent ordering

labels = ["caudate", "putamen", "hippocampus", "thalamus", "accumben", "globus", "icv"]


def construct_labels(labels):
    """
    This function...

    :param labels:
    :return:
    """
    numbers = list(range(1, ((len(labels) * 2) + 1)))
    full_labels = []
    index = 0
    for label in labels:
        full_labels.append("_".join(["left", label]))
        full_labels.append("_".join(["right", label]))
    return full_labels, numbers


def _module_create_labels(labels):
    """
    This function...

    :param labels:
    :return:
    """
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    full_labels, numbers = construct_labels(labels)
    labelMap = zip_longest(full_labels, numbers)
    return OrderedDict(labelMap)  # Use this variable


def format_labels(label):
    """
    Assumes that the label can be split by the '_' character.

    :param label:
    :return:
    """
    side, anatomy = label.split("_")
    if side.lower() in ["l", "left"]:
        side = "left"
    elif side.lower() in ["r", "right"]:
        side = "right"
    else:
        raise ValueError(
            "Label %s is not recognized: cannot determine side %s" % (label, side)
        )
    label = "_".join([side, anatomy])
    return label


def calculate_label_volume(dirname, label):
    """
    This function...

    :param dirname:
    :param label:
    :return:
    """
    labelFile = os.path.join(
        dirname, _config.get("Results", "segmentations"), label + "_seg_seg.nii.gz"
    )
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    nda = sitk.GetArrayFromImage(image)
    maskSum = nda.sum()
    print(maskSum)
    size = image.GetSpacing()
    print(size)
    return maskSum * size[0] * size[1] * size[2]


def calculate_icv(dirname):
    """
    This function...

    :param dirname:
    """
    filename = os.path.join(
        dirname, _config.get("Results", "partials"), "fixed_brainlabels_seg.nii.gz"
    )
    filename = check_file(filename)
    calculate_binary_volume(filename)


def get_volume(args=[], kwds=OrderedDict()):
    """
    This function...

    :param args:
    :param kwds:
    :return:
    """
    dirname = labels = project = subject = session = experimentDir = None
    experimentDir = _config.get("Results", "directory")  # HACK
    for key, value in kwds:
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
            # config needs to be accessible
            dirname = os.path.join(experimentDir, project, subject, session)
        except Exception as err:
            raise err
    volume = 0.0
    for label in labels:
        volume += calculate_label_volume(dirname, label)
    return volume
