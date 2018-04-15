
from itertools import starmap, zip_longest

import os.path
from builtins import map
from builtins import range

import SimpleITK as sitk
import numpy as np

from .partials import calcutateBinaryVolume
from ..common import check_file
from ..config import _config

labels = ['caudate', 'putamen', 'hippocampus', 'thalamus', 'accumben', 'globus', 'icv']


def constructLabels(labels):
    numbers = list(range(1, ((len(labels) * 2) + 1)))
    full_labels = []
    index = 0
    for label in labels:
        full_labels.append('_'.join(['left', label]))
        full_labels.append('_'.join(['right', label]))
    return full_labels, numbers


def _moduleCreateLabels(labels):
    full_labels, numbers = constructLabels(labels)
    labelMap = zip_longest(full_labels, numbers)
    return dict(labelMap)  # Use this variable


def formatLabel(label):
    """
    Assumes that the label can be split by the '_' character.
    """
    side, anatomy = label.split('_')
    if side.lower() in ['l', 'left']:
        side = 'left'
    elif side.lower() in ['r', 'right']:
        side = 'right'
    else:
        raise ValueError('Label %s is not recognized: cannot determine side %s' % (label, side))
    label = '_'.join([side, anatomy])
    return label


def calculateLabelVolume(dirname, label):
    labelFile = os.path.join(dirname, _config.get('Results', 'segmentations'),
                             label + '_seg_seg.nii.gz')
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    nda = sitk.GetArrayFromImage(image)
    maskSum = nda.sum()
    print(maskSum)
    size = image.GetSpacing()
    print(size)
    return maskSum * size[0] * size[1] * size[2]


def calculateICV(dirname):
    filename = os.path.join(dirname, _config.get('Results', 'partials'),
                            'fixed_brainlabels_seg.nii.gz')
    filename = check_file(filename)
    calculateBinaryVolume(filename)


def getVolume(args=[], kwds={}):
    dirname = labels = project = subject = session = experimentDir = None
    experimentDir = _config.get('Results', 'directory')  # HACK
    for key, value in kwds:
        if key == 'dirname':
            dirname = check_file(value)
        elif key == 'labels':
            labels = value
        elif key == 'project':
            project = value
        elif key == 'subject':
            subject = value
        elif key == 'session':
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
        volume += calculateLabelVolume(dirname, label)
    return volume
