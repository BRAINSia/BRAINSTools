import numpy as np
import os.path
from warnings import warn

import SimpleITK as sitk

from ..config import _config

posteriors = ['accumben', 'air', 'caudate', 'crblgm', 'crblwm', 'csf', 'globus', 'hippocampus', 'notcsf', 'notgm', 'notvb', 'notwm', 'putamen', 'surfgm', 'thalamus', 'vb', 'wm']
accumulated = ['background_total', 'gm_total', 'csf_total', 'vb_total', 'globus_total', 'wm_total']

_isAccumulated = False
_tolerance = [0.51, 1.01]
_tolerance_wm = [0.99, 1.01]


def _formatPartialAssertString():
    """
    Returns assertion string for label
    """
    assertString = "Partial label is not recognized: %s\nValid labels are:"
    for p in posteriors + accumulated:
        assertString = '\n'.join([assertString, p.upper() + ','])
    assertString = assertString[:-1]
    return assertString


def _checkLabel(label):
    """
    Verifies string is in the list of valid labels for partial volume resuls
    """
    errorString = _formatPartialAssertString()
    assert label.lower() in posteriors + accumulated, errorString % label


def _setIfAccumulated(label):
    """
    Sets _isAccumulated to True if label is in accumulated list
    """
    if label.lower() in accumulated:
        global _isAccumulated
        _isAccumulated = True


def calculateBinaryVolume(dirname, label):
    label = label.upper()
    _checkLabel(label)
    _setIfAccumulated(label)
    if _isAccumulated:
        fileDir = _config.get('Results', 'accumulated')
    else:
        fileDir = _config.get('Results', 'posteriors')
    labelFile = os.path.join(dirname, fileDir, 'POSTERIOR_'+ label + '.nii.gz')
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    if label == 'WM':
        lowerTol, upperTol = _tolerance_wm
        warn("Using lower threshold of %g for white matter binary calculations." % lowerTol)
    else:
        lowerTol, upperTol = _tolerance
    binary = sitk.BinaryThreshold(image, lowerTol, upperTol, 1, 0)
    nda = sitk.GetArrayFromImage(binary)
    maskSum = nda.sum()
    ## print maskSum
    size = image.GetSpacing()
    ## print size
    return maskSum * size[0] * size[1] * size[2]


def calculatePartialVolume(dirname, label):
    """
    """
    label = label.upper()
    _checkLabel(label)
    _setIfAccumulated(label)
    if _isAccumulated:
        fileDir = _config.get('Results', 'accumulated')
    else:
        fileDir = _config.get('Results', 'posteriors')
    labelFile = os.path.join(dirname, fileDir, 'POSTERIOR_' + label + '.nii.gz')
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    nda = sitk.GetArrayFromImage(image)
    maskSum = nda.sum()
    ## print maskSum
    size = image.GetSpacing()
    ## print size
    return maskSum * size[0] * size[1] * size[2]


def getPartialVolume(args=[], kwds={}):
    """

    """
    dirname = labels = project = subject = session = experimentDir = None
    experimentDir = _config.get('Results', 'directory')
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
    ### DEBUGGING ###
    #print labels
    #print dirname
    #print project
    #print subject
    #print session
    ### END DEBUG ###
    # labels = map(formatLabel, labels) # convert shorthand to human-readable
    if dirname is None:
        try:
            dirname = check_file(os.path.join(experimentDir, project, subject, session))
        except Exception, err:
            raise err
    volume = 0.0
    for label in labels:
        volume += calculateBinaryVolume(dirname, label)
    return volume

