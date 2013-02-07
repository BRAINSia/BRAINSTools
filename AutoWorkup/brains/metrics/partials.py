import numpy as np
import os.path
from warnings import warn

import SimpleITK as sitk

from ..config import _config
from ..common import check_file

partials = ['accumben', 'air', 'caudate', 'crblgm', 'crblwm', 'csf', 'globus', 'hippocampus', 'notcsf', 'notgm', 'notvb', 'notwm', 'putamen', 'surfgm', 'thalamus', 'vb', 'wm']
accumulated = ['background_total', 'gm_total', 'csf_total', 'vb_total', 'globus_total', 'wm_total']

_isAccumulated = False
_tolerance = [0.51, 1.01]


def _formatPartialAssertString():
    """
    Returns assertion string for label
    """
    assertString = "Partial label is not recognized: %s\nValid labels are:"
    for p in partials + accumulated:
        assertString = '\n'.join([assertString, p.upper() + ','])
    assertString = assertString[:-1]
    return assertString


def _checkLabel(label):
    """
    Verifies string is in the list of valid labels for partial volume resuls
    """
    errorString = _formatPartialAssertString()
    if label.lower() == 'icv':
        pass
    else:
        assert label.lower() in partials + accumulated, errorString % label


def _setIfAccumulated(label):
    """
    Sets _isAccumulated to True if label is in accumulated list
    """
    global _isAccumulated
    if label.lower() in accumulated + ['icv']:
        _isAccumulated = True
    else:
        _isAccumulated = False
    # if label == "ICV":
    #     print "_isAccumulated: ", _isAccumulated


def calculateBinaryVolume(dirname, label, _isAccumulated=True, tolerance=_tolerance):
    label = label.upper()

    maskSum = 0.0
    if _isAccumulated:
        fileDir = _config.get('Results', 'accumulated')
    else:
        # print "Not accumulated: ", label
        fileDir = _config.get('Results', 'partials')

    if label == 'ICV':
        for sublabel in accumulated:
            if sublabel == 'background_total':
                continue
            else:
                # print "sublabel: ", sublabel, calculateBinaryVolume(dirname, sublabel, True)
                maskSum += calculateBinaryVolume(dirname, sublabel, True)
        return maskSum

    labelFile = os.path.join(dirname, fileDir, 'POSTERIOR_'+ label + '.nii.gz')
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


def calculatePartialVolume(dirname, label, _isAccumulated=True):
    """
    """
    label = label.upper()

    maskSum = 0.0
    if _isAccumulated:
        fileDir = _config.get('Results', 'accumulated')
    else:
        fileDir = _config.get('Results', 'partials')

    if label == 'ICV':
        for sublabel in accumulated:
            if sublabel == 'background_total':
                continue
            else:
                # print "sublabel: ", sublabel, calculatePartialVolume(dirname, sublabel, True)
                maskSum += calculatePartialVolume(dirname, sublabel, True)
        return maskSum

    labelFile = os.path.join(dirname, fileDir, 'POSTERIOR_' + label + '.nii.gz')
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    nda = sitk.GetArrayFromImage(image)
    maskSum = nda.sum()
    ## print maskSum
    size = image.GetSpacing()
    ## print size
    return maskSum * size[0] * size[1] * size[2]


def getPosteriorVolume(*args, **kwds):
    """

    """
    dirname = labels = project = subject = session = experimentDir = None
    experimentDir = _config.get('Results', 'directory')
    # parse keywords
    for key, value in kwds.items():
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
        elif key == 'binary':
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
                dirname = check_file(os.path.join(experimentDir, project, subject, session))
            except Exception, err:
                raise err
    assert dirname is not None
    volume = 0.0
    if isinstance(labels, str):
        labels = [labels]

    for label in labels:
        label = label.upper()
        _checkLabel(label)
        _setIfAccumulated(label)
        if binary:
            volume += calculateBinaryVolume(dirname, label, _isAccumulated)
        else:
            volume += calculatePartialVolume(dirname, label, _isAccumulated)
    return volume

