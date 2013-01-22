import numpy as np
import os.path
from warnings import warn

import SimpleITK as sitk

from ..config import _config

posteriors = ['accumben', 'air', 'caudate', 'crblgm', 'crblwm', 'csf', 'globus', 'hippocampus', 'notcsf', 'notgm', 'notvb', 'notwm', 'putamen', 'surfgm', 'thalamus', 'vb', 'wm']
_tolerance = [0.99, 1.01]
_tolerance_wm = [0.99, 1.01]


def _formatPosteriorAssertString():
    assertString = "Posterior label is not recognized: %s\nValid labels are:"
    for p in posteriors:
        assertString = '\n'.join([assertString, p.upper() + ','])
    assertString = assertString[:-1]
    return assertString


def _checkPosteriorLabel(label):
    errorString = _formatPosteriorAssertString()
    assert label.lower() in posteriors, errorString % label


def calculateBinaryVolume(dirname, label):
    label = label.upper()
    _checkPosteriorLabel(label)
    labelFile = os.path.join(dirname, 'POSTERIOR_'+ label + '.nii.gz')
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
    label = label.upper()
    _checkPosteriorLabel(label)
    labelFile = os.path.join(dirname, 'POSTERIOR_' + label + '.nii.gz')
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    nda = sitk.GetArrayFromImage(image)
    maskSum = nda.sum()
    ## print maskSum
    size = image.GetSpacing()
    ## print size
    return maskSum * size[0] * size[1] * size[2]


def getPartialVolume(args=[], kwds={}):
    dirname = labels = project = subject = session = experimentDir = None
    experimentDir = _config.get('Results', 'directory')
    partialsDir = _config.get('Results', 'posteriors')
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
            dirname = os.path.join(experimentDir, project, subject, session, partialsDir)
        except Exception, err:
            raise err
    ### DEBUGGING ###
    #print labels
    #print dirname
    #print project
    #print subject
    #print session
    ### END DEBUG ###
    # labels = map(formatLabel, labels) # convert shorthand to human-readable
    volume = 0.0
    for label in labels:
        volume += calculateBinaryVolume(dirname, label)
    return volume

