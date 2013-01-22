import numpy as np
import os.path

import SimpleITK as sitk

from ..config import _config
labels = ['caudate', 'putamen', 'hippocampus', 'thalamus', 'accumben', 'globus']

def constructLabels(labels):
    numbers = range(1,((len(labels)*2) + 1))
    full_labels = []
    index = 0
    for label in labels:
        full_labels.append('_'.join(['left', label]))
        full_labels.append('_'.join(['right', label]))
    return full_labels, numbers

def _moduleCreateLabels(labels):
    full_labels, numbers = constructLabels(labels)
    labelMap = map(None, full_labels, numbers)
    return dict(labelMap) # Use this variable


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
    label = '_'.join([side,anatomy])
    return label


def calculateLabelVolume(dirname, label):
    labelFile = os.path.join(dirname, _config.get('Results', 'segmentations'), label + '_seg.nii.gz')
    assert os.path.exists(labelFile), "File not found: %s" % labelFile
    image = sitk.ReadImage(labelFile)
    nda = sitk.GetArrayFromImage(image)
    maskSum = nda.sum()
    print maskSum
    size = image.GetSpacing()
    print size
    return maskSum * size[0] * size[1] * size[2]


def getVolume(args=[], kwds={}):
    dirname = labels = project = subject = session = experimentDir = None
    experimentDir = _config.get('Results', 'directory') ### HACK
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
        volume += calculateLabelVolume(dirname, label)
    return volume
