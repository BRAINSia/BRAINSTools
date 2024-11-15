"""
train_classifiers.py
======================
Description:

Author:

Usage:

"""

from .training import train_classifier
from .crossvalidation import read_data
import os


def train_classifiers(data_file, cache_dir=os.path.curdir):
    """
    This function...

    :param data_file:
    :param cache_dir:
    """
    data = read_data(data_file)
    for matter in ["WM", "GM"]:
        classifier_file = os.path.join(
            cache_dir, "Classifier", f"{matter}_matter_classifier.pkl"
        )
        if not os.path.exists(os.path.dirname(classifier_file)):
            os.makedirs(os.path.dirname(classifier_file))
        train_classifier(
            data["Features"].values, data["Truth"][matter].values, classifier_file
        )


if __name__ == "__main__":
    train_classifiers(
        "/Shared/sinapse/CACHE/20161105_Davids_CrossValidation/training_data.hdf5",
        "/Shared/sinapse/CACHE/20161105_Davids_CrossValidation",
    )
