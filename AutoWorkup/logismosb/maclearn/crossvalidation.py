import os
import pandas as pd
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.externals import joblib
import cPickle as pickle
import sys


def get_idx_folds(subject_ids, num_folds):
    # Get the split of training and testing indexes
    kfolds = cross_validation.KFold(len(subject_ids), num_folds, shuffle=True)
    return [(x[0], x[1]) for x in kfolds]


def get_subject_id_folds(subject_ids, num_folds):
    splits = get_idx_folds(subject_ids, num_folds)
    return [(subject_ids[train_idx].tolist(), subject_ids[test_idx].tolist()) for train_idx, test_idx in splits]


def get_subject_id_folds_from_data(data, num_folds):
    subject_ids = get_subject_ids(data)
    return get_subject_id_folds(subject_ids, num_folds)


def get_data_with_subject_ids(data, subject_ids):
    return data[data.index.get_level_values(0).isin(subject_ids)]


def get_subject_ids(data):
    return data.index.levels[0].values


def read_data(data_file):
    if ".hdf" in data_file:
        return pd.read_hdf(data_file)
    else:
        sys.exit()


def get_training_and_testing_data(data, fold):
    training_data = get_data_with_subject_ids(data, fold[0])
    testing_data = get_data_with_subject_ids(data, fold[1])
    return training_data, testing_data


def get_truth_from_data(data, matter):
    return data['Truth'][matter].values


def get_features_from_data(data):
    return data['Features'].values


def run_cross_validation_fold(data, fold, output_dir):
    if not os.path.isdir(output_dir):
        print("making output directory: {0}".format(output_dir))
        os.mkdir(output_dir)
    print("splitting training and testing data for fold: {0}".format(fold))
    training_data, testing_data = get_training_and_testing_data(data, fold)
    print("getting data features")
    training_features = get_features_from_data(training_data)
    testing_features = get_features_from_data(testing_data)
    for matter in ["WM", "GM"]:
        clf_file = os.path.join(output_dir, "{0}_classifier.pkl".format(matter))
        print("Training {0} classifier".format(matter))
        clf = train_classifier(training_features, get_truth_from_data(training_data, matter),
                               out_file=clf_file, n_jobs=8)
        print("Getting ROC scores for classifier")
        roc = test_classifier(clf, testing_features, get_truth_from_data(testing_data, matter))

        roc_out_file_name = os.path.join(output_dir, "{0}_roc.pkl".format(matter))
        roc_out_file = open(roc_out_file_name, "wb")
        print("Writing ROC scores to file {0}".format(roc_out_file_name))
        pickle.dump(roc, roc_out_file)


def run_nfold_cross_validation(data_file, nfolds=10, output_dir=os.path.curdir):
    print("reading: {0}".format(data_file))
    data = read_data(data_file)
    print("splitting data into {0} folds".format(nfolds))
    folds = get_subject_id_folds_from_data(data, nfolds)
    for i in range(nfolds):
        print("Fold number: {0}".format(i))
        run_cross_validation_fold(data=data, fold=folds[i], output_dir=os.path.join(output_dir, str(i)))


def train_classifier(train_features, train_targets, n_jobs=-1,
                     clf=RandomForestClassifier(), out_file=None):
    clf.n_jobs = n_jobs
    clf.fit(train_features, train_targets)
    if out_file:
        joblib.dump(clf, out_file)
    return clf


def test_classifier(clf, test_features, test_targets):
    # Predictions
    probas = clf.predict_proba(test_features)
    score_roc = roc_curve(test_targets, probas[:, 1], pos_label=1)
    return score_roc


if __name__ == "__main__":
    run_nfold_cross_validation(data_file="/Shared/sinapse/CACHE/20161025_Davids_CrossValidation/training_data.hdf5",
                               nfolds=10,
                               output_dir="/Shared/sinapse/CACHE/20161025_Davids_CrossValidation/")
