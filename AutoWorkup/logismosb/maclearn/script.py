import preprocess
import training
import testing
import os
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial
import matplotlib.pyplot as plt
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.externals import joblib
import datetime
import cPickle as pickle
from scipy import interp
from training import run_training, get_labeled_region_data

nm_dir = "/Shared/johnsonhj/HDNI/ReferenceData/Neuromorphometrics/20141116_Neuromorphometrics_base_Results/Neuromorphometrics/2012Subscription/"

def plot_feature_importances(list_of_importances, feature_names, out_file=None, title="Feature Importances"):
    std = np.std(list_of_importances, axis=0)
    importances = np.average(list_of_importances, axis=0)
    indices = np.argsort(importances)[::-1]
    num_features = len(feature_names)
    import random
    fig_num = random.randint(1, 10000)

    # Print the feature ranking
    print("Feature ranking:")

    for f in range(num_features):
        print("%d. feature %s (%f)" % (f + 1, feature_names[indices][f], importances[indices[f]]))
    # Plot the feature importances of the forest
    fig = plt.figure(fig_num)
    # We define a fake subplot that is in fact only the plot.
    plot = fig.add_subplot(111)

    # We change the fontsize of minor ticks label
    plot.tick_params(axis='x', which='major', labelsize=15)
    plot.tick_params(axis='y', which='major', labelsize=12)
    plt.title(title, fontsize=14)
    plt.bar(range(num_features), importances[indices],
            color="r", yerr=std[indices], align="center")
    plt.xticks(range(num_features), feature_names[indices])
    plt.xlim([-1, num_features])
    plt.tight_layout()
    if not out_file:
        plt.show()
    else:
        plt.savefig(out_file)


def plot_all_feature_importances(all_importances, feature_names, out_file=None, title=None):
    # get averages for all features
    imp_list = list()
    for matter in all_importances.iterkeys():
        for label in all_importances[matter].iterkeys():
            imp_list.append(np.average(all_importances[matter][label]["Regional"], axis=0))

    plot_feature_importances(imp_list, feature_names, out_file, title)


def split_region_data(train_data, test_data, matter, rg_name, label, n_jobs=-1):
    """
    Splits the training and testing data and returns the region specific training and testing features plus targets
    """
    print("Label Region: {0}".format(label))
    # Training
    train_index = train_data[rg_name][label]
    train_targets = train_data['Targets'][matter][train_index].values
    train_feat = train_data['Features'][train_index].values

    # Testing
    test_index = test_data[rg_name][label]
    test_targets = test_data['Targets'][matter][test_index].values
    test_feat = test_data['Features'][test_index].values

    return train_feat, train_targets, test_feat, test_targets


def train_classifier(train_features, train_targets, n_jobs=-1,
                     clf=RandomForestClassifier()):
    clf.n_jobs = n_jobs
    clf.fit(train_features, train_targets)
    return clf


def test_classifier(clf, test_features, test_targets):
    # Predictions
    probas = clf.predict_proba(test_features)
    score_roc = roc_curve(test_targets, probas[:, 1], pos_label=1)
    return score_roc

cache_dir = "/Shared/sinapse/CACHE/20160510_EdgeDetection"


def runcrossval(idx_split, data_file):
    """
    Runs cross validation after the data splits
    """
    n_jobs = 8
    print("Training Classifiers")

    print("Reading in data")
    # read in data
    data = pd.read_hdf(data_file)

    print("Splitting the indices")
    # get training and testing indices
    train_idx = idx_split[0]
    test_idx = idx_split[1]

    print("Getting the list of ids")
    # get a list of all the ids
    ids = data.index.levels[0].values
    train_ids = ids[train_idx].tolist()
    test_ids = ids[test_idx].tolist()

    print("Creating the fold directory")
    # define the fold directory filepath
    num1, num2 = np.where(np.in1d(ids, test_ids))[0]
    fold_id = "{0}{1}".format(num1, num2)
    fold_dir = os.path.join(cache_dir, fold_id)
    if not os.path.isdir(fold_dir):
        os.makedirs(fold_dir)

    print("Training subjects: {0}".format(train_ids))
    print("Testing subjects: {0}".format(test_ids))

    train_data = data.loc[train_ids]

    classifiers = run_training(train_data, train_base_clf=True, out_dir=fold_dir, n_jobs=n_jobs)

    return test_idx, classifiers


hdf5_file = os.path.join(cache_dir, "alldata.hdf5")
partial_runcrossval = partial(runcrossval, data_file=hdf5_file)


def make_empty_dictionaries(labels):
    import copy
    roc_scores = dict()
    roc_scores_mean = dict()
    roc_auc = dict()
    pixel_counts = dict()

    for matter in ['WM', 'GM']:
        roc_scores[matter] = dict()
        roc_scores_mean[matter] = dict()
        roc_auc[matter] = dict()
        pixel_counts[matter] = dict()

        for label in labels[matter]:
            pixel_counts[matter][label] = list()
            roc_scores[matter][label] = dict()
            roc_scores_mean[matter][label] = dict()
            roc_auc[matter][label] = dict()

            for clf_type in ["Regional", "NonRegional"]:
                roc_scores[matter][label][clf_type] = dict()
                roc_scores_mean[matter][label][clf_type] = dict()
                roc_auc[matter][label][clf_type] = list()

                for pr in ['fpr', 'tpr']:
                    roc_scores[matter][label][clf_type][pr] = list()
                    roc_scores_mean[matter][label][clf_type][pr] = list()

    importances = copy.deepcopy(roc_auc)

    return roc_scores, roc_scores_mean, roc_auc, importances


def get_data(data_file, nm_dir, overwrite=False, out_dir=None):
    # preprocessing
    if not os.path.isfile(data_file) or overwrite:
        if not out_dir:
            out_dir = cache_dir
        csv_file = preprocess.createdatacsv(nm_dir, out_dir, overwrite=overwrite)
        data_samples = training.collectdata(csv_file)
        data = training.combinedata(data_samples)
        preprocess.save_data_frame(data, data_file)
    else:
        if ".hdf" in data_file:
            data = pd.read_hdf(data_file)
        elif ".csv" in data_file:
            data = pd.read_csv(data_file)
        else:
            print("Could not determine reader for data type. Data not read.")
            data = None
    return data


def main():
    overwrite = False
    nm_dir = "/Shared/johnsonhj/HDNI/ReferenceData/Neuromorphometrics/20141116_Neuromorphometrics_base_Results/Neuromorphometrics/2012Subscription/"

    data = get_data(hdf5_file, nm_dir, overwrite=overwrite)

    # Learning

    # Cross Validation
    folds = 10

    # Get the list of subject_ids
    ids = data.index.levels[0].values

    # Get split of training and testing indexes
    kfolds = cross_validation.KFold(len(ids), folds, shuffle=True)
    data_splits = [(x[0], x[1]) for x in kfolds]

    pool = mp.Pool(processes=4)
    results = pool.map(partial_runcrossval, data_splits, 1)
    pool.close()
    pool.join()

    # make empty dictionaries
    labels = dict(WM=data['WMRegions'].columns, GM=data['GMRegions'].columns)
    roc_scores, roc_scores_mean, roc_auc, importances = make_empty_dictionaries(labels)

    # Get ROC scores and feature importance from classifiers
    for test_idx, clf_files in results:
        test_ids = ids[test_idx].tolist()
        test_data = data.loc[test_ids]
        for matter in ['WM', 'GM']:

            # Load classifier for either WM or GM
            base_clf_file = clf_files[matter]['NonRegional']
            base_clf = joblib.load(base_clf_file)

            rg_name = matter + 'Regions'

            for label in data[rg_name].columns:

                # get data for that labeled region
                label_test_features, label_test_targets = get_labeled_region_data(test_data, rg_name, label, matter)

                # load classifier for given region
                regional_clf_file = clf_files[matter]['Regional'][label]
                regional_clf = joblib.load(regional_clf_file)

                for clf_type, clf in [("Regional", regional_clf), ("NonRegional", base_clf)]:

                    # Get TPR and FPR scores for ROC analysis
                    fpr, tpr, _ = test_classifier(clf,
                                                  label_test_features,
                                                  label_test_targets)
                    roc_scores[matter][label][clf_type]['fpr'].append(fpr)
                    roc_scores[matter][label][clf_type]['tpr'].append(tpr)

                    # Get feature importance
                    importances[matter][label][clf_type].append(clf.feature_importances_)

    pickle.dump(roc_scores, open(os.path.join(cache_dir, "roc_scores.pkl"), 'wb'))
    pickle.dump(importances, open(os.path.join(cache_dir, "clf_importances.pkl"), 'wb'))

    mean_fpr = np.linspace(0, 1, 100)
    for j, matter in enumerate(['WM', 'GM']):
        rg_name = matter + 'Regions'
        for i, label in enumerate(data[rg_name].columns):
            for clf_type in ["Regional", "NonRegional"]:

                mean_tpr = 0.0

                for ii, tpr in enumerate(roc_scores[matter][label][clf_type]['tpr']):
                    fpr = roc_scores[matter][label][clf_type]['fpr'][ii]
                    mean_tpr += interp(mean_fpr, fpr, tpr)
                    mean_tpr[0] = 0.0

                mean_tpr /= len(roc_scores[matter][label][clf_type]['tpr'])
                roc_scores_mean[matter][label][clf_type]['tpr'] = mean_tpr
                roc_auc[matter][label][clf_type] = auc(mean_fpr, mean_tpr)
                plt.figure((i + 1) + j * 10, dpi=200)
                plt.plot([0, 1], [0, 1], 'k--')

                # Plot the classifier scores without the regions as baseline
                plt.plot(mean_fpr,
                         roc_scores_mean[matter][label]["NonRegional"]["tpr"],
                         'k', label="NonRegional", linewidth=2)
                plt.plot(mean_fpr,
                         roc_scores_mean[matter][label]["Regional"]["tpr"],
                         label="Regional", linewidth=2)
                fontsize = 16
                plt.xlabel('False Positive Rate', fontsize=fontsize)
                plt.ylabel('True Positive Rate', fontsize=fontsize)
                # plt.title('ROC Curve for {0} Matter Classification in {1}'.format(matter, label))
                plt.legend(loc='best', fontsize=fontsize)
                plt.savefig(os.path.join(cache_dir, '{0}{1}ROC.eps'.format(matter, label)), pad_inches=0)

    pickle.dump(roc_auc, open(os.path.join(cache_dir, "roc_auc_scores.pkl"), 'wb'))

    feat_names_long = np.array([r'T1', r'$\|\nabla\|$', r'$\|\nabla\|\nabla\|\|$',
                                r'$\nabla_x$', r'$\nabla_y$', r'$\nabla_z$',
                                r'$\lambda_1$', r'$\lambda_2$', r'$\lambda_3$'])
    plot_all_feature_importances(importances, feat_names_long,
                                 title="Feature Importance for Region Based Classifiers",
                                 out_file=os.path.join(cache_dir, "RegionalFeatureImportances.eps"))


if __name__ == "__main__":
    print(datetime.datetime.utcnow())
    main()
