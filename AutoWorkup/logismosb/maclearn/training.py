"""
training.py
=================
Description:

Author:

Usage:

"""

import numpy
import os
import csv
import sys
import random
import ast
import SimpleITK as sitk
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import json
from .vesselness import compute_absolute_eigen_values


def get_list_of_features():
    """
    This function...

    :return:
    """
    from ..workflow import get_local_file_location

    _file = open(
        get_local_file_location(os.path.join("maclearn", "data_order.json")), "rb"
    )
    features = json.load(_file)
    _file.close()
    return features


def remove_keys_from_array(array, keys):
    """
    This function...

    :param array:
    :param keys:
    :return:
    """
    for key in keys:
        array.remove(key)
    return array


def mask_with_abc_image(image, abc_image):
    """
    This function...

    :param image:
    :param abc_image:
    :return:
    """
    abc_mask = get_brainmask(abc_image)
    masked_image = sitk.Mask(image, abc_mask)
    return masked_image


def binary_close(image, amount=1):
    """
    This function...

    :param image:
    :param amount:
    :return:
    """
    image = sitk.BinaryErode(image, amount)
    image = sitk.BinaryDilate(image, amount)
    return image


def get_brainmask(abc_image):
    """
    This function...

    :param abc_image:
    :return:
    """
    exclude_image = abc_image < 0
    exclude_codes = [5, 11, 12, 30]
    for code in exclude_codes:
        exclude_image = exclude_image + (abc_image == code)
    exclude_image = binary_close(exclude_image, 2)
    brainmask = abc_image * (exclude_image == 0) > 0

    return brainmask


def masked_image_array(image, mask):
    """
    This function...

    :param image:
    :param mask:
    :return:
    """
    return image_array(sitk.Mask(image, mask))


def mask_array_with_image(array, mask_image):
    """
    This function...

    :param array:
    :param mask_image:
    :return:
    """
    mask_array = image_array(mask_image)
    array[numpy.where(mask_array == 0)] = 0
    return array


def mask_data_with_image(data, mask_image):
    """
    This function...

    :param data:
    :param mask_image:
    :return:
    """
    for i, array in enumerate(data):
        data[i] = mask_array_with_image(array, mask_image)
    return data


def linear_array_from_image_file(image_file):
    """
    This function...

    :param image_file:
    :return:
    """
    image = sitk.ReadImage(image_file)
    return image_array(image)


def image_array(image):
    """
    Returns the 1D array of the numpy matrix
    :param image:
    :return:
    """
    a = sitk.GetArrayFromImage(image)
    a1D = a.reshape(a.size)
    return a1D


def data_by_region(
    data, wmtargets, wmlabelmap, wmlabels, gmtargets, gmlabelmap, gmlabels
):
    """
    Takes in an label map image and devides the data and
    targets into specified regions. Regoins are specified
    by a label list.

    :param data:
    :param wmtargets:
    :param wmlabelmap:
    :param wmlabels:
    :param gmtargets:
    :param gmlabelmap:
    :param gmlabels:
    :return:
    """
    columns = [data]
    keys = ["Features", "WMRegions", "GMRegions", "Targets"]

    wmregions = list()
    for i, label in enumerate(wmlabels):
        wmregions.append(pd.Series(wmlabelmap == label))
    df_wm = pd.concat(wmregions, axis=1, keys=wmlabels)

    gmregions = list()
    for i, label in enumerate(gmlabels):
        gmregions.append(pd.Series(gmlabelmap == label))
    df_gm = pd.concat(gmregions, axis=1, keys=gmlabels)

    df_targets = pd.concat(
        [pd.Series(wmtargets), pd.Series(gmtargets)], axis=1, keys=["WM", "GM"]
    )

    df = pd.concat([data, df_wm, df_gm, df_targets], axis=1, keys=keys)

    return df


def image_data(in_file, modality, abc_file=None, additional_images=None):
    """
    Computes the image features to be used for edge detection. Features are returned as a Pandas DataFrame.

    :param in_file: image file to be read in by SimpleITK
    :param modality: name of the modality
    :param abc_file:
    :param additional_images:
    :return:
    """

    # features can be added or taken out as to optimize the edge detection
    additional_feature_names = get_list_of_features()

    feature_value_arrays = []
    feature_names = []

    # intensity
    image = sitk.ReadImage(in_file, sitk.sitkFloat64)
    feature_value_arrays.append(image_array(image))
    feature_names.append("")

    # gradient magnitude
    feature_value_arrays.append(image_array(sitk.GradientMagnitude(image)))
    feature_names.append("GradMag")

    # second order gradient magnitude
    feature_value_arrays.append(
        image_array(sitk.GradientMagnitude(sitk.GradientMagnitude(image)))
    )
    feature_names.append("GradMag2")

    # Sobel
    feature_names.append("Sobel")
    feature_value_arrays.append(image_array(sitk.SobelEdgeDetection(image)))

    # eigenvalues of hessian
    feature_names.extend([f"Eigen{i}" for i in range(1, 4)])
    feature_value_arrays.extend(
        [eigen.flatten() for eigen in compute_absolute_eigen_values(image, sigma=0)]
    )

    # Laplacian
    feature_names.append("Laplacian")
    feature_value_arrays.append(
        image_array(sitk.Laplacian(image, useImageSpacing=True))
    )

    for sigma in [i * 0.5 for i in range(1, 7)]:
        sigma_str = f"{sigma:.1f}"
        feature_names.extend([f"GaussEigen{i}_{sigma_str}" for i in range(1, 4)])
        feature_value_arrays.extend(
            [
                eigen.flatten()
                for eigen in compute_absolute_eigen_values(image, sigma=sigma)
            ]
        )

        feature_names.append(f"GaussLaplacian_{sigma_str}")
        feature_value_arrays.append(
            image_array(sitk.LaplacianRecursiveGaussian(image, sigma=sigma))
        )

        feature_names.append(f"Gauss_{sigma_str}")
        feature_value_arrays.append(
            image_array(sitk.RecursiveGaussian(image, sigma=sigma))
        )

        feature_value_arrays.append(
            image_array(sitk.GradientMagnitudeRecursiveGaussian(image, sigma=sigma))
        )
        feature_names.append(f"GaussGradMag_{sigma_str}")

    feature_value_series = [pd.Series(array) for array in feature_value_arrays]
    keys = [modality + meas for meas in feature_names]

    for name in additional_feature_names:
        feature_value_series.append(
            pd.Series(image_array(sitk.ReadImage(additional_images[name])))
        )
        keys.append(name)

    data = pd.concat(feature_value_series, keys=keys, axis=1)

    return data


def get_graient_info(t1):
    """
    Takes in an image and computes the gradient, and hessian and returns
    the eigen values of the hessian.

    :param t1:
    :return:
    """
    grad = sitk.Gradient(t1)

    g_array = sitk.GetArrayFromImage(grad)

    gx = sitk.GetImageFromArray(g_array[:, :, :, 0])
    gy = sitk.GetImageFromArray(g_array[:, :, :, 1])
    gz = sitk.GetImageFromArray(g_array[:, :, :, 2])

    for img in [gx, gy, gz]:
        img.SetDirection(t1.GetDirection())
        img.SetOrigin(t1.GetOrigin())
        img.SetSpacing(t1.GetSpacing())

    ggx = sitk.Gradient(gx)
    ggy = sitk.Gradient(gy)
    ggz = sitk.Gradient(gz)

    ggx_array = sitk.GetArrayFromImage(ggx)
    ggy_array = sitk.GetArrayFromImage(ggy)
    ggz_array = sitk.GetArrayFromImage(ggz)

    hessian = numpy.stack((ggx_array, ggy_array, ggz_array), axis=3)
    eigvals = numpy.linalg.eigvals(hessian[:, :, :, :, :])

    gx_array = numpy.abs(numpy.ravel(g_array[:, :, :, 0]))
    gy_array = numpy.abs(numpy.ravel(g_array[:, :, :, 1]))
    gz_array = numpy.abs(numpy.ravel(g_array[:, :, :, 2]))

    eigen1 = numpy.ravel(eigvals[:, :, :, 0])
    eigen2 = numpy.ravel(eigvals[:, :, :, 1])
    eigen3 = numpy.ravel(eigvals[:, :, :, 2])

    return gx_array, gy_array, gz_array, eigen1, eigen2, eigen3


def multimodal_image_data(sample_dict):
    """
    Collects and Combines the image data from multiple modalities

    :param sample_dict:
    :return:
    """
    modals = sample_dict["Modalities"]
    if len(modals) > 1:
        data_list = list()
        for j, modal in enumerate(modals):  # iterate through the modalities
            data_list.append(image_data(sample_dict[modal], modal))
        df = pd.concat(data_list, axis=1)
    else:
        df = image_data(sample_dict[modals[0]], modals[0])
    return df


def collect_data(data_csv):
    """
    Collects the training data from a csv file.
    CSV header format must contain 'Truth', 'Labelmap', 'Labels', and
    'Modalities'.

    :param data_csv:
    :return:
    """

    data_samples = list()
    with open(data_csv, "rb") as csvfile:
        reader = csv.DictReader(csvfile)
        for i, line in enumerate(reader):
            try:

                if i == 0:
                    modalities = ast.literal_eval(line["Modalities"])
                    gmlabels = ast.literal_eval(line["GMLabels"])
                    wmlabels = ast.literal_eval(line["WMLabels"])
                else:
                    # check that the modalities and labels remain constant
                    if not modalities == ast.literal_eval(line["Modalities"]):
                        print(
                            "ERROR: csv line %d - Modalities must be the same for all subjects"
                            % i
                        )
                    elif not gmlabels == ast.literal_eval(line["GMLabels"]):
                        print(
                            "ERROR: csv line %d - GMLabels must be the same for all subjects"
                            % i
                        )
                    elif not wmlabels == ast.literal_eval(line["WMLabels"]):
                        print(
                            "ERROR: csv line %d - WMLabels must be the same for all subjects"
                            % i
                        )

                # replace the string representations with the literal representations
                line["Modalities"] = modalities
                line["GMLabels"] = gmlabels
                line["WMLabels"] = wmlabels

                data_samples.append(line)

            except KeyError as e:
                print(f"ERROR: csv line {i + 1} KeyError: {str(e)}")
                sys.exit()

    return data_samples


def split_data(data_samples, per_testing=0.1):
    """
    Split the data samples into training and testing sets.

    :param data_samples:
    :param per_testing:
    :return:
    """

    if per_testing < 0 or per_testing > 1:
        print("ERROR: Testing percentage must be between 0 and 1")
        sys.exit()

    n = len(data_samples)
    n_test = int(n * per_testing)  # will always round down

    # randomly shuffle the training samples
    train_samples = data_samples
    random.shuffle(data_samples)
    test_samples = list()
    while len(train_samples) > n - n_test:
        test_samples.append(train_samples.pop())

    return train_samples, test_samples


def combine_data(data_samples):
    """
    Takes the given data samples, reads in the images, and combines
    the image data and the targets to be used for classifier
    training.

    :param data_samples:
    :return:
    """

    df_list = list()
    id_list = list()

    for line in data_samples:
        # collect new data
        id_list.append(line["ID"])
        new_data = multimodal_image_data(line)

        # read in the target data
        gm_targets = image_array(sitk.ReadImage(line["GMEdges"]))
        wm_targets = image_array(sitk.ReadImage(line["WMEdges"]))

        # read in the label map
        wmlabelmap = image_array(sitk.ReadImage(line["WMLabelmap"]))
        wmlabels = line["WMLabels"]

        # read in the label map
        gmlabelmap = image_array(sitk.ReadImage(line["GMLabelmap"]))
        gmlabels = line["GMLabels"]

        # split the data by the labeled regions
        df_list.append(
            data_by_region(
                new_data,
                wm_targets,
                wmlabelmap,
                wmlabels,
                gm_targets,
                gmlabelmap,
                gmlabels,
            )
        )

    df_final = pd.concat(df_list, axis=0, keys=id_list)

    return df_final


def get_labeled_region_data(t_data, rg_name, label, matter):
    """
    This function...

    :param t_data:
    :param rg_name:
    :param label:
    :param matter:
    :return:
    """
    # Training
    t_index = t_data[rg_name][label]
    t_targets = t_data["Targets"][matter][t_index].values
    t_feat = t_data["Features"][t_index].values

    return t_feat, t_targets


def train_classifier(
    data, targets, out_file, clf=RandomForestClassifier(), n_jobs=-1, load_clf=True
):
    """
    Trains the classifier and dumps the pickle file

    :param data:
    :param targets:
    :param out_file:
    :param clf:
    :param n_jobs:
    :param load_clf:
    :return:
    """
    if os.path.isfile(out_file):
        print(f"Found classifier {out_file}")
        if not load_clf:
            return
        clf = joblib.load(out_file)
    else:
        print(f"Fitting classifier {out_file}")
        clf.n_jobs = n_jobs
        clf.fit(data, targets)
        joblib.dump(clf, out_file)
    return clf


def run_training(training_data, train_base_clf=False, out_dir=".", n_jobs=-1):
    """
    This function...

    :param training_data:
    :param train_base_clf:
    :param out_dir:
    :param n_jobs:
    :return:
    """
    all_training_features = training_data["Features"].values
    classifiers = dict()
    for matter in ["WM", "GM"]:
        classifiers[matter] = dict()

        print(f"Training {matter}")
        # Get WM training targets
        train_matter_targets = training_data["Targets"][matter].values

        if train_base_clf:
            base_clf_file = os.path.join(out_dir, f"{matter}BaseCLF.pkl")
            train_classifier(
                all_training_features,
                train_matter_targets,
                base_clf_file,
                n_jobs=n_jobs,
                load_clf=False,
            )
            classifiers[matter]["NonRegional"] = base_clf_file

        rg_name = matter + "Regions"

        classifiers[matter]["Regional"] = dict()

        for label in training_data[rg_name].columns:

            # get label specific data
            label_train_features, label_train_targets = get_labeled_region_data(
                training_data, rg_name, label, matter
            )

            # train regional classifier
            regional_clf_file = os.path.join(out_dir, f"{matter}{label}RegionalCLF.pkl")
            train_classifier(
                label_train_features,
                label_train_targets,
                regional_clf_file,
                n_jobs=n_jobs,
                load_clf=False,
            )

            classifiers[matter]["Regional"][label] = regional_clf_file

    return classifiers
