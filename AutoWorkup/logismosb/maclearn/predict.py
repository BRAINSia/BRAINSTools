import SimpleITK as sitk
from sklearn.externals import joblib
import os
import numpy as np
from preprocess import createwatersheds
from training import image_data
import pickle


def get_image_probabilities(clf_dict_file, image_file, aseg_labels_file, filled_labels_file, out_dir=os.getcwd()):
    clf_files = pickle.load(open(clf_dict_file, 'rb'))
    image_data = image_data(image_file, "T1")
    image_features = image_data.values
    wm_labels_image, gm_labels_image, _, _ = createwatersheds(aseg_labels_file, filled_labels_file)
    wm_labels_file = os.path.abspath(os.path.join(out_dir, "wm_regions.nii.gz"))
    sitk.WriteImage(wm_labels_image, wm_labels_file)
    gm_labels_file = os.path.abspath(os.path.join(out_dir, "gm_regions.nii.gz"))
    sitk.WriteImage(gm_labels_image, gm_labels_file)
    wm_proba_file = get_proba_image_file(image_features, clf_files['WM'], image_file, wm_labels_file,
                                         os.path.join(out_dir, "wm_probas.nii.gz"))
    gm_proba_file = get_proba_image_file(image_features, clf_files['GM'], image_file, gm_labels_file,
                                         os.path.join(out_dir, "gm_probas.nii.gz"))
    return wm_proba_file, gm_proba_file, wm_labels_file, gm_labels_file


def get_proba_image_file(features, clf_files, image_file, labels_file, proba_file):
    image = sitk.ReadImage(image_file)
    label_image = sitk.ReadImage(labels_file)
    labels_array = sitk.GetArrayFromImage(label_image)
    prob_array = predict_image_proba(clf_files, features, labels_array).reshape(labels_array.shape)
    prob_image = image_from_array_with_reference_image(prob_array, image)
    sitk.WriteImage(prob_image, proba_file)
    return os.path.abspath(proba_file)


def predict_image_proba(clf_files, features, label_array):
    prob_array = np.zeros(label_array.size)
    for label in clf_files['Regional'].keys():
        idx = label_array.flatten() == label
        clf = joblib.load(clf_files['Regional'][label])
        prob_array[idx] = clf.predict_proba(features[idx])[:, 1]
    return prob_array


def image_from_array(prob_array, origin, spacing, direction):
    prob_image = sitk.GetImageFromArray(prob_array)
    prob_image.SetOrigin(origin)
    prob_image.SetSpacing(spacing)
    prob_image.SetDirection(direction)
    return prob_image


def image_file_from_array_with_reference_image_file(prob_array, reference_image_file, out_file):
    reference_image = sitk.ReadImage(reference_image_file)
    probability_image = image_from_array_with_reference_image(prob_array, reference_image)
    sitk.WriteImage(probability_image, os.path.abspath(out_file))
    return os.path.abspath(out_file)


def image_from_array_with_reference_image(prob_array, reference_image):
    array_3d = shape_array_like_image(prob_array, reference_image)
    prob_image = image_from_array(array_3d, origin=reference_image.GetOrigin(), spacing=reference_image.GetSpacing(),
                                  direction=reference_image.GetDirection())
    return prob_image


def shape_array_like_image(array, image):
    return array.reshape(get_shape_from_image(image))


def get_shape_from_image(image):
    array = sitk.GetArrayFromImage(image)
    return array.shape
