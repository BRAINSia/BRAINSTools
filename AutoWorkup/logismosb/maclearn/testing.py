import SimpleITK as sitk
import os
import numpy
import training

def run_test(clf, sample_dict, label, out_file):
    """
    Test the classifier on the trainign data.
    """
    # read in image data
    data = training.multimodalimagedata(sample_dict)

    # read in the target image
    target_image = sitk.ReadImage(sample_dict["Truth"])
    targets = training.imagearray(target_image)

    # read in the label map
    labelmap = sitk.ReadImage(sample_dict["Labelmap"])

    # split the data by regions
    data_dict, targets_dict, index_dict = training.databyregion(data, targets, labelmap, sample_dict["Labels"])

    # iterate through the labels and make predictions
    proba_file = predict(clf, data_dict[label], target_image, out_file, index_dict[str(label)])
    return proba_file


def predict(clf, data, in_image, out_file, index, neg_proba=False):
    """Make predictions and write predictions to an image file"""
    # TODO: Add scoring to the predictions
    print("Making predictions")
    image_size = sitk.GetArrayFromImage(in_image).size
    image_shape = sitk.GetArrayFromImage(in_image).shape
    data_prob = clf.predict_proba(data)
    predict_prob = numpy.zeros(image_size)
    predict_prob_neg = numpy.zeros(image_size)
    predict_prob[index] = data_prob[:, 1]
    predict_prob_neg[index] = data_prob[:, 0]

    print("Writing probability predictions to file {0}".format(out_file))
    if os.path.isfile(out_file):
        return out_file
    else:
        out_dir = os.path.dirname(out_file)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        if not neg_proba:
            prob_array = predict_prob.reshape(image_shape)
            prob_image = sitk.GetImageFromArray(prob_array)
            prob_image.SetOrigin(in_image.GetOrigin())
            prob_image.SetSpacing(in_image.GetSpacing())
            prob_image.SetDirection(in_image.GetDirection())
            sitk.WriteImage(prob_image, out_file)
        else:
            prob_array_neg = predict_prob_neg.reshape(nm_shape)
            prob_image_neg = sitk.GetImageFromArray(prob_array_neg)
            prob_image_neg.SetOrigin(nm_image.GetOrigin())
            prob_image_neg.SetSpacing(nm_image.GetSpacing())
            prob_image_neg.SetDirection(nm_image.GetDirection())
            sitk.WriteImage(prob_image, out_file)

        return out_file


def combinepredictions(predictions, in_image, out_file):
    """
    Takes in a list of the prediction images that are split by region
    and combines them to make one prediction image.
    """
    out_image = in_image < 0
    for pred in predictions:
        pred_image = sitk.ReadImage(pred)
        out_image = sitk.Cast(out_image, sitk.sitkFloat64) + pred_image
    sitk.WriteImage(out_image, out_file)
