from preprocess import save_data_frame
import os
import pickle
import pandas as pd
from training import linear_array_from_image_file, image_data


def pickle_load(pickled_file):
    _file = open(pickled_file, "rb")
    output = pickle.load(_file)
    _file.close()
    return output


def get_subject_id_from_t1(t1_file):
    return os.path.abspath(t1_file).split("/")[-3]


def collect_training_data(training_files):
    all_data = list()
    for t1_file, additional_files, truth_files in training_files:
        feature_data = image_data(t1_file, "T1", additional_images=additional_files)
        subject_id = get_subject_id_from_t1(t1_file)
        index = pd.MultiIndex.from_tuples([(subject_id, i) for i in feature_data.index])
        feature_data.index = index
        gm_truth_data = pd.Series(linear_array_from_image_file(truth_files["gm"]), name="GM", index=index)
        wm_truth_data = pd.Series(linear_array_from_image_file(truth_files["wm"]), name="WM", index=index)
        truth_data = pd.concat([gm_truth_data, wm_truth_data], axis=1)
        data = pd.concat([feature_data, truth_data], axis=1, keys=["Features", "Truth"])
        all_data.append(data)
    return pd.concat(all_data, axis=0)


def save_training_data():
    cache_dir = "/Shared/sinapse/CACHE/20161105_Davids_CrossValidation"
    training_files = pickle_load(os.path.join(cache_dir, "training_files.pkl"))
    training_data = collect_training_data(training_files)
    training_data_file = os.path.join(cache_dir, "training_data.hdf5")
    save_data_frame(training_data, training_data_file)

if __name__ == "__main__":
    save_training_data()