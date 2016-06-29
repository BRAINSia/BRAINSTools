from script import get_data, nm_dir
from training import run_training
import cPickle as pickle
import os


cache_dir = "/Shared/sinapse/CACHE/20160608_RF_Classifiers_Abs_Gradient"
if not os.path.isdir(cache_dir):
    os.makedirs(cache_dir)
print("Getting Data")
data = get_data(os.path.join(cache_dir, "fs_norm_data.hdf5"), nm_dir, overwrite=True, out_dir=cache_dir)
# data = get_data(os.path.join("/Shared/sinapse/CACHE/20160606_RF_Classifiers", "fs_norm_data.hdf5"), nm_dir,
#                overwrite=False, out_dir=cache_dir)
print("Training Classifiers")
classifiers = run_training(data, train_base_clf=False, out_dir=cache_dir, n_jobs=8)
print("Saving Classifiers")
pickle.dump(classifiers, open(os.path.join(cache_dir, "ClassifierDictionary.pkl"), 'wb'))
