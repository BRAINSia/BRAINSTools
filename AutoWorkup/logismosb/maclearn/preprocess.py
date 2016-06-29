import os
import SimpleITK as sitk
import pandas as pd


def largest_connected_component(image, minSize=1000):
    return sitk.RelabelComponent(sitk.ConnectedComponent(image), minSize) == 1


def fill_mask_holes(image, minSize=0):
    negConnected = largest_connected_component(image == 0, minSize)
    return negConnected == 0


def getedges(nm_file, between_labels=False):
    """Method to extract edges from a neuromorphometrics label map."""
    nm_map = sitk.ReadImage(nm_file, sitk.sitkUInt32)

    # define labels not to include
    # cerebellum and brainstem labels
    background_labels = [0, 4, 11, 32, 33, 34, 35, 38, 39, 40, 41, 46, 61, 62, 71, 72, 73]
    background = sitk.Cast(nm_map * 0, sitk.sitkUInt8)
    for background_label in background_labels:
        background = sitk.Or(nm_map == background_label, background)

    gm_mask = sitk.BinaryThreshold(nm_map, lowerThreshold=99.5, upperThreshold=1000)
    wm_mask = fill_mask_holes(sitk.BinaryThreshold(nm_map, 0.5, 1000) * (gm_mask == 0) * (background == 0))
    # wm_mask = largest_connected_component(wm_mask)

    # get the wm edges/contours
    wm_contours = sitk.LabelContour(wm_mask)

    if between_labels:
        # get the gm edges/contours
        gm_map = sitk.Cast(nm_map, sitk.sitkUInt16) * sitk.Cast(gm_mask, sitk.sitkUInt16)
        gm_label_contours = sitk.LabelContour(gm_map) > 0
        select_gm = gm_label_contours * (gm_label_contours * wm_contours == 0)
        rm_edges = sitk.BinaryDilate(wm_mask, 1)
        gm_contours = select_gm * (rm_edges == 0)
    else:
        gm_contours = (sitk.LabelContour(sitk.Or(gm_mask, wm_mask)) > 0) * (wm_contours == 0)

    return wm_contours, gm_contours


def close_mask(mask, close_param=12):
    mask_dilated = sitk.BinaryDilate(mask, close_param)
    mask_eroded = sitk.BinaryErode(mask_dilated, close_param)
    return mask_eroded, mask_dilated


def get_folds_mask(filled):
    folds_label = 100
    lh_wm_mask = (filled == 255)
    rh_wm_mask = (filled == 127)
    lh_folds, lh_dilated = close_mask(lh_wm_mask)
    rh_folds, rh_dilated = close_mask(rh_wm_mask)
    folds_mask = (lh_folds + rh_folds) > 0
    dilated_mask = (lh_dilated + rh_dilated) > 0
    out_labels = [folds_label, folds_label * 2]
    label_map = ((folds_mask + dilated_mask) * folds_label)
    return label_map, out_labels


def remove_cerebellum(label_map, aseg_label_map, exclude_labels=None):
    dilated_labels = MultilabelDilation(aseg_label_map, radius=4)
    if not exclude_labels:
        exclude_labels = [16, 47, 46, 7, 8]
    cerebellum = aseg_label_map < 0
    for label in exclude_labels:
        cerebellum = cerebellum + (dilated_labels == label)
    return label_map * (cerebellum == 0)


# As advised on the ITK mailing list, label dilation can be implemented via
# distance transforms and watershed transforms. This algorithm is illustrated
# in SimpleITK python code below (courtesy of Bradely Lowekamp)
def MultilabelDilation(img, radius=1, kernel=None):
    if not kernel:
        kernel = sitk.BinaryDilateImageFilter.Ball
    dilatImg = sitk.BinaryDilate(img != 0, radius, kernel)
    wsImg = create_label_watershed(img)
    return sitk.Cast(dilatImg, wsImg.GetPixelID()) * wsImg


def create_label_watershed(labels_image, markWatershedLine=False):
    distImg = sitk.SignedMaurerDistanceMap(labels_image != 0,
                                           insideIsPositive=False,
                                           squaredDistance=False,
                                           useImageSpacing=False)
    wsImg = sitk.MorphologicalWatershedFromMarkers(distImg, labels_image, markWatershedLine=markWatershedLine)
    return wsImg


def cast_to_int(image):
    filt = sitk.StatisticsImageFilter()
    filt.Execute(image)
    image_min = filt.GetMinimum()
    image_max = filt.GetMaximum()
    if image_min >= 0:
        if image_max < 256:
            return sitk.Cast(image, sitk.sitkUInt8)
        elif image_max < 65535:
            return sitk.Cast(image, sitk.sitkUInt16)
    print("Could not cast image")
    print("Max: {0}".format(image_max))
    print("Min: {0}".format(image_min))
    return image


def createwatersheds(aseg_file, filled_file, dilation=6):
    filled = cast_to_int(sitk.ReadImage(filled_file))
    fs = cast_to_int(sitk.ReadImage(aseg_file))
    wm_labels = [10, 11, 12]

    # wm section 11
    # 49 == rh thalamus
    # 53 == rh hippocamus
    # 54 == rh amygdala
    rh_hipthal = ((fs == 49) + (fs == 53) + (fs == 54)) * wm_labels[1] # assigns label 11
    lh_hipthal = ((fs == 10) + (fs == 17) + (fs == 18)) * wm_labels[1] # assigns label 11

    # wm section 2
    # 51 == rh putamen
    # 52 == rh pallidum
    # dilate putamen/pallidum and assign label 12
    rh_pitpul = (
        (sitk.DilateObjectMorphology((fs == 51) + (fs == 52), dilation) * (rh_hipthal == 0)) > 0) * wm_labels[2]
    lh_pitpul = (
        (sitk.DilateObjectMorphology((fs == 12) + (fs == 13), dilation) * (lh_hipthal == 0)) > 0) * wm_labels[2]

    # get the gm map
    gm_map, gm_labels = get_folds_mask(filled)
    gm_map = remove_cerebellum(gm_map, fs)
    gm_mask = gm_map > 0

    wm_regions = rh_hipthal + rh_pitpul + lh_hipthal + lh_pitpul
    wm_map = ((wm_regions == 0) * gm_mask * wm_labels[0]) + wm_regions

    return wm_map, gm_map, wm_labels, gm_labels


def write_image(image, out_file, overwrite=False):
    if not (os.path.isfile(out_file) and not overwrite):
        sitk.WriteImage(image, out_file)


def write_edges(nm_file, out_dir, overwrite=False):
    if overwrite or not check_for_edge_files(out_dir):
        wm_edges, gm_edges = getedges(nm_file)
        for edge, name in [(wm_edges, "wm"), (gm_edges, "gm")]:
            write_image(edge, get_edge_file_name(name, out_dir), overwrite=overwrite)
    return get_edge_file_paths(out_dir)


def get_edge_file_paths(out_dir):
    edge_file_paths = dict()
    for matter in ("wm", "gm"):
        edge_file_paths[matter] = get_edge_file_name(matter, out_dir)
    return edge_file_paths


def check_for_edge_files(directory):
    edge_file_paths = get_edge_file_paths(directory)
    for matter in edge_file_paths:
        if not os.path.isfile(edge_file_paths[matter]):
            return False
    return True


def get_edge_file_name(edge_name, out_dir):
    return os.path.join(out_dir, "{0}_edge.nii.gz".format(edge_name))


def get_nm_fs_scan(nm_file, name='norm.mgz'):
    subject_id = get_nm_subject_id(nm_file)
    nm_fs_t1 = "/Shared/johnsonhj/HDNI/20150206_FS_Neuromorphometric/{0}/mri/{1}".format(subject_id, name)
    return nm_fs_t1


def convert_fs_scan(in_file, out_file, resample_type=None):
    if not os.path.isfile(out_file):
        from nipype.interfaces.freesurfer import MRIConvert
        convert = MRIConvert()
        convert.inputs.in_file = in_file
        convert.inputs.out_file = out_file
        convert.inputs.out_orientation = "LPS"
        convert.inputs.conform = True
        convert.inputs.no_change = True
        if resample_type:
            convert.inputs.resample_type = resample_type
        convert.run()
    return os.path.abspath(out_file)


def get_nm_subject_id(nm_file):
    subject_id = nm_file.split("/")[-3]
    return subject_id


def get_nm_t1(nm_file, fs=False, cache_dir=os.getcwd()):
    if fs:
        in_t1 = get_nm_fs_scan(nm_file)
        out_t1 = convert_fs_scan(in_t1, os.path.join(cache_dir, "t1_{0}.nii.gz".format(get_nm_subject_id(nm_file))))
    else:
        out_t1 = os.path.join(os.path.dirname(nm_file), "t1_average_BRAINSABC.nii.gz")
    return out_t1


def createdatacsv(in_dir, cache_dir, overwrite=False, file_list=None):
    """
    Function that is specific to the neuromorphometrics dataset. This function
    takes in a the neuromorph directory and finds all the necessary files
    and then writes the gm_labels.nii.gz and wm_labels.nii.gz as well as
    a csv that serves as a dictionary for the files needed for this project.
    """
    # Edit this code to change the training data set and to add more
    # types of input images

    import csv
    import glob
    import os
    if not file_list:
        globcard = os.path.join(in_dir, "*", "TissueClassify", "neuro_lbls.nii.gz")
        file_list = glob.glob(globcard)

    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)

    header = ["ID", "WMEdges", "WMLabelmap", "WMLabels",
              "GMEdges", "GMLabelmap", "GMLabels",
              "Modalities", "T1"]
    modalities = ["T1"]

    csv_filename = os.path.join(cache_dir, "project_files.csv")
    if not overwrite and os.path.isfile(csv_filename):
        return csv_filename
    csv_file = open(csv_filename, "wb")
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(header)

    for i, nm_file in enumerate(file_list):
        subject_id = get_nm_subject_id(nm_file)
        t1_file = get_nm_t1(nm_file, fs=True, cache_dir=cache_dir)

        # create subject_dir
        out_dir = os.path.join(cache_dir, subject_id)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        # define output files and check if they exist
        wm_label_file = os.path.join(out_dir, "wm_labels.nii.gz")
        wlf = os.path.isfile(wm_label_file)
        wm_truth_file = os.path.join(out_dir, "wm_edge.nii.gz")
        wtf = os.path.isfile(wm_truth_file)

        gm_label_file = os.path.join(out_dir, "gm_labels.nii.gz")
        glf = os.path.isfile(wm_label_file)
        gm_truth_file = os.path.join(out_dir, "gm_edge.nii.gz")
        gtf = os.path.isfile(gm_truth_file)

        # define label maps
        # these label maps contain the regions that will be used as the
        # separate training regions in the random forest.
        in_fs_file = get_nm_fs_scan(nm_file, name='aseg.mgz')
        out_fs_file = convert_fs_scan(in_fs_file, os.path.join(cache_dir, subject_id, "aseg.nii.gz"),
                                      resample_type='nearest')
        in_filled_file = get_nm_fs_scan(nm_file, name='filled.mgz')
        out_filled_file = convert_fs_scan(in_filled_file, os.path.join(cache_dir, subject_id, "filled.nii.gz"),
                                          resample_type='nearest')
        ws_wm, ws_gm, wm_labels, gm_labels = createwatersheds(out_fs_file, out_filled_file)

        if not (wlf and wtf and glf and gtf) or overwrite:
            # define truth images
            nm_converted_file = convert_fs_scan(nm_file, os.path.join(cache_dir, subject_id, "neuro_labels.nii.gz"),
                                                resample_type='nearest')
            wm_edge, gm_edge = getedges(nm_converted_file)

            # write wm files to subject dir
            write_image(ws_wm, wm_label_file, overwrite)
            write_image(wm_edge, wm_truth_file, overwrite)

            # write gm files to subject dir
            write_image(ws_gm, gm_label_file, overwrite)
            write_image(gm_edge, gm_truth_file, overwrite)

        row = [subject_id,
               wm_truth_file, wm_label_file, wm_labels,
               gm_truth_file, gm_label_file, gm_labels,
               modalities, t1_file]
        csv_writer.writerow(row)

    csv_file.close()
    return csv_filename


def save_data_frame(data, data_file):
    if data_file.endswith(".hdf5"):
        store = pd.HDFStore(data_file, "w", complib=str("zlib"), complevel=5)
        store.put("ImageData", data, data_columns=data.columns)
        store.close()
    elif data_file.endswith(".csv"):
        data.to_csv(data_file)
    else:
        print("Could not determine output type. Data was not written to output.")
