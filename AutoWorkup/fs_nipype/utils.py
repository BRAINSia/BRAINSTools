import nipype
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.io import DataGrabber
from nipype.interfaces.freesurfer import MRIConvert, MRIsConvert
import SimpleITK as sitk
import os

def center_volume(in_file):
    import SimpleITK as sitk
    import os
    img = sitk.ReadImage(in_file)
    size = img.GetSize()
    origin = img.GetOrigin()
    new_origin = [0,0,0]
    for i, xx in enumerate(origin):
        new_origin[i] = float(size[i])/2
        if xx < 0:
            new_origin[i] = -new_origin[i]
    img.SetOrigin(new_origin)
    out_file = os.path.abspath(os.path.basename(in_file))
    sitk.WriteImage(img, out_file)
    return out_file


def recodeLabelMap(in_file, out_file, recode_file):
    """This function has been adapted from BRAINSTools and serves
    as a means to recode a label map based upon an input csv
    file."""
    import SimpleITK as sitk
    import os
    import csv
    import sys

    # Convert csv to RECODE_TABLE
    CSVFile = open(recode_file, 'rb')
    reader = csv.reader(CSVFile)
    header = reader.next()
    n_cols = len(header)
    if n_cols == 4:
        # ignore label names
        label_keys = (0, 2)
    elif n_cols == 2:
        # no label names present
        label_keys = (0, 1)
    else:
        # csv does not match format requirements
        print("ERROR: input csv file for label recoding does meet requirements")
        sys.exit()

    # read csv file
    RECODE_TABLE = list()
    for line in reader:
        RECODE_TABLE.append((int(line[label_keys[0]]), int(line[label_keys[1]])))
        
    def minimizeSizeOfImage(outlabels):
        """This function will find the largest integer value in the labelmap, and
        cast the image to the smallest possible integer size so that no loss of data
        results."""
        measureFilt  = sitk.StatisticsImageFilter()
        measureFilt.Execute(outlabels)
        imgMin=measureFilt.GetMinimum()
        imgMax=measureFilt.GetMaximum()
        if imgMax < (2**8)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt8 )
        elif imgMax < (2**16)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt16 )
        elif imgMax < (2**32)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt32 )
        elif imgMax < (2**64)-1:
            outlabels = sitk.Cast( outlabels, sitk.sitkUInt64 )
        return outlabels
    
    LabelImage=sitk.Cast(sitk.ReadImage(in_file), sitk.sitkUInt32)
    for (old, new) in RECODE_TABLE:
        LabelImage = sitk.Cast((LabelImage == old), sitk.sitkUInt32)*(new - old)+LabelImage
    LabelImage = minimizeSizeOfImage(LabelImage)
    out_file = os.path.abspath(out_file)
    sitk.WriteImage(LabelImage, out_file)
    return out_file


def create_recoding_wf(in_file, out_file=None):
    wf = nipype.Workflow(name="RecodeLabels")

    inputspec = nipype.pipeline.Node(nipype.IdentityInterface(['labelmap',
                                                               'recode_file']), 
                                     name="Inputs")
    inputspec.inputs.recode_file = in_file

    convert_labelmap = nipype.pipeline.Node(MRIConvert(), name="ConvertLabelMap")
    convert_labelmap.inputs.in_type = 'mgz'
    convert_labelmap.inputs.out_type = 'nii'
    convert_labelmap.inputs.out_orientation = 'RAS'
    convert_labelmap.inputs.out_file = 'labelmap.nii'
    wf.connect([(inputspec, convert_labelmap, [('labelmap', 'in_file')])])

    recode = nipype.Node(nipype.Function(['in_file',
                                          'out_file',
                                          'recode_file'],
                                         ['out_file'],
                                         recodeLabelMap), 
                         name = "RecodeLabelMap")
    if out_file == None:
        recode.inputs.out_file = 'recodedlabelmap.nii'
    else:
        recode.inputs.out_file = out_file

    wf.connect([(convert_labelmap, recode, [('out_file', 'in_file')]),
                (inputspec, recode, [('recode_file', 'recode_file')])])

    center_labelmap = nipype.Node(nipype.Function(['in_file'], ['out_file'],
                                                  center_volume),
                                  name="CenterLabelMap")

    wf.connect([(recode, center_labelmap, [('out_file', 'in_file')])])

    outputspec = nipype.Node(nipype.IdentityInterface(['recodedlabelmap']), name="Outputs")

    wf.connect([(center_labelmap, outputspec, [('out_file', 'recodedlabelmap')])])    
    return wf

def createsrcsubj(source_directory):
    """
    Returns a node that acts as the datasource for a source subject such as 
    'fsaverage'
    """
    outfields = ['lh_BA1_exvivo',
                 'lh_BA2_exvivo',
                 'lh_BA3a_exvivo',
                 'lh_BA3b_exvivo',
                 'lh_BA4a_exvivo',
                 'lh_BA4p_exvivo',
                 'lh_BA6_exvivo',
                 'lh_BA44_exvivo',
                 'lh_BA45_exvivo',
                 'lh_V1_exvivo',
                 'lh_V2_exvivo',
                 'lh_MT_exvivo',
                 'lh_entorhinal_exvivo',
                 'lh_perirhinal_exvivo',
                 'lh_BA1_exvivo_thresh',
                 'lh_BA2_exvivo_thresh',
                 'lh_BA3a_exvivo_thresh',
                 'lh_BA3b_exvivo_thresh',
                 'lh_BA4a_exvivo_thresh',
                 'lh_BA4p_exvivo_thresh',
                 'lh_BA6_exvivo_thresh',
                 'lh_BA44_exvivo_thresh',
                 'lh_BA45_exvivo_thresh',
                 'lh_V1_exvivo_thresh',
                 'lh_V2_exvivo_thresh',
                 'lh_MT_exvivo_thresh',
                 'lh_entorhinal_exvivo_thresh',
                 'lh_perirhinal_exvivo_thresh',
                 'rh_BA1_exvivo',
                 'rh_BA2_exvivo',
                 'rh_BA3a_exvivo',
                 'rh_BA3b_exvivo',
                 'rh_BA4a_exvivo',
                 'rh_BA4p_exvivo',
                 'rh_BA6_exvivo',
                 'rh_BA44_exvivo',
                 'rh_BA45_exvivo',
                 'rh_V1_exvivo',
                 'rh_V2_exvivo',
                 'rh_MT_exvivo',
                 'rh_entorhinal_exvivo',
                 'rh_perirhinal_exvivo',
                 'rh_BA1_exvivo_thresh',
                 'rh_BA2_exvivo_thresh',
                 'rh_BA3a_exvivo_thresh',
                 'rh_BA3b_exvivo_thresh',
                 'rh_BA4a_exvivo_thresh',
                 'rh_BA4p_exvivo_thresh',
                 'rh_BA6_exvivo_thresh',
                 'rh_BA44_exvivo_thresh',
                 'rh_BA45_exvivo_thresh',
                 'rh_V1_exvivo_thresh',
                 'rh_V2_exvivo_thresh',
                 'rh_MT_exvivo_thresh',
                 'rh_entorhinal_exvivo_thresh',
                 'rh_perirhinal_exvivo_thresh']
    datasource = pe.Node(nio.DataGrabber(outfields=outfields), name="Source_Subject")
    datasource.inputs.base_directory = source_directory
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(
        lh_BA1_exvivo='label/lh.BA1_exvivo.label',
        lh_BA2_exvivo='label/lh.BA2_exvivo.label',
        lh_BA3a_exvivo='label/lh.BA3a_exvivo.label',
        lh_BA3b_exvivo='label/lh.BA3b_exvivo.label',
        lh_BA4a_exvivo='label/lh.BA4a_exvivo.label',
        lh_BA4p_exvivo='label/lh.BA4p_exvivo.label',
        lh_BA6_exvivo='label/lh.BA6_exvivo.label',
        lh_BA44_exvivo='label/lh.BA44_exvivo.label',
        lh_BA45_exvivo='label/lh.BA45_exvivo.label',
        lh_V1_exvivo='label/lh.V1_exvivo.label',
        lh_V2_exvivo='label/lh.V2_exvivo.label',
        lh_MT_exvivo='label/lh.MT_exvivo.label',
        lh_entorhinal_exvivo='label/lh.entorhinal_exvivo.label',
        lh_perirhinal_exvivo='label/lh.perirhinal_exvivo.label',
        lh_BA1_exvivo_thresh='label/lh.BA1_exvivo.thresh.label',
        lh_BA2_exvivo_thresh='label/lh.BA2_exvivo.thresh.label',
        lh_BA3a_exvivo_thresh='label/lh.BA3a_exvivo.thresh.label',
        lh_BA3b_exvivo_thresh='label/lh.BA3b_exvivo.thresh.label',
        lh_BA4a_exvivo_thresh='label/lh.BA4a_exvivo.thresh.label',
        lh_BA4p_exvivo_thresh='label/lh.BA4p_exvivo.thresh.label',
        lh_BA6_exvivo_thresh='label/lh.BA6_exvivo.thresh.label',
        lh_BA44_exvivo_thresh='label/lh.BA44_exvivo.thresh.label',
        lh_BA45_exvivo_thresh='label/lh.BA45_exvivo.thresh.label',
        lh_V1_exvivo_thresh='label/lh.V1_exvivo.thresh.label',
        lh_V2_exvivo_thresh='label/lh.V2_exvivo.thresh.label',
        lh_MT_exvivo_thresh='label/lh.MT_exvivo.thresh.label',
        lh_entorhinal_exvivo_thresh='label/lh.entorhinal_exvivo.thresh.label',
        lh_perirhinal_exvivo_thresh='label/lh.perirhinal_exvivo.thresh.label',
        rh_BA1_exvivo='label/rh.BA1_exvivo.label',
        rh_BA2_exvivo='label/rh.BA2_exvivo.label',
        rh_BA3a_exvivo='label/rh.BA3a_exvivo.label',
        rh_BA3b_exvivo='label/rh.BA3b_exvivo.label',
        rh_BA4a_exvivo='label/rh.BA4a_exvivo.label',
        rh_BA4p_exvivo='label/rh.BA4p_exvivo.label',
        rh_BA6_exvivo='label/rh.BA6_exvivo.label',
        rh_BA44_exvivo='label/rh.BA44_exvivo.label',
        rh_BA45_exvivo='label/rh.BA45_exvivo.label',
        rh_V1_exvivo='label/rh.V1_exvivo.label',
        rh_V2_exvivo='label/rh.V2_exvivo.label',
        rh_MT_exvivo='label/rh.MT_exvivo.label',
        rh_entorhinal_exvivo='label/rh.entorhinal_exvivo.label',
        rh_perirhinal_exvivo='label/rh.perirhinal_exvivo.label',
        rh_BA1_exvivo_thresh='label/rh.BA1_exvivo.thresh.label',
        rh_BA2_exvivo_thresh='label/rh.BA2_exvivo.thresh.label',
        rh_BA3a_exvivo_thresh='label/rh.BA3a_exvivo.thresh.label',
        rh_BA3b_exvivo_thresh='label/rh.BA3b_exvivo.thresh.label',
        rh_BA4a_exvivo_thresh='label/rh.BA4a_exvivo.thresh.label',
        rh_BA4p_exvivo_thresh='label/rh.BA4p_exvivo.thresh.label',
        rh_BA6_exvivo_thresh='label/rh.BA6_exvivo.thresh.label',
        rh_BA44_exvivo_thresh='label/rh.BA44_exvivo.thresh.label',
        rh_BA45_exvivo_thresh='label/rh.BA45_exvivo.thresh.label',
        rh_V1_exvivo_thresh='label/rh.V1_exvivo.thresh.label',
        rh_V2_exvivo_thresh='label/rh.V2_exvivo.thresh.label',
        rh_MT_exvivo_thresh='label/rh.MT_exvivo.thresh.label',
        rh_entorhinal_exvivo_thresh='label/rh.entorhinal_exvivo.thresh.label',
        rh_perirhinal_exvivo_thresh='label/rh.perirhinal_exvivo.thresh.label')
    return datasource, outfields

