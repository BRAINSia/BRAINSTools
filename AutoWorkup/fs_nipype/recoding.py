import nipype
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
