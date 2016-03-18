import nipype
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.io import DataGrabber
from nipype.interfaces.freesurfer import MRIConvert, MRIsConvert
import SimpleITK as sitk
import os
import sys


def getdefaultconfig():
    config = { 'custom_atlas' : None,
               'cw256' : False,
               'field_strength' : '1.5T',
               'fs_home' : checkenv(),
               'in_T1s' : list(),
               'in_T2' : None,
               'in_FLAIR' : None,
               'longitudinal' : False,
               'long_base' : None,
               'openmp' : None,
               'plugin' : 'Linear',
               'plugin_args' : None,
               'qcache' : False,
               'queue' : None,
               'recoding_file' : None,
               'src_subject_id' : 'fsaverage',
               'subject_id' : None,
               'subjects_dir' : None,
               'timepoints' : list() }
    config['source_subject'] = os.path.join(config['fs_home'], 'subjects',
                                            config['src_subject_id'])
    config['awk_file'] = os.path.join(config['fs_home'], 'bin',
                                      'extract_talairach_avi_QA.awk')
    config['registration_template'] = os.path.join(config['fs_home'], 'average',
                                                   'RB_all_withskull_2014-08-21.gca')
    for hemi in ('lh', 'rh'):
        config['{0}_atlas'.format(hemi)] = os.path.join(
            config['fs_home'], 'average',
            '{0}.average.curvature.filled.buckner40.tif'.format(hemi))
        config['{0}_classifier'.format(hemi)] = os.path.join(
            config['fs_home'], 'average',
            'rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs'.format(hemi))
        config['{0}_classifier2'.format(hemi)] = os.path.join(
            config['fs_home'], 'average',
            '{0}.destrieux.simple.2009-07-29.gcs'.format(hemi))
        config['{0}_classifier3'.format(hemi)] = os.path.join(
            config['fs_home'], 'average',
            '{0}.DKTatlas40.gcs'.format(hemi))
    config['LookUpTable'] = os.path.join(config['fs_home'], 'ASegStatsLUT.txt')
    config['WMLookUpTable'] = os.path.join(config['fs_home'], 'WMParcStatsLUT.txt')
    return config


def checkenv():
    """Check for the necessary FS environment variables"""
    fs_home = os.environ.get('FREESURFER_HOME')
    path = os.environ.get('PATH')
    print("FREESURFER_HOME: {0}".format(fs_home))
    if fs_home == None:
        print("ERROR: please set FREESURFER_HOME before running the workflow")
    elif not os.path.isdir(fs_home):
        print("ERROR: FREESURFER_HOME must be set to a valid directory before " + 
        "running this workflow")
    elif os.path.join(fs_home, 'bin') not in path.replace('//','/'):
        print(path)
        print("ERROR: Could not find necessary executable in path")
        setupscript = os.path.join(fs_home, 'SetUpFreeSurfer.sh')
        if os.path.isfile(setupscript):
            print("Please source the setup script before running the workflow:" +
            "\nsource {0}".format(setupscript))
        else:
            print("Please ensure that FREESURFER_HOME is set to a valid fs " +
            "directory and source the necessary SetUpFreeSurfer.sh script before running " +
            "this workflow")
    else:        
        return fs_home
    sys.exit(2)


def modify_qsub_args(queue, memoryGB, minThreads, maxThreads, stdout='/dev/null', stderr='/dev/null'):
    """
    Code from BRAINSTools:
    https://github.com/BRAINSia/BRAINSTools.git
    BRAINSTools/AutoWorkup/utilities/distributed.py

    Outputs qsub_args string for Nipype nodes
    queue is the string to specify the queue "-q all.q | -q HJ,ICTS,UI"
    memoryGB is a numeric in gigabytes to be given (ie 2.1 will result in "-l mem_free=2.1G")
          if memoryGB = 0, then it is automatically computed.
    minThreads The fewest number of threads to use (if an algorithm has benifits from more than 1 thread)
    maxThreads The max number of threads to use (if an algorithm is not multi-threaded, then just use 1)
    stdout Where to put stdout logs
    stderr Where to put stderr logs

    >>> modify_qsub_args('test', 2, 5, None)
    -S /bin/bash -cwd -pe smp 5 -l mem_free=2G -o /dev/null -e /dev/null test FAIL
    >>> modify_qsub_args('test', 2, 5, -1 )
    -S /bin/bash -cwd -pe smp 5- -l mem_free=2G -o /dev/null -e /dev/null test FAIL
    >>> modify_qsub_args('test', 8, 5, 7)
    -S /bin/bash -cwd -pe smp 5-7 -l mem_free=8G -o /dev/null -e /dev/null test FAIL
    >>> modify_qsub_args('test', 8, 5, 7, -1)
    -S /bin/bash -cwd -pe smp 5-7 -l mem_free=8G -o /dev/null -e /dev/null test FAIL
    >>> modify_qsub_args('test', 1, 5, 7, stdout='/my/path', stderr='/my/error')
    -S /bin/bash -cwd -pe smp 5-7 -l mem_free=1G -o /my/path -e /my/error test FAIL
    """
    import math
    assert memoryGB <= 48 , "Memory must be supplied in GB, so anything more than 24 seems not-useful now."

    ## NOTE: At least 1 thread needs to be requested per 2GB needed
    memoryThreads = int(math.ceil(memoryGB/float(2))) #Ensure that threads are integers
    minThreads = max(minThreads, memoryThreads)
    maxThreads = max(maxThreads, memoryThreads)
    maxThreads=int(maxThreads) # Ensure that threads are integers
    minThreads=int(minThreads) # Ensure that threads are integers

    if maxThreads is None or minThreads == maxThreads:
       threadsRangeString =  '{0}'.format(minThreads)
       maxThreads = minThreads
    elif maxThreads == -1:
       threadsRangeString= '{0}-'.format(minThreads)
       maxThreads = 12345 #HUGE NUMBER!
    else:
       threadsRangeString= "{0}-{1}".format(minThreads,maxThreads)

    if maxThreads < minThreads:
       assert  maxThreads > minThreads, "Must specify maxThreads({0}) > minThreads({1})".format(minThreads,maxThreads)
    format_str = '-q {queue} -S /bin/bash -cwd -pe smp {totalThreads} -o {stdout} -e {stderr}'.format(
                 mint=minThreads, maxt=threadsRangeString,
                 totalThreads=threadsRangeString,
                 mem=memoryGB,
                 stdout=stdout, stderr=stderr, queue=queue)
    return format_str

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

