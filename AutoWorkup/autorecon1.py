import sys
import os
import errno
import nipype
from nipype.interfaces.utility import Function,IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import *


def VerifyInputs(T1sList):
    ##TODO Make this its own node
    ##TODO Convert .mgz files to .nii.gz
    ##TODO Check the FOV
    """Verify size outside of pipeline processing"""
    print "Verifying input T1 size"
    try:
        import SimpleITK as sitk
        T1Length = len(T1sList)
        extension = None
        for i_extension in ['.mgz', '.nii', 'nii.gz']:
            if T1sList[0].endswith(i_extension):
                extension = i_extension
        if T1Length == 0:
            print("ERROR: No T1's Given")
            sys.exit(-1)
        elif T1Length > 1 and extension in ['.nii', '.nii.gz']:
            firstSize = sitk.ReadImage(T1sList[0]).GetSize()
            for otherFilename in T1sList[1:]:
                if firstSize != sitk.ReadImage(otherFilename).GetSize():
                    print("ERROR: T1s not the same size can not process {0} {1} together".format(
                        T1sList[0], otherFilename))
                    sys.exit(-1)
        elif extension == None:
            print "ERROR: Input files must be have '.mgz', '.nii', or '.nii.gz' extension"
            sys.exit(-1)
    except OSError as exc:  # Python >2.5
        print "ERROR: Could not verify input file sizes using SimpleITK"
        print exc
        sys.exit(-1)
    return T1sList

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def outputfilename(subjects_dir, subject_id, filename, subfolder1='mri', subfolder2=''):
    dest_dir = os.path.join(
        subjects_dir, subject_id, subfolder1, subfolder2)
    if not os.path.isdir(dest_dir):
        mkdir_p(dest_dir)
    return os.path.join(dest_dir, filename)

def awk(awk_file, log_file):
    """
    This method uses 'awk' which must be installed prior to running the workflow and is not a
    part of nipype or freesurfer.
    Future work may be done to create a method that achieves the same results using a python
    script.
    """
    import subprocess
    command = 'awk'
    subprocess.call([command, '-f', awk_file, log_file])
    return log_file

def copy_file(in_file, out_file=None):
    """
    Create a function to copy a file that can be modified by a following node without changing the original file
    """
    import os
    import shutil
    if out_file == None:
        out_file = os.path.join(os.getcwd(), os.path.basename(in_file))
    print "copying %s to %s" % (in_file, out_file)
    shutil.copy(in_file, out_file)
    return out_file

def copy_files(in_files, out_files):
    """
    Create a function to copy a file that can be modified by a following node 
    without changing the original file
    """
    import shutil
    if len(in_files) != len(out_files):
        print "ERROR: Length of input files must be identical to the length of \
        outrput files to be copied"
        sys.exit(-1)
    for i, in_file in enumerate(in_files):
        out_file = out_files[i]
        print "copying %s to %s" % (in_file, out_file)
        shutil.copy(in_file, out_file)
    return out_files

def create_preproc_filenames(subjects_dir, subject_id, in_T1s):
    # Create output filenames
    inputvols = list()
    iscaleout = list()
    ltaout = list()
    orig_dir = os.path.join(subjects_dir, subject_id, 'mri', 'orig')

    for i, T1 in enumerate(in_T1s):
        file_num = str(i + 1)
        while len(file_num) < 3:
            file_num = '0' + file_num
        iscaleout.append(os.path.join(orig_dir, file_num + '-iscale.txt'))
        ltaout.append(os.path.join(orig_dir, file_num + '.lta'))
        inputvols.append(os.path.join(orig_dir, file_num + '.mgz'))
    return inputvols, iscaleout, ltaout

def create_AutoRecon1(config):
    # AutoRecon1
    # Workflow
    ar1_wf = pe.Workflow(name='AutoRecon1')

    if not config['longitudinal']:
        ar1_inputs = pe.Node(interface=IdentityInterface(
            fields=['Raw_T1s', 'Raw_T2', 'Raw_FLAIR', 'subject_id', 'subjects_dir']),
                             run_without_submitting=True,
                             name='AutoRecon1_Inputs')
        ar1_inputs.inputs.subject_id = config['current_id']
        ar1_inputs.inputs.subjects_dir = config['subjects_dir']
        ar1_inputs.inputs.Raw_T1s = VerifyInputs(config['in_T1s'])

        inputvols, iscaleout, ltaout = create_preproc_filenames(config['subjects_dir'], config['current_id'], config['in_T1s'])

        # T1 image preparation
        # For all T1's mri_convert ${InputVol} ${out_file}
        T1_image_preparation = pe.MapNode(
            MRIConvert(), iterfield=['in_file', 'out_file'], name="T1_prep")
        T1_image_preparation.inputs.out_file = inputvols

        ar1_wf.connect([(ar1_inputs, T1_image_preparation, [('Raw_T1s', 'in_file')]),
                        ])

        # T2 image preparation
        if config['in_T2'] != None:
            # Create T2raw.mgz
            # mri_convert
            ar1_inputs.inputs.Raw_T2 = config['in_T2']
            T2_convert = pe.Node(MRIConvert(), name="T2_convert")
            T2_convert.inputs.out_file = os.path.join(config['subjects_dir'], config['current_id'], 'mri', 'orig', 'T2raw.mgz')
            T2_convert.inputs.no_scale = True
            ar1_wf.connect([(ar1_inputs, T2_convert, [('Raw_T2', 'in_file')]),
                            ]) 

        if config['in_FLAIR'] != None:
            # FLAIR image preparation
            # Create FLAIRraw.mgz
            # mri_convert
            ar1_inputs.inputs.Raw_FLAIR = config['in_FLAIR']
            FLAIR_convert = pe.Node(MRIConvert(), name="FLAIR_convert")
            FLAIR_convert.inputs.out_file = os.path.join(config['subjects_dir'],
                                                         config['current_id'],
                                                         'mri',
                                                         'orig',
                                                         'FLAIRraw.mgz')
            FLAIR_convert.inputs.no_scale = True
            ar1_wf.connect([(ar1_inputs, FLAIR_convert, [('Raw_FLAIR', 'in_file')]),
                            ])
    else:
        ar1_inputs = pe.Node(interface=IdentityInterface(
            fields=['in_T1s',
                    'iscales',
                    'ltas',
                    'subj_to_template_lta',
                    'template_talairach_xfm',
                    'template_brainmask',
                    'subject_id',
                    'subjects_dir']),
                             run_without_submitting=True,
                             name='AutoRecon1_Inputs')
        ar1_inputs.inputs.subject_id = config['current_id']
        ar1_inputs.inputs.subjects_dir = config['subjects_dir']

        _volnames_, in_iscales, in_ltas = create_preproc_filenames(config['subjects_dir'], config['current_id'], config['in_T1s'])

        copy_ltas = pe.MapNode(Function(['in_file', 'out_file'],
                                        ['out_file'],
                                        copy_file),
                               iterfield=['in_file', 'out_file'],
                               name='Copy_ltas')
        ar1_wf.connect([(ar1_inputs, copy_ltas, [('ltas', 'in_file')])])
        copy_ltas.inputs.out_file = in_ltas

        copy_iscales = pe.MapNode(Function(['in_file', 'out_file'],
                                           ['out_file'],
                                           copy_file),
                                  iterfield=['in_file', 'out_file'],
                                  name='Copy_iscales')
        ar1_wf.connect([(ar1_inputs, copy_iscales, [('iscales', 'in_file')])])
        copy_iscales.inputs.out_file = in_iscales

        concatenate_lta = pe.MapNode(ConcatenateLTA(), iterfield=['in_file'],
                                     name="Concatenate_ltas")
        ar1_wf.connect([(copy_ltas, concatenate_lta, [('out_file', 'in_file')]),
                        (ar1_inputs, concatenate_lta, [('subj_to_template_lta', 'subj_to_base')])])

    
    # Motion Correction
    """
    When there are multiple source volumes, this step will correct for small
    motions between them and then average them together.  The output of the
    motion corrected average is mri/rawavg.mgz which is then conformed to
    255 cubed char images (1mm isotropic voxels) in mri/orig.mgz.
    """

    create_template = pe.Node(RobustTemplate(), name="Robust_Template")
    create_template.inputs.average_metric = 'median'
    create_template.inputs.template_output = outputfilename(config['subjects_dir'], config['current_id'], 'rawavg.mgz')
    create_template.inputs.no_iteration = True
    if config['longitudinal']:
        ar1_wf.connect([(concatenate_lta, create_template, [('out_file', 'initial_transforms')]),
                        (ar1_inputs, create_template, [('in_T1s', 'in_files')]),
                        (copy_iscales, create_template, [('out_file','in_intensity_scales')]),
                        ])
        #TODO: connect mri_concatenate_lta to create_template.inputs.ixforms
        #TODO: connect the copied iscaleout files to the --iscalein input
    else:
        create_template.inputs.fixed_timepoint = True
        create_template.inputs.auto_detect_sensitivity = True
        create_template.inputs.initial_timepoint = 1
        create_template.inputs.scaled_intensity_outputs = iscaleout
        create_template.inputs.transform_outputs = ltaout
        create_template.inputs.subsample_threshold = 200
        create_template.inputs.intensity_scaling = True
        ar1_wf.connect([(T1_image_preparation, create_template, [('out_file', 'in_files')]),
                    ])

    # mri_convert
    conform_template = pe.Node(MRIConvert(), name='Conform_Template')
    conform_template.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], 'orig.mgz')
    if config['longitudinal']:
        conform_template.inputs.out_datatype = 'uchar'
    else:
        conform_template.inputs.conform = True
        conform_template.inputs.cw256 = config['cw256']    
        conform_template.inputs.resample_type = 'cubic'

    ar1_wf.connect(
        [(create_template, conform_template, [('template_output', 'in_file')])])

    add_to_header = pe.Node(AddXFormToHeader(), name="Add_Transform_to_Header")
    add_to_header.inputs.copy_name = True
    add_to_header.inputs.transform = os.path.join(
        config['subjects_dir'], config['current_id'], 'mri', 'transforms', 'talairach.xfm')

    ar1_wf.connect([(conform_template, add_to_header, [('out_file', 'in_file'),
                                                       ('out_file', 'out_file')])])

    # Talairach
    """
    This computes the affine transform from the orig volume to the MNI305 atlas using Avi Snyders 4dfp
    suite of image registration tools, through a FreeSurfer script called talairach_avi.
    Several of the downstream programs use talairach coordinates as seed points.
    """

    bias_correction = pe.Node(MNIBiasCorrection(), name="Bias_correction")
    bias_correction.inputs.iterations = 1
    bias_correction.inputs.protocol_iterations = 1000
    bias_correction.inputs.distance = 50
    bias_correction.inputs.no_rescale = True
    bias_correction.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], 'orig_nu.mgz')

    ar1_wf.connect([(add_to_header, bias_correction, [('out_file', 'in_file')]),
                ])

    if config['longitudinal']:
        copy_template_xfm = pe.Node(Function(['in_file', 'out_file'],
                                             ['out_file'],
                                             copy_file),
                                    name='Copy_Template_Transform')
        copy_template_xfm.inputs.out_file = os.path.join(
            config['subjects_dir'], config['current_id'], 'mri', 'transforms', 'talairach.auto.xfm')

        ar1_wf.connect([(ar1_inputs, copy_template_xfm, [('template_talairach_xfm', 'in_file')])])
    else:
        talairach_avi = pe.Node(TalairachAVI(), name="Compute_Transform")
        talairach_avi.inputs.atlas = '3T18yoSchwartzReactN32_as_orig'
        talairach_avi.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], 
                                                       'talairach.auto.xfm', 'mri', 'transforms')
        
        ar1_wf.connect([(bias_correction, talairach_avi, [('out_file', 'in_file')]),
                    ])
        
    copy_transform = pe.Node(Function(['in_file', 'out_file'],
                                      ['out_file'],
                                      copy_file),
                             name='Copy_Transform')
    copy_transform.inputs.out_file = os.path.join(
        config['subjects_dir'], config['current_id'], 'mri', 'transforms', 'talairach.xfm')

    if config['longitudinal']:
        ar1_wf.connect([(copy_template_xfm, copy_transform, [('out_file', 'in_file')])])
    else:
        ar1_wf.connect([(talairach_avi, copy_transform, [('out_file', 'in_file')])])

    check_alignment = pe.Node(
        CheckTalairachAlignment(), name="Check_Talairach_Alignment")
    check_alignment.inputs.threshold = 0.005
    ar1_wf.connect([(copy_transform, check_alignment, [('out_file', 'in_file')]),
                    ])

    if not config['longitudinal']:
        awk_logfile = pe.Node(Function(['awk_file', 'log_file'],
                                       ['log_file'],
                                       awk),
                              name='Awk')
        awk_logfile.inputs.awk_file = os.path.join(config['FREESURFER_HOME'],
                                                   'bin',
                                                   'extract_talairach_avi_QA.awk')
        ar1_wf.connect([(talairach_avi, awk_logfile, [('out_log', 'log_file')])])

        # TODO datasink the output from TalirachQC...not sure how to do this
        tal_qc = pe.Node(TalairachQC(), name="Detect_Aligment_Failures")
        ar1_wf.connect([(awk_logfile, tal_qc, [('log_file', 'log_file')]),
                    ])

    # Intensity Normalization
    # Performs intensity normalization of the orig volume and places the result in mri/T1.mgz.
    # Attempts to correct for fluctuations in intensity that would otherwise make intensity-based
    # segmentation much more difficult. Intensities for all voxels are scaled so that the mean
    # intensity of the white matter is 110.

    mri_normalize = pe.Node(Normalize(), name="Normalize_T1")
    mri_normalize.inputs.gradient = 1
    mri_normalize.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], 'T1.mgz')
    ar1_wf.connect([(bias_correction, mri_normalize, [('out_file', 'in_file')]),
                    (copy_transform, mri_normalize,
                     [('out_file', 'transform')]),
                    ])
    
    # Skull Strip
    """
    Removes the skull from mri/T1.mgz and stores the result in 
    mri/brainmask.auto.mgz and mri/brainmask.mgz. Runs the mri_watershed program.
    """

    if config['longitudinal']:
        copy_template_brainmask = pe.Node(Function(['in_file', 'out_file'],
                                                   ['out_file'],
                                                   copy_file),
                                          name='Copy_Template_Brainmask')
        copy_template_brainmask.inputs.out_file = os.path.join(
            config['subjects_dir'], config['current_id'], 'mri', 'brainmask_{0}.mgz'.format(config['long_template']))
        ar1_wf.connect([(ar1_inputs, copy_template_brainmask, [('template_brainmask', 'in_file')])])

        mask1 = pe.Node(ApplyMask(), name="ApplyMask1")
        mask1.inputs.keep_mask_deletion_edits = True
        mask1.inputs.out_file = os.path.join(config['subjects_dir'], config['current_id'], 'mri', 'brainmask.auto.mgz')
        
        ar1_wf.connect([(mri_normalize, mask1, [('out_file', 'in_file')]),
                        (copy_template_brainmask, mask1, [('out_file', 'mask_file')])])

        mask2 = pe.Node(ApplyMask(), name="ApplyMask2")
        mask2.inputs.keep_mask_deletion_edits = True
        mask2.inputs.transfer = 255

        ar1_wf.connect([(mask1, mask2, [('out_file', 'in_file'),
                                        ('out_file', 'out_file')]),
                        (copy_template_brainmask, mask2, [('out_file', 'mask_file')])])
                        
    else:    
        mri_em_register = pe.Node(EMRegister(), name="EM_Register")
        mri_em_register.inputs.template = os.path.join(config['FREESURFER_HOME'],
                                                       'average',
                                                       'RB_all_withskull_2014-08-21.gca')
        mri_em_register.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], 
                                                         'talairach_with_skull.lta', 'mri', 'transforms')
        mri_em_register.inputs.skull = True
        if config['openmp'] != None:
            mri_em_register.inputs.num_threads = config['openmp']
        if config['plugin_args'] != None:
            mri_em_register.plugin_args = config['plugin_args']
        ar1_wf.connect([(bias_correction, mri_em_register, [('out_file', 'in_file')]),
                    ])

        watershed_skull_strip = pe.Node(
            WatershedSkullStrip(), name='Watershed_Skull_Strip')
        watershed_skull_strip.inputs.t1 = True
        watershed_skull_strip.inputs.brain_atlas = os.path.join(config['FREESURFER_HOME'],
                                                                'average',
                                                                'RB_all_withskull_2014-08-21.gca')
        watershed_skull_strip.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], 
                                                               'brainmask.auto.mgz')
        ar1_wf.connect([(mri_normalize, watershed_skull_strip, [('out_file', 'in_file')]),
                        (mri_em_register, watershed_skull_strip,
                         [('out_file', 'transform')]),
                    ])

    copy_brainmask = pe.Node(Function(['in_file', 'out_file'],
                                      ['out_file'],
                                      copy_file),
                             name='Copy_Brainmask')
    copy_brainmask.inputs.out_file = os.path.join(config['subjects_dir'],
                                                  config['current_id'],
                                                  'mri',
                                                  'brainmask.mgz')

    if config['longitudinal']:
        ar1_wf.connect([(mask2, copy_brainmask, [('brain_vol', 'in_file')])])
    else:
        ar1_wf.connect([(watershed_skull_strip, copy_brainmask, [('brain_vol', 'in_file')])])
        
    return ar1_wf
