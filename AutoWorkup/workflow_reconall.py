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

def CreateStandardOutFileNames(input_T1s, subject_id, subjects_dir):
    import os
    InputVols = list()
    Iscaleout = list()
    LTAout = list()
    orig_dir = os.path.join(subjects_dir, subject_id, 'mri', 'orig')
    XFMout = os.path.join(
        subjects_dir, subject_id, 'mri', 'transforms', 'talairach.xfm')
    for i, T1 in enumerate(input_T1s):
        file_num = str(i + 1)
        while len(file_num) < 3:
            file_num = '0' + file_num
        Iscaleout.append(os.path.join(orig_dir, file_num + '-iscale.txt'))
        LTAout.append(os.path.join(orig_dir, file_num + '.lta'))
        InputVols.append(os.path.join(orig_dir, file_num + '.mgz'))
    return InputVols, Iscaleout, LTAout, XFMout, subjects_dir

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

def create_reconall(in_T1s, subject_id, in_T2, in_FLAIR, subjects_dir, qcache, cw256, fs_home):

    def outputfilename(filename, subfolder1='mri', subfolder2=''):
        dest_dir = os.path.join(
            subjects_dir, subject_id, subfolder1, subfolder2)
        if not os.path.isdir(dest_dir):
            mkdir_p(dest_dir)
        return os.path.join(dest_dir, filename)

    # Workflow Configurations
    reconall = pe.Workflow(name="recon-all")

    # AutoRecon1
    # Workflow
    ar1_wf = pe.Workflow(name='AutoRecon1')

    ar1_inputs = pe.Node(interface=IdentityInterface(
        fields=['Raw_T1', 'Raw_T2', 'Raw_FLAIR', 'subject_id', 'subjects_dir']),
        run_without_submitting=True,
        name='AutoRecon1_Inputs')

    ar1_inputs.inputs.Raw_T1 = VerifyInputs(in_T1s)
    ar1_inputs.inputs.subject_id = subject_id
    ar1_inputs.inputs.subjects_dir = subjects_dir

    # T1 image preparation
    # For all T1's mri_convert ${InputVol} ${out_file}
    T1_image_preparation = pe.MapNode(
        MRIConvert(), iterfield=['in_file', 'out_file'], name="T1_prep")

    # Create output filenames

    out_fn = pe.Node(Function(['input_T1s', 'subject_id', 'subjects_dir'],
                              ['InputVols', 'Iscaleout', 'LTAout',
                                  'XFMout', 'subjects_dir'],
                              CreateStandardOutFileNames),
                     name="CreateStandardOutFileNames")

    ar1_wf.connect([(ar1_inputs, T1_image_preparation, [('Raw_T1', 'in_file')]),
                    (ar1_inputs, out_fn, [('Raw_T1', 'input_T1s'),
                                          ('subject_id', 'subject_id'),
                                          ('subjects_dir', 'subjects_dir')]),
                    (out_fn, T1_image_preparation,
                     [('InputVols', 'out_file')]),
                    ])


    # Motion Correction
    """
    When there are multiple source volumes, this step will correct for small
    motions between them and then average them together.  The output of the
    motion corrected average is mri/rawavg.mgz which is then conformed to
    255 cubed char images (1mm isotropic voxels) in mri/orig.mgz.
    """

    create_long_template = pe.Node(RobustTemplate(), name="Robust_Template")
    create_long_template.inputs.average_metric = 'median'
    create_long_template.inputs.template_output = outputfilename('rawavg.mgz')
    create_long_template.inputs.auto_detect_sensitivity = True
    create_long_template.inputs.initial_timepoint = 1
    create_long_template.inputs.fixed_timepoint = True
    create_long_template.inputs.no_iteration = True
    create_long_template.inputs.intensity_scaling = True
    create_long_template.inputs.subsample_threshold = 200

    ar1_wf.connect([(T1_image_preparation, create_long_template, [('out_file', 'infiles')]),
                    (out_fn, create_long_template, [('Iscaleout', 'scaled_intensity_outputs'),
                                                    ('LTAout', 'transform_outputs')]),
                    ])

    # mri_convert
    conform_template = pe.Node(MRIConvert(), name='Conform_Template')
    conform_template.inputs.out_file = outputfilename('orig.mgz')
    conform_template.inputs.conform = True
    conform_template.inputs.cw256 = cw256    
    conform_template.inputs.resample_type = 'cubic'

    ar1_wf.connect(
        [(create_long_template, conform_template, [('template_output', 'in_file')])])

    add_to_header = pe.Node(AddXFormToHeader(), name="Add_Transform_to_Header")
    add_to_header.inputs.copy_name = True
    add_to_header.inputs.out_file = outputfilename('orig.mgz')

    ar1_wf.connect([(conform_template, add_to_header, [('out_file', 'in_file')]),
                    (out_fn, add_to_header, [('XFMout', 'transform')]),
                    ])

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
    bias_correction.inputs.out_file = outputfilename('orig_nu.mgz')

    ar1_wf.connect([(add_to_header, bias_correction, [('out_file', 'in_file')]),
                    ])

    talairach_avi = pe.Node(TalairachAVI(), name="Compute_Transform")
    talairach_avi.inputs.atlas = '3T18yoSchwartzReactN32_as_orig'
    talairach_avi.inputs.out_file = outputfilename(
        'talairach.auto.xfm', 'mri', 'transforms')

    ar1_wf.connect([(bias_correction, talairach_avi, [('out_file', 'in_file')]),
                    ])

    copy_transform = pe.Node(Function(['in_file', 'out_file'],
                                      ['out_file'],
                                      copy_file),
                             name='Copy_Transform')

    ar1_wf.connect([(talairach_avi, copy_transform, [('out_file', 'in_file')]),
                    (out_fn, copy_transform, [('XFMout', 'out_file')])
                    ])

    check_alignment = pe.Node(
        CheckTalairachAlignment(), name="Check_Talairach_Alignment")
    check_alignment.inputs.threshold = 0.005
    ar1_wf.connect([(copy_transform, check_alignment, [('out_file', 'in_file')]),
                    ])

    awk_logfile = pe.Node(Function(['awk_file', 'log_file'],
                                   ['log_file'],
                                   awk),
                          name='Awk')
    awk_logfile.inputs.awk_file = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                               'bin',
                                               'extract_talairach_avi_QA.awk')
    ar1_wf.connect([(talairach_avi, awk_logfile, [('out_log', 'log_file')]),
                    ])

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
    mri_normalize.inputs.out_file = outputfilename('T1.mgz')
    ar1_wf.connect([(bias_correction, mri_normalize, [('out_file', 'in_file')]),
                    (copy_transform, mri_normalize,
                     [('out_file', 'transform')]),
                    ])

    # Skull Strip
    """
    Removes the skull from mri/T1.mgz and stores the result in 
    mri/brainmask.auto.mgz and mri/brainmask.mgz. Runs the mri_watershed program.
    """

    mri_em_register = pe.Node(EMRegister(), name="EM_Register")
    mri_em_register.inputs.template = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                                   'average',
                                                   'RB_all_withskull_2014-08-21.gca')
    mri_em_register.inputs.out_file = outputfilename(
        'talairach_with_skull.lta', 'mri', 'transforms')
    mri_em_register.inputs.skull = True
    ar1_wf.connect([(bias_correction, mri_em_register, [('out_file', 'in_file')]),
                    ])

    watershed_skull_strip = pe.Node(
        WatershedSkullStrip(), name='Watershed_Skull_Strip')
    watershed_skull_strip.inputs.t1 = True
    watershed_skull_strip.inputs.brain_atlas = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                                            'average',
                                                            'RB_all_withskull_2014-08-21.gca')
    watershed_skull_strip.inputs.out_file = outputfilename(
        'brainmask.auto.mgz')
    ar1_wf.connect([(mri_normalize, watershed_skull_strip, [('out_file', 'in_file')]),
                    (mri_em_register, watershed_skull_strip,
                     [('out_file', 'transform')]),
                    ])

    copy_brainmask = pe.Node(Function(['in_file', 'out_file'],
                                      ['out_file'],
                                      copy_file),
                             name='Copy_Brainmask')
    copy_brainmask.inputs.out_file = outputfilename('brainmask.mgz')

    ar1_wf.connect([(watershed_skull_strip, copy_brainmask, [('brain_vol', 'in_file')]),
                    ])

    # AutoRecon2
    # Workflow
    ar2_wf = pe.Workflow(name="AutoRecon2")

    # Input node
    ar2_inputs = pe.Node(IdentityInterface(fields=['orig',
                                                   'brainmask',
                                                   'transform',
                                                   'subject_id',
                                                   'lh',
                                                   'rh',
                                                   'subjects_dir',
                                                   'freesurfer_home']),
                         run_without_submitting=True,
                         name='AutoRecon2_Inputs')
    ar2_inputs.inputs.lh = 'lh'
    ar2_inputs.inputs.rh = 'rh'
    ar2_inputs.inputs.subjects_dir = subjects_dir

    # NU Intensity Correction
    """
    Non-parametric Non-uniform intensity Normalization (N3), corrects for 
    intensity non-uniformity in MR data, making relatively few assumptions about
    the data. This runs the MINC tool 'nu_correct'.
    """
    intensity_correction = pe.Node(
        MNIBiasCorrection(), name="Intensity_Correction")
    intensity_correction.inputs.iterations = 1
    intensity_correction.inputs.protocol_iterations = 1000
    intensity_correction.inputs.out_file = outputfilename('nu.mgz')
    ar2_wf.connect([(ar2_inputs, intensity_correction, [('orig', 'in_file'),
                                                        ('brainmask', 'mask'),
                                                        ('transform',
                                                         'transform')
                                                        ])
                    ])

    add_to_header_nu = pe.Node(AddXFormToHeader(), name="Add_XForm_to_NU")
    add_to_header_nu.inputs.copy_name = True
    ar2_wf.connect([(intensity_correction, add_to_header_nu, [('out_file', 'in_file'),
                                                              ('out_file',
                                                               'out_file')
                                                              ]),
                    (ar2_inputs, add_to_header_nu,
                     [('transform', 'transform')])
                    ])

    # EM Registration
    """
    Computes the transform to align the mri/nu.mgz volume to the default GCA 
    atlas found in FREESURFER_HOME/average (see -gca flag for more info).
    """
    align_transform = pe.Node(EMRegister(), name="Align_Transform")
    align_transform.inputs.template = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                                   'average',
                                                   'RB_all_2014-08-21.gca')
    align_transform.inputs.out_file = outputfilename(
        'talairach.lta', 'mri', 'transforms')
    align_transform.inputs.nbrspacing = 3
    ar2_wf.connect([(ar2_inputs, align_transform, [('brainmask', 'mask')]),
                    (add_to_header_nu, align_transform,
                     [('out_file', 'in_file')])
                    ])

    # CA Normalize
    """
    Further normalization, based on GCA model. The normalization is based on an
    estimate of the most certain segmentation voxels, which it then uses to 
    estimate the bias field/scalings. Creates mri/norm.mgz.
    """
    ca_normalize = pe.Node(CANormalize(), name='CA_Normalize')
    ca_normalize.inputs.out_file = outputfilename('norm.mgz')
    ca_normalize.inputs.control_points = outputfilename('ctrl_pts.mgz')
    ca_normalize.inputs.atlas = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                             'average',
                                             'RB_all_2014-08-21.gca')
    ar2_wf.connect([(align_transform, ca_normalize, [('out_file', 'transform')]),
                    (ar2_inputs, ca_normalize, [('brainmask', 'mask')]),
                    (add_to_header_nu, ca_normalize, [('out_file', 'in_file')])])

    # CA Register
    # Computes a nonlinear transform to align with GCA atlas.
    ca_register = pe.Node(CARegister(), name='CA_Register')
    ca_register.inputs.align = 'after'
    ca_register.inputs.no_big_ventricles = True
    ca_register.inputs.template = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                               'average',
                                               'RB_all_2014-08-21.gca')
    ca_register.inputs.out_file = outputfilename(
        'talairach.m3z', 'mri', 'transforms')
    ar2_wf.connect([(ca_normalize, ca_register, [('out_file', 'in_file')]),
                    (ar2_inputs, ca_register, [('brainmask', 'mask')]),
                    (align_transform, ca_register, [('out_file', 'transform')])
                    ])

    # Remove Neck
    """
    The neck region is removed from the NU-corrected volume mri/nu.mgz. Makes use
    of transform computed from prior CA Register stage.
    """
    remove_neck = pe.Node(RemoveNeck(), name='Remove_Neck')
    remove_neck.inputs.radius = 25
    remove_neck.inputs.out_file = outputfilename('nu_noneck.mgz')
    remove_neck.inputs.template = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                               'average',
                                               'RB_all_2014-08-21.gca')
    ar2_wf.connect([(ca_register, remove_neck, [('out_file', 'transform')]),
                    (add_to_header_nu, remove_neck, [('out_file', 'in_file')])
                    ])

    # EM Registration, with Skull
    # Computes transform to align volume mri/nu_noneck.mgz with GCA volume
    # possessing the skull.
    em_reg_withskull = pe.Node(EMRegister(), name='EM_Register_withSkull')
    em_reg_withskull.inputs.skull = True
    em_reg_withskull.inputs.out_file = outputfilename(
        'talairach_with_skull_2.lta', 'mri', 'transforms')
    em_reg_withskull.inputs.template = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                                    'average',
                                                    'RB_all_withskull_2014-08-21.gca')

    ar2_wf.connect([(align_transform, em_reg_withskull, [('out_file', 'transform')]),
                    (remove_neck, em_reg_withskull, [('out_file', 'in_file')])
                    ])

    # CA Label
    # Labels subcortical structures, based in GCA model.
    ca_label = pe.Node(CALabel(), name='CA_Label')
    ca_label.inputs.relabel_unlikely = (9, .3)
    ca_label.inputs.prior = 0.5
    ca_label.inputs.align = True
    ca_label.inputs.out_file = outputfilename('aseg.auto_noCCseg.mgz')
    ca_label.inputs.template = os.path.join(os.path.abspath(os.environ.get('FREESURFER_HOME')),
                                            'average',
                                            'RB_all_2014-08-21.gca')
    ar2_wf.connect([(ca_normalize, ca_label, [('out_file', 'in_file')]),
                    (ca_register, ca_label, [('out_file', 'transform')])
                    ])

    # mri_cc - segments the corpus callosum into five separate labels in the
    # subcortical segmentation volume 'aseg.mgz'
    segment_cc = pe.Node(SegmentCC(), name="Segment_CorpusCallosum")
    segment_cc.inputs.out_rotation = outputfilename(
        'cc_up.lta', 'mri', 'transforms')
    segment_cc.inputs.out_file = outputfilename('aseg.auto.mgz')
    ar2_wf.connect([(ar2_inputs, segment_cc, [('subject_id', 'subject_id'),
                                              ('subjects_dir', 'subjects_dir')]),
                    (ca_label, segment_cc, [('out_file', 'in_file')]),
                    (ca_normalize, segment_cc, [('out_file', 'in_norm')]),
                    ])

    copy_cc = pe.Node(Function(['in_file', 'out_file'],
                               ['out_file'],
                               copy_file),
                      name='Copy_CCSegmentation')
    copy_cc.inputs.out_file = outputfilename('aseg.presurf.mgz')

    ar2_wf.connect([(segment_cc, copy_cc, [('out_file', 'in_file')])
                    ])

    # Normalization2
    """
    Performs a second (major) intensity correction using only the brain volume a
    s the input (so that it has to be done after the skull strip). Intensity 
    normalization works better when the skull has been removed. Creates a new 
    brain.mgz volume. The -autorecon2-cp stage begins here.
    """
    normalization2 = pe.Node(Normalize(), name="Normalization2")
    normalization2.inputs.out_file = outputfilename('brain.mgz')
    ar2_wf.connect([(copy_cc, normalization2, [('out_file', 'segmentation')]),
                    (ar2_inputs, normalization2, [('brainmask', 'mask')]),
                    (ca_normalize, normalization2, [('out_file', 'in_file')])
                    ])

    # Mask Brain Final Surface

    # Applies brainmask.mgz to brain.mgz to create brain.finalsurfs.mgz.
    mri_mask = pe.Node(ApplyMask(), name="Mask_Brain_Final_Surface")
    mri_mask.inputs.mask_thresh = 5
    mri_mask.inputs.out_file = outputfilename('brain.finalsurfs.mgz')

    ar2_wf.connect([(normalization2, mri_mask, [('out_file', 'in_file')]),
                    (ar2_inputs, mri_mask, [('brainmask', 'mask_file')])
                    ])

    # WM Segmentation
    """
    Attempts to separate white matter from everything else. The input is 
    mri/brain.mgz, and the output is mri/wm.mgz. Uses intensity, neighborhood,
    and smoothness constraints. This is the volume that is edited when manually
    fixing defects. Calls mri_segment, mri_edit_wm_with_aseg, and mri_pretess.
    """

    wm_seg = pe.Node(SegmentWM(), name="Segment_WM")
    wm_seg.inputs.out_file = outputfilename('wm.seg.mgz')
    ar2_wf.connect([(normalization2, wm_seg, [('out_file', 'in_file')])
                    ])

    edit_wm = pe.Node(EditWMwithAseg(), name='Edit_WhiteMatter')
    edit_wm.inputs.out_file = outputfilename('wm.asegedit.mgz')
    edit_wm.inputs.keep_in = True
    ar2_wf.connect([(wm_seg, edit_wm, [('out_file', 'in_file')]),
                    (copy_cc, edit_wm, [('out_file', 'seg_file')]),
                    (normalization2, edit_wm, [('out_file', 'brain_file')])
                    ])

    pretess = pe.Node(MRIPretess(), name="MRI_Pretess")
    pretess.inputs.out_file = outputfilename('wm.mgz')
    pretess.inputs.label = 'wm'
    ar2_wf.connect([(edit_wm, pretess, [('out_file', 'in_filled')]),
                    (ca_normalize, pretess, [('out_file', 'in_norm')])
                    ])

    # Cut/Fill
    """ This creates the subcortical mass from which the orig surface is created. 
    The mid brain is cut from the cerebrum, and the hemispheres are cut from each 
    other. The left hemisphere is binarized to 255. The right hemisphere is binarized 
    to 127. The input is mri/wm.mgz and the output is mri/filled.mgz. Calls mri_fill.
    """

    cut_and_fill = pe.Node(MRIFill(), name="Cut_and_Fill")
    cut_and_fill.inputs.log_file = outputfilename('ponscc.cut.log', 'scripts')
    cut_and_fill.inputs.out_file = outputfilename('filled.mgz')
    ar2_wf.connect([(pretess, cut_and_fill, [('out_file', 'in_file')]),
                    (align_transform, cut_and_fill,
                     [('out_file', 'transform')]),
                    (ca_label, cut_and_fill, [('out_file', 'segmentation')]),
                    ])

    # Split by Hemisphere
    # fuction to define the filenames that are unique to each hemisphere
    def hemisphere_names(hemisphere, subject_id, subjects_dir):
        import os
        if hemisphere == 'lh':
            label = 255
        else:
            label = 127
        dest_dir = os.path.join(subjects_dir, subject_id, 'surf')
        if not os.path.isdir(dest_dir):
            mkdir_p(dest_dir)
        orig = os.path.join(dest_dir, hemisphere + '.orig')
        orig_nofix = orig + '.nofix'
        smoothwm = os.path.join(dest_dir, hemisphere + '.smoothwm')
        smoothwm_nofix = smoothwm + '.nofix'
        inflated = os.path.join(dest_dir, hemisphere + '.inflated')
        inflated_nofix = inflated + '.nofix'
        qsphere_nofix = os.path.join(dest_dir, hemisphere + '.qsphere.nofix')
        sulc = os.path.join(dest_dir, hemisphere + '.sulc')
        stats = os.path.join(
            subjects_dir, subject_id, 'stats', hemisphere + '.curv.stats')

        return label, hemisphere, orig, smoothwm, inflated, orig_nofix, smoothwm_nofix, inflated_nofix, qsphere_nofix, sulc, stats, subjects_dir

    hemispheres = pe.Node(Function(['hemisphere', 'subject_id', 'subjects_dir'],
                                   ['label',
                                    'hemisphere',
                                    'orig',
                                    'smoothwm',
                                    'inflated',
                                    'orig_nofix',
                                    'smoothwm_nofix',
                                    'inflated_nofix',
                                    'qsphere_nofix',
                                    'sulc',
                                    'stats',
                                    'subjects_dir'],
                                   hemisphere_names),
                          name='Hemispheres')
    ar2_lh = pe.Workflow("AutoRecon2_Left")

    # Tessellation
    """
    This is the step where the orig surface (ie, surf/?h.orig.nofix) is created.
    The surface is created by covering the filled hemisphere with triangles. 
    Runs mri_pretess to create a connected WM volume (neighboring voxels must 
    have faces in common) and then mri_tessellate to create the surface. The 
    places where the points of the triangles meet are called vertices. Creates
    the file surf/?h.orig.nofix Note: the topology fixer will create the surface
    ?h.orig. Finally mris_extract_main_component will remove small surface 
    components, not connected to the main body.
    """

    pretess2 = pe.Node(MRIPretess(), name='Pretess_by_Hemisphere')
    pretess2.inputs.out_file = 'filled-pretess.mgz'
    ar2_lh.connect([(hemispheres, pretess2, [('label', 'label')]),
                    ])

    tesselate = pe.Node(MRITessellate(), name="Tesselation")
    ar2_lh.connect([(hemispheres, tesselate, [('orig_nofix', 'out_file'),
                                              ('label', 'label_value')
                                              ]),
                    (pretess2, tesselate, [('out_file', 'in_file')]),
                    ])

    extract_main_component = pe.Node(
        ExtractMainComponent(), name="Extract_Main_Component")
    ar2_lh.connect([(tesselate, extract_main_component, [('surface', 'in_file'),
                                                         ('surface',
                                                          'out_file')
                                                         ]),
                    ])

    copy_orig = pe.Node(Function(['in_file', 'out_file'],
                                 ['out_file'],
                                 copy_file),
                        name='Copy_Orig')
    ar2_lh.connect([(extract_main_component, copy_orig, [('out_file', 'in_file')]),
                    (hemispheres, copy_orig, [('orig', 'out_file')])
                    ])

    # Orig Surface Smoothing 1
    """
    After tesselation, the orig surface is very jagged because each triangle is
    on the edge of a voxel face and so are at right angles to each other. The 
    vertex positions are adjusted slightly here to reduce the angle. This is 
    only necessary for the inflation processes. Creates surf/?h.smoothwm(.nofix).
    Calls mris_smooth. Smooth1 is the step just after tessellation.
    """

    smooth1 = pe.Node(SmoothTessellation(), name="Smooth1")
    smooth1.inputs.disable_estimates = True
    smooth1.inputs.seed = 1234

    ar2_lh.connect([(hemispheres, smooth1, [('smoothwm_nofix', 'out_file')]),
                    (extract_main_component, smooth1,
                     [('out_file', 'in_file')])
                    ])

    # Inflation 1
    """
    Inflation of the surf/?h.smoothwm(.nofix) surface to create surf/?h.inflated.
    The inflation attempts to minimize metric distortion so that distances and
    areas are preserved (ie, the surface is not stretched). In this sense, it is
    like inflating a paper bag and not a balloon. Inflate1 is the step just after
    tessellation.
    """

    inflate1 = pe.Node(MRIsInflate(), name="inflate1")
    inflate1.inputs.no_save_sulc = True

    copy_inflate1 = pe.Node(Function(['in_file', 'out_file'],
                                     ['out_file'],
                                     copy_file),
                            name='Copy_Inflate1')

    ar2_lh.connect([(smooth1, inflate1, [('surface', 'in_file')]),
                    (hemispheres, inflate1, [('inflated_nofix', 'out_file')]),
                    (inflate1, copy_inflate1, [('out_file', 'in_file')]),
                    (hemispheres, copy_inflate1, [('inflated', 'out_file')])
                    ])

    # Sphere
    """
    This is the initial step of automatic topology fixing. It is a 
    quasi-homeomorphic spherical transformation of the inflated surface designed
    to localize topological defects for the subsequent automatic topology fixer. 
    Calls mris_sphere.
    """

    qsphere = pe.Node(Sphere(), name="Sphere")
    qsphere.inputs.seed = 1234
    qsphere.inputs.magic = True

    ar2_lh.connect([(inflate1, qsphere, [('out_file', 'in_file')]),
                    (hemispheres, qsphere, [('qsphere_nofix', 'out_file')]),
                    ])

    # Automatic Topology Fixer
    """
    Finds topological defects (ie, holes in a filled hemisphere) using 
    surf/?h.qsphere.nofix, and changes the orig surface (surf/?h.orig.nofix) to 
    remove the defects. Changes the number of vertices. All the defects will be
    removed, but the user should check the orig surface in the volume to make 
    sure that it looks appropriate.
    """

    # This mris_fix_topology does not take in the {lh,rh}.orig file, but instead takes in the
    # subject ID and hemisphere and tries to find it from the subjects
    # directory

    fix_topology = pe.Node(FixTopology(), name="Fix_Topology")
    fix_topology.inputs.mgz = True
    fix_topology.inputs.ga = True
    fix_topology.inputs.seed = 1234

    ar2_lh.connect([(copy_orig, fix_topology, [('out_file', 'in_orig')]),
                    (copy_inflate1, fix_topology,
                     [('out_file', 'in_inflated')]),
                    (hemispheres, fix_topology,
                     [('hemisphere', 'hemisphere')]),
                    (qsphere, fix_topology, [('out_file', 'sphere')]),
                    ])

    euler_number = pe.Node(EulerNumber(), name="Euler_Number")

    ar2_lh.connect([(fix_topology, euler_number, [('out_file', 'in_file')]),
                    ])

    remove_intersection = pe.Node(
        RemoveIntersection(), name="Remove_Intersection")

    ar2_lh.connect([(euler_number, remove_intersection, [('out_file', 'in_file'),
                                                         ('out_file',
                                                          'out_file')
                                                         ]),
                    ])

    # The inflated file is removed in AutoRecon2 after it is used for the Fix
    # Topology step
    def rmfile(in_file, dependent):
        import os
        os.remove(in_file)
        out_file = in_file
        return out_file

    remove_inflate1 = pe.Node(Function(['in_file', 'dependent'],
                                       ['out_file'],
                                       rmfile),
                              name="Remove_Inflate1")

    ar2_lh.connect([(copy_inflate1, remove_inflate1, [('out_file', 'in_file')]),
                    (remove_intersection, remove_inflate1,
                     [('out_file', 'dependent')])
                    ])

    # White

    # This function implicitly calls other inputs based on the subject_id
    # need to make sure files are data sinked to the correct folders before
    # calling
    make_surfaces = pe.Node(MakeSurfaces(), name="Make_Surfaces")
    make_surfaces.inputs.noaparc = True
    make_surfaces.inputs.mgz = True
    make_surfaces.inputs.white_only = True

    ar2_lh.connect([(remove_intersection, make_surfaces, [('out_file', 'in_orig')]),
                    (hemispheres, make_surfaces,
                     [('hemisphere', 'hemisphere')])
                    ])

    # Orig Surface Smoothing 2

    # After tesselation, the orig surface is very jagged because each triangle is on
    # the edge of a voxel face and so are at right angles to each other. The vertex
    # positions are adjusted slightly here to reduce the angle. This is only necessary
    # for the inflation processes. Smooth2 is the step just after topology
    # fixing.

    smooth2 = pe.Node(SmoothTessellation(), name="Smooth2")
    smooth2.inputs.disable_estimates = True
    smooth2.inputs.smoothing_iterations = 3
    smooth2.inputs.seed = 1234

    ar2_lh.connect([(hemispheres, smooth2, [('smoothwm', 'out_file')]),
                    (make_surfaces, smooth2, [('out_white', 'in_file')])
                    ])

    # Inflation 2

    # Inflation of the surf/?h.smoothwm(.nofix) surface to create surf/?h.inflated.
    # The inflation attempts to minimize metric distortion so that distances and areas
    # are preserved (ie, the surface is not stretched). In this sense, it is like
    # inflating a paper bag and not a balloon. Inflate2 is the step just after
    # topology fixing.
    inflate2 = pe.Node(MRIsInflate(), name="inflate2")
    ar2_lh.connect([(smooth2, inflate2, [('surface', 'in_file')]),
                    (remove_inflate1, inflate2, [('out_file', 'out_file')]),
                    (hemispheres, inflate2, [('sulc', 'out_sulc')]),
                    ])

    # Compute Curvature
    # ?

    curvature1 = pe.Node(Curvature(), name="Curvature1")
    curvature1.inputs.save = True
    ar2_lh.connect([(make_surfaces, curvature1, [('out_white', 'in_file')]),
                    ])

    curvature2 = pe.Node(Curvature(), name="Curvature2")
    curvature2.inputs.threshold = .999
    curvature2.inputs.n = True
    curvature2.inputs.averages = 5
    curvature2.inputs.save = True
    curvature2.inputs.distances = (10, 10)

    ar2_lh.connect([(inflate2, curvature2, [('out_file', 'in_file')]),
                    ])

    curvature_stats = pe.Node(CurvatureStats(), name="Curvature_Stats")
    curvature_stats.inputs.min_max = True
    curvature_stats.inputs.write = True
    curvature_stats.inputs.values = True

    ar2_lh.connect([(smooth2, curvature_stats, [('surface', 'surface')]),
                    (make_surfaces, curvature_stats,
                     [('out_curv', 'in_curv')]),
                    (inflate2, curvature_stats, [('out_sulc', 'in_sulc')]),
                    (hemispheres, curvature_stats, [('hemisphere', 'hemisphere'),
                                                    ('stats', 'out_file')]),
                    ])

    ar2_rh = ar2_lh.clone(name="AutoRecon2_Right")

    # Connect inputs for the hemisphere workflows
    ar2_wf.connect([(ar2_inputs, ar2_lh, [('lh', 'Hemispheres.hemisphere'),
                                          ('subject_id',
                                           'Hemispheres.subject_id'),
                                          ('subjects_dir', 'Hemispheres.subjects_dir')]),
                    (ar2_inputs, ar2_rh, [('rh', 'Hemispheres.hemisphere'),
                                          ('subject_id',
                                           'Hemispheres.subject_id'),
                                          ('subjects_dir', 'Hemispheres.subjects_dir')]),
                    (ca_normalize, ar2_lh, [
                     ('out_file', 'Pretess_by_Hemisphere.in_norm')]),
                    (cut_and_fill, ar2_lh, [
                     ('out_file', 'Pretess_by_Hemisphere.in_filled')]),
                    (ar2_inputs, ar2_lh, [('subject_id', 'Fix_Topology.subject_id'),
                                          ('subjects_dir', 'Fix_Topology.subjects_dir')]),
                    (ar2_inputs, ar2_lh, [('subject_id', 'Make_Surfaces.subject_id'),
                                          ('subjects_dir', 'Make_Surfaces.subjects_dir')]),
                    (ar2_inputs, ar2_lh, [('subject_id', 'Curvature_Stats.subject_id'),
                                          ('subjects_dir', 'Curvature_Stats.subjects_dir')]),
                    (copy_cc, ar2_lh, [('out_file', 'Make_Surfaces.in_aseg')]),
                    (mri_mask, ar2_lh, [('out_file', 'Make_Surfaces.in_T1')]),
                    (cut_and_fill, ar2_lh, [
                     ('out_file', 'Make_Surfaces.in_filled')]),
                    (pretess, ar2_lh, [('out_file', 'Make_Surfaces.in_wm')]),

                    (ca_normalize, ar2_rh, [
                     ('out_file', 'Pretess_by_Hemisphere.in_norm')]),
                    (cut_and_fill, ar2_rh, [
                     ('out_file', 'Pretess_by_Hemisphere.in_filled')]),
                    (ar2_inputs, ar2_rh, [('subject_id', 'Fix_Topology.subject_id'),
                                          ('subjects_dir', 'Fix_Topology.subjects_dir')]),
                    (ar2_inputs, ar2_rh, [('subject_id', 'Make_Surfaces.subject_id'),
                                          ('subjects_dir', 'Make_Surfaces.subjects_dir')]),
                    (ar2_inputs, ar2_rh, [('subject_id', 'Curvature_Stats.subject_id'),
                                          ('subjects_dir', 'Curvature_Stats.subjects_dir')]),
                    (copy_cc, ar2_rh, [('out_file', 'Make_Surfaces.in_aseg')]),
                    (mri_mask, ar2_rh, [('out_file', 'Make_Surfaces.in_T1')]),
                    (cut_and_fill, ar2_rh, [
                     ('out_file', 'Make_Surfaces.in_filled')]),
                    (pretess, ar2_rh, [('out_file', 'Make_Surfaces.in_wm')]),
                    ])

    # AutoRecon3
    # Workflow
    ar3_wf = pe.Workflow(name="AutoRecon3")
    ar3_lh_wf = pe.Workflow(name="AutoRecon3_Left")

    # Input Node
    ar3_inputs = pe.Node(IdentityInterface(fields=['lh',
                                                   'rh',
                                                   'subject_id',
                                                   'subjects_dir',
                                                   'lh_inflated',
                                                   'rh_inflated',
                                                   'lh_smoothwm',
                                                   'rh_smoothwm',
                                                   'lh_white',
                                                   'rh_white',
                                                   'lh_white_H',
                                                   'rh_white_H',
                                                   'lh_white_K',
                                                   'rh_white_K',
                                                   'lh_cortex_label',
                                                   'rh_cortex_label',
                                                   'lh_orig',
                                                   'rh_orig',
                                                   'lh_sulc',
                                                   'rh_sulc',
                                                   'lh_area',
                                                   'rh_area',
                                                   'lh_curv',
                                                   'rh_curv',
                                                   'lh_atlas',
                                                   'rh_atlas',
                                                   'lh_classifier',
                                                   'rh_classifier',
                                                   'lh_orig_nofix',
                                                   'rh_orig_nofix',
                                                   'aseg_presurf',
                                                   'brain_finalsurfs',
                                                   'wm',
                                                   'filled',
                                                   'T2raw',
                                                   'brainmask',
                                                   'transform',
                                                   'orig_mgz',
                                                   'rawavg',
                                                   'norm']),
                         name='AutoRecon3_Inputs')

    ar3_inputs.inputs.lh = 'lh'
    ar3_inputs.inputs.rh = 'rh'

    ar3_inputs.inputs.lh_atlas = os.path.join(
        fs_home, 'average/lh.average.curvature.filled.buckner40.tif')
    ar3_inputs.inputs.lh_classifier = os.path.join(
        fs_home, 'average/lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs')
    ar3_inputs.inputs.rh_atlas = os.path.join(
        fs_home, 'average/rh.average.curvature.filled.buckner40.tif')
    ar3_inputs.inputs.rh_classifier = os.path.join(
        fs_home, 'average/rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs')

    ar3_lh_inputs = pe.Node(IdentityInterface(fields=['hemisphere',
                                                      'subject_id',
                                                      'subjects_dir',
                                                      'inflated',
                                                      'smoothwm',
                                                      'white',
                                                      'cortex_label',
                                                      'orig',
                                                      'aseg_presurf',
                                                      'brain_finalsurfs',
                                                      'wm',
                                                      'filled',
                                                      'T2raw',
                                                      'sphere',
                                                      'sulc',
                                                      'area',
                                                      'curv',
                                                      'classifier',
                                                      'atlas']),
                            name="Inputs")

    # Spherical Inflation

    # Inflates the orig surface into a sphere while minimizing metric distortion.
    # This step is necessary in order to register the surface to the spherical
    # atlas (also known as the spherical morph). Calls mris_sphere. Creates
    # surf/?h.sphere. The -autorecon3 stage begins here.

    ar3_sphere = pe.Node(Sphere(), name="Spherical_Inflation")
    ar3_sphere.inputs.seed = 1234

    ar3_lh_wf.connect([(ar3_lh_inputs, ar3_sphere, [('inflated', 'in_file'),
                                                    ('smoothwm',
                                                     'in_smoothwm'),
                                                    ('sphere', 'out_file'),
                                                    ])])

    # Ipsilateral Surface Registation (Spherical Morph)

    # Registers the orig surface to the spherical atlas through surf/?h.sphere.
    # The surfaces are first coarsely registered by aligning the large scale
    # folding patterns found in ?h.sulc and then fine tuned using the small-scale
    # patterns as in ?h.curv. Calls mris_register. Creates surf/?h.sphere.reg.

    ar3_surfreg = pe.Node(Register(), name="Surface_Registration")

    ar3_lh_wf.connect([(ar3_sphere, ar3_surfreg, [('out_file', 'in_surf')]),
                       (ar3_lh_inputs, ar3_surfreg, [('smoothwm', 'in_smoothwm'),
                                                     ('smoothwm', 'curv'),
                                                     ('sulc', 'in_sulc'),
                                                     ('atlas', 'target')])
                       ])

    # Jacobian

    # Computes how much the white surface was distorted in order to register to
    # the spherical atlas during the -surfreg step.

    ar3_jacobian = pe.Node(Jacobian(), name="Jacobian")

    ar3_lh_wf.connect([(ar3_lh_inputs, ar3_jacobian, [('white', 'in_origsurf')]),
                       (ar3_surfreg, ar3_jacobian,
                        [('out_file', 'in_mappedsurf')])
                       ])

    # Average Curvature

    # Resamples the average curvature from the atlas to that of the subject.
    # Allows the user to display activity on the surface of an individual
    # with the folding pattern (ie, anatomy) of a group.

    ar3_paint = pe.Node(Paint(), name="Average_Curvature")
    ar3_paint.inputs.averages = 5
    ar3_paint.inputs.template_param = 6

    ar3_lh_wf.connect([(ar3_surfreg, ar3_paint, [('out_file', 'in_surf')]),
                       (ar3_lh_inputs, ar3_paint, [('atlas', 'template')])])

    # Cortical Parcellation

    # Assigns a neuroanatomical label to each location on the cortical
    # surface. Incorporates both geometric information derived from the
    # cortical model (sulcus and curvature), and neuroanatomical convention.

    ar3_parcellation = pe.Node(MRIsCALabel(), "Cortical_Parcellation")
    ar3_parcellation.inputs.seed = 1234

    ar3_lh_wf.connect([(ar3_lh_inputs, ar3_parcellation, [('smoothwm', 'smoothwm'),
                                                          ('cortex_label',
                                                           'label'),
                                                          ('aseg_presurf',
                                                           'aseg'),
                                                          ('hemisphere',
                                                           'hemisphere'),
                                                          ('subject_id',
                                                           'subject_id'),
                                                          ('subjects_dir',
                                                           'subjects_dir'),
                                                          ('classifier', 'classifier')]),
                       (ar3_surfreg, ar3_parcellation,
                        [('out_file', 'canonsurf')])
                       ])

    # Pial Surface

    ar3_pial = pe.Node(MakeSurfaces(), name="Make_Pial_Surface")
    ar3_pial.inputs.fix_mtl = True
    ar3_pial.inputs.no_white = True
    ar3_pial.inputs.mgz = True

    ar3_lh_wf.connect([(ar3_lh_inputs, ar3_pial, [('wm', 'in_wm'),
                                                  ('orig', 'in_orig'),
                                                  ('filled', 'in_filled'),
                                                  ('white', 'orig_pial'),
                                                  ('white', 'orig_white'),
                                                  ('brain_finalsurfs',
                                                   'in_T1'),
                                                  ('aseg_presurf', 'in_aseg'),
                                                  ('hemisphere', 'hemisphere'),
                                                  ('subject_id', 'subject_id'),
                                                  ('subjects_dir', 'subjects_dir')]),
                       (ar3_parcellation, ar3_pial, [('out_file', 'in_label')])
                       ])

    # Surface Volume
    """
    Creates the ?h.volume file by first creating the ?h.mid.area file by
    adding ?h.area(.white) to ?h.area.pial, then dividing by two. Then ?h.volume
    is created by multiplying ?.mid.area with ?h.thickness.
    """

    ar3_add = pe.Node(MRIsCalc(), name="Add_Pial_Area")
    ar3_add.inputs.action = "add"
    ar3_add.inputs.out_file = 'area.mid'
    ar3_lh_wf.connect([(ar3_pial, ar3_add, [('out_area', 'in_file2')]),
                       (ar3_lh_inputs, ar3_add, [('area', 'in_file1')]),
                       ])

    ar3_divide = pe.Node(MRIsCalc(), name="Mid_Pial")
    ar3_divide.inputs.action = "div"
    ar3_divide.inputs.in_int = 2
    ar3_divide.inputs.out_file = 'area.mid'
    ar3_lh_wf.connect([(ar3_add, ar3_divide, [('out_file', 'in_file1')]),
                       ])

    ar3_volume = pe.Node(MRIsCalc(), name="Calculate_Volume")
    ar3_volume.inputs.action = "mul"
    ar3_volume.inputs.out_file = 'volume'
    ar3_lh_wf.connect([(ar3_divide, ar3_volume, [('out_file', 'in_file1')]),
                       (ar3_pial, ar3_volume, [('out_thickness', 'in_file2')]),
                       ])

    # Workflow1 Outputs
    outputs1 = pe.Node(IdentityInterface(fields=['pial']),
                       name="Outputs")
    ar3_lh_wf.connect([(ar3_pial, outputs1, [('out_pial', 'pial')]),

                       ])

    ar3_rh_wf = ar3_lh_wf.clone(name="AutoRecon3_Right")

    # Cortical Ribbon Mask
    """
    Creates binary volume masks of the cortical ribbon
    ie, each voxel is either a 1 or 0 depending upon whether it falls in the ribbon or not.
    """
    volume_mask = pe.Node(VolumeMask(), name="Mask_Ribbon")
    volume_mask.inputs.left_whitelabel = 2
    volume_mask.inputs.left_ribbonlabel = 3
    volume_mask.inputs.right_whitelabel = 41
    volume_mask.inputs.right_ribbonlabel = 42
    volume_mask.inputs.save_ribbon = True

    ar3_wf.connect([(ar3_inputs, volume_mask, [('aseg_presurf', 'in_aseg'),
                                               ('subject_id', 'subject_id'),
                                               ('subjects_dir',
                                                'subjects_dir'),
                                               ('lh_white', 'lh_white'),
                                               ('rh_white', 'rh_white'),
                                               ]),
                    (ar3_lh_wf, volume_mask, [
                     ('Make_Pial_Surface.out_pial', 'lh_pial')]),
                    (ar3_rh_wf, volume_mask, [
                     ('Make_Pial_Surface.out_pial', 'rh_pial')]),
                    ])


    for hemisphere, workflow in [('lh', ar3_lh_wf),  ('rh', ar3_rh_wf)]:
        if hemisphere == 'lh':
            opp_hemi = 'rh'
            opp_wf = ar3_rh_wf
        else:
            opp_hemi = 'lh'
            opp_wf = ar3_lh_wf
        # Connect Workflows
        ar3_wf.connect([(ar3_inputs, workflow, [('{0}_inflated'.format(hemisphere), 'Inputs.inflated'),
                                                ('{0}_smoothwm'.format(hemisphere),
                                                 'Inputs.smoothwm'),
                                                ('{0}'.format(hemisphere), 'Inputs.hemisphere'),
                                                ('{0}_white'.format(hemisphere), 'Inputs.white'),
                                                ('{0}_cortex_label'.format(hemisphere),
                                                 'Inputs.cortex_label'),
                                                ('{0}_orig'.format(hemisphere), 'Inputs.orig'),
                                                ('{0}_sulc'.format(hemisphere), 'Inputs.sulc'),
                                                ('{0}_area'.format(hemisphere), 'Inputs.area'),
                                                ('{0}_curv'.format(hemisphere), 'Inputs.curv'),
                                                ('subject_id',
                                                 'Inputs.subject_id'),
                                                ('subjects_dir',
                                                 'Inputs.subjects_dir'),
                                                ('aseg_presurf',
                                                 'Inputs.aseg_presurf'),
                                                ('brain_finalsurfs',
                                                 'Inputs.brain_finalsurfs'),
                                                ('wm', 'Inputs.wm'),
                                                ('filled', 'Inputs.filled'),
                                                ('{0}_atlas'.format(hemisphere), 'Inputs.atlas'),
                                                ('{0}_classifier'.format(hemisphere),
                                                 'Inputs.classifier')
                                                ])
                        ])


        # Parcellation Statistics
        """
        Runs mris_anatomical_stats to create a summary table of cortical parcellation statistics for each structure, including
        structure name
        number of vertices
        total surface area (mm^2)
        total gray matter volume (mm^3)
        average cortical thickness (mm)
        standard error of cortical thicknessr (mm)
        integrated rectified mean curvature
        integrated rectified Gaussian curvature
        folding index
        intrinsic curvature index.
        """
        parcellation_stats_white = pe.Node(
            ParcellationStats(), name="Parcellation_Stats_{0}_White".format(hemisphere) )
        parcellation_stats_white.inputs.mgz = True
        parcellation_stats_white.inputs.tabular_output = True
        parcellation_stats_white.inputs.surface = 'white'
        parcellation_stats_white.inputs.out_color = outputfilename('aparc.annot.ctab', 'label')
        parcellation_stats_white.inputs.out_table = outputfilename('{0}.aparc.stats'.format(hemisphere), 'stats')

        ar3_wf.connect([(ar3_inputs, parcellation_stats_white, [('subject_id', 'subject_id'),
                                                                ('subjects_dir',
                                                                 'subjects_dir'),
                                                                ('wm', 'wm'),
                                                                ('lh_white',
                                                                 'lh_white'),
                                                                ('rh_white',
                                                                 'rh_white'),
                                                                ('transform',
                                                                 'transform'),
                                                                ('brainmask',
                                                                 'brainmask'),
                                                                ('aseg_presurf',
                                                                 'aseg'),
                                                                ('{0}_cortex_label'.format(hemisphere),
                                                                 'in_cortex'),
                                                                ('{0}'.format(hemisphere),
                                                                 'hemisphere'),
                                                                ]),
                        (workflow, parcellation_stats_white, [('Make_Pial_Surface.out_pial', '{0}_pial'.format(hemisphere)),
                                                              ('Make_Pial_Surface.out_thickness',
                                                               'thickness'),
                                                              ('Cortical_Parcellation.out_file',
                                                               'in_annotation')
                                                              ]),
                        (opp_wf, parcellation_stats_white, [('Make_Pial_Surface.out_pial', '{0}_pial'.format(opp_hemi)),
                                                           ]),
                        (volume_mask, parcellation_stats_white, [('out_ribbon', 'ribbon')
                                                             ]),
                    ])
        
        parcellation_stats_pial = pe.Node(
            ParcellationStats(), name="Parcellation_Stats_{0}_Pial".format(hemisphere) )
        parcellation_stats_pial.inputs.mgz = True
        parcellation_stats_pial.inputs.tabular_output = True
        parcellation_stats_pial.inputs.surface = 'pial'
        parcellation_stats_pial.inputs.out_color = outputfilename('aparc.annot.ctab', 'label')
        parcellation_stats_pial.inputs.out_table = outputfilename('{0}.aparc.pial.stats'.format(hemisphere), 'stats')

        ar3_wf.connect([(ar3_inputs, parcellation_stats_pial, [('subject_id', 'subject_id'),
                                                             ('subjects_dir',
                                                              'subjects_dir'),
                                                             ('wm', 'wm'),
                                                             ('{0}_white'.format(hemisphere),
                                                              '{0}_white'.format(hemisphere)),
                                                             ('{0}_white'.format(opp_hemi),
                                                              '{0}_white'.format(opp_hemi)),
                                                             ('transform',
                                                              'transform'),
                                                             ('brainmask',
                                                              'brainmask'),
                                                             ('aseg_presurf',
                                                              'aseg'),
                                                             ('{0}_cortex_label'.format(hemisphere),
                                                              'in_cortex'),
                                                             ('{0}'.format(hemisphere),
                                                              'hemisphere'),
                                                             ]),
                    (workflow, parcellation_stats_pial, [('Make_Pial_Surface.out_pial', '{0}_pial'.format(hemisphere)),
                                                            ('Make_Pial_Surface.out_thickness',
                                                             'thickness'),
                                                            ('Cortical_Parcellation.out_file',
                                                             'in_annotation')
                                                            ]),
                    (opp_wf, parcellation_stats_pial, [('Make_Pial_Surface.out_pial', '{0}_pial'.format(opp_hemi)),
                                                            ]),
                    (volume_mask, parcellation_stats_pial, [('out_ribbon', 'ribbon')
                                                              ]),
                    ])

        # Cortical Parcellation 2
        cortical_parcellation_2 = pe.Node(MRIsCALabel(), name="Cortical_Parcellation_{0}_2".format(hemisphere))
        if hemisphere == 'lh':
            lh_cort_parc2 = cortical_parcellation_2
        else:
            rh_cort_parc2 = cortical_parcellation_2
            
        cortical_parcellation_2.inputs.classifier = os.path.join(
            fs_home, 'average', '{0}.destrieux.simple.2009-07-29.gcs'.format(hemisphere))
        cortical_parcellation_2.inputs.out_file = outputfilename('{0}.aparc.a2009s.annot'.format(hemisphere), 'label')
        cortical_parcellation_2.inputs.seed = 1234

        ar3_wf.connect([(ar3_inputs, cortical_parcellation_2, [('subject_id', 'subject_id'),
                                                     ('subjects_dir',
                                                      'subjects_dir'),
                                                     ('{0}'.format(hemisphere), 'hemisphere'),
                                                     ('{0}_smoothwm'.format(hemisphere), 'smoothwm'),
                                                     ('aseg_presurf', 'aseg'),
                                                     ('{0}_cortex_label'.format(hemisphere), 'label')]),
                        (workflow, cortical_parcellation_2, [
                            ('Surface_Registration.out_file', 'canonsurf')])
                        ])

        # Parcellation Statistics 2
        parcellation_stats_white_2 = parcellation_stats_white.clone(
            name="Parcellation_Statistics_{0}_2".format(hemisphere))
        parcellation_stats_white_2.inputs.out_color = outputfilename('aparc.annot.a2009s.ctab', 'label')
        parcellation_stats_white_2.inputs.out_table = outputfilename('{0}.aparc.a2009s.stats'.format(hemisphere), 'stats')
        ar3_wf.connect([(ar3_inputs, parcellation_stats_white_2, [('subject_id', 'subject_id'),
                                                                  ('subjects_dir',
                                                                   'subjects_dir'),
                                                                  ('wm', 'wm'),
                                                                  ('lh_white',
                                                                   'lh_white'),
                                                                  ('rh_white',
                                                                   'rh_white'),
                                                                  ('transform',
                                                                   'transform'),
                                                                  ('brainmask',
                                                                   'brainmask'),
                                                                  ('aseg_presurf',
                                                                   'aseg'),
                                                                  ('{0}_cortex_label'.format(hemisphere),
                                                                    'in_cortex'),
                                                                  ('{0}'.format(hemisphere),
                                                                    'hemisphere'),
                                                              ]),
                        (workflow, parcellation_stats_white_2, [('Make_Pial_Surface.out_pial', '{0}_pial'.format(hemisphere)),
                                                                ('Make_Pial_Surface.out_thickness',
                                                                 'thickness'),
                                                            ]),
                        (opp_wf, parcellation_stats_white_2, [('Make_Pial_Surface.out_pial', '{0}_pial'.format(opp_hemi)),
                                                          ]),
                        (volume_mask, parcellation_stats_white_2, [('out_ribbon', 'ribbon')
                                                               ]),
                        (cortical_parcellation_2, parcellation_stats_white_2,
                         [('out_file', 'in_annotation')])
                    ])

        # Cortical Parcellation 3
        cortical_parcellation_3 = pe.Node(MRIsCALabel(), name="Cortical_Parcellation_{0}_3".format(hemisphere))
        cortical_parcellation_3.inputs.classifier = os.path.join(
            fs_home, 'average', '{0}.DKTatlas40.gcs'.format(hemisphere))
        cortical_parcellation_3.inputs.out_file = outputfilename('{0}.aparc.DKTatlas40.annot'.format(hemisphere), 'label')
        cortical_parcellation_3.inputs.seed = 1234
        ar3_wf.connect([(ar3_inputs, cortical_parcellation_3, [('subject_id', 'subject_id'),
                                                               ('subjects_dir',
                                                                'subjects_dir'),
                                                               ('{0}'.format(hemisphere), 'hemisphere'),
                                                               ('{0}_smoothwm'.format(hemisphere), 'smoothwm'),
                                                               ('aseg_presurf', 'aseg'),
                                                               ('{0}_cortex_label'.format(hemisphere), 'label')]),
                        (workflow, cortical_parcellation_3, [
                            ('Surface_Registration.out_file', 'canonsurf')])
                    ])

        # Parcellation Statistics 3
        parcellation_stats_white_3 = parcellation_stats_white.clone(
            name="Parcellation_Statistics_{0}_3".format(hemisphere))
        parcellation_stats_white_3.inputs.out_color = outputfilename('aparc.annot.DKTatlas40.ctab', 'label')
        parcellation_stats_white_3.inputs.out_table = outputfilename('{0}.aparc.DKTatlas40.stats'.format(hemisphere), 'stats')

        ar3_wf.connect([(ar3_inputs, parcellation_stats_white_3, [('subject_id', 'subject_id'),
                                                                   ('subjects_dir',
                                                                    'subjects_dir'),
                                                                   ('wm', 'wm'),
                                                                   ('lh_white',
                                                                    'lh_white'),
                                                                   ('rh_white',
                                                                    'rh_white'),
                                                                   ('transform',
                                                                    'transform'),
                                                                   ('brainmask',
                                                                    'brainmask'),
                                                                   ('aseg_presurf',
                                                                    'aseg'),
                                                                   ('{0}_cortex_label'.format(hemisphere),
                                                                    'in_cortex'),
                                                                   ('{0}'.format(hemisphere),
                                                                    'hemisphere'),
                                                               ]),
                        (workflow, parcellation_stats_white_3, [('Make_Pial_Surface.out_pial', '{0}_pial'.format(hemisphere)),
                                                                  ('Make_Pial_Surface.out_thickness',
                                                                   'thickness'),
                                                              ]),
                        (opp_wf, parcellation_stats_white_3, [('Make_Pial_Surface.out_pial', '{0}_pial'.format(opp_hemi)),
                                                              ]),
                        (volume_mask, parcellation_stats_white_3, [('out_ribbon', 'ribbon')
                                                                ]),
                        (cortical_parcellation_3, parcellation_stats_white_3,
                         [('out_file', 'in_annotation')])
                    ])

        # WM/GM Contrast
        contrast = pe.Node(Contrast(), name="WM_GM_Contrast_{0}".format(hemisphere))
        if hemisphere == 'lh':
            lh_contrast = contrast
        else:
            rh_contrast = contrast

        ar3_wf.connect([(ar3_inputs, contrast, [('orig_mgz', 'orig'),
                                               ('rawavg', 'rawavg'),
                                               ('subject_id', 'subject_id'),
                                               ('subjects_dir',
                                                'subjects_dir'),
                                               ('{0}_white'.format(hemisphere), 'white'),
                                               ('{0}_cortex_label'.format(hemisphere), 'cortex'),
                                               ('{0}'.format(hemisphere), 'hemisphere')]),
                    (workflow, contrast, [('Make_Pial_Surface.out_thickness', 'thickness'),
                                              ('Cortical_Parcellation.out_file', 'annotation')])
                    ])
        #End for


    # Relabel Hypointensities
    relabel_hypos = pe.Node(
        RelabelHypointensities(), name="Relabel_Hypointensities")
    ar3_wf.connect([(ar3_inputs, relabel_hypos, [('aseg_presurf', 'aseg'),
                                                 ('lh_white', 'lh_white'),
                                                 ('rh_white', 'rh_white'),
                                                 ])])

    # APARC to ASEG
    # Adds information from the ribbon into the aseg.mgz (volume parcellation).
    aparc_2_aseg = pe.Node(Aparc2Aseg(), name="Aparc2Aseg")
    aparc_2_aseg.inputs.volmask = True
    ar3_wf.connect([(ar3_inputs, aparc_2_aseg, [('lh_white', 'lh_white'),
                                                ('rh_white', 'rh_white'),
                                                ('subject_id', 'subject_id'),
                                                ('subjects_dir',
                                                 'subjects_dir'),
                                                ]),
                    (ar3_lh_wf, aparc_2_aseg, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                               ('Cortical_Parcellation.out_file',
                                                'lh_annotation'),
                                               ]),
                    (ar3_rh_wf, aparc_2_aseg, [('Make_Pial_Surface.out_pial', 'rh_pial'),
                                               ('Cortical_Parcellation.out_file',
                                                'rh_annotation'),
                                               ]),
                    (volume_mask, aparc_2_aseg, [('rh_ribbon', 'rh_ribbon'),
                                                 ('lh_ribbon', 'lh_ribbon'),
                                                 ('out_ribbon', 'ribbon'),
                                                 ]),
                    (relabel_hypos, aparc_2_aseg, [('out_file', 'aseg')])
                    ])

    aparc_2_aseg_2009 = pe.Node(Aparc2Aseg(), name="Aparc2Aseg_2009")
    aparc_2_aseg_2009.inputs.volmask = True
    ar3_wf.connect([(ar3_inputs, aparc_2_aseg_2009, [('lh_white', 'lh_white'),
                                                     ('rh_white', 'rh_white'),
                                                     ('subject_id',
                                                      'subject_id'),
                                                     ('subjects_dir',
                                                      'subjects_dir'),
                                                     ]),
                    (ar3_lh_wf, aparc_2_aseg_2009, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                                    ]),
                    (lh_cort_parc2, aparc_2_aseg_2009,
                     [('out_file', 'lh_annotation')]),
                    (rh_cort_parc2, aparc_2_aseg_2009,
                     [('out_file', 'rh_annotation')]),
                    (ar3_rh_wf, aparc_2_aseg_2009, [('Make_Pial_Surface.out_pial', 'rh_pial'),
                                                    ]),
                    (volume_mask, aparc_2_aseg_2009, [('rh_ribbon', 'rh_ribbon'),
                                                      ('lh_ribbon',
                                                       'lh_ribbon'),
                                                      ('out_ribbon', 'ribbon'),
                                                      ]),
                    (relabel_hypos, aparc_2_aseg_2009, [('out_file', 'aseg')])
                    ])

    apas_2_aseg = pe.Node(Apas2Aseg(), name="Apas_2_Aseg")
    ar3_wf.connect([(aparc_2_aseg, apas_2_aseg, [('out_file', 'in_file')])])

    # Segmentation Stats
    """
    Computes statistics on the segmented subcortical structures found in 
    mri/aseg.mgz. Writes output to file stats/aseg.stats.
    """

    segstats = pe.Node(SegStatsReconAll(), name="Segmentation_Statistics")
    segstats.inputs.color_table_file = os.path.join(
        fs_home, 'ASegStatsLUT.txt')
    segstats.inputs.empty = True
    segstats.inputs.brain_vol = 'brain-vol-from-seg'
    segstats.inputs.exclude_ctx_gm_wm = True
    segstats.inputs.supratent = True
    segstats.inputs.subcort_gm = True
    segstats.inputs.etiv = True
    segstats.inputs.wm_vol_from_surf = True
    segstats.inputs.cortex_vol_from_surf = True
    segstats.inputs.total_gray = True
    segstats.inputs.euler = True
    segstats.inputs.exclude_id = 0
    segstats.inputs.intensity_units = "MR"

    ar3_wf.connect([(apas_2_aseg, segstats, [('out_file', 'segmentation_file')]),
                    (ar3_inputs, segstats, [('subject_id', 'subject_id'),
                                            ('subjects_dir', 'subjects_dir'),
                                            ('lh_white', 'lh_white'),
                                            ('rh_white', 'rh_white'),
                                            ('aseg_presurf', 'presurf_seg'),
                                            ('transform', 'transform'),
                                            ('norm', 'in_intensity'),
                                            ('norm', 'partial_volume_file'),
                                            ('brainmask', 'brainmask_file'),
                                            ('lh_orig_nofix', 'lh_orig_nofix'),
                                            ('rh_orig_nofix', 'rh_orig_nofix'),
                                            ]),
                    (volume_mask, segstats, [('out_ribbon', 'ribbon')]),
                    (ar3_lh_wf, segstats, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                           ]),
                    (ar3_rh_wf, segstats, [('Make_Pial_Surface.out_pial', 'rh_pial'),
                                           ]),
                    ])

    # White Matter Parcellation

    # Adds WM Parcellation info into the aseg and computes stat.

    wm_parcellation = pe.Node(Aparc2Aseg(), name="WM_Parcellation")
    wm_parcellation.inputs.volmask = True
    wm_parcellation.inputs.label_wm = True
    wm_parcellation.inputs.hypo_wm = True
    wm_parcellation.inputs.rip_unknown = True

    ar3_wf.connect([(ar3_inputs, wm_parcellation, [('lh_white', 'lh_white'),
                                                   ('rh_white', 'rh_white'),
                                                   ('subject_id',
                                                    'subject_id'),
                                                   ('subjects_dir',
                                                    'subjects_dir'),
                                                   ]),
                    (ar3_lh_wf, wm_parcellation, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                                  ('Cortical_Parcellation.out_file',
                                                   'lh_annotation'),
                                                  ]),
                    (ar3_rh_wf, wm_parcellation, [('Make_Pial_Surface.out_pial', 'rh_pial'),
                                                  ('Cortical_Parcellation.out_file',
                                                   'rh_annotation'),
                                                  ]),
                    (volume_mask, wm_parcellation, [('rh_ribbon', 'rh_ribbon'),
                                                    ('lh_ribbon', 'lh_ribbon'),
                                                    ('out_ribbon', 'ribbon'),
                                                    ]),
                    (apas_2_aseg, wm_parcellation, [('out_file', 'aseg')]),
                    (aparc_2_aseg, wm_parcellation, [('out_file', 'ctxseg')])
                    ])

    # White Matter Segmentation Stats

    wm_segstats = pe.Node(
        SegStatsReconAll(), name="WM_Segmentation_Statistics")
    wm_segstats.inputs.color_table_file = os.path.join(
        fs_home, 'WMParcStatsLUT.txt')
    wm_segstats.inputs.intensity_units = "MR"
    wm_segstats.inputs.wm_vol_from_surf = True
    wm_segstats.inputs.etiv = True
    wm_segstats.inputs.exclude_id = 0

    ar3_wf.connect([(wm_parcellation, wm_segstats, [('out_file', 'segmentation_file')]),
                    (ar3_inputs, wm_segstats, [('subject_id', 'subject_id'),
                                               ('subjects_dir',
                                                'subjects_dir'),
                                               ('lh_white', 'lh_white'),
                                               ('rh_white', 'rh_white'),
                                               ('aseg_presurf', 'presurf_seg'),
                                               ('transform', 'transform'),
                                               ('norm', 'in_intensity'),
                                               ('norm', 'partial_volume_file'),
                                               ('brainmask', 'brainmask_file'),
                                               ('lh_orig_nofix',
                                                'lh_orig_nofix'),
                                               ('rh_orig_nofix',
                                                'rh_orig_nofix'),
                                               ]),
                    (volume_mask, wm_segstats, [('out_ribbon', 'ribbon')]),
                    (ar3_lh_wf, wm_segstats, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                              ]),
                    (ar3_rh_wf, wm_segstats, [('Make_Pial_Surface.out_pial', 'rh_pial'),
                                              ]),
                    ])

    # Brodmann Area Maps (BA Maps) and Hinds V1 Atlas

    ba_inputs = pe.Node(IdentityInterface(fields=['lh_sphere_reg',
                                                  'rh_sphere_reg',
                                                  'lh_white',
                                                  'rh_white',
                                                  'lh_pial',
                                                  'rh_pial',
                                                  'transform',
                                                  'lh_thickness',
                                                  'rh_thickness',
                                                  'brainmask',
                                                  'aseg',
                                                  'ribbon',
                                                  'wm',
                                                  'subject_id',
                                                  'subjects_dir'
                                                  ]),
                        name="BA_Maps_Inputs")

    ba_WF = pe.Workflow(name="Brodmann_Area_Maps")
    labels = ["BA1", "BA2", "BA3a", "BA3b", "BA4a", "BA4p", "BA6",
              "BA44", "BA45", "V1", "V2", "MT", "entorhinal", "perirhinal"]
    for hemisphere in ['lh', 'rh']:
        for threshold in [True, False]:
            if threshold:
                node_name = 'BA_Maps_' + hemisphere + '_Tresh'
            else:
                node_name = 'BA_Maps_' + hemisphere

            node = pe.MapNode(
                Label2Label(), name=node_name, iterfield=['label'])
            node.inputs.hemisphere = hemisphere
            node.inputs.label = labels
            node.inputs.threshold = threshold

            label2annot = pe.Node(Label2Annot(), name=node_name + '_2_Annot')
            label2annot.inputs.hemisphere = hemisphere
            label2annot.inputs.color_table = os.path.join(
                fs_home, 'average', 'colortable_BA.txt')
            label2annot.inputs.verbose_off = True
            label2annot.inputs.keep_max = True
            if threshold:
                label2annot.inputs.out_annot = "BA_exvivo.thresh"
            else:
                label2annot.inputs.out_annot = "BA_exvivo"

            stats_node = pe.Node(
                ParcellationStats(), name=node_name + '_Stats')
            stats_node.inputs.hemisphere = hemisphere
            stats_node.inputs.mgz = True
            stats_node.inputs.surface = 'white'
            stats_node.inputs.tabular_output = True

            ba_WF.connect([(ba_inputs, node, [('{0}_sphere_reg'.format(hemisphere), 'sphere_reg'),
                                              ('{0}_white'.format(
                                                  hemisphere), 'white'),
                                              ('subject_id', 'subject_id'),
                                              ('subjects_dir', 'subjects_dir'),
                                              ]),
                           (node, label2annot, [('out_file', 'in_labels')]),
                           (ba_inputs, label2annot, [('subject_id', 'subject_id'),
                                                     ('subjects_dir',
                                                      'subjects_dir'),
                                                     ]),
                           (label2annot, stats_node,
                            [('out_file', 'in_annotation')]),
                           (ba_inputs, stats_node, [('{0}_thickness'.format(hemisphere), 'thickness'),
                                                    ('subject_id',
                                                     'subject_id'),
                                                    ('subjects_dir',
                                                     'subjects_dir'),
                                                    ('lh_white', 'lh_white'),
                                                    ('rh_white', 'rh_white'),
                                                    ('lh_pial', 'lh_pial'),
                                                    ('rh_pial', 'rh_pial'),
                                                    ('transform', 'transform'),
                                                    ('brainmask', 'brainmask'),
                                                    ('aseg', 'aseg'),
                                                    ('wm', 'wm'),
                                                    ('ribbon', 'ribbon')])])

    ar3_wf.connect([(ar3_lh_wf, ba_WF, [('Surface_Registration.out_file', 'BA_Maps_Inputs.lh_sphere_reg'),
                                        ('Make_Pial_Surface.out_thickness',
                                         'BA_Maps_Inputs.lh_thickness'),
                                        ('Make_Pial_Surface.out_pial',
                                         'BA_Maps_Inputs.lh_pial'),
                                        ]),
                    (ar3_rh_wf, ba_WF, [('Surface_Registration.out_file', 'BA_Maps_Inputs.rh_sphere_reg'),
                                        ('Make_Pial_Surface.out_thickness',
                                         'BA_Maps_Inputs.rh_thickness'),
                                        ('Make_Pial_Surface.out_pial',
                                         'BA_Maps_Inputs.rh_pial'),
                                        ]),
                    (ar3_inputs, ba_WF, [('subject_id', 'BA_Maps_Inputs.subject_id'),
                                         ('subjects_dir',
                                          'BA_Maps_Inputs.subjects_dir'),
                                         ('lh_white',
                                          'BA_Maps_Inputs.lh_white'),
                                         ('rh_white',
                                          'BA_Maps_Inputs.rh_white'),
                                         ('transform',
                                          'BA_Maps_Inputs.transform'),
                                         ('aseg_presurf',
                                          'BA_Maps_Inputs.aseg'),
                                         ('brainmask',
                                          'BA_Maps_Inputs.brainmask'),
                                         ('wm', 'BA_Maps_Inputs.wm')]),
                    (volume_mask, ba_WF, [
                        ('out_ribbon', 'BA_Maps_Inputs.ribbon')])
                    ])

    reconall.connect([(ar1_wf, ar3_wf, [('AutoRecon1_Inputs.subject_id', 'AutoRecon3_Inputs.subject_id'),
                                        ('AutoRecon1_Inputs.subjects_dir',
                                         'AutoRecon3_Inputs.subjects_dir'),
                                        ('Copy_Brainmask.out_file',
                                         'AutoRecon3_Inputs.brainmask'),
                                        ('Copy_Transform.out_file',
                                         'AutoRecon3_Inputs.transform'),
                                        ('Add_Transform_to_Header.out_file',
                                         'AutoRecon3_Inputs.orig_mgz'),
                                        ('Robust_Template.template_output',
                                         'AutoRecon3_Inputs.rawavg'),
                                        ]),
                      (ar1_wf, ar2_wf, [('Copy_Brainmask.out_file', 'AutoRecon2_Inputs.brainmask'),
                                        ('Copy_Transform.out_file',
                                         'AutoRecon2_Inputs.transform'),
                                        ('Add_Transform_to_Header.out_file',
                                         'AutoRecon2_Inputs.orig'),
                                        ('AutoRecon1_Inputs.subject_id',
                                         'AutoRecon2_Inputs.subject_id'),
                                        ('AutoRecon1_Inputs.subjects_dir',
                                         'AutoRecon2_Inputs.subjects_dir'),
                                        ]),
                      (ar2_lh, ar3_wf, [('inflate2.out_file', 'AutoRecon3_Inputs.lh_inflated'),
                                        ('Smooth2.surface',
                                         'AutoRecon3_Inputs.lh_smoothwm'),
                                        ('Make_Surfaces.out_white',
                                         'AutoRecon3_Inputs.lh_white'),
                                        ('Make_Surfaces.out_cortex',
                                         'AutoRecon3_Inputs.lh_cortex_label'),
                                        ('Make_Surfaces.out_area',
                                         'AutoRecon3_Inputs.lh_area'),
                                        ('Make_Surfaces.out_curv',
                                         'AutoRecon3_Inputs.lh_curv'),
                                        ('inflate2.out_sulc',
                                         'AutoRecon3_Inputs.lh_sulc'),
                                        ('Extract_Main_Component.out_file',
                                         'AutoRecon3_Inputs.lh_orig_nofix'),
                                        ('Remove_Intersection.out_file',
                                         'AutoRecon3_Inputs.lh_orig'),
                                        ('Curvature1.out_mean',
                                         'AutoRecon3_Inputs.lh_white_H'),
                                        ('Curvature1.out_gauss',
                                         'AutoRecon3_Inputs.lh_white_K'),
                                        ]),
                      (ar2_rh, ar3_wf, [('inflate2.out_file', 'AutoRecon3_Inputs.rh_inflated'),
                                        ('Smooth2.surface',
                                         'AutoRecon3_Inputs.rh_smoothwm'),
                                        ('Make_Surfaces.out_white',
                                         'AutoRecon3_Inputs.rh_white'),
                                        ('Make_Surfaces.out_cortex',
                                         'AutoRecon3_Inputs.rh_cortex_label'),
                                        ('Make_Surfaces.out_area',
                                         'AutoRecon3_Inputs.rh_area'),
                                        ('Make_Surfaces.out_curv',
                                         'AutoRecon3_Inputs.rh_curv'),
                                        ('inflate2.out_sulc',
                                         'AutoRecon3_Inputs.rh_sulc'),
                                        ('Extract_Main_Component.out_file',
                                         'AutoRecon3_Inputs.rh_orig_nofix'),
                                        ('Remove_Intersection.out_file',
                                         'AutoRecon3_Inputs.rh_orig'),
                                        ('Curvature1.out_mean',
                                         'AutoRecon3_Inputs.rh_white_H'),
                                        ('Curvature1.out_gauss',
                                         'AutoRecon3_Inputs.rh_white_K'),
                                        ]),
                      (ar2_wf, ar3_wf, [('Copy_CCSegmentation.out_file', 'AutoRecon3_Inputs.aseg_presurf'),
                                        ('Mask_Brain_Final_Surface.out_file',
                                         'AutoRecon3_Inputs.brain_finalsurfs'),
                                        ('MRI_Pretess.out_file',
                                         'AutoRecon3_Inputs.wm'),
                                        ('Cut_and_Fill.out_file',
                                         'AutoRecon3_Inputs.filled'),
                                        ('CA_Normalize.out_file',
                                         'AutoRecon3_Inputs.norm'),
                                        ]),
                      ])

    if in_T2 != None:
        # T2 image preparation
        # Create T2raw.mgz
        # mri_convert
        ar1_inputs.inputs.Raw_T2 = in_T2
        T2_convert = pe.Node(MRIConvert(), name="T2_convert")
        T2_convert.inputs.out_file = outputfilename('T2raw.mgz', 'mri', 'orig')
        T2_convert.inputs.no_scale = True
        ar1_wf.connect([(ar1_inputs, T2_convert, [('Raw_T2', 'in_file')]),
                        ]) 
        ar3_wf.connect([(ar3_inputs, workflow, [('T2raw', 'Inputs.T2raw')])])
        reconall.connect([(ar1_wf, ar3_wf, [('T2_convert.out_file',
                                             'AutoRecon3_Inputs.T2raw')])])

    if in_FLAIR != None:
        # FLAIR image preparation
        # Create FLAIRraw.mgz
        # mri_convert
        ar1_inputs.inputs.Raw_FLAIR = in_FLAIR
        FLAIR_convert = pe.Node(MRIConvert(), name="FLAIR_convert")
        FLAIR_convert.inputs.out_file = outputfilename(
            'FLAIRraw.mgz', 'mri', 'orig')
        FLAIR_convert.inputs.no_scale = True
        ar1_wf.connect([(ar1_inputs, FLAIR_convert, [('Raw_FLAIR', 'in_file')]),
                        ])


    if qcache:
        for hemisphere in ['lh', 'rh']:
            if hemisphere == 'lh':
                hemi_wf = ar3_lh_wf
                hemi_contrast = lh_contrast
            else:
                hemi_wf = ar3_rh_wf
                hemi_contrast = rh_contrast
            for connection, meas_file, meas_name in [(hemi_wf, 'Make_Pial_Surface.out_thickness', 'thickness'),
                                                     (ar3_inputs, '{0}_area'.format(
                                                         hemisphere), 'area'),
                                                     (hemi_wf, 'Make_Pial_Surface.out_area',
                                                      'area.pial'),
                                                     (hemi_wf, 'Calculate_Volume.out_file',
                                                      'volume'),
                                                     (ar3_inputs, '{0}_curv'.format(
                                                         hemisphere), 'curv'),
                                                     (ar3_inputs, '{0}_sulc'.format(
                                                         hemisphere), 'sulc'),
                                                     (ar3_inputs, '{0}_white_K'.format(
                                                         hemisphere), 'white.K'),
                                                     (ar3_inputs, '{0}_white_H'.format(
                                                         hemisphere), 'white.H'),
                                                     (hemi_wf, 'Jacobian.out_file',
                                                      'jacobian_white'),
                                                     (hemi_contrast, 'out_contrast', 'w-g.pct.mgh')]:
                preprocess = pe.Node(MRISPreprocReconAll(), name="QCache_Preproc_{0}_{1}".format(
                    hemisphere, meas_name.replace('.', '_')))
                target_id = 'fsaverage'
                preprocess.inputs.out_file = outputfilename(
                    '{0}.{1}.{2}.mgh'.format(hemisphere, meas_name, target_id), 'surf')
                target_dir = os.path.join(subjects_dir, target_id)
                if not os.path.isdir(target_dir):
                    # link fsaverage if it doesn't exist
                    target_home = os.path.join(fs_home, 'subjects', target_id)
                    # Create a symlink
                    os.symlink(target_home, target_dir)
                preprocess.inputs.target = target_id
                preprocess.inputs.hemi = hemisphere
                ar3_wf.connect([(ar3_inputs, preprocess, [('subject_id', 'subject_id'),
                                                          ('subjects_dir',
                                                           'subjects_dir'),
                                                          ])])
                ar3_wf.connect([(hemi_wf, preprocess, [('Surface_Registration.out_file', 'surfreg_file')]),
                                (connection, preprocess,
                                 [(meas_file, 'surf_measure_file')])
                                ])

                for value in range(0, 26, 5):
                    surf2surf = pe.Node(SurfaceSmooth(), name="Qcache_{0}_{1}_fwhm{2}".format(
                        hemisphere, meas_name.replace('.', '_'), value))
                    surf2surf.inputs.fwhm = value
                    surf2surf.inputs.cortex = True
                    surf2surf.inputs.subject_id = target_id
                    surf2surf.inputs.hemi = hemisphere
                    tval_file = "{0}.{1}.fwhm{2}.fsaverage.mgh".format(
                        hemisphere, meas_name, value)
                    surf2surf.inputs.out_file = outputfilename(
                        tval_file, 'surf')
                    ar3_lh_wf.connect([(preprocess, surf2surf, [('out_file', 'in_file')]),
                                       (ar3_inputs, surf2surf,
                                        [('subjects_dir', 'subjects_dir')]),
                                       ])

    return reconall

