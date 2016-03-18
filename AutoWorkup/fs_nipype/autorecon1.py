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


def awk(awk_file, log_file):
    """
    This method uses 'awk' which must be installed prior to running the workflow and is not a
    part of nipype or freesurfer.
    Future work may be done to create a method that achieves the same results using a python
    script.
    """
    import subprocess
    import os
    command = ['awk', '-f', awk_file, log_file]
    print(''.join(command))
    subprocess.call(command)
    log_file = os.path.abspath(log_file)
    return log_file

def copy_file(in_file, out_file=None):
    """
    Create a function to copy a file that can be modified by a following node without changing the original file
    """
    import os
    import shutil
    if out_file == None:
        out_file = os.path.join(os.getcwd(), os.path.basename(in_file))
    if type(in_file) is list and len(in_file) == 1:
        in_file = in_file[0]
    out_file = os.path.abspath(out_file)
    in_file = os.path.abspath(in_file)
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


def create_orig_dir(in_orig, in_tal):
    """This function creates a temporary directory that includes the talairach.
    This way the talairach can be added to the header of the orig image."""
    import os
    import shutil

    out_orig = os.path.abspath(os.path.basename(in_orig))
    shutil.copy(in_orig, out_orig)
    
    tal_base = os.path.basename(in_tal)
    out_tal = os.path.abspath(os.path.join('transforms', tal_base))
    shutil.copy(in_tal, out_tal)

    return out_orig, out_tal


def create_preproc_filenames(in_T1s):
    # Create output filenames
    inputvols = list()
    iscaleout = list()
    ltaout = list()

    for i, T1 in enumerate(in_T1s):
        file_num = str(i + 1)
        while len(file_num) < 3:
            file_num = '0' + file_num
        iscaleout.append(file_num + '-iscale.txt')
        ltaout.append(file_num + '.lta')
        inputvols.append(file_num + '.mgz')
    return inputvols, iscaleout, ltaout

def create_AutoRecon1(config):
    # AutoRecon1
    # Workflow
    ar1_wf = pe.Workflow(name='AutoRecon1')

    if not config['longitudinal']:
        # single session processing
        inputSpec = pe.Node(interface=IdentityInterface(
            fields=['in_t1s', 'in_t2', 'in_flair']),
                             run_without_submitting=True,
                             name='Inputs')
        inputSpec.inputs.in_t1s = VerifyInputs(config['in_T1s'])

        origvols, iscaleout, ltaout = create_preproc_filenames(config['in_T1s'])

        # T1 image preparation
        # For all T1's mri_convert ${InputVol} ${out_file}
        T1_image_preparation = pe.MapNode(MRIConvert(),
                                          iterfield=['in_file', 'out_file'],
                                          name="T1_prep")
        T1_image_preparation.inputs.out_file = origvols

        ar1_wf.connect([(inputSpec, T1_image_preparation, [('in_t1s', 'in_file')]),
                        ])

        # T2 image preparation
        if config['in_T2'] != None:
            # Create T2raw.mgz
            # mri_convert
            inputSpec.inputs.in_t2 = config['in_T2']
            T2_convert = pe.Node(MRIConvert(), name="T2_convert")
            T2_convert.inputs.out_file = 'T2raw.mgz'
            T2_convert.inputs.no_scale = True
            ar1_wf.connect([(inputSpec, T2_convert, [('in_t2', 'in_file')]),
                            ]) 

        # FLAIR image preparation
        if config['in_FLAIR'] != None:
            # Create FLAIRraw.mgz
            # mri_convert
            inputSpec.inputs.in_flair = config['in_FLAIR']
            FLAIR_convert = pe.Node(MRIConvert(), name="FLAIR_convert")
            FLAIR_convert.inputs.out_file = 'FLAIRraw.mgz'
            FLAIR_convert.inputs.no_scale = True
            ar1_wf.connect([(inputSpec, FLAIR_convert, [('in_flair', 'in_file')]),
                            ])
    else:
        # longitudinal inputs
        inputSpec = pe.Node(interface=IdentityInterface(
            fields=['in_t1s',
                    'iscales',
                    'ltas',
                    'subj_to_template_lta',
                    'template_talairach_xfm',
                    'template_brainmask']),
                             run_without_submitting=True,
                             name='Inputs')

        _volnames_, in_iscales, in_ltas = create_preproc_filenames(config['in_T1s'])

        copy_ltas = pe.MapNode(Function(['in_file', 'out_file'],
                                        ['out_file'],
                                        copy_file),
                               iterfield=['in_file', 'out_file'],
                               name='Copy_ltas')
        ar1_wf.connect([(inputSpec, copy_ltas, [('ltas', 'in_file')])])
        copy_ltas.inputs.out_file = in_ltas

        copy_iscales = pe.MapNode(Function(['in_file', 'out_file'],
                                           ['out_file'],
                                           copy_file),
                                  iterfield=['in_file', 'out_file'],
                                  name='Copy_iscales')
        ar1_wf.connect([(inputSpec, copy_iscales, [('iscales', 'in_file')])])
        copy_iscales.inputs.out_file = in_iscales

        concatenate_lta = pe.MapNode(ConcatenateLTA(), iterfield=['in_file'],
                                     name="Concatenate_ltas")
        ar1_wf.connect([(copy_ltas, concatenate_lta, [('out_file', 'in_file')]),
                        (inputSpec, concatenate_lta, [('subj_to_template_lta', 'subj_to_base')])])

    
    # Motion Correction
    """
    When there are multiple source volumes, this step will correct for small
    motions between them and then average them together.  The output of the
    motion corrected average is mri/rawavg.mgz which is then conformed to
    255 cubed char images (1mm isotropic voxels) in mri/orig.mgz.
    """

    if not config['longitudinal'] and len(config['in_T1s']) == 1:
        # if only one T1 scan is input, just copy the converted scan as the rawavg.mgz
        create_template = pe.Node(Function(['in_file', 'out_file'],
                                           ['out_file'],
                                            copy_file),
                                  name="Robust_Template")
        create_template.inputs.out_file = 'rawavg.mgz'
        
        ar1_wf.connect([(T1_image_preparation, create_template, [('out_file', 'in_file')]),
                    ])
        
        print("WARNING: only one run found. This is OK, but motion correction" +
              "cannot be performed on one run, so I'll copy the run to rawavg" +
              "and continue.")
        
    else:
        # if multiple T1 scans are given
        create_template = pe.Node(RobustTemplate(), name="Robust_Template")
        create_template.inputs.average_metric = 'median'
        create_template.inputs.out_file = 'rawavg.mgz'
        create_template.inputs.no_iteration = True
            
        if config['longitudinal']:
            # if running longitudinally
            ar1_wf.connect([(concatenate_lta, create_template, [('out_file', 'initial_transforms')]),
                            (inputSpec, create_template, [('in_t1s', 'in_files')]),
                            (copy_iscales, create_template, [('out_file','in_intensity_scales')]),
                        ])
        else:
            # if running single session
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
    conform_template.inputs.out_file = 'orig.mgz'
    if config['longitudinal']:
        conform_template.inputs.out_datatype = 'uchar'
    else:
        conform_template.inputs.conform = True
        if len(config['in_T1s']) != 1:
            conform_template.inputs.cw256 = config['cw256']    
            conform_template.inputs.resample_type = 'cubic'
            
    ar1_wf.connect(
        [(create_template, conform_template, [('out_file', 'in_file')])])

    # Talairach
    """
    This computes the affine transform from the orig volume to the MNI305 atlas using Avi Snyders 4dfp
    suite of image registration tools, through a FreeSurfer script called talairach_avi.
    Several of the downstream programs use talairach coordinates as seed points.
    """

    bias_correction = pe.Node(MNIBiasCorrection(), name="Bias_correction")
    bias_correction.inputs.iterations = 1
    bias_correction.inputs.protocol_iterations = 1000
    if config['field_strength'] == '3T':
        # 3T params from Zheng, Chee, Zagorodnov 2009 NeuroImage paper
        # "Improvement of brain segmentation accuracy by optimizing
        # non-uniformity correction using N3"
        # namely specifying iterations, proto-iters and distance: 
        bias_correction.inputs.distance = 50
    else:
        # 1.5T default
        bias_correction.inputs.distance = 200
    # per c.larsen, decrease convergence threshold (default is 0.001)
    bias_correction.inputs.stop = 0.0001
    # per c.larsen, decrease shrink parameter: finer sampling (default is 4)
    bias_correction.inputs.shrink =  2
    # add the mask, as per c.larsen, bias-field correction is known to work
    # much better when the brain area is properly masked, in this case by
    # brainmask.mgz.
    bias_correction.inputs.no_rescale = True
    bias_correction.inputs.out_file = 'orig_nu.mgz'

    ar1_wf.connect([(conform_template, bias_correction, [('out_file', 'in_file')]),
                ])

    if config['longitudinal']:
        # longitudinal processing
        # Just copy the template xfm
        talairach_avi = pe.Node(Function(['in_file', 'out_file'],
                                             ['out_file'],
                                             copy_file),
                                    name='Copy_Template_Transform')
        talairach_avi.inputs.out_file = 'talairach.auto.xfm'

        ar1_wf.connect([(inputSpec, talairach_avi, [('template_talairach_xfm', 'in_file')])])
    else:
        # single session processing
        talairach_avi = pe.Node(TalairachAVI(), name="Compute_Transform")
        
        if config['custom_atlas'] != None:
            # allows to specify a custom atlas
            talairach_avi.inputs.atlas = config['custom_atlas']
            
        talairach_avi.inputs.out_file = 'talairach.auto.xfm'
        
        ar1_wf.connect([(bias_correction, talairach_avi, [('out_file', 'in_file')]),
                    ])
        
    copy_transform = pe.Node(Function(['in_file', 'out_file'],
                                      ['out_file'],
                                      copy_file),
                             name='Copy_Transform')
    copy_transform.inputs.out_file = 'talairach.xfm'

    ar1_wf.connect([(talairach_avi, copy_transform, [('out_file', 'in_file')])])


    # In recon-all the talairach.xfm is added to orig, even though
    # it does not exist yet. This is a compromise to keep from
    # having to change the time stamp of the orig volume after talairaching.
    # Here we are going to add xfm to the header after the xfm has been created.
    # This may mess up the timestamp.

    add_xform_to_orig = pe.Node(AddXFormToHeader(), name="Add_Transform_to_Header")
    add_xform_to_orig.inputs.copy_name = True
    add_xform_to_orig.inputs.out_file = conform_template.inputs.out_file

    ar1_wf.connect([(conform_template, add_xform_to_orig, [('out_file', 'in_file')]),
                    (copy_transform, add_xform_to_orig, [('out_file', 'transform')])])


        
    # check the alignment of the talairach
    # TODO: Figure out how to read output from this node.
    check_alignment = pe.Node(CheckTalairachAlignment(),
                              name="Check_Talairach_Alignment")
    check_alignment.inputs.threshold = 0.005
    ar1_wf.connect([(copy_transform, check_alignment, [('out_file', 'in_file')]),
                    ])

    if not config['longitudinal']:
        awk_logfile = pe.Node(Function(['awk_file', 'log_file'],
                                       ['log_file'],
                                       awk),
                              name='Awk')
        awk_logfile.inputs.awk_file = config['awk_file']
                                       
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
    mri_normalize.inputs.out_file = 'T1.mgz'
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
        copy_template_brainmask.inputs.out_file = 'brainmask_{0}.mgz'.format(config['long_template'])
        
        ar1_wf.connect([(inputSpec, copy_template_brainmask, [('template_brainmask', 'in_file')])])

        mask1 = pe.Node(ApplyMask(), name="ApplyMask1")
        mask1.inputs.keep_mask_deletion_edits = True
        mask1.inputs.out_file = 'brainmask.auto.mgz'
        
        ar1_wf.connect([(mri_normalize, mask1, [('out_file', 'in_file')]),
                        (copy_template_brainmask, mask1, [('out_file', 'mask_file')])])

        brainmask = pe.Node(ApplyMask(), name="ApplyMask2")
        brainmask.inputs.keep_mask_deletion_edits = True
        brainmask.inputs.transfer = 255
        brainmask.inputs.out_file = mask1.inputs.out_file

        ar1_wf.connect([(mask1, brainmask, [('out_file', 'in_file')]),
                        (copy_template_brainmask, brainmask, [('out_file', 'mask_file')])])
                        
    else:    
        mri_em_register = pe.Node(EMRegister(), name="EM_Register")
        mri_em_register.inputs.template = config['registration_template']
        mri_em_register.inputs.out_file = 'talairach_with_skull.lta'
        mri_em_register.inputs.skull = True
        
        if config['openmp'] != None:
            mri_em_register.inputs.num_threads = config['openmp']
        if config['plugin_args'] != None:
            mri_em_register.plugin_args = config['plugin_args']
            
        ar1_wf.connect([(bias_correction, mri_em_register, [('out_file', 'in_file')]),
                    ])

        brainmask = pe.Node(
            WatershedSkullStrip(), name='Watershed_Skull_Strip')
        brainmask.inputs.t1 = True
        brainmask.inputs.brain_atlas = config['registration_template']
        brainmask.inputs.out_file = 'brainmask.auto.mgz'
        ar1_wf.connect([(mri_normalize, brainmask, [('out_file', 'in_file')]),
                        (mri_em_register, brainmask,
                         [('out_file', 'transform')]),
                    ])

    copy_brainmask = pe.Node(Function(['in_file', 'out_file'],
                                      ['out_file'],
                                      copy_file),
                             name='Copy_Brainmask')
    copy_brainmask.inputs.out_file = 'brainmask.mgz'

    ar1_wf.connect([(brainmask, copy_brainmask, [('out_file', 'in_file')])])

    outputs = ['origvols',
               't2_raw',
               'flair',
               'rawavg',
               'orig_nu',
               'orig',
               'talairach_auto',
               'talairach',
               't1',
               'talskull',
               'brainmask_auto',
               'brainmask',
               'braintemplate']
    outputspec = pe.Node(IdentityInterface(fields=outputs),
                         name="Outputs")
    
    ar1_wf.connect([(T1_image_preparation, outputspec, [('out_file', 'origvols')]),
                    (T2_convert, outputspec, [('out_file', 't2_raw')]),
                    (FLAIR_convert, outputspec, [('out_file', 'flair')]),
                    (create_template, outputspec, [('out_file', 'rawavg')]),
                    (add_xform_to_orig, outputspec, [('out_file', 'orig')]),
                    (bias_correction, outputspec, [('out_file', 'oirg_nu')]),
                    (talairach_avi, outputspec, [('out_file', 'talairach_auto')]),
                    (copy_transform, outputspec, [('out_file', 'talairach')]),
                    (mri_normalize, outputspec, [('out_file', 't1')]),
                    (brainmask, outputspec, [('out_file', 'brainmask_auto')]),
                    (copy_brainmask, outputspec, [('out_file', 'brainmask')]),                    
                    ])

    
    if config['longitudinal']:
        ar1_wf.connect([(copy_template_brainmask, outputspec, [('out_file', 'braintemplate')]),
                        ])
    else:
        ar1_wf.connect([(mri_em_register, outputspec, [('out_file', 'talskull')]),
                        ])

    
    return ar1_wf, outputs
