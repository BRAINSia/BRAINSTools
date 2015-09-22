import nipype
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.io import DataGrabber, FreeSurferSource
from autorecon1 import mkdir_p, create_AutoRecon1
from autorecon2 import create_AutoRecon2
from autorecon3 import create_AutoRecon3

def create_reconall(in_T1s, subject_id, in_T2, in_FLAIR, subjects_dir, qcache, cw256, fs_home, longitudinal, long_base):
    ar1_wf = create_AutoRecon1(subjects_dir, subject_id, fs_home, in_T1s, in_T2, in_FLAIR, cw256, longitudinal, long_base)
    ar2_wf, ar2_lh, ar2_rh = create_AutoRecon2(subjects_dir, subject_id, fs_home)
    ar3_wf = create_AutoRecon3(subjects_dir, subject_id, fs_home, qcache)

    # Connect workflows 
    reconall = pe.Workflow(name="recon-all")
    if longitudinal:
        # grab files from the initial single session run
        grab_inittp_files = pe.Node(DataGrabber(), name="Grab_Initial_Files".format(subject_id),
                                    infields=['subject_id'],
                                    outfileds=['inputvols', 'iscales', 'ltas'])
        grab_inittp_files.inputs.template = '*'
        grab_inittp_files.inputs.base_directory = subjects_dir
        grab_inittp_files.inputs.field_template = dict(inputvols='%s/mri/orig/0*.mgz',
                                                       iscales='%s/mri/orig/0*-iscale.txt',
                                                       ltas='%s/mri/orig/0*.lta')
        
        grab_inittp_files.inputs.template_args = dict(inputvols=[['subject_id']],
                                                      iscales=[['subject_id']],
                                                      ltas=[['subject_id']])

        reconall.connect([(grab_inittp_files, ar1_wf, [('inputvols', 'AutoRecon1_Inputs.in_T1s'),
                                                       ('iscales', 'AutoRecon1_Inputs.iscales'),
                                                       ('ltas', 'AutoRecon1_Inputs.ltas')])])

        # datasource files from the template run
        ds_template_files = pe.Node(FreeSurferSource(), name="Datasource_Template_Files")
        ds_template_files.inputs.subject_id = subject_id
        ds_template_files.inputs.subjects_dir = subjects_dir

        reconall.connect([(ds_template_files, ar1_wf, [('brainmask', 'AutoRecon1_Inputs.brainmask')])])

        # grab files from template run
        grab_template_files = pe.Node(DataGrabber(), name="Grab_Template_Files",
                                      infields=['subject_id', 'long_base'],
                                      outfields=['subj_to_base',
                                                 'talairach_xfm',
                                                 'talairach_lta',
                                                 'label_intensities'])
        grab_template_files.inputs.template = '*'
        grab_template_files.inputs.base_directory = subjects_dir
        grab_template_files.inputs.subject_id = subject_id
        grab_template_files.inputs.long_base = long_base
        grab_template_files.inputs.field_template = dict(
            subj_to_base='%s/mri/transforms/%s_to_%s.lta',
            talairach_xfm='%s/mri/transfroms/talairach.xfm',
            talairach_lta='%s/mri/transfroms/talairach.lta',
            talairach_m3z='%s/mri/transfroms/talairach.m3z',
            label_intensities='%s/mri/aseg.auto_noCCseg.label_intensities.txt')
        
        grab_template_files.inputs.template_args = dict(subj_to_base=[['long_base', 'subject_id', 'long_base']],
                                                        talairach_xfm=[['long_base']],
                                                        talairach_lta=[['long_base']],
                                                        talairach_m3z=[['long_base']],
                                                        label_intensities=[['long_base']])
        reconall.connect([(grab_template_files, ar1_wf, [('subj_to_base', 'AutoRecon1_Inputs.subj_to_base'),
                                                         ('talairach_xfm', 'AutoRecon1_Inputs.talairach')])])
        # end longitudinal data collection
        
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


    return reconall

