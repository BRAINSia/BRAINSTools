import nipype
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.io import DataGrabber, FreeSurferSource
from nipype.interfaces.utility import Merge
from autorecon1 import mkdir_p, create_AutoRecon1
from autorecon2 import create_AutoRecon2
from autorecon3 import create_AutoRecon3

def create_reconall(config):
    ar1_wf = create_AutoRecon1(config)
    ar2_wf, ar2_lh, ar2_rh = create_AutoRecon2(config)
    ar3_wf = create_AutoRecon3(config)

    # Connect workflows 
    reconall = pe.Workflow(name="recon-all")
    if config['longitudinal']:
        # grab files from the initial single session run
        grab_inittp_files = pe.Node(DataGrabber(), name="Grab_Initial_Files",
                                    infields=['subject_id'],
                                    outfileds=['inputvols', 'iscales', 'ltas'])
        grab_inittp_files.inputs.template = '*'
        grab_inittp_files.inputs.base_directory = config['subjects_dir']
        grab_inittp_files.inputs.field_template = dict(inputvols='%s/mri/orig/0*.mgz',
                                                       iscales='%s/mri/orig/0*-iscale.txt',
                                                       ltas='%s/mri/orig/0*.lta')
        
        grab_inittp_files.inputs.template_args = dict(inputvols=[['subject_id']],
                                                      iscales=[['subject_id']],
                                                      ltas=[['subject_id']])

        reconall.connect([(grab_inittp_files, ar1_wf, [('inputvols', 'AutoRecon1_Inputs.in_T1s'),
                                                       ('iscales', 'AutoRecon1_Inputs.iscales'),
                                                       ('ltas', 'AutoRecon1_Inputs.ltas')])])

        merge_norms = pe.Node(Merge(len(config['timepoints'])), name="Merge_Norms")
        merge_segs = pe.Node(Merge(len(config['timepoints'])), name="Merge_Segmentations")
        merge_segs_noCC = pe.Node(Merge(len(config['timepoints'])), name="Merge_Segmentations_noCC")
        merge_template_ltas = pe.Node(Merge(len(config['timepoints'])), name="Merge_Template_ltas")

        for i, tp in enumerate(config['timepoints']):
            # datasource timepoint files
            tp_data_source = pe.Node(FreeSurferSource(), name="{0}_DataSource".format(tp))
            tp_data_source.inputs.subject_id = tp
            tp_data_source.inputs.subjects_dir = config['subjects_dir']

            tp_data_grabber = pe.Node(DataGrabber(), name="{0}_DataGrabber".format(tp),
                                      infields=['tp', 'long_tempate'],
                                      outfileds=['subj_to_template_lta', 'seg_noCC', 'seg_presurf'])
            tp_data_grabber.inputs.template = '*'
            tp_data_grabber.inputs.base_directory = config['subjects_dir']
            tp_data_grabber.inputs.field_template = dict(
                subj_to_template_lta='%s/mri/transforms/%s_to_%s.lta',
                seg_noCC='%s/mri/aseg.auto_noCCseg.mgz',
                seg_presurf='%s/mri/aseg.presurf.mgz',)

            tp_data_grabber.inputs.template_args = dict(
                subj_to_template_lta=[['long_template', 'tp', 'long_template']],
                seg_noCC=[['tp']],
                seg_presurf=[['tp']])
                        
            reconall.connect([(tp_data_source, merge_norms, [('norm', 'in{0}'.format(i))]),
                              (tp_data_grabber, merge_segs, [('seg_presurf', 'in{0}'.format(i))]),
                              (tp_data_grabber, merge_segs_noCC, [('seg_noCC', 'in{0}'.format(i))]),
                              (tp_data_grabber, merge_template_ltas, [('subj_to_template_lta', 'in{0}'.format(i))])])

            if tp == config['subject_id']:
                reconall.connect([(tp_data_source, ar2_wf, [('wm', 'AutoRecon2_Inputs.init_wm')]),
                                  (tp_data_grabber, ar2_wf, [('subj_to_template_lta', 'AutoRecon2_Inputs.subj_to_template_lta')]),
                                  (tp_data_grabber, ar2_wf, [('subj_to_template_lta', 'AutoRecon1_Inputs.subj_to_template_lta')])])

        reconall.connect([(merge_norms, ar2_wf, [('out', 'AutoRecon2_Inputs.alltps_norms')]),
                          (merge_segs, ar2_wf, [('out', 'AutoRecon2_Inputs.alltps_segs')]),
                          (merge_template_ltas, ar2_wf, [('out', 'AutoRecon2_Inputs.alltps_to_template_ltas')]),
                          (merge_segs_noCC, ar2_wf, [('out', 'AutoRecon2_Inputs.alltps_segs_noCC')])])

                        

        # datasource files from the template run
        ds_template_files = pe.Node(FreeSurferSource(), name="Datasource_Template_Files")
        ds_template_files.inputs.subject_id = config['subject_id']
        ds_template_files.inputs.subjects_dir = config['subjects_dir']

        reconall.connect([(ds_template_files, ar1_wf, [('brainmask', 'AutoRecon1_Inputs.template_brainmask')]),
                          (ds_template_files, ar2_wf, [('aseg', 'AutoRecon2_Inputs.template_aseg')])])

        # grab files from template run
        grab_template_files = pe.Node(DataGrabber(), name="Grab_Template_Files",
                                      infields=['subject_id', 'long_template'],
                                      outfields=['template_talairach_xfm',
                                                 'template_talairach_lta',
                                                 'template_talairach_m3z',
                                                 'template_label_intensities',
                                                 'template_lh_white',
                                                 'template_rh_white',
                                                 'template_lh_pial',
                                                 'template_rh_pial'])
        grab_template_files.inputs.template = '*'
        grab_template_files.inputs.base_directory = config['subjects_dir']
        grab_template_files.inputs.subject_id = config['subject_id']
        grab_template_files.inputs.long_template = config['long_template']
        grab_template_files.inputs.field_template = dict(
            template_talairach_xfm='%s/mri/transfroms/talairach.xfm',
            template_talairach_lta='%s/mri/transfroms/talairach.lta',
            template_talairach_m3z='%s/mri/transfroms/talairach.m3z',
            template_label_intensities='%s/mri/aseg.auto_noCCseg.label_intensities.txt',
            template_lh_white='%s/surf/lh.white',
            template_rh_white='%s/surf/rh.white',
            template_lh_pial='%s/surf/lh.pial',
            template_rh_pial='%s/surf/rh.pial')
        
        grab_template_files.inputs.template_args = dict(
            template_talairach_xfm=[['long_template']],
            template_talairach_lta=[['long_template']],
            template_talairach_m3z=[['long_template']],
            template_lh_white=[['long_template']],
            template_rh_white=[['long_template']],
            template_lh_pial=[['long_template']],
            template_rh_pial=[['long_template']])
        reconall.connect([(grab_template_files, ar1_wf, [('template_talairach_xfm', 'AutoRecon1_Inputs.template_talairach_xfm')]),
                          (grab_template_files, ar2_wf, [('template_talairach_lta', 'AutoRecon2_Inputs.template_talairach_lta'),
                                                         ('template_talairach_m3z', 'AutoRecon2_Inputs.template_talairach_m3z'),
                                                         ('template_label_intensities', 'AutoRecon2_Inputs.template_label_intensities'),
                                                         ('template_lh_white', 'AutoRecon2_Inputs.template_lh_white'),
                                                         ('template_rh_white', 'AutoRecon2_Inputs.template_rh_white'),
                                                         ('template_lh_pial', 'AutoRecon2_Inputs.template_lh_pial'),
                                                         ('template_rh_pial', 'AutoRecon2_Inputs.template_rh_pial'),])
                          ])
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
                                        ('Fill.out_file',
                                         'AutoRecon3_Inputs.filled'),
                                        ('CA_Normalize.out_file',
                                         'AutoRecon3_Inputs.norm'),
                                        ]),
                      ])

    return reconall

