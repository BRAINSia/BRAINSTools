import nipype
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.io import DataGrabber, FreeSurferSource, DataSink
from nipype.interfaces.utility import Merge, IdentityInterface
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

        reconall.connect([(grab_inittp_files, ar1_wf, [('inputvols', 'Inputs.in_T1s'),
                                                       ('iscales', 'Inputs.iscales'),
                                                       ('ltas', 'Inputs.ltas')])])

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
                reconall.connect([(tp_data_source, ar2_wf, [('wm', 'Inputs.init_wm')]),
                                  (tp_data_grabber, ar2_wf, [('subj_to_template_lta', 'Inputs.subj_to_template_lta')]),
                                  (tp_data_grabber, ar2_wf, [('subj_to_template_lta', 'Inputs.subj_to_template_lta')])])

        reconall.connect([(merge_norms, ar2_wf, [('out', 'Inputs.alltps_norms')]),
                          (merge_segs, ar2_wf, [('out', 'Inputs.alltps_segs')]),
                          (merge_template_ltas, ar2_wf, [('out', 'Inputs.alltps_to_template_ltas')]),
                          (merge_segs_noCC, ar2_wf, [('out', 'Inputs.alltps_segs_noCC')])])

                        

        # datasource files from the template run
        ds_template_files = pe.Node(FreeSurferSource(), name="Datasource_Template_Files")
        ds_template_files.inputs.subject_id = config['subject_id']
        ds_template_files.inputs.subjects_dir = config['subjects_dir']

        reconall.connect([(ds_template_files, ar1_wf, [('brainmask', 'Inputs.template_brainmask')]),
                          (ds_template_files, ar2_wf, [('aseg', 'Inputs.template_aseg')])])

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
        reconall.connect([(grab_template_files, ar1_wf, [('template_talairach_xfm', 'Inputs.template_talairach_xfm')]),
                          (grab_template_files, ar2_wf, [('template_talairach_lta', 'Inputs.template_talairach_lta'),
                                                         ('template_talairach_m3z', 'Inputs.template_talairach_m3z'),
                                                         ('template_label_intensities', 'Inputs.template_label_intensities'),
                                                         ('template_lh_white', 'Inputs.template_lh_white'),
                                                         ('template_rh_white', 'Inputs.template_rh_white'),
                                                         ('template_lh_pial', 'Inputs.template_lh_pial'),
                                                         ('template_rh_pial', 'Inputs.template_rh_pial'),])
                          ])
        # end longitudinal data collection

    # connect autorecon 1 - 3 
    reconall.connect([(ar1_wf, ar3_wf, [('Inputs.subject_id', 'Inputs.subject_id'),
                                        ('Inputs.subjects_dir', 'Inputs.subjects_dir'),
                                        ('Outputs.brainmask', 'Inputs.brainmask'),
                                        ('Outputs.talairach', 'Inputs.transform'),
                                        ('Outputs.orig', 'Inputs.orig_mgz'),
                                        ('Outputs.rawavg', 'Inputs.rawavg'),
                                        ]),
                      (ar1_wf, ar2_wf, [('Outputs.brainmask', 'Inputs.brainmask'),
                                        ('Outputs.talairach', 'Inputs.transform'),
                                        ('Outputs.orig', 'Inputs.orig'),
                                        ('Inputs.subject_id', 'Inputs.subject_id'),
                                        ('Inputs.subjects_dir', 'Inputs.subjects_dir'),
                                        ]),
                      (ar2_lh, ar3_wf, [('inflate2.out_file', 'Inputs.lh_inflated'),
                                        ('Smooth2.surface',
                                         'Inputs.lh_smoothwm'),
                                        ('Make_Surfaces.out_white',
                                         'Inputs.lh_white'),
                                        ('Make_Surfaces.out_cortex',
                                         'Inputs.lh_cortex_label'),
                                        ('Make_Surfaces.out_area',
                                         'Inputs.lh_area'),
                                        ('Make_Surfaces.out_curv',
                                         'Inputs.lh_curv'),
                                        ('inflate2.out_sulc',
                                         'Inputs.lh_sulc'),
                                        ('Extract_Main_Component.out_file',
                                         'Inputs.lh_orig_nofix'),
                                        ('Remove_Intersection.out_file',
                                         'Inputs.lh_orig'),
                                        ('Curvature1.out_mean',
                                         'Inputs.lh_white_H'),
                                        ('Curvature1.out_gauss',
                                         'Inputs.lh_white_K'),
                                        ]),
                      (ar2_rh, ar3_wf, [('inflate2.out_file', 'Inputs.rh_inflated'),
                                        ('Smooth2.surface', 'Inputs.rh_smoothwm'),
                                        ('Make_Surfaces.out_white', 'Inputs.rh_white'),
                                        ('Make_Surfaces.out_cortex', 'Inputs.rh_cortex_label'),
                                        ('Make_Surfaces.out_area', 'Inputs.rh_area'),
                                        ('Make_Surfaces.out_curv', 'Inputs.rh_curv'),
                                        ('inflate2.out_sulc', 'Inputs.rh_sulc'),
                                        ('Extract_Main_Component.out_file', 'Inputs.rh_orig_nofix'),
                                        ('Remove_Intersection.out_file', 'Inputs.rh_orig'),
                                        ('Curvature1.out_mean', 'Inputs.rh_white_H'),
                                        ('Curvature1.out_gauss', 'Inputs.rh_white_K'),
                                        ]),
                      (ar2_wf, ar3_wf, [('Outputs.aseg_presurf', 'Inputs.aseg_presurf'),
                                        ('Outputs.brain_finalsurfs', 'Inputs.brain_finalsurfs'),
                                        ('Outputs.wm', 'Inputs.wm'),
                                        ('Outputs.filled', 'Inputs.filled'),
                                        ('Outputs.norm', 'Inputs.norm'),
                                        ]),
                      ])

    #TODO: add more outputs to outputspec
    outputspec = pe.Node(IdentityInterface(fields=['t2_raw',
                                                   'flair',
                                                   'rawavg',
                                                   'orig',
                                                   'orig_nu',
                                                   'talairach_auto',
                                                   'talairach',
                                                   't1',
                                                   'brainmask_auto',
                                                   'brainmask',
                                                   'braintemplate',
                                                   'nu',
                                                   'tal_lta',
                                                   'norm',
                                                   'ctrl_pts',
                                                   'tal_m3z',
                                                   'nu_noneck',
                                                   'talskull2',
                                                   'aseg_noCC',
                                                   'cc_up',
                                                   'aseg_auto',
                                                   'aseg_presurf',
                                                   'brain',
                                                   'brain_finalsurfs',
                                                   'wm_seg',
                                                   'wm_aseg',
                                                   'wm',
                                                   'filled',
                                                   'ponscc_log',
                                                   'aseg',
                                                   'recoded_labelmap']),
                         name="Outputs")
    
    reconall.connect([(ar1_wf, outputspec, [('Outputs.t2_raw', 't2_raw'),
                                            ('Outputs.flair', 'flair'),
                                            ('Outputs.rawavg', 'rawavg'),
                                            ('Outputs.orig', 'orig'),
                                            ('Outputs.orig_nu', 'orig_nu'),
                                            ('Outputs.talairach_auto', 'talairach_auto'),
                                            ('Outputs.talairach', 'talairach'),
                                            ('Outputs.t1', 't1'),
                                            ('Outputs.brainmask_auto', 'brainmask_auto'),
                                            ('Outputs.brainmask', 'brainmask'),
                                            ('Outputs.braintemplate', 'braintemplate'),
                                        ]),
                      (ar2_wf, outputspec, [('Outputs.nu', 'nu'),
                                            ('Outputs.tal_lta', 'tal_lta'),
                                            ('Outputs.norm', 'norm'),
                                            ('Outputs.ctrl_pts', 'ctrl_pts'),
                                            ('Outputs.tal_m3z', 'tal_m3z'),
                                            ('Outputs.nu_noneck', 'nu_noneck'),
                                            ('Outputs.talskull2', 'talskull2'),
                                            ('Outputs.aseg_noCC', 'aseg_noCC'),
                                            ('Outputs.cc_up', 'cc_up'),
                                            ('Outputs.aseg_auto', 'aseg_auto'),
                                            ('Outputs.aseg_presurf', 'aseg_presurf'),
                                            ('Outputs.brain', 'brain'),
                                            ('Outputs.brain_finalsurfs', 'brain_finalsurfs'),
                                            ('Outputs.wm_seg', 'wm_seg'),
                                            ('Outputs.wm_aseg', 'wm_aseg'),
                                            ('Outputs.wm', 'wm'),
                                            ('Outputs.filled', 'filled'),
                                            ('Outputs.ponscc_log', 'ponscc_log'),
                                        ]),
                      (ar3_wf, outputspec, [('Outputs.aseg', 'aseg'),
                                        ]),
                      ])

    #TODO: Datasink outputs
    datasink = pe.Node(DataSink(), name="DataSink")
    datasink.inputs.container = config['current_id']
    datasink.inputs.base_directory = config['subjects_dir']
    
    reconall.connect([(outputspec, datasink, [('t2_raw', 'mri.orig.@t2raw'),
                                              ('flair', 'mri.orig.@flair'),
                                              ('rawavg', 'mri.@rawavg'),
                                              ('orig', 'mri.@orig'),
                                              ('orig_nu', 'mri.@orig_nu'),
                                              ('talairach_auto', 'mri.transforms.@tal_auto'),
                                              ('talairach', 'mri.transforms.@tal'),
                                              ('t1', 'mri.@t1'),
                                              ('brainmask_auto', 'mri.@brainmask_auto'),
                                              ('brainmask', 'mri.@brainmask'),
                                              ('braintemplate', 'mri.@braintemplate'),
                                              ('nu', 'mri.@nu'),
                                              ('tal_lta', 'mri.transforms.@tal_lta'),
                                              ('norm', 'mri.@norm'),
                                              ('ctrl_pts', 'mri.@ctrl_pts'),
                                              ('tal_m3z', 'mri.transforms.@tal_m3z'),
                                              ('nu_noneck', 'mri.@nu_noneck'),
                                              ('talskull2', 'mri.transforms.@talskull2'),
                                              ('aseg_noCC', 'mri.@aseg_noCC'),
                                              ('cc_up', 'mri.transforms.@cc_up'),
                                              ('aseg_auto', 'mri.@aseg_auto'),
                                              ('aseg_presurf', 'mri.@aseg_presuf'),
                                              ('brain', 'mri.@brain'),
                                              ('brain_finalsurfs', 'mri.@brain_finalsurfs'),
                                              ('wm_seg', 'mri.@wm_seg'),
                                              ('wm_aseg', 'mri.@wm_aseg'),
                                              ('wm', 'mri.@wm'),
                                              ('filled', 'mri.@filled'),
                                              ('ponscc_log', 'mri.@ponscc_log'),
                                          ]),
                      ])

    #TODO: DataSink from T1_image_prep
    
    
    #### Workflow additions go here
    if config['recoding_file'] != None:
        from recoding import create_recoding_wf
        recode = create_recoding_wf(config['recoding_file'])
        reconall.connect([(ar3_wf, recode, [('Outputs.aseg', 'Inputs.labelmap')]),
                          (recode, outputspec, [('Outputs.recodedlabelmap', 'recoded_labelmap')])])
        

    return reconall

