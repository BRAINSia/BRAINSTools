import nipype
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.io import DataGrabber, FreeSurferSource, DataSink
from nipype.interfaces.utility import Merge, IdentityInterface
from autorecon1 import mkdir_p, create_AutoRecon1
from autorecon2 import create_AutoRecon2
from autorecon3 import create_AutoRecon3

def create_reconall(config):
    ar1_wf, ar1_outputs = create_AutoRecon1(config)
    ar2_wf, ar2_outputs = create_AutoRecon2(config)
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
                        
            reconall.connect([(tp_data_source, merge_norms, [('norm',
                                                              'in{0}'.format(i))]),
                              (tp_data_grabber, merge_segs, [('seg_presurf',
                                                              'in{0}'.format(i))]),
                              (tp_data_grabber, merge_segs_noCC, [('seg_noCC',
                                                                   'in{0}'.format(i))]),
                              (tp_data_grabber, merge_template_ltas, [('subj_to_template_lta',
                                                                       'in{0}'.format(i))])])

            if tp == config['subject_id']:
                reconall.connect([(tp_data_source, ar2_wf, [('wm', 'Inputs.init_wm')]),
                                  (tp_data_grabber, ar2_wf, [('subj_to_template_lta',
                                                              'Inputs.subj_to_template_lta')]),
                                  (tp_data_grabber, ar2_wf, [('subj_to_template_lta',
                                                              'Inputs.subj_to_template_lta')])])

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
        reconall.connect([(grab_template_files, ar1_wf, [('template_talairach_xfm',
                                                          'Inputs.template_talairach_xfm')]),
                          (grab_template_files, ar2_wf, [('template_talairach_lta',
                                                          'Inputs.template_talairach_lta'),
                                                         ('template_talairach_m3z',
                                                          'Inputs.template_talairach_m3z'),
                                                         ('template_label_intensities',
                                                          'Inputs.template_label_intensities'),
                                                         ('template_lh_white', 'Inputs.template_lh_white'),
                                                         ('template_rh_white', 'Inputs.template_rh_white'),
                                                         ('template_lh_pial', 'Inputs.template_lh_pial'),
                                                         ('template_rh_pial', 'Inputs.template_rh_pial'),])
                          ])
        # end longitudinal data collection

    # connect autorecon 1 - 3 
    reconall.connect([(ar1_wf, ar3_wf, [('Outputs.brainmask', 'Inputs.brainmask'),
                                        ('Outputs.talairach', 'Inputs.transform'),
                                        ('Outputs.orig', 'Inputs.orig_mgz'),
                                        ('Outputs.rawavg', 'Inputs.rawavg'),
                                        ]),
                      (ar1_wf, ar2_wf, [('Outputs.brainmask', 'Inputs.brainmask'),
                                        ('Outputs.talairach', 'Inputs.transform'),
                                        ('Outputs.orig', 'Inputs.orig'),
                                        ]),
                      (ar2_wf, ar3_wf, [('Outputs.aseg_presurf', 'Inputs.aseg_presurf'),
                                        ('Outputs.brain_finalsurfs', 'Inputs.brain_finalsurfs'),
                                        ('Outputs.wm', 'Inputs.wm'),
                                        ('Outputs.filled', 'Inputs.filled'),
                                        ('Outputs.norm', 'Inputs.norm'),
                                        ]),
                      ])
    for hemi in ('lh', 'rh'):
        reconall.connect([(ar2_wf, ar3_wf, [('Outputs.{0}_inflated'.format(hemi),
                                             'Inputs.{0}_inflated'.format(hemi)),
                                            ('Outputs.{0}_smoothwm'.format(hemi),
                                             'Inputs.{0}_smoothwm'.format(hemi)),
                                            ('Outputs.{0}_white'.format(hemi),
                                             'Inputs.{0}_white'.format(hemi)),
                                            ('Outputs.{0}_cortex'.format(hemi),
                                             'Inputs.{0}_cortex_label'.format(hemi)),
                                            ('Outputs.{0}_area'.format(hemi),
                                             'Inputs.{0}_area'.format(hemi)),
                                            ('Outputs.{0}_curv'.format(hemi),
                                             'Inputs.{0}_curv'.format(hemi)),
                                            ('Outputs.{0}_sulc'.format(hemi),
                                             'Inputs.{0}_sulc'.format(hemi)),
                                            ('Outputs.{0}_orig_nofix'.format(hemi),
                                             'Inputs.{0}_orig_nofix'.format(hemi)),
                                            ('Outputs.{0}_orig'.format(hemi),
                                             'Inputs.{0}_orig'.format(hemi)),
                                            ('Outputs.{0}_white_H'.format(hemi),
                                             'Inputs.{0}_white_H'.format(hemi)),
                                            ('Outputs.{0}_white_K'.format(hemi),
                                             'Inputs.{0}_white_K'.format(hemi))])])


    #TODO: add more outputs to outputspec
    outputs = ar1_outputs + ar2_outputs
    outputspec = pe.Node(IdentityInterface(fields=outputs),
                         name="Outputs")

    for outfields, wf in [(ar1_outputs, ar1_wf), (ar2_outputs, ar2_wf)]:
        for field in outfields:
            reconall.connect([(wf, outputspec, [('Outputs.' + field, field)])])

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
                                              ('lh_orig_nofix', 'surf.@lh_orig_nofix'),
                                              ('lh_orig', 'surf.@lh_orig'),
                                              ('lh_smoothwm_nofix', 'surf.@lh_smoothwm_nofix'),
                                              ('lh_inflated_nofix', 'surf.@lh_inflated_nofix'),
                                              ('lh_qsphere_nofix', 'surf.@lh_qsphere_nofix'),
                                              ('lh_white', 'surf.@lh_white'),
                                              ('lh_curv', 'surf.@lh_curv'),
                                              ('lh_area', 'surf.@lh_area'),
                                              ('lh_cortex', 'surf.@lh_cortex'),
                                              ('lh_pial', 'surf.@lh_pial'),
                                              ('lh_thickness', 'surf.@lh_thickness'),
                                              ('lh_smoothwm', 'surf.@lh_smoothwm'),
                                              ('lh_sulc', 'surf.@lh_sulc'),
                                              ('lh_inflated', 'surf.@lh_inflated'),
                                              ('lh_white_H', 'surf.@lh_white_H'),
                                              ('lh_white_K', 'surf.@lh_white_K'),
                                              ('lh_inflated_H', 'surf.@lh_inflated_H'),
                                              ('lh_inflated_K', 'surf.@lh_inflated_K'),
                                              ('lh_curv_stats', 'surf.@lh_curv_stats'),
                                              ('rh_orig_nofix', 'surf.@rh_orig_nofix'),
                                              ('rh_orig', 'surf.@rh_orig'),
                                              ('rh_smoothwm_nofix', 'surf.@rh_smoothwm_nofix'),
                                              ('rh_inflated_nofix', 'surf.@rh_inflated_nofix'),
                                              ('rh_qsphere_nofix', 'surf.@rh_qsphere_nofix'),
                                              ('rh_white', 'surf.@rh_white'),
                                              ('rh_curv', 'surf.@rh_curv'),
                                              ('rh_area', 'surf.@rh_area'),
                                              ('rh_cortex', 'surf.@rh_cortex'),
                                              ('rh_pial', 'surf.@rh_pial'),
                                              ('rh_thickness', 'surf.@rh_thickness'),
                                              ('rh_smoothwm', 'surf.@rh_smoothwm'),
                                              ('rh_sulc', 'surf.@rh_sulc'),
                                              ('rh_inflated', 'surf.@rh_inflated'),
                                              ('rh_white_H', 'surf.@rh_white_H'),
                                              ('rh_white_K', 'surf.@rh_white_K'),
                                              ('rh_inflated_H', 'surf.@rh_inflated_H'),
                                              ('rh_inflated_K', 'surf.@rh_inflated_K'),
                                              ('rh_curv_stats', 'surf.@rh_curv_stats').
                                              ('lh_aparc_annot_ctab', 'label.@aparc_annot_ctab'),
                                              ('lh_aparc_stats', 'stats.@lh_aparc_stats'),
                                              ('rh_aparc_stats', 'stats.@rh_aparc_stats'),
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

