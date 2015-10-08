import os
import nipype
from nipype.interfaces.utility import Function,IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import *
from autorecon1 import outputfilename

def create_AutoRecon3(config['subjects_dir'], config['current_id'], config['FREESURFER_HOME'], config['qcache']):
    
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
                                                   'brainmask',
                                                   'transform',
                                                   'orig_mgz',
                                                   'rawavg',
                                                   'norm']),
                         name='AutoRecon3_Inputs')

    ar3_inputs.inputs.lh = 'lh'
    ar3_inputs.inputs.rh = 'rh'

    ar3_inputs.inputs.lh_atlas = os.path.join(
        config['FREESURFER_HOME'], 'average/lh.average.curvature.filled.buckner40.tif')
    ar3_inputs.inputs.lh_classifier = os.path.join(
        config['FREESURFER_HOME'], 'average/lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs')
    ar3_inputs.inputs.rh_atlas = os.path.join(
        config['FREESURFER_HOME'], 'average/rh.average.curvature.filled.buckner40.tif')
    ar3_inputs.inputs.rh_classifier = os.path.join(
        config['FREESURFER_HOME'], 'average/rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs')

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
    if config['plugin_args'] != None:
        ar3_sphere.plugin_args = config['plugin_args']
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
    if config['plugin_args'] != None:
        ar3_parcellation.plugin_args = config['plugin_args']
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
        parcellation_stats_white.inputs.out_color = outputfilename(config['subjects_dir'], config['current_id'], 'aparc.annot.ctab', 'label')
        parcellation_stats_white.inputs.out_table = outputfilename(config['subjects_dir'], config['current_id'], '{0}.aparc.stats'.format(hemisphere), 'stats')

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
        parcellation_stats_pial.inputs.out_color = outputfilename(config['subjects_dir'], config['current_id'], 'aparc.annot.ctab', 'label')
        parcellation_stats_pial.inputs.out_table = outputfilename(config['subjects_dir'], config['current_id'], '{0}.aparc.pial.stats'.format(hemisphere), 'stats')

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
            config['FREESURFER_HOME'], 'average', '{0}.destrieux.simple.2009-07-29.gcs'.format(hemisphere))
        cortical_parcellation_2.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], '{0}.aparc.a2009s.annot'.format(hemisphere), 'label')
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
        parcellation_stats_white_2.inputs.out_color = outputfilename(config['subjects_dir'], config['current_id'], 'aparc.annot.a2009s.ctab', 'label')
        parcellation_stats_white_2.inputs.out_table = outputfilename(config['subjects_dir'], config['current_id'], '{0}.aparc.a2009s.stats'.format(hemisphere), 'stats')
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
            config['FREESURFER_HOME'], 'average', '{0}.DKTatlas40.gcs'.format(hemisphere))
        cortical_parcellation_3.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], '{0}.aparc.DKTatlas40.annot'.format(hemisphere), 'label')
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
        parcellation_stats_white_3.inputs.out_color = outputfilename(config['subjects_dir'], config['current_id'], 'aparc.annot.DKTatlas40.ctab', 'label')
        parcellation_stats_white_3.inputs.out_table = outputfilename(config['subjects_dir'], config['current_id'], '{0}.aparc.DKTatlas40.stats'.format(hemisphere), 'stats')

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
        config['FREESURFER_HOME'], 'ASegStatsLUT.txt')
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
        config['FREESURFER_HOME'], 'WMParcStatsLUT.txt')
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
                config['FREESURFER_HOME'], 'average', 'colortable_BA.txt')
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

    if config['qcache']:
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
                preprocess.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], 
                    '{0}.{1}.{2}.mgh'.format(hemisphere, meas_name, target_id), 'surf')
                target_dir = os.path.join(config['subjects_dir'], target_id)
                if not os.path.isdir(target_dir):
                    # link fsaverage if it doesn't exist
                    target_home = os.path.join(config['FREESURFER_HOME'], 'subjects', target_id)
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
                    surf2surf.inputs.out_file = outputfilename(config['subjects_dir'], config['current_id'], 
                        tval_file, 'surf')
                    ar3_lh_wf.connect([(preprocess, surf2surf, [('out_file', 'in_file')]),
                                       (ar3_inputs, surf2surf,
                                        [('subjects_dir', 'subjects_dir')]),
                                       ])

    return ar3_wf
