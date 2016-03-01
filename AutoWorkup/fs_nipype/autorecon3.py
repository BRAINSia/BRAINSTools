import os
import nipype
from nipype.interfaces.utility import Function,IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import *
from ba_maps import create_ba_maps_wf

def create_AutoRecon3(config):
    
    # AutoRecon3
    # Workflow
    ar3_wf = pe.Workflow(name="AutoRecon3")

    # Input Node
    inputspec = pe.Node(IdentityInterface(fields=['lh',
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
                         name='Inputs')

    inputspec.inputs.lh_atlas = os.path.join(
        config['FREESURFER_HOME'], 'average', 
        'lh.average.curvature.filled.buckner40.tif')
    inputspec.inputs.lh_classifier = os.path.join(
        config['FREESURFER_HOME'], 'average', 
        'lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs')
    inputspec.inputs.rh_atlas = os.path.join(
        config['FREESURFER_HOME'], 'average',
        'rh.average.curvature.filled.buckner40.tif')
    inputspec.inputs.rh_classifier = os.path.join(
        config['FREESURFER_HOME'], 'average',
        'rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs')

    ar3_lh_wf1 = pe.Workflow(name="AutoRecon3_Left_1")
    ar3_rh_wf1 = pe.Workflow(name="AutoRecon3_Right_1")
    for hemisphere, hemi_wf in [('lh', ar3_lh_wf1), ('rh', ar3_rh_wf1)]:
        hemi_inputspec1 = pe.Node(IdentityInterface(fields=['subject_id',
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
        ar3_sphere.inputs.out_file = '{0}.sphere'.format(hemisphere)
        if config['openmp'] != None:
            ar3_sphere.inputs.num_threads = config['openmp']
            if config['plugin_args'] != None:
                ar3_sphere.plugin_args = config['plugin_args']
        hemi_wf.connect([(hemi_inputspec1, ar3_sphere, [('inflated', 'in_file'),
                                                        ('smoothwm', 'in_smoothwm')])])

        # Ipsilateral Surface Registation (Spherical Morph)

        # Registers the orig surface to the spherical atlas through surf/?h.sphere.
        # The surfaces are first coarsely registered by aligning the large scale
        # folding patterns found in ?h.sulc and then fine tuned using the small-scale
        # patterns as in ?h.curv. Calls mris_register. Creates surf/?h.sphere.reg.

        ar3_surfreg = pe.Node(Register(), name="Surface_Registration")
        ar3_surfreg.inputs.out_file = '{0}.sphere.reg'.format(hemisphere)
        ar3_surfreg.inputs.curv = True
        hemi_wf.connect([(ar3_sphere, ar3_surfreg, [('out_file', 'in_surf')]),
                         (hemi_inputspec1, ar3_surfreg, [('smoothwm', 'in_smoothwm'),
                                                         ('sulc', 'in_sulc'),
                                                         ('atlas', 'target')])
                       ])

        # Jacobian

        # Computes how much the white surface was distorted in order to register to
        # the spherical atlas during the -surfreg step.

        ar3_jacobian = pe.Node(Jacobian(), name="Jacobian")
        ar3_jacobian.inputs.out_file = '{0}.jacobian_white'.format(hemisphere)
        hemi_wf.connect([(hemi_inputspec1, ar3_jacobian, [('white', 'in_origsurf')]),
                         (ar3_surfreg, ar3_jacobian, [('out_file', 'in_mappedsurf')])
                       ])

        # Average Curvature

        # Resamples the average curvature from the atlas to that of the subject.
        # Allows the user to display activity on the surface of an individual
        # with the folding pattern (ie, anatomy) of a group.

        ar3_paint = pe.Node(Paint(), name="Average_Curvature")
        ar3_paint.inputs.averages = 5
        ar3_paint.inputs.template_param = 6
        ar3_paint.inputs.out_file = "{0}.avg_curv".format(hemisphere)
        hemi_wf.connect([(ar3_surfreg, ar3_paint, [('out_file', 'in_surf')]),
                           (hemi_inputspec1, ar3_paint, [('atlas', 'template')])])

        # Cortical Parcellation

        # Assigns a neuroanatomical label to each location on the cortical
        # surface. Incorporates both geometric information derived from the
        # cortical model (sulcus and curvature), and neuroanatomical convention.

        ar3_parcellation = pe.Node(MRIsCALabel(), "Cortical_Parcellation")
        ar3_parcellation.inputs.seed = 1234
        ar3_parcellation.inputs.hemisphere = hemisphere
        ar3_parcellation.inputs.copy_smoothwm = True
        ar3_parcellation.inputs.out_file = "{0}.aparc.annot".format(hemisphere)        

        if config['openmp'] != None:
            ar3_parcellation.inputs.num_threads = config['openmp']
            if config['plugin_args'] != None:
                ar3_parcellation.plugin_args = config['plugin_args']
        hemi_wf.connect([(hemi_inputspec1, ar3_parcellation, [('smoothwm', 'smoothwm'),
                                                              ('cortex_label', 'label'),
                                                              ('aseg_presurf',  'aseg'),
                                                              ('classifier', 'classifier')]),
                         (ar3_surfreg, ar3_parcellation, [('out_file', 'canonsurf')])
                               ])

        # Pial Surface

        ar3_pial = pe.Node(MakeSurfaces(), name="Make_Pial_Surface")
        ar3_pial.inputs.no_white = True
        ar3_pial.inputs.mgz = True
        ar3_pial.inputs.hemisphere = hemisphere

        hemi_wf.connect([(hemi_inputspec1, ar3_pial, [('wm', 'in_wm'),
                                                      ('orig', 'in_orig'),
                                                      ('filled', 'in_filled'),
                                                      ('white', 'orig_pial'),
                                                      ('white', 'orig_white'),
                                                      ('brain_finalsurfs', 'in_T1'),
                                                      ('aseg_presurf', 'in_aseg')]),
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
        ar3_add.inputs.out_file = '{0}.area.mid'.format(hemisphere)
        hemi_wf.connect([(ar3_pial, ar3_add, [('out_area', 'in_file2')]),
                         (hemi_inputspec1, ar3_add, [('area', 'in_file1')]),
                       ])

        ar3_divide = pe.Node(MRIsCalc(), name="Mid_Pial")
        ar3_divide.inputs.action = "div"
        ar3_divide.inputs.in_int = 2
        ar3_divide.inputs.out_file = '{0}.area.mid'.format(hemisphere)
        hemi_wf.connect([(ar3_add, ar3_divide, [('out_file', 'in_file1')]),
                       ])

        ar3_volume = pe.Node(MRIsCalc(), name="Calculate_Volume")
        ar3_volume.inputs.action = "mul"
        ar3_volume.inputs.out_file = '{0}.volume'.format(hemisphere)
        hemi_wf.connect([(ar3_divide, ar3_volume, [('out_file', 'in_file1')]),
                         (ar3_pial, ar3_volume, [('out_thickness', 'in_file2')]),
                       ])

        # Connect the inputs
        ar3_wf.connect([(inputspec, hemi_wf, [('{0}_inflated'.format(hemisphere), 'Inputs.inflated'),
                                                ('{0}_smoothwm'.format(hemisphere),
                                                 'Inputs.smoothwm'),
                                                ('{0}_white'.format(hemisphere), 'Inputs.white'),
                                                ('{0}_cortex_label'.format(hemisphere),
                                                 'Inputs.cortex_label'),
                                                ('{0}_orig'.format(hemisphere), 'Inputs.orig'),
                                                ('{0}_sulc'.format(hemisphere), 'Inputs.sulc'),
                                                ('{0}_area'.format(hemisphere), 'Inputs.area'),
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

        # Workflow1 Outputs
        hemi_outputs1 = ['sphere',
                         'sphere_reg',
                         'jacobian_white',
                         'avg_curv',
                         'aparc_annot',
                         'area',
                         'curv',
                         'pial',
                         'thickness',
                         'area_mid',
                         'volume']
        hemi_outputspec1 = pe.Node(IdentityInterface(fields=hemi_outputs1),
                              name="Outputs")
        hemi_wf.connect([(ar3_pial, hemi_outputspec1, [('out_pial', 'pial'),
                                                       ('out_curv', 'curv'),
                                                       ('out_area', 'area'),
                                                       ('out_thickness', 'thickness')]),
                         (ar3_divide, hemi_outputspec1, [('out_file', 'area_mid')]),
                         (ar3_volume, hemi_outputspec1, [('out_file', 'volume')]),
                         (ar3_parcellation, hemi_outputspec1, [('out_file', 'aparc_annot')]),
                         (ar3_paint, hemi_outputspec1, [('out_file', 'jacobian_white')]),
                         (ar3_surfreg, hemi_outputspec1, [('out_file', 'sphere_reg')]),
                         (ar3_sphere, hemi_outputspec1, [('out_file', 'sphere')])
                           ])


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
    volume_mask.inputs.copy_inputs = True

    
    ar3_wf.connect([(inputspec, volume_mask, [('aseg_presurf', 'in_aseg'),
                                              ('lh_white', 'lh_white'),
                                              ('rh_white', 'rh_white'),
                                               ]),
                    (ar3_lh_wf1, volume_mask, [('Outputs.pial', 'lh_pial')]),
                    (ar3_rh_wf1, volume_mask, [('Outputs.pial', 'rh_pial')]),
                    ])


    for hemisphere, workflow in [('lh', ar3_lh_wf1),  ('rh', ar3_rh_wf1)]:
        if hemisphere == 'lh':
            opp_hemi = 'rh'
            opp_wf = ar3_rh_wf1
        else:
            opp_hemi = 'lh'
            opp_wf = ar3_lh_wf1

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
        parcellation_stats_white.inputs.out_color = os.path.join(
            config['subjects_dir'], config['current_id'], 
            'label', 'aparc.annot.ctab')
        parcellation_stats_white.inputs.out_table = os.path.join(
            config['subjects_dir'], config['current_id'], 
            'stats', '{0}.aparc.stats'.format(hemisphere))

        ar3_wf.connect([(inputspec, parcellation_stats_white, [('subject_id', 'subject_id'),
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
                        (workflow, parcellation_stats_white, [('Make_Pial_Surface.out_pial',
                                                               '{0}_pial'.format(hemisphere)),
                                                              ('Make_Pial_Surface.out_thickness',
                                                               'thickness'),
                                                              ('Cortical_Parcellation.out_file',
                                                               'in_annotation')
                                                              ]),
                        (opp_wf, parcellation_stats_white, [('Make_Pial_Surface.out_pial',
                                                             '{0}_pial'.format(opp_hemi)),
                                                           ]),
                        (volume_mask, parcellation_stats_white, [('out_ribbon', 'ribbon')
                                                             ]),
                    ])
        
        parcellation_stats_pial = pe.Node(
            ParcellationStats(), name="Parcellation_Stats_{0}_Pial".format(hemisphere) )
        parcellation_stats_pial.inputs.mgz = True
        parcellation_stats_pial.inputs.tabular_output = True
        parcellation_stats_pial.inputs.surface = 'pial'
        parcellation_stats_pial.inputs.out_color = os.path.join(
            config['subjects_dir'], config['current_id'], 
            'label', 'aparc.annot.ctab')
        parcellation_stats_pial.inputs.out_table = os.path.join(
            config['subjects_dir'], config['current_id'], 
            'stats', '{0}.aparc.pial.stats'.format(hemisphere))

        ar3_wf.connect([(inputspec, parcellation_stats_pial, [('subject_id', 'subject_id'),
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
                    (workflow, parcellation_stats_pial, [('Make_Pial_Surface.out_pial',
                                                          '{0}_pial'.format(hemisphere)),
                                                            ('Make_Pial_Surface.out_thickness',
                                                             'thickness'),
                                                            ('Cortical_Parcellation.out_file',
                                                             'in_annotation')
                                                            ]),
                    (opp_wf, parcellation_stats_pial, [('Make_Pial_Surface.out_pial',
                                                        '{0}_pial'.format(opp_hemi)),
                                                            ]),
                    (volume_mask, parcellation_stats_pial, [('out_ribbon', 'ribbon')
                                                              ]),
                    ])

        # Cortical Parcellation 2
        cortical_parcellation_2 = pe.Node(MRIsCALabel(),
                                          name="Cortical_Parcellation_{0}_2".format(hemisphere))
        if hemisphere == 'lh':
            lh_cort_parc2 = cortical_parcellation_2
        else:
            rh_cort_parc2 = cortical_parcellation_2
            
        cortical_parcellation_2.inputs.classifier = os.path.join(
            config['FREESURFER_HOME'], 'average', 
            '{0}.destrieux.simple.2009-07-29.gcs'.format(hemisphere))
        cortical_parcellation_2.inputs.out_file = os.path.join(
            config['subjects_dir'], config['current_id'], 
            'label', '{0}.aparc.a2009s.annot'.format(hemisphere))
        cortical_parcellation_2.inputs.seed = 1234

        ar3_wf.connect([(inputspec, cortical_parcellation_2, [('subject_id', 'subject_id'),
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
        parcellation_stats_white_2.inputs.out_color = os.path.join(
            config['subjects_dir'], config['current_id'],
            'label', 'aparc.annot.a2009s.ctab')
        parcellation_stats_white_2.inputs.out_table = os.path.join(
            config['subjects_dir'], config['current_id'], 
            'stats', '{0}.aparc.a2009s.stats'.format(hemisphere))
        ar3_wf.connect([(inputspec, parcellation_stats_white_2, [('subject_id', 'subject_id'),
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
                        (workflow, parcellation_stats_white_2, [('Make_Pial_Surface.out_pial',
                                                                 '{0}_pial'.format(hemisphere)),
                                                                ('Make_Pial_Surface.out_thickness',
                                                                 'thickness'),
                                                            ]),
                        (opp_wf, parcellation_stats_white_2, [('Make_Pial_Surface.out_pial',
                                                               '{0}_pial'.format(opp_hemi)),
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
        cortical_parcellation_3.inputs.out_file = os.path.join(
            config['subjects_dir'], config['current_id'], 
            'label', '{0}.aparc.DKTatlas40.annot'.format(hemisphere))
        cortical_parcellation_3.inputs.seed = 1234
        ar3_wf.connect([(inputspec, cortical_parcellation_3, [('subject_id', 'subject_id'),
                                                               ('subjects_dir',
                                                                'subjects_dir'),
                                                               ('{0}'.format(hemisphere), 'hemisphere'),
                                                               ('{0}_smoothwm'.format(hemisphere),
                                                                'smoothwm'),
                                                               ('aseg_presurf', 'aseg'),
                                                               ('{0}_cortex_label'.format(hemisphere),
                                                                'label')]),
                        (workflow, cortical_parcellation_3, [
                            ('Surface_Registration.out_file', 'canonsurf')])
                    ])

        # Parcellation Statistics 3
        parcellation_stats_white_3 = parcellation_stats_white.clone(
            name="Parcellation_Statistics_{0}_3".format(hemisphere))
        parcellation_stats_white_3.inputs.out_color = os.path.join(
            config['subjects_dir'], config['current_id'],
            'label', 'aparc.annot.DKTatlas40.ctab')
        parcellation_stats_white_3.inputs.out_table = os.path.join(
            config['subjects_dir'], config['current_id'], 
            'stats', '{0}.aparc.DKTatlas40.stats'.format(hemisphere))

        ar3_wf.connect([(inputspec, parcellation_stats_white_3, [('subject_id', 'subject_id'),
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
                        (workflow, parcellation_stats_white_3, [('Make_Pial_Surface.out_pial',
                                                                 '{0}_pial'.format(hemisphere)),
                                                                  ('Make_Pial_Surface.out_thickness',
                                                                   'thickness'),
                                                              ]),
                        (opp_wf, parcellation_stats_white_3, [('Make_Pial_Surface.out_pial',
                                                               '{0}_pial'.format(opp_hemi)),
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

        ar3_wf.connect([(inputspec, contrast, [('orig_mgz', 'orig'),
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
    ar3_wf.connect([(inputspec, relabel_hypos, [('aseg_presurf', 'aseg'),
                                                 ('lh_white', 'lh_white'),
                                                 ('rh_white', 'rh_white'),
                                                 ])])

    # APARC to ASEG
    # Adds information from the ribbon into the aseg.mgz (volume parcellation).
    aparc_2_aseg = pe.Node(Aparc2Aseg(), name="Aparc2Aseg")
    aparc_2_aseg.inputs.volmask = True
    ar3_wf.connect([(inputspec, aparc_2_aseg, [('lh_white', 'lh_white'),
                                                ('rh_white', 'rh_white'),
                                                ('subject_id', 'subject_id'),
                                                ('subjects_dir',
                                                 'subjects_dir'),
                                                ]),
                    (ar3_lh_wf1, aparc_2_aseg, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                               ('Cortical_Parcellation.out_file',
                                                'lh_annotation'),
                                               ]),
                    (ar3_rh_wf1, aparc_2_aseg, [('Make_Pial_Surface.out_pial', 'rh_pial'),
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
    ar3_wf.connect([(inputspec, aparc_2_aseg_2009, [('lh_white', 'lh_white'),
                                                     ('rh_white', 'rh_white'),
                                                     ('subject_id',
                                                      'subject_id'),
                                                     ('subjects_dir',
                                                      'subjects_dir'),
                                                     ]),
                    (ar3_lh_wf1, aparc_2_aseg_2009, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                                    ]),
                    (lh_cort_parc2, aparc_2_aseg_2009,
                     [('out_file', 'lh_annotation')]),
                    (rh_cort_parc2, aparc_2_aseg_2009,
                     [('out_file', 'rh_annotation')]),
                    (ar3_rh_wf1, aparc_2_aseg_2009, [('Make_Pial_Surface.out_pial', 'rh_pial'),
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
                    (inputspec, segstats, [('subject_id', 'subject_id'),
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
                    (ar3_lh_wf1, segstats, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                           ]),
                    (ar3_rh_wf1, segstats, [('Make_Pial_Surface.out_pial', 'rh_pial'),
                                           ]),
                    ])

    # White Matter Parcellation

    # Adds WM Parcellation info into the aseg and computes stat.

    wm_parcellation = pe.Node(Aparc2Aseg(), name="WM_Parcellation")
    wm_parcellation.inputs.volmask = True
    wm_parcellation.inputs.label_wm = True
    wm_parcellation.inputs.hypo_wm = True
    wm_parcellation.inputs.rip_unknown = True

    ar3_wf.connect([(inputspec, wm_parcellation, [('lh_white', 'lh_white'),
                                                   ('rh_white', 'rh_white'),
                                                   ('subject_id',
                                                    'subject_id'),
                                                   ('subjects_dir',
                                                    'subjects_dir'),
                                                   ]),
                    (ar3_lh_wf1, wm_parcellation, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                                  ('Cortical_Parcellation.out_file',
                                                   'lh_annotation'),
                                                  ]),
                    (ar3_rh_wf1, wm_parcellation, [('Make_Pial_Surface.out_pial', 'rh_pial'),
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
                    (inputspec, wm_segstats, [('subject_id', 'subject_id'),
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
                    (ar3_lh_wf1, wm_segstats, [('Make_Pial_Surface.out_pial', 'lh_pial'),
                                              ]),
                    (ar3_rh_wf1, wm_segstats, [('Make_Pial_Surface.out_pial', 'rh_pial'),
                                              ]),
                    ])

    # add brodman area maps to the workflow
    ba_WF = create_ba_maps_wf(config)

    ar3_wf.connect([(ar3_lh_wf1, ba_WF, [('Surface_Registration.out_file', 'Inputs.lh_sphere_reg'),
                                        ('Make_Pial_Surface.out_thickness',
                                         'Inputs.lh_thickness'),
                                        ('Make_Pial_Surface.out_pial',
                                         'Inputs.lh_pial'),
                                        ]),
                    (ar3_rh_wf1, ba_WF, [('Surface_Registration.out_file', 'Inputs.rh_sphere_reg'),
                                        ('Make_Pial_Surface.out_thickness',
                                         'Inputs.rh_thickness'),
                                        ('Make_Pial_Surface.out_pial',
                                         'Inputs.rh_pial'),
                                        ]),
                    (inputspec, ba_WF, [('subject_id', 'Inputs.subject_id'),
                                         ('subjects_dir',
                                          'Inputs.subjects_dir'),
                                         ('lh_white',
                                          'Inputs.lh_white'),
                                         ('rh_white',
                                          'Inputs.rh_white'),
                                         ('transform',
                                          'Inputs.transform'),
                                         ('aseg_presurf',
                                          'Inputs.aseg'),
                                         ('brainmask',
                                          'Inputs.brainmask'),
                                         ('wm', 'Inputs.wm')]),
                    (volume_mask, ba_WF, [
                        ('out_ribbon', 'Inputs.ribbon')])
                    ])

    if config['qcache']:
        for hemisphere in ['lh', 'rh']:
            if hemisphere == 'lh':
                hemi_wf = ar3_lh_wf1
                hemi_contrast = lh_contrast
            else:
                hemi_wf = ar3_rh_wf1
                hemi_contrast = rh_contrast
            for connection, meas_file, meas_name in [(hemi_wf,
                                                      'Make_Pial_Surface.out_thickness',
                                                      'thickness'),
                                                     (inputspec, '{0}_area'.format(
                                                         hemisphere), 'area'),
                                                     (hemi_wf, 'Make_Pial_Surface.out_area',
                                                      'area.pial'),
                                                     (hemi_wf, 'Calculate_Volume.out_file',
                                                      'volume'),
                                                     (inputspec, '{0}_curv'.format(
                                                         hemisphere), 'curv'),
                                                     (inputspec, '{0}_sulc'.format(
                                                         hemisphere), 'sulc'),
                                                     (inputspec, '{0}_white_K'.format(
                                                         hemisphere), 'white.K'),
                                                     (inputspec, '{0}_white_H'.format(
                                                         hemisphere), 'white.H'),
                                                     (hemi_wf, 'Jacobian.out_file',
                                                      'jacobian_white'),
                                                     (hemi_contrast, 'out_contrast', 'w-g.pct.mgh')]:
                preprocess = pe.Node(MRISPreprocReconAll(), name="QCache_Preproc_{0}_{1}".format(
                    hemisphere, meas_name.replace('.', '_')))
                target_id = 'fsaverage'
                preprocess.inputs.out_file = os.path.join(
                    config['subjects_dir'], config['current_id'], 
                    'surf', '{0}.{1}.{2}.mgh'.format(hemisphere, meas_name, target_id))
                target_dir = os.path.join(config['subjects_dir'], target_id)
                if not os.path.isdir(target_dir):
                    # link fsaverage if it doesn't exist
                    target_home = os.path.join(config['FREESURFER_HOME'], 'subjects', target_id)
                    # Create a symlink
                    os.symlink(target_home, target_dir)
                preprocess.inputs.target = target_id
                preprocess.inputs.hemi = hemisphere
                ar3_wf.connect([(inputspec, preprocess, [('subject_id', 'subject_id'),
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
                    surf2surf.inputs.out_file = os.path.join(
                        config['subjects_dir'], config['current_id'], 
                        'surf', tval_file)
                    ar3_lh_wf1.connect([(preprocess, surf2surf, [('out_file', 'in_file')]),
                                       (inputspec, surf2surf,
                                        [('subjects_dir', 'subjects_dir')]),
                                       ])

    #TODO: Add outputs to outputspec
    outputspec = pe.Node(IdentityInterface(fields=['aseg']),
                         name="Outputs")
    ar3_wf.connect([(apas_2_aseg, outputspec, [('out_file', 'aseg')])])

    return ar3_wf
