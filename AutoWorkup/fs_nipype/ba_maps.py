import os
import nipype
from nipype.interfaces.utility import Function,IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import *

def create_ba_maps_wf(config):
    # Brodmann Area Maps (BA Maps) and Hinds V1 Atlas
    inputs = ['lh_sphere_reg',
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
              'wm']
    
    inputSpec = pe.Node(IdentityInterface(fields=inputs),
                        name="Inputs")

    ba_WF = pe.Workflow(name="Brodmann_Area_Maps")
        
    outputSpec = pe.Node(IdentityInterface(fields=['lh_table',
                                                   'lh_color',
                                                   'lh_thresh_table',
                                                   'lh_thresh_color',
                                                   'rh_table',
                                                   'rh_color',
                                                   'rh_thresh_table',
                                                   'rh_thresh_color']),
                         name="Outputs")

    labels = ["BA1", "BA2", "BA3a", "BA3b", "BA4a", "BA4p", "BA6",
              "BA44", "BA45", "V1", "V2", "MT", "entorhinal", "perirhinal"]
    for hemisphere in ['lh', 'rh']:
        for threshold in [True, False]:
            ## TODO: Create sources subject and merge outputs to a list
            source_subject = pe.Node(nio.DataGrabber(outfields=labels+['sphere_reg']),
                                     name="Source_Subject")
            source_subject.inputs.base_directory = config['source_subject']
            source_subject.inputs.template = '*'
            field_template = dict(sphere_reg='surf/{0}.sphere.reg'.format(hemisphere))
            if threshold:
                for label in labels:
                    field_tempalte[label] = 'label/{0}.{1}_exvivo.thresh.label'.format(hemisphere, label)
                node_name = 'BA_Maps_' + hemisphere + '_Tresh'
            else:
                for label in labels:
                    field_tempalte[label] = 'label/{0}.{1}_exvivo.label'.format(hemisphere, label)
                node_name = 'BA_Maps_' + hemisphere
            source_subject.inputs.field_template = field_template

            merge = pe.Node(merge(len(labels)), name="Merge_SRCLabels")

            for i,label in enumerate(labels):
                ba_WF.connect([(source_subject, merge, [(label, 'in{0}'.format(i))])])

            node = pe.MapNode(
                Label2Label(), name=node_name, iterfield=['source_label'])
            node.inputs.hemisphere = hemisphere
            node.inputs.threshold = threshold

            ba_WF.connect([(merge, node, [('out', 'source_label')]),
                           (source_subject, node, [('sphere_reg', 'source_sphere_reg')])])
            
            label2annot = pe.Node(Label2Annot(), name=node_name + '_2_Annot')
            label2annot.inputs.hemisphere = hemisphere
            label2annot.inputs.color_table = os.path.join(
                config['FREESURFER_HOME'], 'average', 'colortable_BA.txt')
            label2annot.inputs.verbose_off = True
            label2annot.inputs.keep_max = True

            stats_node = pe.Node(
                ParcellationStats(), name=node_name + '_Stats')
            stats_node.inputs.hemisphere = hemisphere
            stats_node.inputs.mgz = True
            stats_node.inputs.surface = 'white'
            stats_node.inputs.tabular_output = True

            if threshold:
                label2annot.inputs.out_annot = "BA_exvivo.thresh"
                ba_WF.connect([(stats_node, outputSpec, [('out_color', 
                                                          '{0}_thresh_color'.format(hemisphere)),
                                                         ('out_table', 
                                                          '{0}_thresh_table'.format(hemisphere))])])
            else:
                label2annot.inputs.out_annot = "BA_exvivo"
                ba_WF.connect([(stats_node, outputSpec, [('out_color', 
                                                          '{0}_color'.format(hemisphere)),
                                                         ('out_table', 
                                                          '{0}_table'.format(hemisphere))])])


            ba_WF.connect([(inputSpec, node, [('{0}_sphere_reg'.format(hemisphere),
                                               'sphere_reg'),
                                              ('{0}_white'.format(
                                                  hemisphere), 'white'),
                                              ]),
                           (node, label2annot, [('out_file', 'in_labels')]),
                           (label2annot, stats_node,
                            [('out_file', 'in_annotation')]),
                           (inputSpec, stats_node, [('{0}_thickness'.format(hemisphere), 
                                                     'thickness'),
                                                    ('lh_white', 'lh_white'),
                                                    ('rh_white', 'rh_white'),
                                                    ('lh_pial', 'lh_pial'),
                                                    ('rh_pial', 'rh_pial'),
                                                    ('transform', 'transform'),
                                                    ('brainmask', 'brainmask'),
                                                    ('aseg', 'aseg'),
                                                    ('wm', 'wm'),
                                                    ('ribbon', 'ribbon')])])


    return ba_WF
