from nipype.interfaces.utility import IdentityInterface, Function
from nipype.pipeline import Node, Workflow
from nipype.interfaces.freesurfer import MRIConvert, MRIsConvert
from nipype.interfaces.semtools import BRAINSResample
from interfaces import *
from freesurfer_utils import SplitLabels, SurfaceMask, recode_labelmap, create_ones_image, MultiLabelDilation
import json


def get_local_file_location(relative_file_name):
    # credit: http://stackoverflow.com/questions/7165749/open-file-in-a-relative-location-in-python
    import os
    script_dir = os.path.abspath(os.path.dirname(__file__))  # <-- absolute dir the script is in
    return os.path.join(script_dir, relative_file_name)


def read_json_config(relative_file_name):
    # HACK: read the json file in this folder
    config_file_path = get_local_file_location(relative_file_name)
    with open(config_file_path, "rb") as config_file:
        config = json.load(config_file)
    return config


def create_logb_workflow(name="LOGISMOSB_WF", master_config=None, plugin_args=None):
    logb_wf = Workflow(name=name)

    config = read_json_config("config.json")
    config['atlas_info'] = get_local_file_location(config['atlas_info'])

    inputs_node = Node(
        IdentityInterface(
            fields=['t1_file',
                    't2_file',
                    'posterior_files',
                    'joint_fusion_file',
                    'brainlabels_file',
                    'hncma_atlas']), name="inputspec")
    inputs_node.run_without_submitting = True

    # ensure that t1 and t2 are in the same voxel lattice
    input_t2 = Node(BRAINSResample(), "ResampleInputT2Volume")
    input_t2.inputs.outputVolume = "t2_resampled.nii.gz"
    input_t2.inputs.pixelType = 'ushort'
    input_t2.inputs.interpolationMode = "Linear"

    logb_wf.connect([(inputs_node, input_t2, [('t1_file', 'referenceVolume'),
                                              ('t2_file', 'inputVolume')])])

    white_matter_masking_node = Node(interface=WMMasking(), name="WMMasking")
    white_matter_masking_node.inputs.dilation = config['WMMasking']['dilation']
    white_matter_masking_node.inputs.csf_threshold = config['WMMasking']['csf_threshold']
    if master_config and master_config['labelmap_colorlookup_table']:
        white_matter_masking_node.inputs.atlas_info = master_config['labelmap_colorlookup_table']
    else:
        white_matter_masking_node.inputs.atlas_info = config['atlas_info']

    logb_wf.connect([(inputs_node, white_matter_masking_node, [("posterior_files", "posterior_files"),
                                                               ("joint_fusion_file", "atlas_file"),
                                                               ("brainlabels_file", "brainlabels_file"),
                                                               ("hncma_atlas", "hncma_file")])])

    gm_labels = Node(interface=CreateGMLabelMap(), name="GM_Labelmap")
    gm_labels.inputs.atlas_info = config['atlas_info']
    logb_wf.connect([(inputs_node, gm_labels, [('joint_fusion_file', 'atlas_file')])])

    logismosb_output_node = create_output_spec(["wmsurface_file", "gmsurface_file"], config["hemisphere_names"],
                                               name="outputspec")

    for hemisphere in config["hemisphere_names"]:
        genus_zero_filter = Node(interface=GenusZeroImageFilter(), name="{0}_GenusZeroImageFilter".format(hemisphere))
        genus_zero_filter.inputs.connectivity = config['GenusZeroImageFilter']['connectivity']
        genus_zero_filter.inputs.biggestComponent = config['GenusZeroImageFilter']['biggestComponent']
        genus_zero_filter.inputs.connectedComponent = config['GenusZeroImageFilter']['connectedComponent']
        genus_zero_filter.inputs.out_mask = "{0}_genus_zero_white_matter.nii.gz".format(hemisphere)

        logb_wf.connect([(white_matter_masking_node, genus_zero_filter, [('{0}_wm'.format(hemisphere), 'in_file')])])

        surface_generation = Node(interface=BRAINSSurfaceGeneration(),
                                  name="{0}_BRAINSSurfaceGeneration".format(hemisphere))
        surface_generation.inputs.smoothSurface = config['BRAINSSurfaceGeneration']['smoothSurface']
        surface_generation.inputs.numIterations = config['BRAINSSurfaceGeneration']['numIterations']
        surface_generation.inputs.out_file = "{0}_white_matter_surface.vtk".format(hemisphere)

        logb_wf.connect([(genus_zero_filter, surface_generation, [('out_file', 'in_file')])])

        logismosb = Node(interface=LOGISMOSB(), name="{0}_LOGISMOSB".format(hemisphere))
        logismosb.inputs.smoothnessConstraint = config['LOGISMOSB']['smoothnessConstraint']
        logismosb.inputs.nColumns = config['LOGISMOSB']['nColumns']
        logismosb.inputs.columnChoice = config['LOGISMOSB']['columnChoice']
        logismosb.inputs.columnHeight = config['LOGISMOSB']['columnHeight']
        logismosb.inputs.nodeSpacing = config['LOGISMOSB']['nodeSpacing']
        logismosb.inputs.w = config['LOGISMOSB']['w']
        logismosb.inputs.a = config['LOGISMOSB']['a']
        logismosb.inputs.nPropagate = config['LOGISMOSB']['nPropagate']
        logismosb.inputs.basename = hemisphere
        if config['LOGISMOSB']['thickRegions']:
            logismosb.inputs.thick_regions = config['LOGISMOSB']['thickRegions']
        else:
            logismosb.inputs.useHNCMALabels = True

        if plugin_args:
            logismosb.plugin_args = plugin_args

        logb_wf.connect([(inputs_node, logismosb, [("t1_file", "t1_file"),
                                                   ('hncma_atlas', 'atlas_file')]),
                         (input_t2, logismosb, [("outputVolume", "t2_file")]),
                         (genus_zero_filter, logismosb, [("out_file", "wm_file")]),
                         (surface_generation, logismosb, [("out_file", "mesh_file")]),
                         (white_matter_masking_node, logismosb, [('{0}_boundary'.format(hemisphere),
                                                                  'brainlabels_file')]),
                         (logismosb, logismosb_output_node, [("gmsurface_file",
                                                              "{0}_gmsurface_file".format(hemisphere)),
                                                             ("wmsurface_file",
                                                              "{0}_wmsurface_file".format(hemisphere))])])

    return logb_wf


def create_output_spec(outputs, hemisphere_names, name):
    final_output_names = list()
    for output in outputs:
        for hemisphere in hemisphere_names:
            final_output_names.append("{0}_".format(hemisphere) + output)
    return Node(IdentityInterface(final_output_names), name)


def create_fs_compatible_logb_workflow(name="LOGISMOSB", plugin_args=None, config=None):
    """Create a workflow to run LOGISMOS-B from FreeSurfer Inputs"""

    if not config:
        config = read_json_config("fs_logb_config.json")

    wf = Workflow(name)

    inputspec = Node(IdentityInterface(['t1_file', 't2_file', 'white', 'aseg', 'hemi', 'recoding_file', 'gm_proba',
                                        'wm_proba', 'lut_file', 'hncma_atlas']), name="inputspec")

    # convert the white mesh to a vtk file with scanner coordinates
    to_vtk = Node(MRIsConvert(), name="WhiteVTK")
    to_vtk.inputs.out_file = "white.vtk"
    to_vtk.inputs.to_scanner = True

    wf.connect(inputspec, 'white', to_vtk, 'in_file')

    # convert brainslabels to nifti
    aseg_to_nifti = Node(MRIConvert(), "ABCtoNIFTI")
    aseg_to_nifti.inputs.out_file = "aseg.nii.gz"
    aseg_to_nifti.inputs.out_orientation = "LPS"
    wf.connect(inputspec, 'aseg', aseg_to_nifti, 'in_file')

    # create brainslabels from aseg
    aseg2brains = Node(Function(['in_file', 'recode_file', 'out_file'],
                                ['out_file'],
                                recode_labelmap), name="ConvertAseg2BRAINSLabels")
    aseg2brains.inputs.out_file = "brainslabels.nii.gz"

    wf.connect([(inputspec, aseg2brains, [('recoding_file', 'recode_file')]),
                (aseg_to_nifti, aseg2brains, [('out_file', 'in_file')])])

    t1_to_nifti = Node(MRIConvert(), "T1toNIFTI")
    t1_to_nifti.inputs.out_file = "t1.nii.gz"
    t1_to_nifti.inputs.out_orientation = "LPS"
    wf.connect(inputspec, 't1_file', t1_to_nifti, 'in_file')

    def t2_convert(in_file=None, reference_file=None, out_file=None):
        import os
        from nipype.interfaces.freesurfer import MRIConvert
        from nipype.interfaces.traits_extension import Undefined
        from nipype import Node
        if in_file:
            t2_to_nifti = Node(MRIConvert(), "T2toNIFTI")
            t2_to_nifti.inputs.in_file = in_file
            t2_to_nifti.inputs.out_file = os.path.abspath(out_file)
            t2_to_nifti.inputs.out_orientation = "LPS"
            if reference_file:
                t2_to_nifti.inputs.reslice_like = reference_file
            result = t2_to_nifti.run()
            out_file = os.path.abspath(result.outputs.out_file)
        else:
            out_file = Undefined
        return out_file

    t2_node = Node(Function(['in_file', 'reference_file', 'out_file'], ['out_file'], t2_convert), name="T2Convert")
    t2_node.inputs.out_file = "t2.nii.gz"
    wf.connect(inputspec, 't2_file', t2_node, 'in_file')
    wf.connect(t1_to_nifti, 'out_file', t2_node, 'reference_file')

    # convert raw t1 to lia
    t1_to_ras = Node(MRIConvert(), "T1toRAS")
    t1_to_ras.inputs.out_orientation = "LIA"
    t1_to_ras.inputs.out_file = "t1_lia.mgz"
    wf.connect(inputspec, 't1_file', t1_to_ras, 'in_file')

    # Create ones image for use when masking the white matter
    ones = Node(Function(['in_volume', 'out_file'],
                         ['out_file'],
                         create_ones_image),
                name="Ones_Image")
    ones.inputs.out_file = "ones.mgz"

    wf.connect(t1_to_ras, 'out_file', ones, 'in_volume')

    # use the ones image to obtain a white matter mask
    surfmask = Node(SurfaceMask(), name="WhiteMask")
    surfmask.inputs.out_file = "white_ras.mgz"

    wf.connect(ones, 'out_file', surfmask, 'in_volume')
    wf.connect(inputspec, 'white', surfmask, 'in_surface')

    surfmask_to_nifti = Node(MRIConvert(), "MasktoNIFTI")
    surfmask_to_nifti.inputs.out_file = "white.nii.gz"
    surfmask_to_nifti.inputs.out_orientation = "LPS"

    wf.connect(surfmask, 'out_file', surfmask_to_nifti, 'in_file')

    # create hemi masks

    split = Node(SplitLabels(), name="SplitLabelMask")
    split.inputs.out_file = "HemiBrainLabels.nii.gz"
    wf.connect([(aseg2brains, split, [('out_file', 'in_file')]),
                (inputspec, split, [('lut_file', 'lookup_table')]),
                (aseg_to_nifti, split, [('out_file', 'labels_file')]),
                (inputspec, split, [('hemi', 'hemi')])])

    dilate = Node(MultiLabelDilation(), "DilateLabels")
    dilate.inputs.out_file = "DilatedBrainLabels.nii.gz"
    dilate.inputs.radius = 1
    wf.connect(split, 'out_file', dilate, 'in_file')

    convert_label_map = Node(MRIConvert(), "ConvertLabelMapToMatchT1")
    convert_label_map.inputs.resample_type = "nearest"
    convert_label_map.inputs.out_file = "BrainLabelsFromAsegInT1Space.nii.gz"
    wf.connect(t1_to_nifti, 'out_file', convert_label_map, 'reslice_like')
    wf.connect(dilate, 'out_file', convert_label_map, 'in_file')

    logb = Node(LOGISMOSB(), name="LOGISMOS-B")
    logb.inputs.smoothnessConstraint = config['LOGISMOSB']['smoothnessConstraint']
    logb.inputs.nColumns = config['LOGISMOSB']['nColumns']
    logb.inputs.columnChoice = config['LOGISMOSB']['columnChoice']
    logb.inputs.columnHeight = config['LOGISMOSB']['columnHeight']
    logb.inputs.nodeSpacing = config['LOGISMOSB']['nodeSpacing']
    logb.inputs.w = config['LOGISMOSB']['w']
    logb.inputs.a = config['LOGISMOSB']['a']
    logb.inputs.nPropagate = config['LOGISMOSB']['nPropagate']

    if plugin_args:
        logb.plugin_args = plugin_args

    wf.connect([(t1_to_nifti, logb, [('out_file', 't1_file')]),
                (t2_node, logb, [('out_file', 't2_file')]),
                (inputspec, logb, [('hemi', 'basename'),
                                   ('hncma_atlas', 'atlas_file'),
                                   ('wm_proba', 'wm_proba_file'),
                                   ('gm_proba', 'gm_proba_file')]),
                (to_vtk, logb, [('converted', 'mesh_file')]),
                (surfmask_to_nifti, logb, [('out_file', 'wm_file')]),
                (convert_label_map, logb, [('out_file', 'brainlabels_file')])])

    outputspec = Node(IdentityInterface(['gmsurface_file', 'wmsurface_file']), name="outputspec")

    wf.connect([(logb, outputspec, [('gmsurface_file', 'gmsurface_file'),
                                    ('wmsurface_file', 'wmsurface_file')])])

    return wf


def create_fs_logb_workflow_for_both_hemispheres(name="FSLOGB", plugin_args=None, ml=False, config=None):
    """Creates a workflow that connects FreeSurfer with LOGISMOS-B"""

    fslogb_wf = Workflow(name=name)

    inputspec = Node(IdentityInterface(['recoding_file', 'lut_file', 'aseg_presurf', 'rawavg', 't2_raw', 'lh_white',
                                        'rh_white', 'hncma_atlas']),
                     name="inputspec")

    inputspec.inputs.recoding_file = get_local_file_location("abc_fs_equivelants.json")
    inputspec.inputs.lut_file = get_local_file_location("FreeSurferColorLUT.csv")

    # create outputspec with gm and wm surfaces
    outputs = ['lh_gm_surf_file', 'lh_wm_surf_file', 'rh_gm_surf_file', 'rh_wm_surf_file']
    outputspec = Node(IdentityInterface(outputs), name="outputspec")

    for hemi in ('lh', 'rh'):
        hemi_logb_wf = create_fs_compatible_logb_workflow("{0}_LOGBWF".format(hemi), plugin_args=plugin_args,
                                                          config=config)
        hemi_logb_wf.inputs.inputspec.hemi = hemi
        fslogb_wf.connect([(inputspec, hemi_logb_wf, [('aseg_presurf', 'inputspec.aseg'),
                                                      ('rawavg', 'inputspec.t1_file'),
                                                      ('t2_raw', 'inputspec.t2_file'),
                                                      ('hncma_atlas', 'inputspec.hncma_atlas'),
                                                      ('{0}_white'.format(hemi), 'inputspec.white')]),
                           (inputspec, hemi_logb_wf, [('recoding_file', 'inputspec.recoding_file'),
                                                      ('lut_file', 'inputspec.lut_file')])])

        # move the outputs from logb to the outputspec
        fslogb_wf.connect([(hemi_logb_wf, outputspec, [('outputspec.gmsurface_file', '{0}_gm_surf_file'.format(hemi)),
                                                       ('outputspec.wmsurface_file',
                                                        '{0}_wm_surf_file'.format(hemi))])])

    return fslogb_wf
