"""
ba_maps.py
=================
Description:

Author:

Usage:

"""

import os
from nipype.interfaces.utility import IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import *


def create_ba_maps_wf(config):
    """
    This function...

    :param config:
    :return:
    """
    # Brodmann Area Maps (BA Maps) and Hinds V1 Atlas
    inputSpec = pe.Node(
        IdentityInterface(
            fields=[
                "lh_sphere_reg",
                "rh_sphere_reg",
                "lh_white",
                "rh_white",
                "lh_pial",
                "rh_pial",
                "transform",
                "lh_thickness",
                "rh_thickness",
                "brainmask",
                "aseg",
                "ribbon",
                "wm",
                "subject_id",
                "subjects_dir",
            ]
        ),
        name="Inputs",
    )

    outputSpec = pe.Node(
        IdentityInterface(
            fields=[
                "lh_table",
                "lh_color",
                "lh_thresh_table",
                "lh_thresh_color",
                "rh_table",
                "rh_color",
                "rh_thresh_table",
                "rh_thresh_color",
            ]
        ),
        name="Outputs",
    )

    ba_WF = pe.Workflow(name="Brodmann_Area_Maps")
    labels = [
        "BA1",
        "BA2",
        "BA3a",
        "BA3b",
        "BA4a",
        "BA4p",
        "BA6",
        "BA44",
        "BA45",
        "V1",
        "V2",
        "MT",
        "entorhinal",
        "perirhinal",
    ]
    for hemisphere in ["lh", "rh"]:
        for threshold in [True, False]:
            if threshold:
                node_name = "BA_Maps_" + hemisphere + "_Tresh"
            else:
                node_name = "BA_Maps_" + hemisphere

            node = pe.MapNode(Label2Label(), name=node_name, iterfield=["label"])
            node.inputs.hemisphere = hemisphere
            node.inputs.label = labels
            node.inputs.threshold = threshold

            label2annot = pe.Node(Label2Annot(), name=node_name + "_2_Annot")
            label2annot.inputs.hemisphere = hemisphere
            label2annot.inputs.color_table = os.path.join(
                config["FREESURFER_HOME"], "average", "colortable_BA.txt"
            )
            label2annot.inputs.verbose_off = True
            label2annot.inputs.keep_max = True

            stats_node = pe.Node(ParcellationStats(), name=node_name + "_Stats")
            stats_node.inputs.hemisphere = hemisphere
            stats_node.inputs.mgz = True
            stats_node.inputs.surface = "white"
            stats_node.inputs.tabular_output = True

            if threshold:
                label2annot.inputs.out_annot = "BA_exvivo.thresh"
                ba_WF.connect(
                    [
                        (
                            stats_node,
                            outputSpec,
                            [
                                ("out_color", f"{hemisphere}_thresh_color"),
                                ("out_table", f"{hemisphere}_thresh_table"),
                            ],
                        )
                    ]
                )
            else:
                label2annot.inputs.out_annot = "BA_exvivo"
                ba_WF.connect(
                    [
                        (
                            stats_node,
                            outputSpec,
                            [
                                ("out_color", f"{hemisphere}_color"),
                                ("out_table", f"{hemisphere}_table"),
                            ],
                        )
                    ]
                )

            ba_WF.connect(
                [
                    (
                        inputSpec,
                        node,
                        [
                            (f"{hemisphere}_sphere_reg", "sphere_reg"),
                            (f"{hemisphere}_white", "white"),
                            ("subject_id", "subject_id"),
                            ("subjects_dir", "subjects_dir"),
                        ],
                    ),
                    (node, label2annot, [("out_file", "in_labels")]),
                    (
                        inputSpec,
                        label2annot,
                        [
                            ("subject_id", "subject_id"),
                            ("subjects_dir", "subjects_dir"),
                        ],
                    ),
                    (label2annot, stats_node, [("out_file", "in_annotation")]),
                    (
                        inputSpec,
                        stats_node,
                        [
                            (f"{hemisphere}_thickness", "thickness"),
                            ("subject_id", "subject_id"),
                            ("subjects_dir", "subjects_dir"),
                            ("lh_white", "lh_white"),
                            ("rh_white", "rh_white"),
                            ("lh_pial", "lh_pial"),
                            ("rh_pial", "rh_pial"),
                            ("transform", "transform"),
                            ("brainmask", "brainmask"),
                            ("aseg", "aseg"),
                            ("wm", "wm"),
                            ("ribbon", "ribbon"),
                        ],
                    ),
                ]
            )

    return ba_WF
