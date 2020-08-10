#! /usr/bin/env python
"""
template.py
================
This program is used to generate the subject- and session-specific workflows for BRAINSTool processing

Usage:
  template.py [--rewrite-datasinks] [--wfrun PLUGIN] [--use-sentinal] [--use-shuffle] [--dotfilename PFILE] --workphase WORKPHASE --pe ENV --ExperimentConfig FILE SUBJECTS...
  template.py -v | --version
  template.py -h | --help

Arguments:
  SUBJECTS              List of subject IDs to process

Options:
  -h, --help            Show this help and exit
  -v, --version         Print the version and exit
  --dotfilename=PFILE   Turn on printing pipeline to file PFILE
  --rewrite-datasinks   Turn on the Nipype option to overwrite all files in the 'results' directory
  --use-sentinal        Use the t1_average file as a marker to determine if session needs to be run
  --use-shuffle         Shuffle the subjects randomly to minimize multiple runs collision probability
  --pe=ENV              The processing environment to use from configuration file
  --wfrun=PLUGIN        The name of the workflow plugin option (default: 'local')
  --workphase WORKPHASE The type of processing to be done only VALID is ['subject-template-generation']
  --ExperimentConfig=FILE   The configuration file

Examples:
  $ template.py --pe OSX --ExperimentConfig my_baw.config all
  $ template.py --wfrun helium_all.q --pe OSX --ExperimentConfig my_baw.config 1058 1059
  $ template.py --rewrite-datasinks --pe OSX --ExperimentConfig my_baw.config 2001

"""

import glob
import os
import sys
import traceback
from builtins import range

from .baw_exp import open_subject_database


def get_processed_subjects(resultdir, input_subjects_list):
    """
  This function...

  :param resultdir:
  :param input_subjects_list:
  :return:
  """
    import glob

    required_files = [
        "AvgTemplateHeadregion.nii.gz",
        "AvgHncmaAtlas.nii.gz",
        "AvgLAccumbenProbabilityMap.nii.gz",
        "AvgLCaudateProbabilityMap.nii.gz",
        "AvgLGlobusProbabilityMap.nii.gz",
        "AvgLHippocampusProbabilityMap.nii.gz",
        "AvgLMKS.fcsv",
        "AvgLPutamenProbabilityMap.nii.gz",
        "AvgLThalamusProbabilityMap.nii.gz",
        "AvgPhi.nii.gz",
        "AvgRAccumbenProbabilityMap.nii.gz",
        "AvgRCaudateProbabilityMap.nii.gz",
        "AvgRGlobusProbabilityMap.nii.gz",
        "AvgRHippocampusProbabilityMap.nii.gz",
        "AvgRho.nii.gz",
        "AvgRPutamenProbabilityMap.nii.gz",
        "AvgRThalamus_ProbabilityMap.nii.gz",
        "AvgT1.nii.gz",
        "AvgT2.nii.gz",
        "AvgTemplateHeadregion.nii.gz",
        "AvgTemplateLeftHemisphere.nii.gz",
        "AvgTemplateNacLabels.nii.gz",
        "AvgTemplateRightHemisphere.nii.gz",
        "AvgTemplateVentricles.nii.gz",
        "AvgTemplateWMPM2Labels.nii.gz",
        "AvgTheta.nii.gz",
        "ClippedAvgAIR.nii.gz",
        "ClippedAvgBASAL.nii.gz",
        "ClippedAvgBRAINMASK.nii.gz",
        "ClippedAvgCRBLGM.nii.gz",
        "ClippedAvgCRBLWM.nii.gz",
        "ClippedAvgCSF.nii.gz",
        "ClippedAvgGLOBUS.nii.gz",
        "ClippedAvgHIPPOCAMPUS.nii.gz",
        "ClippedAvgNOTCSF.nii.gz",
        "ClippedAvgNOTGM.nii.gz",
        "ClippedAvgNOTVB.nii.gz",
        "ClippedAvgNOTWM.nii.gz",
        "ClippedAvgSURFGM.nii.gz",
        "ClippedAvgTHALAMUS.nii.gz",
        "ClippedAvgVB.nii.gz",
        "ClippedAvgWM.nii.gz",
    ]
    # resultdir/subject_dir/Atlas/AVG_T1.nii.gz
    for subject in input_subjects_list:
        sential_file_pattern = "*{0}/Atlas/Avg_template_rightHemisphere.nii.gz".format(
            subject
        )
        glob_search_pattern = os.path.join(resultdir, sential_file_pattern)
        processedSubjectsPaths = glob.glob(glob_search_pattern)
        processedSubjects = list()
        partial_done = list()
        processedDirectories = [os.path.dirname(s) for s in processedSubjectsPaths]
        for testDirectory in processedDirectories:
            all_files_exists = True
            for testFile in required_files:
                if not os.path.exists(os.path.join(testDirectory, testFile)):
                    all_files_exists = False
                    print(
                        (
                            "MISSING FILE: {0}: ".format(
                                os.path.join(testDirectory, testFile)
                            )
                        )
                    )
            if all_files_exists:
                processedSubjects.append(
                    os.path.basename(os.path.dirname(testDirectory))
                )
            else:
                partial_done.append(os.path.basename(os.path.dirname(testDirectory)))
        print("SKIPPING COMPLETED SUBJECTS: {0}".format(processedSubjects))
        print("Finishing Incomplete SUBJECTS: {0}".format(partial_done))
        return processedSubjects


def get_subjects_sessions_dictionary(
    input_subjects, cache, resultdir, prefix, dbfile, useSentinal, shuffle=False
):
    """
  This function...

  :param input_subjects:
  :param cache:
  :param resultdir:
  :param prefix:
  :param dbfile:
  :param useSentinal:
  :param shuffle:
  :return:
  """
    import random
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    _temp = open_subject_database(cache, ["all"], prefix, dbfile)
    if "all" in input_subjects:
        input_subjects = _temp.get_all_subjects()
    if useSentinal:
        print("=" * 80)
        print("Using Sentinal Files to Limit Jobs Run")
        _all_subjects = set(input_subjects)
        _processed_subjects = set(get_processed_subjects(resultdir, input_subjects))
        subjects = list(
            _all_subjects - _processed_subjects
        )  # NOTE - in set operation notation removes values
    else:
        subjects = input_subjects

    if shuffle:
        random.shuffle(subjects)  # randomly shuffle to get max cluster efficiency
    subject_sessions_dictionary = OrderedDict()
    for subject in subjects:
        subject_sessions_dictionary[subject] = _temp.get_sessions_from_subject(subject)
    return subjects, subject_sessions_dictionary


## Merge the different subjects together
def merge_by_extended_list_elements(
    t1s, t2s, pds, fls, labels, posteriors, passive_intensities, passive_masks
):
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    """
  *** NOTE:  All input lists MUST have the same number of elements (even if they are null) ***

  output = [{'T1':        os.path.join(mydatadir, '01_T1_half.nii.gz'),
             'INV_T1':    os.path.join(mydatadir, '01_T1_inv_half.nii.gz'),
             'LABEL_MAP': os.path.join(mydatadir, '01_T1_inv_half.nii.gz')
            },
            {'T1':        os.path.join(mydatadir, '02_T1_half.nii.gz'),
             'INV_T1':    os.path.join(mydatadir, '02_T1_inv_half.nii.gz'),
             'LABEL_MAP': os.path.join(mydatadir, '02_T1_inv_half.nii.gz')
            },
            {'T1':        os.path.join(mydatadir, '03_T1_half.nii.gz'),
             'INV_T1':    os.path.join(mydatadir, '03_T1_inv_half.nii.gz'),
             'LABEL_MAP': os.path.join(mydatadir, '03_T1_inv_half.nii.gz')
            }
           ]
  #          SUBJECT_01                    SUBJECT_02                        SUBJECT_03
  labels = ['brain_label_seg.nii.gz',      'brain_label_seg.nii.gz',          ...      ]
  pds    = [None,                          None,                              ...      ]
  t1s    = ['t1_average_BRAINSABC.nii.gz', 't1_average_BRAINSABC.nii.gz',     ...      ]
  t2s    = ['t2_average_BRAINSABC.nii.gz', 't2_average_BRAINSABC.nii.gz',     ...      ]

  :param t1s:
  :param t2s:
  :param pds:
  :param fls:
  :param labels:
  :param posteriors:
  :param passive_intensities:
  :param passive_masks:
  :return:
  """
    # print "t1s", t1s
    # print "t2s", t2s
    # print "pds", pds
    # print "fls", fls
    # print "labels", labels
    # print "$$$$$$$$$$$$$$$$$$$$$$$"
    # print "posteriors", posteriors
    ListOfImagesDictionaries = [
        OrderedDict() for i in t1s
    ]  # Initial list with empty dictionaries
    # HACK:  Need to make it so that AVG_AIR.nii.gz has a background value of 1
    registrationImageTypes = ["T1"]  # ['T1','T2'] someday.
    DefaultContinuousInterpolationType = (
        "Linear"
    )  # or 'LanczosWindowedSinc' ('Linear' for speed)
    interpolationMapping = {
        "T1": DefaultContinuousInterpolationType,
        "T2": DefaultContinuousInterpolationType,
        "PD": DefaultContinuousInterpolationType,
        "FL": DefaultContinuousInterpolationType,
        "BRAINMASK": "MultiLabel",
    }
    for subject_index in range(len(t1s)):
        if t1s[subject_index] is not None:
            ListOfImagesDictionaries[subject_index]["T1"] = t1s[subject_index]
        if isinstance(t2s, list) and t2s[subject_index] is not None:
            ListOfImagesDictionaries[subject_index]["T2"] = t2s[subject_index]
        if isinstance(pds, list) and pds[subject_index] is not None:
            ListOfImagesDictionaries[subject_index]["PD"] = pds[subject_index]
        if isinstance(fls, list) and fls[subject_index] is not None:
            ListOfImagesDictionaries[subject_index]["FL"] = fls[subject_index]
        if labels[subject_index] is not None:
            ListOfImagesDictionaries[subject_index]["BRAINMASK"] = labels[subject_index]
        print(ListOfImagesDictionaries[subject_index])
        for key, value in list(posteriors.items()):
            # print "key:", key, " -> value:", value
            ListOfImagesDictionaries[subject_index][key] = value[subject_index]
            interpolationMapping[key] = DefaultContinuousInterpolationType
        for key, value in list(passive_intensities.items()):
            # print "key:", key, " -> value:", value
            ListOfImagesDictionaries[subject_index][key] = value[subject_index]
            interpolationMapping[key] = DefaultContinuousInterpolationType
        for key, value in list(passive_masks.items()):
            # print "key:", key, " -> value:", value
            ListOfImagesDictionaries[subject_index][key] = value[subject_index]
            interpolationMapping[key] = "MultiLabel"

    # print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    # print "ListOfImagesDictionaries", ListOfImagesDictionaries
    # print "registrationImageTypes", registrationImageTypes
    # print "interpolationMapping", interpolationMapping
    return ListOfImagesDictionaries, registrationImageTypes, interpolationMapping


def xml_filename(subject):
    """
  This function...

  :param subject:
  :return:
  """
    return "AtlasDefinition_{0}.xml".format(subject)


def get_session_from_subject_dictionary(subject_session_dictionary, subject):
    """
  This function...

  :param subject_session_dictionary:
  :param subject:
  :return:
  """
    print("#" + subject + "#" * 80 + "\n")
    print(subject_session_dictionary[subject])
    if len(subject_session_dictionary[subject]) == 0:
        import sys

        print(subject_session_dictionary)
        print("ERROR:  No sessions for subject {0}".format(subject))
        sys.exit(-1)
    return subject_session_dictionary[subject]


def _template_runner(argv, environment, experiment, pipeline_options, cluster):
    """
  This function...

  :param argv:
  :param environment:
  :param experiment:
  :param pipeline_options:
  :param cluster:
  :return:
  """
    print("Getting subjects from database...")
    # subjects = argv["--subjects"].split(',')
    subjects, subjects_sessions_dictionary = get_subjects_sessions_dictionary(
        argv["SUBJECTS"],
        experiment["cachedir"],
        experiment["resultdir"],
        environment["prefix"],
        experiment["dbfile"],
        argv["--use-sentinal"],
        argv["--use-shuffle"],
    )  # Build database before parallel section

    # unused variable?
    # useSentinal = argv["--use-sentinal"]

    # Quick preliminary sanity check
    for thisSubject in subjects:
        if len(subjects_sessions_dictionary[thisSubject]) == 0:
            print(
                (
                    "ERROR: subject {0} has no sessions found.  Did you supply a valid subject id on the command line?".format(
                        thisSubject
                    )
                )
            )
            sys.exit(-1)

    for thisSubject in subjects:
        print("Processing atlas generation for this subject: {0}".format(thisSubject))
        print("=" * 80)
        print("Copying Atlas directory and determining appropriate Nipype options...")
        subj_pipeline_options = nipype_options(
            argv, pipeline_options, cluster, experiment, environment
        )  # Generate Nipype options
        print("Dispatching jobs to the system...")
        ######
        # Now start workflow construction
        ######
        # Set universal pipeline options
        nipype_config.update_config(subj_pipeline_options)

        ready_for_template_building = True
        for thisSession in subjects_sessions_dictionary[thisSubject]:
            path_test = os.path.join(
                experiment["previousresult"],
                "*/{0}/{1}/TissueClassify/t1_average_BRAINSABC.nii.gz".format(
                    thisSubject, thisSession
                ),
            )
            t1_file_result = glob.glob(path_test)
            if len(t1_file_result) != 1:
                print(
                    (
                        "Incorrect number of t1 images found for data grabber {0}".format(
                            t1_file_result
                        )
                    )
                )
                print("     at path {0}".format(path_test))
                ready_for_template_building = False
        if not ready_for_template_building:
            print("TEMPORARY SKIPPING:  Not ready to process {0}".format(thisSubject))
            continue

        base_output_directory = os.path.join(
            subj_pipeline_options["logging"]["log_directory"], thisSubject
        )
        template = pe.Workflow(name="SubjectAtlas_Template_" + thisSubject)
        template.base_dir = base_output_directory

        subjectNode = pe.Node(
            interface=IdentityInterface(fields=["subject"]),
            run_without_submitting=True,
            name="99_subjectIterator",
        )
        subjectNode.inputs.subject = thisSubject

        sessionsExtractorNode = pe.Node(
            Function(
                function=get_session_from_subject_dictionary,
                input_names=["subject_session_dictionary", "subject"],
                output_names=["sessions"],
            ),
            run_without_submitting=True,
            name="99_sessionsExtractor",
        )
        sessionsExtractorNode.inputs.subject_session_dictionary = (
            subjects_sessions_dictionary
        )

        baselineOptionalDG = pe.MapNode(
            nio.DataGrabber(
                infields=["subject", "session"],
                outfields=["t2_average", "pd_average", "fl_average"],
                run_without_submitting=True,
            ),
            run_without_submitting=True,
            iterfield=["session"],
            name="BaselineOptional_DG",
        )

        baselineOptionalDG.inputs.base_directory = experiment["previousresult"]
        baselineOptionalDG.inputs.sort_filelist = True
        baselineOptionalDG.inputs.raise_on_empty = False
        baselineOptionalDG.inputs.template = "*"

        baselineOptionalDG.inputs.field_template = {
            "t2_average": "*/%s/%s/TissueClassify/t2_average_BRAINSABC.nii.gz",
            "pd_average": "*/%s/%s/TissueClassify/pd_average_BRAINSABC.nii.gz",
            "fl_average": "*/%s/%s/TissueClassify/fl_average_BRAINSABC.nii.gz",
        }
        baselineOptionalDG.inputs.template_args = {
            "t2_average": [["subject", "session"]],
            "pd_average": [["subject", "session"]],
            "fl_average": [["subject", "session"]],
        }

        baselineRequiredDG = pe.MapNode(
            nio.DataGrabber(
                infields=["subject", "session"],
                outfields=[
                    "t1_average",
                    "brainMaskLabels",
                    "posteriorImages",
                    "passive_intensities",
                    "passive_masks",
                    "BCD_ACPC_Landmarks_fcsv",
                ],
                run_without_submitting=True,
            ),
            run_without_submitting=True,
            iterfield=["session"],
            name="Baseline_DG",
        )

        baselineRequiredDG.inputs.base_directory = experiment["previousresult"]
        baselineRequiredDG.inputs.sort_filelist = True
        baselineRequiredDG.inputs.raise_on_empty = True
        baselineRequiredDG.inputs.template = "*"
        posterior_files = [
            "AIR",
            "BASAL",
            "CRBLGM",
            "CRBLWM",
            "CSF",
            "GLOBUS",
            "HIPPOCAMPUS",
            "NOTCSF",
            "NOTGM",
            "NOTVB",
            "NOTWM",
            "SURFGM",
            "THALAMUS",
            "VB",
            "WM",
        ]
        passive_intensities_files = [
            "rho.nii.gz",
            "phi.nii.gz",
            "theta.nii.gz",
            "lThalamusProbabilityMap.nii.gz",
            "rAccumbenProbabilityMap.nii.gz",
            "lGlobusProbabilityMap.nii.gz",
            "lAccumbenProbabilityMap.nii.gz",
            "lCaudateProbabilityMap.nii.gz",
            "lPutamenProbabilityMap.nii.gz",
            "rThalamusProbabilityMap.nii.gz",
            "rPutamenProbabilityMap.nii.gz",
            "rCaudateProbabilityMap.nii.gz",
            "rHippocampusProbabilityMap.nii.gz",
            "rGlobusProbabilityMap.nii.gz",
            "lHippocampusProbabilityMap.nii.gz",
        ]
        passive_mask_files = [
            "templateWMPM2Labels.nii.gz",
            "hncmaAtlas.nii.gz",
            "templateNacLabels.nii.gz",
            "templateLeftHemisphere.nii.gz",
            "templateRightHemisphere.nii.gz",
            "templateVentricles.nii.gz",
            "templateHeadregion.nii.gz",
        ]

        baselineRequiredDG.inputs.field_template = {
            "t1_average": "*/%s/%s/TissueClassify/t1AverageBRAINSABC.nii.gz",
            "brainMaskLabels": "*/%s/%s/TissueClassify/completeBrainlabelsSeg.nii.gz",
            "BCD_ACPC_Landmarks_fcsv": "*/%s/%s/ACPCAlign/BCDACPCLandmarks.fcsv",
            "posteriorImages": "*/%s/%s/TissueClassify/POSTERIOR%s.nii.gz",
            "passive_intensities": "*/%s/%s/WarpedAtlas2Subject/%s",
            "passive_masks": "*/%s/%s/WarpedAtlas2Subject/%s",
        }
        baselineRequiredDG.inputs.template_args = {
            "t1_average": [["subject", "session"]],
            "brainMaskLabels": [["subject", "session"]],
            "BCD_ACPC_Landmarks_fcsv": [["subject", "session"]],
            "posteriorImages": [["subject", "session", posterior_files]],
            "passive_intensities": [["subject", "session", passive_intensities_files]],
            "passive_masks": [["subject", "session", passive_mask_files]],
        }

        MergeByExtendListElementsNode = pe.Node(
            Function(
                function=merge_by_extended_list_elements,
                input_names=[
                    "t1s",
                    "t2s",
                    "pds",
                    "fls",
                    "labels",
                    "posteriors",
                    "passive_intensities",
                    "passive_masks",
                ],
                output_names=[
                    "ListOfImagesDictionaries",
                    "registrationImageTypes",
                    "interpolationMapping",
                ],
            ),
            run_without_submitting=True,
            name="99_MergeByExtendListElements",
        )

        template.connect(
            [
                (subjectNode, baselineRequiredDG, [("subject", "subject")]),
                (subjectNode, baselineOptionalDG, [("subject", "subject")]),
                (subjectNode, sessionsExtractorNode, [("subject", "subject")]),
                (sessionsExtractorNode, baselineRequiredDG, [("sessions", "session")]),
                (sessionsExtractorNode, baselineOptionalDG, [("sessions", "session")]),
                (
                    baselineRequiredDG,
                    MergeByExtendListElementsNode,
                    [
                        ("t1_average", "t1s"),
                        ("brainMaskLabels", "labels"),
                        (
                            (
                                "posteriorImages",
                                convert_sessions_list_of_posterior_list_to_dictionary_of_session_lists,
                            ),
                            "posteriors",
                        ),
                    ],
                ),
                (
                    baselineOptionalDG,
                    MergeByExtendListElementsNode,
                    [
                        ("t2_average", "t2s"),
                        ("pd_average", "pds"),
                        ("fl_average", "fls"),
                    ],
                ),
                (
                    baselineRequiredDG,
                    MergeByExtendListElementsNode,
                    [
                        (
                            (
                                "passive_intensities",
                                convert_sessions_list_of_posterior_list_to_dictionary_of_session_lists,
                            ),
                            "passive_intensities",
                        )
                    ],
                ),
                (
                    baselineRequiredDG,
                    MergeByExtendListElementsNode,
                    [
                        (
                            (
                                "passive_masks",
                                convert_sessions_list_of_posterior_list_to_dictionary_of_session_lists,
                            ),
                            "passive_masks",
                        )
                    ],
                ),
            ]
        )

        myInitAvgWF = pe.Node(
            interface=ants.AverageImages(), name="Atlas_antsSimpleAverage"
        )  # was 'Phase1_antsSimpleAverage'
        myInitAvgWF.inputs.dimension = 3
        myInitAvgWF.inputs.normalize = True
        myInitAvgWF.inputs.num_threads = -1
        template.connect(baselineRequiredDG, "t1_average", myInitAvgWF, "images")
        ####################################################################################################
        # TEMPLATE_BUILD_RUN_MODE = 'MULTI_IMAGE'
        # if numSessions == 1:
        #     TEMPLATE_BUILD_RUN_MODE = 'SINGLE_IMAGE'
        ####################################################################################################
        CLUSTER_QUEUE = cluster["queue"]
        CLUSTER_QUEUE_LONG = cluster["long_q"]
        buildTemplateIteration1 = baw_ants_registration_template_build_single_iteration_wf(
            "iteration01", CLUSTER_QUEUE, CLUSTER_QUEUE_LONG
        )
        # buildTemplateIteration2 = buildTemplateIteration1.clone(name='buildTemplateIteration2')
        buildTemplateIteration2 = baw_ants_registration_template_build_single_iteration_wf(
            "Iteration02", CLUSTER_QUEUE, CLUSTER_QUEUE_LONG
        )

        CreateAtlasXMLAndCleanedDeformedAveragesNode = pe.Node(
            interface=Function(
                function=create_atlas_xml_and_cleaned_deformed_averages,
                input_names=[
                    "t1_image",
                    "deformed_list",
                    "AtlasTemplate",
                    "outDefinition",
                ],
                output_names=["outAtlasFullPath", "clean_deformed_list"],
            ),
            # This is a lot of work, so submit it run_without_submitting=True,
            run_without_submitting=True,  # HACK:  THIS NODE REALLY SHOULD RUN ON THE CLUSTER!
            name="99_CreateAtlasXMLAndCleanedDeformedAverages",
        )

        if subj_pipeline_options["plugin_name"].startswith(
            "SGE"
        ):  # for some nodes, the qsub call needs to be modified on the cluster

            CreateAtlasXMLAndCleanedDeformedAveragesNode.plugin_args = {
                "template": subj_pipeline_options["plugin_args"]["template"],
                "qsub_args": modify_qsub_args(cluster["queue"], 1, 1, 1),
                "overwrite": True,
            }
            for bt in [buildTemplateIteration1, buildTemplateIteration2]:
                BeginANTS = bt.get_node("BeginANTS")
                BeginANTS.plugin_args = {
                    "template": subj_pipeline_options["plugin_args"]["template"],
                    "overwrite": True,
                    "qsub_args": modify_qsub_args(cluster["queue"], 7, 4, 16),
                }
                wimtdeformed = bt.get_node("wimtdeformed")
                wimtdeformed.plugin_args = {
                    "template": subj_pipeline_options["plugin_args"]["template"],
                    "overwrite": True,
                    "qsub_args": modify_qsub_args(cluster["queue"], 2, 2, 2),
                }

                # AvgAffineTransform = bt.get_node("AvgAffineTransform")
                # AvgAffineTransform.plugin_args = {'template': subj_pipeline_options['plugin_args']['template'], 'overwrite': True,
                #                                  'qsub_args': modify_qsub_args(cluster['queue'], 2, 1, 1)}

                wimtPassivedeformed = bt.get_node("wimtPassivedeformed")
                wimtPassivedeformed.plugin_args = {
                    "template": subj_pipeline_options["plugin_args"]["template"],
                    "overwrite": True,
                    "qsub_args": modify_qsub_args(cluster["queue"], 2, 2, 4),
                }

        # Running off previous baseline experiment
        NACCommonAtlas = make_atlas_node(
            experiment["atlascache"],
            "NACCommonAtlas_{0}".format("subject"),
            ["S_BRAINSABCSupport"],
        )  # HACK : replace 'subject' with subject id once this is a loop rather than an iterable.
        template.connect(
            [
                (
                    myInitAvgWF,
                    buildTemplateIteration1,
                    [("output_average_image", "inputspec.fixed_image")],
                ),
                (
                    MergeByExtendListElementsNode,
                    buildTemplateIteration1,
                    [
                        (
                            "ListOfImagesDictionaries",
                            "inputspec.ListOfImagesDictionaries",
                        ),
                        ("registrationImageTypes", "inputspec.registrationImageTypes"),
                        ("interpolationMapping", "inputspec.interpolationMapping"),
                    ],
                ),
                (
                    buildTemplateIteration1,
                    buildTemplateIteration2,
                    [("outputspec.template", "inputspec.fixed_image")],
                ),
                (
                    MergeByExtendListElementsNode,
                    buildTemplateIteration2,
                    [
                        (
                            "ListOfImagesDictionaries",
                            "inputspec.ListOfImagesDictionaries",
                        ),
                        ("registrationImageTypes", "inputspec.registrationImageTypes"),
                        ("interpolationMapping", "inputspec.interpolationMapping"),
                    ],
                ),
                (
                    subjectNode,
                    CreateAtlasXMLAndCleanedDeformedAveragesNode,
                    [(("subject", xml_filename), "outDefinition")],
                ),
                (
                    NACCommonAtlas,
                    CreateAtlasXMLAndCleanedDeformedAveragesNode,
                    [("ExtendedAtlasDefinition_xml_in", "AtlasTemplate")],
                ),
                (
                    buildTemplateIteration2,
                    CreateAtlasXMLAndCleanedDeformedAveragesNode,
                    [
                        ("outputspec.template", "t1_image"),
                        ("outputspec.passive_deformed_templates", "deformed_list"),
                    ],
                ),
            ]
        )

        ## Genearate an average lmks file.
        myAverageLmk = pe.Node(interface=GenerateAverageLmkFile(), name="myAverageLmk")
        myAverageLmk.inputs.outputLandmarkFile = "AVG_LMKS.fcsv"
        template.connect(
            baselineRequiredDG,
            "BCD_ACPC_Landmarks_fcsv",
            myAverageLmk,
            "inputLandmarkFiles",
        )

        # Create DataSinks
        SubjectAtlas_DataSink = pe.Node(nio.DataSink(), name="Subject_DS")
        SubjectAtlas_DataSink.overwrite = subj_pipeline_options["ds_overwrite"]
        SubjectAtlas_DataSink.inputs.base_directory = experiment["resultdir"]

        template.connect(
            [
                (subjectNode, SubjectAtlas_DataSink, [("subject", "container")]),
                (
                    CreateAtlasXMLAndCleanedDeformedAveragesNode,
                    SubjectAtlas_DataSink,
                    [("outAtlasFullPath", "Atlas.@definitions")],
                ),
                (
                    CreateAtlasXMLAndCleanedDeformedAveragesNode,
                    SubjectAtlas_DataSink,
                    [("clean_deformed_list", "Atlas.@passive_deformed_templates")],
                ),
                (
                    subjectNode,
                    SubjectAtlas_DataSink,
                    [(("subject", outputPattern), "regexp_substitutions")],
                ),
                (
                    buildTemplateIteration2,
                    SubjectAtlas_DataSink,
                    [("outputspec.template", "Atlas.@template")],
                ),
                (
                    myAverageLmk,
                    SubjectAtlas_DataSink,
                    [("outputLandmarkFile", "Atlas.@outputLandmarkFile")],
                ),
            ]
        )

        dotfilename = argv["--dotfilename"]
        if dotfilename is not None:
            print("WARNING: Printing workflow, but not running pipeline")
            print_workflow(
                template,
                plugin=subj_pipeline_options["plugin_name"],
                dotfilename=dotfilename,
            )
        else:
            run_workflow(
                template,
                plugin=subj_pipeline_options["plugin_name"],
                plugin_args=subj_pipeline_options["plugin_args"],
            )


if __name__ == "__main__":
    from .AutoWorkup import setup_environment

    from docopt import docopt

    argv = docopt(__doc__, version="1.1")
    print(argv)
    if argv["--workphase"] != "subject-template-generation":
        print(
            "ERROR: Only --workphase subject-template-generation supported for template building"
        )
        sys.exit(-1)
    print("=" * 100)
    environment, experiment, pipeline, cluster = setup_environment(argv)
    from nipype import config as nipype_config
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from nipype.interfaces.utility import IdentityInterface, Function
    import nipype.interfaces.ants as ants

    from .PipeLineFunctionHelpers import (
        convert_sessions_list_of_posterior_list_to_dictionary_of_session_lists,
    )
    from .workflows.atlasNode import (
        make_atlas_node,
        create_atlas_xml_and_cleaned_deformed_averages,
    )
    from .utilities.misc import generate_subject_output_pattern as outputPattern
    from .utilities.distributed import modify_qsub_args
    from .workflows.utils import run_workflow, print_workflow
    from .BAWantsRegistrationBuildTemplate import (
        baw_ants_registration_template_build_single_iteration_wf,
    )
    from .utilities.configFileParser import nipype_options
    from nipype.interfaces.semtools.testing.generateaveragelmkfile import (
        GenerateAverageLmkFile,
    )

    exit_status = _template_runner(argv, environment, experiment, pipeline, cluster)
    sys.exit(exit_status)
