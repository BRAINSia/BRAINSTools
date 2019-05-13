"""
BAWantsRegistrationTemplateBuildTemplate.py
==============================================
Description:

Author:

Usage:

"""


from builtins import range
from builtins import str

import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe
from nipype.interfaces.ants import (
    Registration,
    ApplyTransforms,
    AverageImages,
    MultiplyImages,
    AverageAffineTransform,
)
from nipype.interfaces.utility import Function

from .utilities.distributed import modify_qsub_args
from .utilities.misc import common_ants_registration_settings


def make_list_of_one_element(inputFile):
    """
    This functions...

    :param inputFile:
    :return:
    """
    outputList = [inputFile]
    return outputList


def get_first_list_element(this_list):
    """
    This function..

    :param this_list:
    :return:
    """
    return this_list[0]


def make_transform_list_with_gradient_warps(averageAffineTranform, gradientStepWarp):
    """
    This function...

    :param averageAffineTranform:
    :param gradientStepWarp:
    :return:
    """
    return [
        averageAffineTranform,
        gradientStepWarp,
        gradientStepWarp,
        gradientStepWarp,
        gradientStepWarp,
    ]


def renest_deformed_passive_images(
    deformedPassiveImages, flattened_image_nametypes, interpolationMapping
):
    import os
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    """ Now make a list of lists of images where the outter list is per image type,
    and the inner list is the same size as the number of subjects to be averaged.
    In this case, the first element will be a list of all the deformed T2's, and
    the second element will be a list of all deformed POSTERIOR_AIR,  etc..

    :param deformedPassiveImages:
    :param flattened_image_nametypes:
    :param interpolationMapping:
    :return:
    """
    all_images_size = len(deformedPassiveImages)
    image_dictionary_of_lists = OrderedDict()
    nested_imagetype_list = list()
    outputAverageImageName_list = list()
    image_type_list = list()
    nested_interpolation_type = list()
    ## make empty_list, this is not efficient, but it works
    for name in flattened_image_nametypes:
        image_dictionary_of_lists[name] = list()
    for index in range(0, all_images_size):
        curr_name = flattened_image_nametypes[index]
        curr_file = deformedPassiveImages[index]
        image_dictionary_of_lists[curr_name].append(curr_file)
    for image_type, image_list in list(image_dictionary_of_lists.items()):
        nested_imagetype_list.append(image_list)
        outputAverageImageName_list.append("AVG_" + image_type + ".nii.gz")
        image_type_list.append("WARP_AVG_" + image_type)
        if image_type in interpolationMapping:
            nested_interpolation_type.append(interpolationMapping[image_type])
        else:
            nested_interpolation_type.append("Linear")  # Linear is the default.
    print(("\n" * 10))
    print(("HACK: ", nested_imagetype_list))
    print(("HACK: ", outputAverageImageName_list))
    print(("HACK: ", image_type_list))
    print(("HACK: ", nested_interpolation_type))
    return (
        nested_imagetype_list,
        outputAverageImageName_list,
        image_type_list,
        nested_interpolation_type,
    )


def split_composite_to_component_transform(transformFilename):
    """
    This function...

    :param transformFilename:
    :return:
    """
    ### Nota bene: The outputs will include the initial_moving_transform from Registration (which depends on what
    ###            the invert_initial_moving_transform is set to)
    import os
    import subprocess
    import sys

    if transformFilename.endswith(".h5"):
        decomposedOutputPrefix = "decomposedTransform"
        commandLine = (
            "CompositeTransformUtil  --disassemble "
            + transformFilename
            + " "
            + decomposedOutputPrefix
        )

        affineTransformName = "00_" + decomposedOutputPrefix + "_AffineTransform.mat"
        warpTransformName = (
            "01_" + decomposedOutputPrefix + "_DisplacementFieldTransform.nii.gz"
        )

        script_name = "decomposedTransform" + "_script.sh"
        print(script_name)
        script = open(script_name, "w")
        script.write("#!/bin/bash\n")
        script.write("\npushd " + " " + os.path.abspath(".") + "\n")
        script.write(commandLine)
        script.write("\npopd \n")
        script.close()
        os.chmod(script_name, 0o777)
        script_name = os.path.abspath(script_name)
        print(("XX" * 40))
        print(("Starting CompositeTransformUtil script: {0}".format(script_name)))
        scriptStatus = subprocess.check_call([script_name], shell=True)
        if scriptStatus != 0:
            sys.exit(scriptStatus)
        print("Ending CompositeTransformUtil")
        if os.path.exists(affineTransformName) and os.path.exists(warpTransformName):
            affine_component_list = os.path.abspath(affineTransformName)
            warp_component_list = os.path.abspath(warpTransformName)
        else:
            print("There is no decomposed trasforms generated")
            print(affineTransformName)
            print(warpTransformName)
            sys.exit(-1)
    else:
        print("There is no composite transform to split")
        print(transformFilename)
        sys.exit(-1)
    return affine_component_list, warp_component_list


## Flatten and return equal length transform and images lists.


def flatten_transform_and_images_list(
    ListOfPassiveImagesDictionaries,
    transforms,
    interpolationMapping,
    invert_transform_flags=None,
):
    """
    This function...

    :param ListOfPassiveImagesDictionaries:
    :param transforms:
    :param interpolationMapping:
    :param invert_transform_flags: None
    :return:
    """
    import sys

    print(
        (
            "HACK:  DEBUG: ListOfPassiveImagesDictionaries\n{lpi}\n".format(
                lpi=ListOfPassiveImagesDictionaries
            )
        )
    )
    subjCount = len(ListOfPassiveImagesDictionaries)
    tranCount = len(transforms)
    if subjCount != tranCount:
        print(
            (
                "ERROR:  subjCount must equal tranCount {0} != {1}".format(
                    subjCount, tranCount
                )
            )
        )
        sys.exit(-1)
    if invert_transform_flags is None:
        invert_transform_flags = [False] * subjCount
    invertTfmsFlagsCount = len(invert_transform_flags)
    if subjCount != invertTfmsFlagsCount:
        print(
            (
                "ERROR:  subjCount must equal invertTfmsFlags {0} != {1}".format(
                    subjCount, invertTfmsFlagsCount
                )
            )
        )
        sys.exit(-1)
    flattened_images = list()
    flattened_image_nametypes = list()
    flattened_transforms = list()
    flattened_invert_transform_flags = list()
    flattened_interpolation_type = list()
    passiveImagesCount = len(ListOfPassiveImagesDictionaries[0])
    for subjIndex in range(0, subjCount):
        # if passiveImagesCount != len(ListOfPassiveImagesDictionaries[subjIndex]):
        #    print "ERROR:  all image lengths must be equal {0} != {1}".format(passiveImagesCount,len(ListOfPassiveImagesDictionaries[subjIndex]))
        #    sys.exit(-1)
        subjImgDictionary = ListOfPassiveImagesDictionaries[subjIndex]
        subjToAtlasTransform = transforms[subjIndex]
        subjToAtlasInvertFlags = invert_transform_flags[subjIndex]
        for imgname, img in list(subjImgDictionary.items()):
            flattened_images.append(img)
            flattened_image_nametypes.append(imgname)
            flattened_transforms.append(subjToAtlasTransform)
            flattened_invert_transform_flags.append(subjToAtlasInvertFlags)
            if imgname in interpolationMapping:
                flattened_interpolation_type.append(interpolationMapping[imgname])
            else:
                flattened_interpolation_type.append("Linear")  # Linear is the default.
    print(("HACK: flattened images    {0}\n".format(flattened_images)))
    print(("HACK: flattened nametypes {0}\n".format(flattened_image_nametypes)))
    print(("HACK: flattened txfms     {0}\n".format(flattened_transforms)))
    print(("HACK: flattened txfmsFlags{0}\n".format(flattened_invert_transform_flags)))
    return (
        flattened_images,
        flattened_transforms,
        flattened_invert_transform_flags,
        flattened_image_nametypes,
        flattened_interpolation_type,
    )


def get_moving_images(
    ListOfImagesDictionaries, registrationImageTypes, interpolationMapping
):
    """ This currently ONLY works when registrationImageTypes has
        length of exactly 1.  When the new multi-variate registration
        is introduced, it will be expanded.

        :param ListOfImagesDictionaries:
        :param registrationImageTypes:
        :param interpolationMapping:
        :return:
    """
    if len(registrationImageTypes) != 1:
        print("ERROR:  Multivariate imageing not supported yet!")
        return []
    moving_images = [
        mdict[registrationImageTypes[0]] for mdict in ListOfImagesDictionaries
    ]
    moving_interpolation_type = interpolationMapping[registrationImageTypes[0]]
    return moving_images, moving_interpolation_type


def get_passive_images(ListOfImagesDictionaries, registrationImageTypes):
    """
    This function...

    :param ListOfImagesDictionaries:
    :param registrationImageTypes:
    :return:
    """
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    if len(registrationImageTypes) != 1:
        print("ERROR:  Multivariate imageing not supported yet!")
        return [OrderedDict()]
    passive_images = list()
    for mdict in ListOfImagesDictionaries:
        ThisSubjectPassiveImages = OrderedDict()
        for key, value in list(mdict.items()):
            if key not in registrationImageTypes:
                ThisSubjectPassiveImages[key] = value
        passive_images.append(ThisSubjectPassiveImages)
    return passive_images


##
## NOTE:  The modes can be either 'SINGLE_IMAGE' or 'MULTI'
##        'SINGLE_IMAGE' is quick shorthand when you are building an atlas with a single subject, then registration can
##                    be short-circuted
##        any other string indicates the normal mode that you would expect and replicates the shell script build_template_parallel.sh


def baw_ants_registration_template_build_single_iteration_wf(
    iterationPhasePrefix, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG
):
    """

    Inputs::

           inputspec.images :
           inputspec.fixed_image :
           inputspec.ListOfPassiveImagesDictionaries :
           inputspec.interpolationMapping :

    Outputs::

           outputspec.template :
           outputspec.transforms_list :
           outputspec.passive_deformed_templates :
    """
    TemplateBuildSingleIterationWF = pe.Workflow(
        name="antsRegistrationTemplateBuildSingleIterationWF_"
        + str(iterationPhasePrefix)
    )

    inputSpec = pe.Node(
        interface=util.IdentityInterface(
            fields=[
                "ListOfImagesDictionaries",
                "registrationImageTypes",
                # 'maskRegistrationImageType',
                "interpolationMapping",
                "fixed_image",
            ]
        ),
        run_without_submitting=True,
        name="inputspec",
    )
    ## HACK: TODO: We need to have the AVG_AIR.nii.gz be warped with a default voxel value of 1.0
    ## HACK: TODO: Need to move all local functions to a common untility file, or at the top of the file so that
    ##             they do not change due to re-indenting.  Otherwise re-indenting for flow control will trigger
    ##             their hash to change.
    ## HACK: TODO: REMOVE 'transforms_list' it is not used.  That will change all the hashes
    ## HACK: TODO: Need to run all python files through the code beutifiers.  It has gotten pretty ugly.
    outputSpec = pe.Node(
        interface=util.IdentityInterface(
            fields=["template", "transforms_list", "passive_deformed_templates"]
        ),
        run_without_submitting=True,
        name="outputspec",
    )

    ### NOTE MAP NODE! warp each of the original images to the provided fixed_image as the template
    BeginANTS = pe.MapNode(
        interface=Registration(), name="BeginANTS", iterfield=["moving_image"]
    )
    # SEE template.py many_cpu_BeginANTS_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,4,2,8), 'overwrite': True}
    ## This is set in the template.py file BeginANTS.plugin_args = BeginANTS_cpu_sge_options_dictionary
    common_ants_registration_settings(
        antsRegistrationNode=BeginANTS,
        registrationTypeDescription="SixStageAntsRegistrationT1Only",
        output_transform_prefix=str(iterationPhasePrefix) + "_tfm",
        output_warped_image="atlas2subject.nii.gz",
        output_inverse_warped_image="subject2atlas.nii.gz",
        save_state="SavedantsRegistrationNodeSyNState.h5",
        invert_initial_moving_transform=False,
        initial_moving_transform=None,
    )

    GetMovingImagesNode = pe.Node(
        interface=util.Function(
            function=get_moving_images,
            input_names=[
                "ListOfImagesDictionaries",
                "registrationImageTypes",
                "interpolationMapping",
            ],
            output_names=["moving_images", "moving_interpolation_type"],
        ),
        run_without_submitting=True,
        name="99_GetMovingImagesNode",
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec,
        "ListOfImagesDictionaries",
        GetMovingImagesNode,
        "ListOfImagesDictionaries",
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec,
        "registrationImageTypes",
        GetMovingImagesNode,
        "registrationImageTypes",
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec, "interpolationMapping", GetMovingImagesNode, "interpolationMapping"
    )

    TemplateBuildSingleIterationWF.connect(
        GetMovingImagesNode, "moving_images", BeginANTS, "moving_image"
    )
    TemplateBuildSingleIterationWF.connect(
        GetMovingImagesNode, "moving_interpolation_type", BeginANTS, "interpolation"
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec, "fixed_image", BeginANTS, "fixed_image"
    )

    ## Now warp all the input_images images
    wimtdeformed = pe.MapNode(
        interface=ApplyTransforms(),
        iterfield=["transforms", "input_image"],
        # iterfield=['transforms', 'invert_transform_flags', 'input_image'],
        name="wimtdeformed",
    )
    wimtdeformed.inputs.interpolation = "Linear"
    wimtdeformed.default_value = 0
    # HACK: Should try using forward_composite_transform
    ##PREVIOUS TemplateBuildSingleIterationWF.connect(BeginANTS, 'forward_transform', wimtdeformed, 'transforms')
    TemplateBuildSingleIterationWF.connect(
        BeginANTS, "composite_transform", wimtdeformed, "transforms"
    )
    ##PREVIOUS TemplateBuildSingleIterationWF.connect(BeginANTS, 'forward_invert_flags', wimtdeformed, 'invert_transform_flags')
    ## NOTE: forward_invert_flags:: List of flags corresponding to the forward transforms
    # wimtdeformed.inputs.invert_transform_flags = [False,False,False,False,False]
    TemplateBuildSingleIterationWF.connect(
        GetMovingImagesNode, "moving_images", wimtdeformed, "input_image"
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec, "fixed_image", wimtdeformed, "reference_image"
    )

    ##  Shape Update Next =====
    ## Now  Average All input_images deformed images together to create an updated template average
    AvgDeformedImages = pe.Node(interface=AverageImages(), name="AvgDeformedImages")
    AvgDeformedImages.inputs.dimension = 3
    AvgDeformedImages.inputs.output_average_image = (
        str(iterationPhasePrefix) + ".nii.gz"
    )
    AvgDeformedImages.inputs.normalize = True
    TemplateBuildSingleIterationWF.connect(
        wimtdeformed, "output_image", AvgDeformedImages, "images"
    )

    ## Now average all affine transforms together
    AvgAffineTransform = pe.Node(
        interface=AverageAffineTransform(), name="AvgAffineTransform"
    )
    AvgAffineTransform.inputs.dimension = 3
    AvgAffineTransform.inputs.output_affine_transform = (
        "Avererage_" + str(iterationPhasePrefix) + "_Affine.h5"
    )

    SplitCompositeTransform = pe.MapNode(
        interface=util.Function(
            function=split_composite_to_component_transform,
            input_names=["transformFilename"],
            output_names=["affine_component_list", "warp_component_list"],
        ),
        iterfield=["transformFilename"],
        run_without_submitting=True,
        name="99_SplitCompositeTransform",
    )
    TemplateBuildSingleIterationWF.connect(
        BeginANTS, "composite_transform", SplitCompositeTransform, "transformFilename"
    )
    ## PREVIOUS TemplateBuildSingleIterationWF.connect(BeginANTS, 'forward_transforms', SplitCompositeTransform, 'transformFilename')
    TemplateBuildSingleIterationWF.connect(
        SplitCompositeTransform,
        "affine_component_list",
        AvgAffineTransform,
        "transforms",
    )

    ## Now average the warp fields togther
    AvgWarpImages = pe.Node(interface=AverageImages(), name="AvgWarpImages")
    AvgWarpImages.inputs.dimension = 3
    AvgWarpImages.inputs.output_average_image = (
        str(iterationPhasePrefix) + "warp.nii.gz"
    )
    AvgWarpImages.inputs.normalize = True
    TemplateBuildSingleIterationWF.connect(
        SplitCompositeTransform, "warp_component_list", AvgWarpImages, "images"
    )

    ## Now average the images together
    ## TODO:  For now GradientStep is set to 0.25 as a hard coded default value.
    GradientStep = 0.25
    GradientStepWarpImage = pe.Node(
        interface=MultiplyImages(), name="GradientStepWarpImage"
    )
    GradientStepWarpImage.inputs.dimension = 3
    GradientStepWarpImage.inputs.second_input = -1.0 * GradientStep
    GradientStepWarpImage.inputs.output_product_image = (
        "GradientStep0.25_" + str(iterationPhasePrefix) + "_warp.nii.gz"
    )
    TemplateBuildSingleIterationWF.connect(
        AvgWarpImages, "output_average_image", GradientStepWarpImage, "first_input"
    )

    ## Now create the new template shape based on the average of all deformed images
    UpdateTemplateShape = pe.Node(
        interface=ApplyTransforms(), name="UpdateTemplateShape"
    )
    UpdateTemplateShape.inputs.invert_transform_flags = [True]
    UpdateTemplateShape.inputs.interpolation = "Linear"
    UpdateTemplateShape.default_value = 0

    TemplateBuildSingleIterationWF.connect(
        AvgDeformedImages,
        "output_average_image",
        UpdateTemplateShape,
        "reference_image",
    )
    TemplateBuildSingleIterationWF.connect(
        [
            (
                AvgAffineTransform,
                UpdateTemplateShape,
                [(("affine_transform", make_list_of_one_element), "transforms")],
            )
        ]
    )
    TemplateBuildSingleIterationWF.connect(
        GradientStepWarpImage,
        "output_product_image",
        UpdateTemplateShape,
        "input_image",
    )

    ApplyInvAverageAndFourTimesGradientStepWarpImage = pe.Node(
        interface=util.Function(
            function=make_transform_list_with_gradient_warps,
            input_names=["averageAffineTranform", "gradientStepWarp"],
            output_names=["TransformListWithGradientWarps"],
        ),
        run_without_submitting=True,
        name="99_MakeTransformListWithGradientWarps",
    )
    # ApplyInvAverageAndFourTimesGradientStepWarpImage.inputs.ignore_exception = True

    TemplateBuildSingleIterationWF.connect(
        AvgAffineTransform,
        "affine_transform",
        ApplyInvAverageAndFourTimesGradientStepWarpImage,
        "averageAffineTranform",
    )
    TemplateBuildSingleIterationWF.connect(
        UpdateTemplateShape,
        "output_image",
        ApplyInvAverageAndFourTimesGradientStepWarpImage,
        "gradientStepWarp",
    )

    ReshapeAverageImageWithShapeUpdate = pe.Node(
        interface=ApplyTransforms(), name="ReshapeAverageImageWithShapeUpdate"
    )
    ReshapeAverageImageWithShapeUpdate.inputs.invert_transform_flags = [
        True,
        False,
        False,
        False,
        False,
    ]
    ReshapeAverageImageWithShapeUpdate.inputs.interpolation = "Linear"
    ReshapeAverageImageWithShapeUpdate.default_value = 0
    ReshapeAverageImageWithShapeUpdate.inputs.output_image = (
        "ReshapeAverageImageWithShapeUpdate.nii.gz"
    )
    TemplateBuildSingleIterationWF.connect(
        AvgDeformedImages,
        "output_average_image",
        ReshapeAverageImageWithShapeUpdate,
        "input_image",
    )
    TemplateBuildSingleIterationWF.connect(
        AvgDeformedImages,
        "output_average_image",
        ReshapeAverageImageWithShapeUpdate,
        "reference_image",
    )
    TemplateBuildSingleIterationWF.connect(
        ApplyInvAverageAndFourTimesGradientStepWarpImage,
        "TransformListWithGradientWarps",
        ReshapeAverageImageWithShapeUpdate,
        "transforms",
    )
    TemplateBuildSingleIterationWF.connect(
        ReshapeAverageImageWithShapeUpdate, "output_image", outputSpec, "template"
    )

    ######
    ######
    ######  Process all the passive deformed images in a way similar to the main image used for registration
    ######
    ######
    ######
    ##############################################
    ## Now warp all the ListOfPassiveImagesDictionaries images
    FlattenTransformAndImagesListNode = pe.Node(
        Function(
            function=flatten_transform_and_images_list,
            input_names=[
                "ListOfPassiveImagesDictionaries",
                "transforms",
                "interpolationMapping",
                "invert_transform_flags",
            ],
            output_names=[
                "flattened_images",
                "flattened_transforms",
                "flattened_invert_transform_flags",
                "flattened_image_nametypes",
                "flattened_interpolation_type",
            ],
        ),
        run_without_submitting=True,
        name="99_FlattenTransformAndImagesList",
    )

    GetPassiveImagesNode = pe.Node(
        interface=util.Function(
            function=get_passive_images,
            input_names=["ListOfImagesDictionaries", "registrationImageTypes"],
            output_names=["ListOfPassiveImagesDictionaries"],
        ),
        run_without_submitting=True,
        name="99_GetPassiveImagesNode",
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec,
        "ListOfImagesDictionaries",
        GetPassiveImagesNode,
        "ListOfImagesDictionaries",
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec,
        "registrationImageTypes",
        GetPassiveImagesNode,
        "registrationImageTypes",
    )

    TemplateBuildSingleIterationWF.connect(
        GetPassiveImagesNode,
        "ListOfPassiveImagesDictionaries",
        FlattenTransformAndImagesListNode,
        "ListOfPassiveImagesDictionaries",
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec,
        "interpolationMapping",
        FlattenTransformAndImagesListNode,
        "interpolationMapping",
    )
    TemplateBuildSingleIterationWF.connect(
        BeginANTS,
        "composite_transform",
        FlattenTransformAndImagesListNode,
        "transforms",
    )
    ## FlattenTransformAndImagesListNode.inputs.invert_transform_flags = [False,False,False,False,False,False]
    ## TODO: Please check of invert_transform_flags has a fixed number.
    ## PREVIOUS TemplateBuildSingleIterationWF.connect(BeginANTS, 'forward_invert_flags', FlattenTransformAndImagesListNode, 'invert_transform_flags')
    wimtPassivedeformed = pe.MapNode(
        interface=ApplyTransforms(),
        iterfield=[
            "transforms",
            "invert_transform_flags",
            "input_image",
            "interpolation",
        ],
        name="wimtPassivedeformed",
    )
    wimtPassivedeformed.default_value = 0
    TemplateBuildSingleIterationWF.connect(
        AvgDeformedImages,
        "output_average_image",
        wimtPassivedeformed,
        "reference_image",
    )
    TemplateBuildSingleIterationWF.connect(
        FlattenTransformAndImagesListNode,
        "flattened_interpolation_type",
        wimtPassivedeformed,
        "interpolation",
    )
    TemplateBuildSingleIterationWF.connect(
        FlattenTransformAndImagesListNode,
        "flattened_images",
        wimtPassivedeformed,
        "input_image",
    )
    TemplateBuildSingleIterationWF.connect(
        FlattenTransformAndImagesListNode,
        "flattened_transforms",
        wimtPassivedeformed,
        "transforms",
    )
    TemplateBuildSingleIterationWF.connect(
        FlattenTransformAndImagesListNode,
        "flattened_invert_transform_flags",
        wimtPassivedeformed,
        "invert_transform_flags",
    )

    RenestDeformedPassiveImagesNode = pe.Node(
        Function(
            function=renest_deformed_passive_images,
            input_names=[
                "deformedPassiveImages",
                "flattened_image_nametypes",
                "interpolationMapping",
            ],
            output_names=[
                "nested_imagetype_list",
                "outputAverageImageName_list",
                "image_type_list",
                "nested_interpolation_type",
            ],
        ),
        run_without_submitting=True,
        name="99_RenestDeformedPassiveImages",
    )
    TemplateBuildSingleIterationWF.connect(
        inputSpec,
        "interpolationMapping",
        RenestDeformedPassiveImagesNode,
        "interpolationMapping",
    )
    TemplateBuildSingleIterationWF.connect(
        wimtPassivedeformed,
        "output_image",
        RenestDeformedPassiveImagesNode,
        "deformedPassiveImages",
    )
    TemplateBuildSingleIterationWF.connect(
        FlattenTransformAndImagesListNode,
        "flattened_image_nametypes",
        RenestDeformedPassiveImagesNode,
        "flattened_image_nametypes",
    )
    ## Now  Average All passive input_images deformed images together to create an updated template average
    AvgDeformedPassiveImages = pe.MapNode(
        interface=AverageImages(),
        iterfield=["images", "output_average_image"],
        name="AvgDeformedPassiveImages",
    )
    AvgDeformedPassiveImages.inputs.dimension = 3
    AvgDeformedPassiveImages.inputs.normalize = False
    TemplateBuildSingleIterationWF.connect(
        RenestDeformedPassiveImagesNode,
        "nested_imagetype_list",
        AvgDeformedPassiveImages,
        "images",
    )
    TemplateBuildSingleIterationWF.connect(
        RenestDeformedPassiveImagesNode,
        "outputAverageImageName_list",
        AvgDeformedPassiveImages,
        "output_average_image",
    )

    ## -- TODO:  Now neeed to reshape all the passive images as well
    ReshapeAveragePassiveImageWithShapeUpdate = pe.MapNode(
        interface=ApplyTransforms(),
        iterfield=["input_image", "reference_image", "output_image", "interpolation"],
        name="ReshapeAveragePassiveImageWithShapeUpdate",
    )
    ReshapeAveragePassiveImageWithShapeUpdate.inputs.invert_transform_flags = [
        True,
        False,
        False,
        False,
        False,
    ]
    ReshapeAveragePassiveImageWithShapeUpdate.default_value = 0
    TemplateBuildSingleIterationWF.connect(
        RenestDeformedPassiveImagesNode,
        "nested_interpolation_type",
        ReshapeAveragePassiveImageWithShapeUpdate,
        "interpolation",
    )
    TemplateBuildSingleIterationWF.connect(
        RenestDeformedPassiveImagesNode,
        "outputAverageImageName_list",
        ReshapeAveragePassiveImageWithShapeUpdate,
        "output_image",
    )
    TemplateBuildSingleIterationWF.connect(
        AvgDeformedPassiveImages,
        "output_average_image",
        ReshapeAveragePassiveImageWithShapeUpdate,
        "input_image",
    )
    TemplateBuildSingleIterationWF.connect(
        AvgDeformedPassiveImages,
        "output_average_image",
        ReshapeAveragePassiveImageWithShapeUpdate,
        "reference_image",
    )
    TemplateBuildSingleIterationWF.connect(
        ApplyInvAverageAndFourTimesGradientStepWarpImage,
        "TransformListWithGradientWarps",
        ReshapeAveragePassiveImageWithShapeUpdate,
        "transforms",
    )
    TemplateBuildSingleIterationWF.connect(
        ReshapeAveragePassiveImageWithShapeUpdate,
        "output_image",
        outputSpec,
        "passive_deformed_templates",
    )

    return TemplateBuildSingleIterationWF
