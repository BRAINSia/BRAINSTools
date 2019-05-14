def image_autounwrap(wrapped_inputfn, unwrapped_outputbasefn):
    """ Find optimal image roll in each direction
    to roll the image with circular boundaries such
    that the resulting head is not split across the
    image boundaries

    :param wrapped_inputfn:
    :param unwrapped_outputfn:
    :return:
    """
    import SimpleITK as sitk
    import numpy as np
    from scipy.signal import savgol_filter

    def flip_permute_to_identity(sitkImageIn):
        """
        This function...

        :param sitkImageIn:
        :return:
        """
        dc = np.array(sitkImageIn.GetDirection())
        dc = dc.reshape(3, 3)
        permute_values = [7, 7, 7]
        for i in range(0, 3):
            permute_values[i] = np.argmax(np.abs(dc[i, :]))
        permuted_image = sitk.PermuteAxes(sitkImageIn, [int(x) for x in permute_values])

        dc = np.array(permuted_image.GetDirection())
        dc = dc.reshape(3, 3)
        flip_values = [False, False, False]
        for i in range(0, 3):
            if dc[i, i] < 0:
                flip_values[i] = True
        flipped_permuted_image = sitk.Flip(permuted_image, flip_values)

        return flipped_permuted_image

    # ensure that normal strings are used here
    # via typecasting.  ReadImage requires types
    # to be strings
    wrapped_inputfn = [str(ii) for ii in wrapped_inputfn]
    unwrapped_outputbasefn = [str(ii) for ii in unwrapped_outputbasefn]

    def one_axis_unwrap(wrapped_image, axis):
        """
        This function...

        :param wrapped_image:
        :param axis:
        :return:
        """
        slice_values = list()
        sitkAxis = wrapped_image.GetDimension() - 1 - axis

        last_slice = wrapped_image.GetSize()[sitkAxis]
        mask = 1.0 - sitk.OtsuThreshold(wrapped_image)
        mask = sitk.BinaryClosingByReconstruction(mask, 6)  ## Fill some small holes

        image_as_np = sitk.GetArrayFromImage(
            wrapped_image * sitk.Cast(mask, wrapped_image.GetPixelIDValue())
        )
        for ii in range(0, last_slice):
            next_index = (ii + 1) % last_slice
            if axis == 0:
                curr_slice = image_as_np[ii, :, :].flatten()
                next_slice = image_as_np[next_index, :, :].flatten()
            elif axis == 1:
                curr_slice = image_as_np[:, ii, :].flatten()
                next_slice = image_as_np[:, next_index, :].flatten()
            elif axis == 2:
                curr_slice = image_as_np[:, :, ii].flatten()
                next_slice = image_as_np[:, :, next_index].flatten()
            else:
                curr_slice = 0
                next_slice = 0
                metric_value = 0
                print("FATAL ERROR")
            diff = curr_slice - next_slice
            diff = diff * diff
            metric_value = np.sum(diff)
            if ii == 0:
                ref_slice_limit = 5 * metric_value
            if metric_value > ref_slice_limit:
                metric_value = ref_slice_limit
            slice_values.append(metric_value)
        del image_as_np
        ## Call smoothing function to remove small noise
        # return slice_values,slice_values
        window_length = 3  # 2*(208/2)+1
        polyorder = 1
        slice_values = savgol_filter(
            np.array(slice_values), window_length, polyorder, deriv=1, mode="wrap"
        )

        min_slice = np.argmax(slice_values)

        axis_max = wrapped_image.GetSize()[sitkAxis] - 1
        if min_slice > axis_max / 2:
            zRoll = min_slice - axis_max
        else:
            zRoll = min_slice
        orig_image_as_np = sitk.GetArrayFromImage(wrapped_image)
        unwrapped_image_as_np = np.roll(orig_image_as_np, zRoll, axis)
        outim = sitk.GetImageFromArray(unwrapped_image_as_np)
        outim.CopyInformation(wrapped_image)
        return outim, zRoll, slice_values

    unwrapped_outputfn = []
    for index in range(0, len(wrapped_inputfn)):
        ii = wrapped_inputfn[index]
        wrapped_image = sitk.ReadImage(str(ii))
        identdc_wrapped_image = flip_permute_to_identity(wrapped_image)
        del wrapped_image
        if 0 == 1:  # THIS DOES NOT WORK ROBUSTLY YET
            unwrapped_image, rotationZ, zslicevalues = one_axis_unwrap(
                identdc_wrapped_image, 0
            )
            unwrapped_image, rotationY, yslicevalues = one_axis_unwrap(
                unwrapped_image, 1
            )
            unwrapped_image, rotationX, xslicevalues = one_axis_unwrap(
                unwrapped_image, 2
            )

            new_origin = identdc_wrapped_image.TransformContinuousIndexToPhysicalPoint(
                (-rotationX, -rotationY, -rotationZ)
            )
            del identdc_wrapped_image
            unwrapped_image.SetOrigin(new_origin)
        else:
            unwrapped_image = identdc_wrapped_image
        import os

        unwrapped_outputfn1 = os.path.realpath(unwrapped_outputbasefn[index])
        sitk.WriteImage(unwrapped_image, unwrapped_outputfn1)
        unwrapped_outputfn.append(unwrapped_outputfn1)

    return unwrapped_outputfn
