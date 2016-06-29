import numpy as np
import SimpleITK as sitk
import sys


def divide_nonzero(array1, array2):
    """
    Divides two arrays. Returns zero when dividing by zero.
    """
    nonzero_idx = array2 != 0
    output = np.zeros_like(array1)
    output[nonzero_idx] = np.divide(array1[nonzero_idx], array2[nonzero_idx])
    return output


def create_image_like(data, image):
    new_image = sitk.GetImageFromArray(data)
    new_image.CopyInformation(image)
    return new_image


def get_directional_images(gradient_image):
    """
    Warning: Directions are defined according to the voxel latice and not the RAS space.
    """
    gradient_array = sitk.GetArrayFromImage(gradient_image)
    x_gradient = create_image_like(gradient_array[:, :, :, 0], gradient_image)
    y_gradient = create_image_like(gradient_array[:, :, :, 1], gradient_image)
    z_gradient = create_image_like(gradient_array[:, :, :, 2], gradient_image)
    return x_gradient, y_gradient, z_gradient


def compute_directional_gradient_images(image, sigma):
    gradient_image = compute_gradient(image, sigma=sigma)
    return get_directional_images(gradient_image)


def compute_gradient(image, sigma=0):
    if sigma > 0:
        return sitk.GradientRecursiveGaussian(image, sigma=sigma, useImageDirection=True)
    else:
        return sitk.Gradient(image, useImageSpacing=True, useImageDirection=True)


def compute_hessian_matrix(image, sigma):
    """
    Computes the hessian matrix for an image. This can be used to detect vesselness as well as other image features.

    First Derivative:
    [ gx,  gy,  gz ]

    Hessian Matrix:
    [ gxx, gxy, gxz]
    [ gyx, gyy, gyz]
    [ gzx, gzy, gzz]
    """
    x_gradient, y_gradient, z_gradient = compute_directional_gradient_images(image, sigma=sigma)
    x2_gradient, y2_gradient, z2_gradient = [compute_gradient(grad, sigma=sigma) for grad in [x_gradient,
                                                                                              y_gradient,
                                                                                              z_gradient]]
    x2_array, y2_array, z2_array = [sitk.GetArrayFromImage(grad2) for grad2 in [x2_gradient,
                                                                                y2_gradient,
                                                                                z2_gradient]]
    hessian = np.stack((x2_array, y2_array, z2_array), axis=-1)
    if sigma > 0:
        return hessian * np.square(sigma)
    return hessian


def separate_eigen_values(eigen_values):
    return [eigen_values[:, :, :, i] for i in range(eigen_values.shape[-1])]


def compute_eigen_values_from_hessian(hessian):
    eigen_values = np.linalg.eigvalsh(hessian[:, :, :, :, :])
    eigen_values_list = separate_eigen_values(eigen_values)
    if check_eigen_values(*eigen_values_list):
        return eigen_values_list
    else:
        sorted_eigen_values = sortbyabs(eigen_values, axis=-1)
        sorted_eigen_values_list = separate_eigen_values(sorted_eigen_values)
        if check_eigen_values(*sorted_eigen_values_list):
            return sorted_eigen_values_list
        else:
            print("Could not sort eigen values")
            sys.exit(1)


def sortbyabs(a, axis=0):
    """Sort array along a given axis by the absolute value
    modified from: http://stackoverflow.com/a/11253931/4067734
    """
    import numpy as np
    index = list(np.ix_(*[np.arange(i) for i in a.shape]))
    index[axis] = np.abs(a).argsort(axis)
    return a[index]


def compute_measures(eigen1, eigen2, eigen3):
    """
    RA - plate-like structures
    RB - blobl-like structures
    S - background
    """
    Ra = divide_nonzero(np.abs(eigen1), np.abs(eigen2))
    Rb = divide_nonzero(np.abs(eigen1), np.sqrt(np.abs(np.multiply(eigen2, eigen3))))
    S = np.sqrt(np.square(eigen1) + np.square(eigen2) + np.square(eigen3))
    return Ra, Rb, S


def compute_plate_like_factor(Ra, alpha):
    return 1 - np.exp(np.negative(np.square(Ra)) / (2 * np.square(alpha)))


def compute_blob_like_factor(Rb, beta):
    return np.exp(np.negative(np.square(Rb) / (2 * np.square(beta))))


def compute_background_factor(S, c):
    return 1 - np.exp(np.negative(np.square(S)) / (2 * np.square(c)))


def compute_vesselness(eigen1, eigen2, eigen3, alpha, beta, c, black_white):
    Ra, Rb, S = compute_measures(eigen1, eigen2, eigen3)
    plate = compute_plate_like_factor(Ra, alpha)
    blob = compute_blob_like_factor(Rb, beta)
    background = compute_background_factor(S, c)
    return filter_out_background(plate * blob * background, black_white, eigen2, eigen3)


def compute_vesselness_image(image, sigma=0, alpha=0.5, beta=0.5, frangic=500, black_white=False):
    eigen1, eigen2, eigen3 = compute_eigen_values(image, sigma=sigma)
    voxel_data = compute_vesselness(eigen1=eigen1, eigen2=eigen2, eigen3=eigen3, black_white=black_white,
                                    alpha=alpha, beta=beta, c=frangic)
    return create_image_like(voxel_data, image)


def compute_eigen_values(image, sigma):
    hessian = compute_hessian_matrix(image, sigma=sigma)
    eigen1, eigen2, eigen3 = compute_eigen_values_from_hessian(hessian)
    return eigen1, eigen2, eigen3


def compute_absolute_eigen_values(image, sigma):
    return [np.abs(eigen) for eigen in compute_eigen_values(image, sigma=sigma)]


def check_eigen_values(eigen1, eigen2, eigen3):
    """
    Check that |eigen1| <= |eigen2| <= |eigen3|
    """
    if np.all(np.abs(eigen1) <= np.abs(eigen2)) and np.all(np.abs(eigen2) <= np.abs(eigen3)):
        return True
    else:
        return False


def filter_out_background(voxel_data, black_white, eigen2, eigen3):
    """
    Set black_white to true if vessels are darker than the background and to false if
    vessels are brighter than the background.
    """
    if black_white:
        voxel_data[eigen2 < 0] = 0
        voxel_data[eigen3 < 0] = 0
    else:
        voxel_data[eigen2 > 0] = 0
        voxel_data[eigen3 > 0] = 0
    voxel_data[np.isnan(voxel_data)] = 0
    return voxel_data


def compute_vesselness_image_with_smoothing(image, sigmas, alpha=0.5, beta=0.5, frangic=500, black_white=True):
    output_vesselness_image = None
    for sigma in sigmas:
        str_sigma = "{0:.1f}".format(sigma).replace(".", "_")
        print("--- sigma = {sigma} ---".format(sigma=str_sigma))

        vessel_image = compute_vesselness_image(image, sigma=sigma, alpha=alpha, beta=beta, frangic=frangic,
                                                black_white=black_white)

        if output_vesselness_image:
            output_vesselness_image = sitk.Maximum(output_vesselness_image, vessel_image)
        else:
            output_vesselness_image = vessel_image
    return output_vesselness_image