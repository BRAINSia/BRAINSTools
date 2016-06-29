import vtk, math
import sys, os
import SimpleITK as sitk
from fs_thickness_measurements import read_poly_data
from nipype.interfaces.base import BaseInterface, TraitedSpec, BaseInterfaceInputSpec, traits


class Mesh2MaskInputSpec(BaseInterfaceInputSpec):
    input_mesh = traits.File(mandatory=True, exists=True)
    output_image = traits.File(mandatory=True)
    input_image = traits.File(exists=True, mandatory=False)


class Mesh2MaskOutputSpec(TraitedSpec):
    output_image = traits.File()


class Mesh2Mask(BaseInterface):
    input_spec = Mesh2MaskInputSpec
    output_spec = Mesh2MaskOutputSpec

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["output_image"] = os.path.abspath(self.inputs.output_image)
        return outputs

    def _run_interface(self, runtime):
        mesh2mask(inputMesh=self.inputs.input_mesh, outputImage=self.inputs.output_image,
                  inputImage=self.inputs.input_image)
        return runtime


def patch(mesh):
        edges = vtk.vtkFeatureEdges()
        edges.SetInputData(mesh)
        edges.FeatureEdgesOff()
        edges.NonManifoldEdgesOn()
        edges.ManifoldEdgesOff()
        edges.BoundaryEdgesOn()
        edges.Update()
        print("{0} cells and {1} points".format(edges.GetOutput().GetNumberOfCells(),
                                                edges.GetOutput().GetNumberOfPoints()))
        if edges.GetOutput().GetNumberOfPoints() == 0:
            print("No defects")
            return mesh
        
        patchSkel = vtk.vtkStripper()
        patchSkel.SetInputData(edges.GetOutput())
        patchSkel.Update()
        
        patchPoly = vtk.vtkPolyData()
        patchPoly.Initialize()
        patchPoly.SetPoints(patchSkel.GetOutput().GetPoints())
        patchPoly.SetPolys(patchSkel.GetOutput().GetLines())
        
        patchTri = vtk.vtkTriangleFilter()
        patchTri.SetInputData(patchPoly);
        patchTri.Update()
        
        meshAppend = vtk.vtkAppendPolyData()
        meshAppend.AddInputData(patchTri.GetOutput())
        meshAppend.AddInputData(mesh)
        meshAppend.Update()
        
        poly = vtk.vtkCleanPolyData()
        poly.SetInputData(meshAppend.GetOutput())
        poly.Update()
        
        return poly.GetOutput()


def preprocess_mesh(mesh, num_patches=10):

    # run patching
    for i in range(num_patches):
        mesh = patch(mesh)

    edges2 = vtk.vtkFeatureEdges()
    edges2.SetInputData(mesh)
    edges2.FeatureEdgesOff()
    edges2.NonManifoldEdgesOn()
    edges2.ManifoldEdgesOff()
    edges2.BoundaryEdgesOn()
    edges2.Update()
    edgesPd = edges2.GetOutput()
    print("{0} cells and {1} points".format(edgesPd.GetNumberOfCells(),
                                            edgesPd.GetNumberOfPoints()))
    
    toDelete = vtk.vtkFloatArray()
    toDelete.SetNumberOfComponents(1)
    toDelete.SetNumberOfValues(mesh.GetNumberOfPoints())
    
    for i in range(mesh.GetNumberOfPoints()):
        toDelete.SetValue(i, 0)
        
    for i in range(edgesPd.GetNumberOfPoints()):
        edgePnt = edgesPd.GetPoint(i)
        meshVertId = mesh.FindPoint(edgePnt)
        toDelete.SetValue(meshVertId, 1)
        
    mesh.GetPointData().SetScalars(toDelete)
    
    clipper = vtk.vtkClipPolyData()
    clipper.SetValue(0.5)
    clipper.SetInputData(mesh)
    clipper.InsideOutOn()
    clipper.Update()

    return clipper.GetOutput()


def read_vtk_image(inputImage):
    # check image extension
    img_ext = str(inputImage).split(".", 1)[-1]
    if "nii" in img_ext:
            imageReader = vtk.vtkNIFTIImageReader()
    else:
        print("ERROR: Input must be NIFTI image file. \n\tUnrecognized extension: {0}".format(img_ext))
        return sys.exit(os.EX_IOERR)
    # read image in
    imageReader.SetFileName(inputImage)
    imageReader.Update()
    return imageReader.GetOutput()


def mesh2mask(inputMesh, outputImage, inputImage=None, superRes=False, spacing=(1.0, 1.0, 1.0)):
    """
    This program takes in a closed 3D surface, vtkPolyData, and converts it into volume 
    representation (vtkImageData) where the foreground voxels are 1 and the background voxels
    are 0. Internally vtkPolyDataToImageStencil is utilized. The resultant image is saved to disk
    in NIFTIimage file format. SimpleITK is used to convert images to standard orientation used 
    for 3D medical images.

    Inputs:
    inputMesh == a vtkPolyData file of a 3D surface
    outputImage == output file path for NIFTI image
    inputImage == reference image to get desired spacing, origin, and direction.
    """
    
    VTK_MAJOR_VERSION = str(vtk.vtkVersion().GetVTKVersion()).split(".")[0]

    outputImage = str(outputImage)
    if inputImage:
        inputImage = str(inputImage)

    # check output image extension
    out_ext = outputImage.split(".", 1)[-1]
    if "nii" in out_ext:
        writer = vtk.vtkNIFTIImageWriter()
    else:
        print("ERROR: Output must be NIFTI image file. \n\tUnrecognized extension: {0}".format(out_ext))
        return sys.exit(os.EX_IOERR)

    if inputImage:
        image = read_vtk_image(inputImage)
    else:
        image = None


    # read the mesh in
    pd = read_poly_data(inputMesh)
    pd = preprocess_mesh(pd)

    # allocate whiteImage
    whiteImage = vtk.vtkImageData()

    # polygonal data -. image stencil:
    pol2stenc = vtk.vtkPolyDataToImageStencil()
    if VTK_MAJOR_VERSION <= 5:
        pol2stenc.SetInput(pd)
    else:
        pol2stenc.SetInputData(pd)
        
    # get the bounds
    bounds = pd.GetBounds()

    if image:
        spacing = image.GetSpacing()
        dim = image.GetDimensions()
        # set VTK direction to RAS
        vtk_direction = (-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0)

        # vtk does not get the correct origin
        #origin = image.GetOrigin ()

        # use SimpleITK instead
        image_sitk = sitk.ReadImage(inputImage)
        origin = image_sitk.GetOrigin()
        direction = image_sitk.GetDirection()

        print direction
        print origin
        print spacing
        
        # superRes slows down the script quite a bit
        if superRes:
            """Creates an image with pixels a fourth the size of the original
            This helps allivaite some of the partial voluming effect that
            can take place."""
            denom = 2
            spacing = (spacing[0] / float(denom), spacing[1] / float(denom), spacing[2] / float(denom))
            dim = (dim[0] * denom, dim[1] * denom, dim[2] * denom)

        # VTKImages seem to always have the same direction
        origin = (origin[0]*vtk_direction[0],
                  origin[1]*vtk_direction[4],
                  origin[2]*vtk_direction[8])
        if direction != vtk_direction:
            if direction == (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
                origin = (origin[0] - spacing[0]*(dim[0] - 1),
                          origin[1] - spacing[1]*(dim[1] - 1),
                          origin[2])
            else:
                print("ERROR: Not sure how to handle input image direction")
                sys.exit()
        print origin
    else:
        if superRes:
            spacing = (spacing[0] / float(2), spacing[1] / float(2), spacing[2] / float(2))

        # compute dimensions
        dim = [0, 0, 0]
        for i in range(3):
            dim[i] = int(math.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]))
        dim = tuple(dim)

        # get origin
        origin = [0, 0, 0]
        origin[0] = bounds[0] + spacing[0] / float(2)
        origin[1] = bounds[2] + spacing[1] / float(2)
        origin[2] = bounds[4] + spacing[2] / float(2)
        origin = tuple(origin)

    whiteImage.SetSpacing(spacing)
    whiteImage.SetOrigin(origin)

    pol2stenc.SetOutputOrigin(origin)
    pol2stenc.SetOutputSpacing(spacing)
    
    # set dimensions
    whiteImage.SetDimensions(dim)
    whiteImage.SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1)

    if VTK_MAJOR_VERSION <= 5:
        whiteImage.SetScalarTypeToUnsignedChar()
        whiteImage.AllocateScalars()
    else:
        whiteImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,1)

    # fill the image with foreground voxels:
    inval = 1
    outval = 0
    count = whiteImage.GetNumberOfPoints()
    for i in range(count):
        whiteImage.GetPointData().GetScalars().SetTuple1(i, inval)

    # update pol2stenc
    pol2stenc.SetOutputWholeExtent(whiteImage.GetExtent())
    pol2stenc.Update()

    # cut the corresponding white image and set the background:
    imgstenc = vtk.vtkImageStencil()
    if VTK_MAJOR_VERSION <= 5:
        imgstenc.SetInput(whiteImage)
        imgstenc.SetStencil(pol2stenc.GetOutput())
    else:
        imgstenc.SetInputData(whiteImage)
        imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())

    imgstenc.ReverseStencilOff()
    imgstenc.SetBackgroundValue(outval)
    imgstenc.Update()

    # write image to file
    
    writer.SetFileName( outputImage )

    if inputImage != None and direction != vtk_direction:
        flipFilter = vtk.vtkImageFlip()
        flipFilter.SetFilteredAxis(1) # flip y axis
        flipFilter.SetInputData( imgstenc.GetOutput() )
        flipFilter.SetFlipAboutOrigin( 1 )
        flipFilter.Update()

        flipFilter2 = vtk.vtkImageFlip()
        flipFilter2.SetFilteredAxis(0) # flip x axis
        flipFilter2.SetInputData( flipFilter.GetOutput() )
        flipFilter2.SetFlipAboutOrigin( 1 )
        flipFilter2.Update()
    
        if VTK_MAJOR_VERSION <= 5:
            writer.SetInput(flipFilter2.GetOutput())
        else:
            writer.SetInputData(flipFilter2.GetOutput())
        writer.Write()

        itk_image = sitk.ReadImage(inputImage)
        out_image = sitk.ReadImage(outputImage)
        out_image.SetDirection(itk_image.GetDirection())
        out_image.SetOrigin(itk_image.GetOrigin())
        sitk.WriteImage(out_image, outputImage)
    else:
        if VTK_MAJOR_VERSION <= 5:
            writer.SetInput(imgstenc.GetOutput())
        else:
            writer.SetInputData(imgstenc.GetOutput())
        writer.Write()    
    
    return os.path.abspath(outputImage)


def usage():
    print(mesh2mask.__doc__)


if __name__ == "__main__":
    import getopt, sys
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:r:", ["help", "inMesh=", "outImage=", "refImage="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    outputImage = None
    inputMesh = None
    refImage = None
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--outImage"):
            outputImage = a
        elif o in ("-i", "--inMesh"):
            inputMesh = a
            if not os.path.isfile(inputMesh):
                print("ERROR: {0} does not exist.".format(inputMesh))

        elif o in ("-r", "--refImage"):
            refImage = a
            if not os.path.isfile(refImage):
                    print("ERROR: {0} does not exist.".format(refImage))
        else:
            assert False, "unhandled option"

    # check inputs
    if inputMesh != None:
        print("Input mesh: {0}".format(inputMesh))
    else:
        print("Error: input required for --inMesh")
        usage()
        sys.exit()
    if outputImage != None:
        print("Output image: {0}".format(outputImage))
    else:
        print("Error: input required for --outImage")
        usage()
        sys.exit()

    print("Input image: {0}".format(refImage))
    

    print("Converting mesh to mask")
    print("Mask created: {0}".format( mesh2mask(inputMesh, outputImage, refImage) ) )

