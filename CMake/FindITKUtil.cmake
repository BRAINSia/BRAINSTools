#
# encapsulate calling find_package(ITK)

macro(FindITKUtil)
  if(Slicer_BUILD_BRAINSTOOLS) ## Slicer has it's own internal MGHIO that conflicts with ITK
    set(FindITK_MGHIO "")
  else()
    set(FindITK_MGHIO MGHIO )
  endif()

  # ITK_FOUND needs to be reset, or it won't redo
  # setting up the include directories
  set(ITK_FOUND OFF)
  find_package(ITK COMPONENTS
    # Everything needs ITKCommon
    ITKCommon
    # Common depends on thes modules
    ITKVNLInstantiation
    ITKKWSys
    ITKDoubleConversion
    ITKVNLInstantiation
    ITKVNL
    # IO Components
    ITKIOImageBase
    ITKIOBMP
    ITKIOBioRad
    ITKIOGDCM
    ITKIOGIPL
    ITKIOJPEG
    ITKIOLSM
    ITKIOMeta
    ITKIONIFTI
    ITKIONRRD
    ITKIOPNG
    ITKIOStimulate
    ITKIOTIFF
    ITKIOTransformInsightLegacy
    ITKIOVTK
    ITKIOSpatialObjects
    ITKIOTransformBase
    ITKIOHDF5
    ITKIOTransformMatlab
    ITKIOTransformHDF5
    ITKIOGE
    ${FindITK_MGHIO}
    # other modules specific to the current directory
    ${ARGN}
    REQUIRED)

  if(Slicer_BUILD_BRAINSTOOLS)
    set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
  endif()
  include(${ITK_USE_FILE})

endmacro()
