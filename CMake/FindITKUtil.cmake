#
# encapsulate calling find_package(ITK)

## ITK_VAR_PREFIX is a required parameter to this macro
## "ITK_VAR_PREFIX" is the new prefix variable to be used instead
## of "ITK" for [_LIBRARIES|_INCLUDE_DIRS|_LIBRARY_DIRS] when you want to link
## build or link against ITK
## Instead of using ${ITK_INCLUDE_DIRS}, you now use ${${ITK_VAR_PREFIX}_INCLUDE_DIRS}
## Instead of using ${ITK_LIBRARIES}, you now use ${${ITK_VAR_PREFIX}_LIBRARIES}
## Instead of using ${ITK_LIBRARY_DIRS}, you now use ${${ITK_VAR_PREFIX}_LIBRARY_DIRS}
## After ITK_VAR_PREFIX the user of the FindITKUtil macro should provide
## a list of optional modules that are needed and will be available in the ${ARGN}
## cmake variable.
macro(FindITKUtil ITK_VAR_PREFIX )

  # ITK_FOUND needs to be reset, or it won't redo
  # setting up the include directories
  unset(ITK_FIND_COMPONENTS)
  unset(ITK_MODULES_REQUESTED)
  unset(ITK_FOUND)
  set( ITK_INCLUSION_MODULES
    # Everything needs ITKCommon
    ITKCommon
    # Common depends on these modules
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
    ITKIOMRC
    ITKIOSpatialObjects
    ITKIOTransformBase
    ITKIOHDF5
    ITKIOTransformMatlab
    ITKIOTransformHDF5
    ITKIOGE
    MGHIO
    # other modules specific to the current directory
    ${ARGN}
    )
  list(REMOVE_DUPLICATES ITK_INCLUSION_MODULES)
  find_package(ITK 4.5 COMPONENTS ${ITK_INCLUSION_MODULES} REQUIRED)
  #itk_module_config(${ITK_VAR_PREFIX} ${ITK_INCLUSION_MODULES})

  include(${ITK_USE_FILE})
  set(${ITK_VAR_PREFIX}_INCLUDE_DIRS ${ITK_INCLUDE_DIRS})
  set(${ITK_VAR_PREFIX}_LIBRARY_DIRS ${ITK_LIBRARY_DIRS})
  set(${ITK_VAR_PREFIX}_LIBRARIES ${ITK_LIBRARIES})

  unset(ITK_INCLUDE_DIRS)
  unset(ITK_LIBRARY_DIRS)
  unset(ITK_LIBRARIES)
  unset(ITK_FIND_COMPONENTS)
  unset(ITK_MODULES_REQUESTED)
  unset(ITK_FOUND)

  ## Now include directories associated with the requested modules
  include_directories(${${ITK_VAR_PREFIX}_INCLUDE_DIRS})
endmacro()
