macro( _set_if_not_empty var value )
  if( NOT "${value}" STREQUAL "" )
    set( ${var} "${value}" )
  endif()
endmacro()

set( _GIT_VERSION_MAJOR "5" )
set( _GIT_VERSION_MINOR "4" )
_set_if_not_empty( _GIT_VERSION_PATCH "0" )
_set_if_not_empty( _GIT_VERSION_TWEAK "" )
_set_if_not_empty( _GIT_VERSION_RC "" )
## DO NOT SET THE CONSTANTLY CHANGING VERSION_POST for DISTRIBUTION RELEASES
_set_if_not_empty( _GIT_VERSION_POST "POSTTAG" )
## DO NOT SET THE CONSTANTLY CHANGING VERSION_DEV for DISTRIBUTION RELEASES
_set_if_not_empty( _GIT_VERSION_DEV "DEV" )
## DO NOT SET THE CONSTANTLY CHANGING HASH FOR DISTRIBUTION RELEASES
_set_if_not_empty( _GIT_VERSION_HASH "DIST" )
