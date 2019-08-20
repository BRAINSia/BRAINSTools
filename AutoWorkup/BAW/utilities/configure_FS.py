"""
configure_FS.py
============================
Description:
    The purpose of this is to...
Author:

Usage:

By using
         execfile('path/to/configure_FS', OrderedDict(env={os.environ-like dictionary})
    you will set the FreeSurfer environment driven by the configuration file.
"""


try:
    list(env.keys())
    FS_VARS
except NameError as AttributeError:
    raise AssertionError(
        "Run this file like: execfile('path/to/configure_FS', OrderedDict(env={}) \
        where 'env' is set to an os.environ-like dictionary"
    )
import os
from . import misc

#####################################################################################
#  FreeSurfer is extraordinarly finicky and is easily confused and incorrect.

# INFO: Remove all FREESURFER vars in subsequent scripts!!!
#  Force that all the FREESURFER env vars are set in subsequent scripts by
#  ensuring that rough versions of these environmental variables are not
#  set internal to this script.

for ENVVAR_TO_CHECK in FS_VARS:
    if ENVVAR_TO_CHECK in os.environ:
        if ENVVAR_TO_CHECK in env:
            os.environ[ENVVAR_TO_CHECK] = env[ENVVAR_TO_CHECK]
        else:
            raise EnvironmentError(
                "Freesurfer variable {0}={1} exists! \
                Please unset before continuing.".format(
                    ENVVAR_TO_CHECK, os.environ[ENVVAR_TO_CHECK]
                )
            )
    else:
        pass
