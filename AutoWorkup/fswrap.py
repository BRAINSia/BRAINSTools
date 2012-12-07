# Standard library imports
import os

# Local imports
from nipype.interfaces.base import (TraitedSpec, File, traits, InputMultiPath, OutputMultiPath, isdefined)
from nipype.utils.filemanip import split_filename
from nipype.interfaces.base import Command, CommandInputSpec


class FSScriptInputSpec(CommandInputSpec):
    T1_image = File(argstr='--T1image %s', exists=True, desc='Original T1 image')
    brainmask = File(argstr='--Brainmask %s', exists=True,
                     desc='The normalized T1 image with the skull removed. Normalized 0-110 where white matter=110.')
    wm_prob = File(argstr='--WMProb %s', exists=True, desc='')
    fs_home = Directory(argstr='--FSHomeDir %s',
                         desc='Location of FreeSurfer (differs for Mac and Linux environments')
    fs_subj_dir = Directory(argstr='--FSSubjDir %s', desc='FreeSurfer subjects directory')
    fs_source = traits.Str(argstr='--FSSource %s', default='${FREESURFER_HOME}/FreeSurferEnv.csh',
                           desc='')
    subj_id = traits.Str(argstr='--SubjID %s', desc='Subject_Session')

class FSScriptOutputSpec(TraitedSpec):
    pass

class FSScript(Command):
    """
    Examples
    --------
    """

    _cmd = 'fsscript'
    input_spec = FSScriptInputSpec
    output_spec = FSScriptOutputSpec

    def _list_outputs(self):
        outputs = self._outputs().get()
        return outputs
