import os
import sys

from nipype.interfaces.base import (TraitedSpec, Directory, File, traits, InputMultiPath, OutputMultiPath, isdefined)
from nipype.utils.filemanip import split_filename
from nipype.interfaces.base import CommandLine, CommandLineInputSpec


class FSScriptInputSpec(CommandLineInputSpec):
## CROSS first cross sectional analysis
## BASE  second generate a subject specific base reference (template building)
## LONG  third use the BASE, and CROSS to generate a new better informed cross sectional result.
## Universal used
    subcommand = traits.Str('autorecon', argstr='%s', position=0, usedefault=True, desc='Define which subcommand to run: options ="autorecon", "template", "longitudinal"')
    subjects_dir = Directory(argstr='--subjects_dir %s', desc='FreeSurfer subjects directory')

## auto-recon CROSS flags
    T1_files = File(argstr='--T1_files %s', exists=True, desc='Original T1 image')
    brainmask = File(argstr='--brainmask %s', exists=True,
                     desc='The normalized T1 image with the skull removed. Normalized 0-110 where white matter=110.')

## CROSS and LONG flags
    subj_session_id = traits.Str(argstr='--subj_session_id %s', desc='Subject_Session used for "-subjid <> in cross sectional and used in -long <> for longitudinal')

## BASE/Template building flags
    list_all_subj_session_ids = traits.List(traits.Str(), argstr='--list_all_subj_session_ids %s', desc='List of sessions for a subject template')

## LONG and BASE flags
    base_template_id = traits.Str(argstr='--base_template_id %s', desc='The name of the result subdirectory (not full path) for the base/template processing to occur')


class FSScriptOutputSpec(TraitedSpec):
    T1_out = File(exist=True, desc='brain.mgz')
    label1_out = File(exist=True, desc='aparc+aseg.nii.gz')
    label2_out = File(exist=True, desc='aparc.a2009+aseg.nii.gz')
    processed_output_name = traits.Str(desc='The name of the subdirectory (not a full path) for this processing stage')
    outDir = Directory(exist=True, desc='Full path to the output directory for this stage of processing')
    ## HACK: TEST
    subj_session_id = traits.Str(desc='Subject_Session used for "-subjid <> in cross sectional and used in -long <> for longitudinal')


class FSScript(CommandLine):
    """
    Examples
    --------
    """
    import inspect
    import os
    this_file = inspect.getfile(inspect.currentframe())  # script filename (usually with path)
    this_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))  # script)
    _cmd = '{0} {1}'.format(sys.executable, os.path.join(this_path, 'fsscript.py'))
    # _cmd = 'fsscript.py'
    input_spec = FSScriptInputSpec
    output_spec = FSScriptOutputSpec

    def _format_arg(self, opt, spec, val):
        if opt == 'subprocess':
            if not val in ['autorecon', 'template', 'longitudinal']:
                raise ValueException("{0} is not a valid value for 'subprocess'".format(val))
        return super(FSScript, self)._format_arg(opt, spec, val)

    def _list_outputs(self):
        outputs = self._outputs().get()
        ## HACK: TEST
        outputs['subj_session_id'] = self.inputs.subj_session_id
        if self.inputs.subcommand == 'autorecon':
            outputs['T1_out'] = os.path.join(self.inputs.subjects_dir, 'mri', 'brain.mgz')
            outputs['label1_out'] = os.path.join(self.inputs.subjects_dir, 'mri_nifti', 'aparc+aseg.nii.gz')
            outputs['label2_out'] = os.path.join(self.inputs.subjects_dir, 'mri_nifti', 'aparc.a2009+aseg.nii.gz')
            outputs['processed_output_name'] = self.inputs.subj_session_id
            outputs['outDir'] = os.path.join(self.inputs.subjects_dir, self.inputs.subj_session_id)
        elif self.inputs.subcommand == 'template':
            outputs['processed_output_name'] = self.inputs.base_template_id
            outputs['outDir'] = os.path.join(self.inputs.subjects_dir, self.inputs.base_template_id)
        elif self.inputs.subcommand == 'longitudinal':
            longitudinal_processed_output_name = self.inputs.subj_session_id + ".long." + self.inputs.base_template_id
            outputs['processed_output_name'] = longitudinal_processed_output_name
            outputs['outDir'] = os.path.join(self.inputs.subjects_dir, longitudinal_processed_output_name)
        return outputs
