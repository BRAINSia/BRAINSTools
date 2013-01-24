import os
import sys

from nipype.interfaces.base import (TraitedSpec, Directory, File, traits, InputMultiPath, OutputMultiPath, isdefined)
from nipype.utils.filemanip import split_filename
from nipype.interfaces.base import CommandLine, CommandLineInputSpec


class FSScriptInputSpec(CommandLineInputSpec):
    subject_id = traits.Str(argstr='--subject_id %s', desc='Subject_Session')
    subjects_dir = Directory(argstr='--subjects_dir %s', desc='FreeSurfer subjects directory')
    subcommand = traits.Str('autorecon', argstr='%s', position=0, usedefault=True, desc='Define which subcommand to run: options ="autorecon", "template", "longitudinal"')
    T1_files = File(argstr='--T1_files %s', exists=True, desc='Original T1 image')
    brainmask = File(argstr='--brainmask %s', exists=True,
                     desc='The normalized T1 image with the skull removed. Normalized 0-110 where white matter=110.')
    subjectTemplate_id = traits.Str(argstr='--subjectTemplate_id %s', desc='Subject_template')
    session_ids = traits.List(traits.Str(), argstr='--session_ids %s', desc='List of sessions for a subject template')
    session_id = traits.Str(argstr='--session_id %s', desc='Session for a subject longitudinal analysis')
    template_id = traits.Str(desc='Template ID used in longitudinal processing')
    # TODO: fs_env_script = traits.Str(argstr='--FSSource %s', default='${FREESURFER_HOME}/FreeSurferEnv.sh', desc='')
    # TODO: fs_home = Directory(argstr='--FSHomeDir %s', desc='Location of FreeSurfer (differs for Mac and Linux environments')


class FSScriptOutputSpec(TraitedSpec):
    T1_out = File(exist=True, desc='brain.mgz')
    label1_out = File(exist=True, desc='aparc+aseg.nii.gz')
    label2_out = File(exist=True, desc='aparc.a2009+aseg.nii.gz')
    outDir = Directory(exist=True, desc='Template directory for subject')

class FSScript(CommandLine):
    """
    Examples
    --------
    """
    import inspect, os
    this_file = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
    this_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script)
    _cmd = '{0} {1}'.format(sys.executable, os.path.join(this_path,'fsscript.py'))
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
        if self.inputs.subcommand == 'autorecon':
            outputs['T1_out'] = os.path.join(os.getcwd(), 'mri', 'brain.mgz')
            outputs['label1_out'] = os.path.join(os.getcwd(), 'mri_nifti', 'aparc+aseg.nii.gz')
            outputs['label2_out'] = os.path.join(os.getcwd(), 'mri_nifti', 'aparc.a2009+aseg.nii.gz')
            outputs['outDir'] = self.inputs.subjects_dir
        elif self.inputs.subcommand == 'template':
            outputs['outDir'] = os.path.join(os.getcwd(), self.inputs.subjectTemplate_id)
        elif self.inputs.subcommand == 'longitudinal':
            templateFile = self.inputs.template_id + "_" + self.inputs.session_id + ".long"
            outputs['outDir'] = os.path.join(os.getcwd(), templateFile)
        return outputs
