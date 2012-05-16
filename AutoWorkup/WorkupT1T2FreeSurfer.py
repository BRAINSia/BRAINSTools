#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import ReconAll


"""
    from WorkupT1T2FreeSurfer import CreateFreeSurferWorkflow
    myLocalFSWF= CreateFreeSurferWorkflow("HansFSTest")
    baw200.connect(uidSource,'uid',myLocalFSWF,'InputSpec.subject_id')
    baw200.connect(SplitAvgBABC,'avgBABCT1',myLocalFSWF,'InputSpec.T1_files')
"""

def CreateFreeSurferWorkflow(WFname):
    freesurferWF= pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['subject_id','T1_files','T2_files']), name='InputSpec' )
    
    if 0 == 1:     
        print("""Run Freesurfer ReconAll at""")
        fs_reconall = pe.Node(interface=ReconAll(),name="40_FS510")
        freesurfer_sge_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 1 -l mem_free=3100M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
        fs_reconall.plugin_args=freesurfer_sge_options_dictionary
        fs_reconall.inputs.directive = 'all'
        freesurferWF.connect(inputsSpec,'subject_id',fs_reconall,'subject_id')
        freesurferWF.connect(inputsSpec,'T1_files',  fs_reconall,'T1_files')

    if 0 == 1:
        print("""Run Freesurfer ReconAll at""")
        fs_autorecon1 = pe.Node(interface=ReconAll(),name="41_FS510_autorecon1")
        freesurfer_sge_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 1 -l mem_free=3100M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
        fs_autorecon1.plugin_args=freesurfer_sge_options_dictionary
        fs_autorecon1.inputs.directive = 'autorecon1'
        freesurferWF.connect(inputsSpec,'subject_id',fs_autorecon1,'subject_id')
        freesurferWF.connect(inputsSpec,'T1_files',  fs_autorecon1,'T1_files')

        fs_autorecon2_cp = pe.Node(interface=ReconAll(),name="42cp_FS510_autorecon2_cp")
        freesurfer_sge_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 1 -l mem_free=3100M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
        fs_autorecon2_cp.plugin_args=freesurfer_sge_options_dictionary
        fs_autorecon2_cp.inputs.directive = 'autorecon2-cp'
        freesurferWF.connect(fs_autorecon1,'subject_id',fs_autorecon2_cp,'subject_id')
        freesurferWF.connect(fs_autorecon1,'subjects_dir',fs_autorecon2_cp,'subjects_dir')

        fs_autorecon2_wm = pe.Node(interface=ReconAll(),name="42cp_FS510_autorecon2_wm")
        freesurfer_sge_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 1 -l mem_free=3100M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
        fs_autorecon2_wm.plugin_args=freesurfer_sge_options_dictionary
        fs_autorecon2_wm.inputs.directive = 'autorecon2-wm'
        freesurferWF.connect(fs_autorecon2_cp,'subject_id',fs_autorecon2_wm,'subject_id')
        freesurferWF.connect(fs_autorecon2_cp,'subjects_dir',fs_autorecon2_wm,'subjects_dir')

        fs_autorecon3 = pe.Node(interface=ReconAll(),name="43_FS510_autorecon3")
        freesurfer_sge_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 1 -l mem_free=3100M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
        fs_autorecon3.plugin_args=freesurfer_sge_options_dictionary
        fs_autorecon3.inputs.directive = 'autorecon2-wm'
        freesurferWF.connect(fs_autorecon2_wm,'subject_id',fs_autorecon3,'subject_id')
        freesurferWF.connect(fs_autorecon2_wm,'subjects_dir',fs_autorecon3,'subjects_dir')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['subject_id','subjects_dir']), name='OutputSpec' )

    return freesurferWF
