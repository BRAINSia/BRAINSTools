from __future__ import print_function
__author__ = 'johnsonhj'

######################################################################################
###### Now ensure that all the required packages can be read in from this custom path
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# print sys.path
from nipype import config  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
config.enable_debug_mode()  # NOTE:  This needs to occur AFTER the PYTHON_AUX_PATHS has been modified
#config.enable_provenance()

##############################################################################
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine
import nipype.interfaces.io as nio   # Data i/o
from nipype.interfaces.freesurfer import ReconAll

from nipype.utils.misc import package_check
# package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
package_check('numpy', '1.3', 'tutorial1')
package_check('scipy', '0.7', 'tutorial1')
package_check('networkx', '1.0', 'tutorial1')
package_check('IPython', '0.10', 'tutorial1')

import os
## Check to ensure that SimpleITK can be found
#import SimpleITK as sitk

SLICER_REFERENCE_DIR='/scratch/DWI_DATA'
SLICER_RESULTS_DIR='/scratch/DWI_DATA/SlicerResults'
ExperimentBaseDirectoryCache = SLICER_RESULTS_DIR

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#####################################################################################
#     Prepend the shell environment search paths
PROGRAM_PATHS = '/Users/johnsonhj/src/BT-build/bin:/usr/local/bin'
PROGRAM_PATHS = PROGRAM_PATHS.split(':')
PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))
os.environ['PATH'] = ':'.join(PROGRAM_PATHS)

from nipype.interfaces.semtools import *

print("Building Pipeline")
########### PIPELINE INITIALIZATION #############
DWI_AutoProcess = pe.Workflow(name="DWI_20130515")  # HACK: This needs to be specified in the config file.
DWI_AutoProcess.config['execution'] = {
    'plugin': 'Linear',
    #'stop_on_first_crash':'true',
    #'stop_on_first_rerun': 'true',
    'stop_on_first_crash': 'false',
    'stop_on_first_rerun': 'false',  # This stops at first attempt to rerun, before running, and before deleting previous results.
    'hash_method': 'timestamp',
    'single_thread_matlab': 'true',  # Multi-core 2011a  multi-core for matrix multiplication.
    'remove_unnecessary_outputs': 'true', #remove any interface outputs not needed by the workflow
    'use_relative_paths': 'false',  # relative paths should be on, require hash update when changed.
    'remove_node_directories': 'false',  # Experimental
    'local_hash_check': 'true',
    'job_finished_timeout': 45
}
DWI_AutoProcess.config['logging'] = {
    'workflow_level': 'DEBUG',
    'filemanip_level': 'DEBUG',
    'interface_level': 'DEBUG',
    'log_directory': ExperimentBaseDirectoryCache
}
DWI_AutoProcess.base_dir = ExperimentBaseDirectoryCache




DWI_DG = pe.Node(interface=nio.DataGrabber(sort_filelist = False,outfields=['DWIQCed','ACPCT1','ACPCT2','FST1','FSWMParc']),
    run_without_submitting=True,
    name='ImageInputs'
    )
DWI_DG.inputs.base_directory = SLICER_REFERENCE_DIR
DWI_DG.inputs.template = '*'
DWI_DG.inputs.field_template = dict(
    DWIQCed='20121019_DTIPrep/HDNI_001/068044003/068044003_20120522_30/DTIPrepOutput/068044003_068044003_20120522_30_DWI_CONCAT_QCed.nrrd',
    ACPCT1='20130109_TrackOn_Results/HDNI_001/068044003/068044003_20120522_30/TissueClassify/t1_average_BRAINSABC.nii.gz',
    ACPCT2='20130109_TrackOn_Results/HDNI_001/068044003/068044003_20120522_30/TissueClassify/t2_average_BRAINSABC.nii.gz',
    FST1='20130109_TrackOn_Results/FREESURFER52_SUBJECTS/068044003_068044003_20120522_30/mri/T1.mgz',
    FSWMParc='20130109_TrackOn_Results/FREESURFER52_SUBJECTS/068044003_068044003_20120522_30/mri/wmparc.mgz'
)
#print os.path.join(SLICER_REFERENCE_DIR,'20130109_TrackOn_Results/HDNI_*/068044003/068044003_20120522_30/TissueClassify/t1_average_BRAINSABC.nii.gz')
#print DWI_DG.outputs

##  First Generate a transform from DWIQCed space to ACPC Space
BRAINSFitDWI2ACPC=pe.Node(interface=BRAINSFit(),name='DWI2ACPC_TransformEstimation')
BRAINSFitDWI2ACPC.inputs.useRigid = True
BRAINSFitDWI2ACPC.inputs.initializeTransformMode = 'useMomentsAlign'
BRAINSFitDWI2ACPC.inputs.maskProcessingMode = 'ROIAUTO'
BRAINSFitDWI2ACPC.inputs.ROIAutoDilateSize = 7
BRAINSFitDWI2ACPC.inputs.outputTransform = 'DWI_to_ACPC.h5'
BRAINSFitDWI2ACPC.inputs.writeOutputTransformInFloat = True
#print BF.cmdline

DWI_AutoProcess.connect(DWI_DG,'DWIQCed',BRAINSFitDWI2ACPC,'movingVolume')
DWI_AutoProcess.connect(DWI_DG,'ACPCT2',BRAINSFitDWI2ACPC,'fixedVolume')

DWI_ResampleInPlace=pe.Node(interface=gtractResampleDWIInPlace(),name='DWI_ResampleInPlace')
DWI_ResampleInPlace.inputs.outputVolume = 'ACPC_DWI.nrrd'

DWI_AutoProcess.connect(DWI_DG,'DWIQCed',DWI_ResampleInPlace,'inputVolume')
DWI_AutoProcess.connect(BRAINSFitDWI2ACPC,'outputTransform',DWI_ResampleInPlace,'inputTransform')

##  Second Generate a transform from FreeSurfer space to ACPC Space

BRAINSFitFS2ACPC=pe.Node(interface=BRAINSFit(),name='FS2ACPC_TransformEstimation')
BRAINSFitFS2ACPC.inputs.useRigid = True
BRAINSFitFS2ACPC.inputs.initializeTransformMode = 'useMomentsAlign'
BRAINSFitFS2ACPC.inputs.maskProcessingMode = 'ROIAUTO'
BRAINSFitFS2ACPC.inputs.ROIAutoDilateSize = 12
BRAINSFitFS2ACPC.inputs.outputTransform = 'FS_to_ACPC.h5'
BRAINSFitFS2ACPC.inputs.writeOutputTransformInFloat = True
#print BF.cmdline

DWI_AutoProcess.connect(DWI_DG,'FST1',BRAINSFitFS2ACPC,'movingVolume')
DWI_AutoProcess.connect(DWI_DG,'ACPCT1',BRAINSFitFS2ACPC,'fixedVolume')

## Resample from FS space to ACPC space with DWI_ACPC voxel grid
ResampleWMParc2DWIACPC = pe.Node(interface=BRAINSResample(),name='ResampleWMParc2DWIACPC')
ResampleWMParc2DWIACPC.inputs.interpolationMode = 'NearestNeighbor'
ResampleWMParc2DWIACPC.inputs.outputVolume = 'FSWMParcInDWIACPC.nrrd'
DWI_AutoProcess.connect(BRAINSFitFS2ACPC,'outputTransform',ResampleWMParc2DWIACPC,'warpTransform')
DWI_AutoProcess.connect(DWI_DG,'FSWMParc',ResampleWMParc2DWIACPC,'inputVolume')
DWI_AutoProcess.connect(DWI_ResampleInPlace,'outputVolume',ResampleWMParc2DWIACPC,'referenceVolume')

ResampleWMParc2ACPC = pe.Node(interface=BRAINSResample(),name='ResampleWMParcToACPC')
ResampleWMParc2ACPC.inputs.interpolationMode = 'NearestNeighbor'
ResampleWMParc2ACPC.inputs.outputVolume = 'FSWMParcInACPC.nrrd'

DWI_AutoProcess.connect(BRAINSFitFS2ACPC,'outputTransform',ResampleWMParc2ACPC,'warpTransform')
DWI_AutoProcess.connect(DWI_DG,'FSWMParc',ResampleWMParc2ACPC,'inputVolume')
DWI_AutoProcess.connect(DWI_DG,'ACPCT1',ResampleWMParc2ACPC,'referenceVolume')


DWI_AutoProcess.write_graph()
#DWI_AutoProcess._write_report_info()
DWI_AutoProcess.run()

x=BRAINSResample()
print(x.inputs)
