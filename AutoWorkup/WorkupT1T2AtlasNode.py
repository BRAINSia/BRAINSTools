
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine
import os

#######################  HACK:  Needed to make some global variables for quick
#######################         processing needs
#Generate by running a file system list "ls -1 $AtlasDir *.nii.gz *.xml *.fcsv *.wgts"
#atlas_file_list="AtlasPVDefinition.xml ALLPVAIR.nii.gz ALLPVBASALTISSUE.nii.gz ALLPVCRBLGM.nii.gz ALLPVCRBLWM.nii.gz ALLPVCSF.nii.gz ALLPVNOTCSF.nii.gz ALLPVNOTGM.nii.gz ALLPVNOTVB.nii.gz ALLPVNOTWM.nii.gz ALLPVSURFGM.nii.gz ALLPVVB.nii.gz ALLPVWM.nii.gz avg_t1.nii.gz avg_t2.nii.gz tempNOTVBBOX.nii.gz template_ABC_lables.nii.gz template_WMPM2_labels.nii.gz template_WMPM2_labels.txt template_brain.nii.gz template_cerebellum.nii.gz template_class.nii.gz template_headregion.nii.gz template_leftHemisphere.nii.gz template_nac_lables.nii.gz template_nac_lables.txt template_rightHemisphere.nii.gz template_t1.nii.gz template_t1_clipped.nii.gz template_t2.nii.gz template_t2_clipped.nii.gz template_ventricles.nii.gz probabilityMaps/l_caudate_ProbabilityMap.nii.gz probabilityMaps/r_caudate_ProbabilityMap.nii.gz probabilityMaps/l_hippocampus_ProbabilityMap.nii.gz probabilityMaps/r_hippocampus_ProbabilityMap.nii.gz probabilityMaps/l_putamen_ProbabilityMap.nii.gz probabilityMaps/r_putamen_ProbabilityMap.nii.gz probabilityMaps/l_thalamus_ProbabilityMap.nii.gz probabilityMaps/r_thalamus_ProbabilityMap.nii.gz spatialImages/phi.nii.gz spatialImages/rho.nii.gz spatialImages/theta.nii.gz"
#atlas_file_names=atlas_file_list.split(' ')
## HACK
atlas_file_names=["AtlasPVDefinition.xml","AtlasPVDefinition.xml.in","ALLPVAIR.nii.gz",
                      "ALLPVBASALTISSUE.nii.gz","ALLPVCRBLGM.nii.gz",
                      "ALLPVCRBLWM.nii.gz","ALLPVCSF.nii.gz","ALLPVNOTCSF.nii.gz",
                      "ALLPVNOTGM.nii.gz","ALLPVNOTVB.nii.gz","ALLPVNOTWM.nii.gz",
                      "ALLPVSURFGM.nii.gz","ALLPVVB.nii.gz","ALLPVWM.nii.gz",
                      "avg_t1.nii.gz","avg_t2.nii.gz","tempNOTVBBOX.nii.gz",
                      "template_ABC_lables.nii.gz","template_WMPM2_labels.nii.gz",
                      "template_WMPM2_labels.txt","template_brain.nii.gz",
                      "template_cerebellum.nii.gz","template_class.nii.gz",
                      "template_headregion.nii.gz","template_leftHemisphere.nii.gz",
                      "template_nac_lables.nii.gz","template_nac_lables.txt",
                      "template_rightHemisphere.nii.gz","template_t1.nii.gz",
                      "template_t1_clipped.nii.gz","template_t2.nii.gz",
                      "template_t2_clipped.nii.gz","template_ventricles.nii.gz",
                      "template_landmarks.fcsv","template_landmark_weights.csv",
                      "template_landmarks_31.fcsv","template_landmark_weights_31.csv",

                   "probabilityMaps/l_accumben_ProbabilityMap.nii.gz",
                   "probabilityMaps/r_accumben_ProbabilityMap.nii.gz",
                   "probabilityMaps/l_caudate_ProbabilityMap.nii.gz",
                   "probabilityMaps/r_caudate_ProbabilityMap.nii.gz",
                   "probabilityMaps/l_globus_ProbabilityMap.nii.gz",
                   "probabilityMaps/r_globus_ProbabilityMap.nii.gz",
                   "probabilityMaps/l_hippocampus_ProbabilityMap.nii.gz",
                   "probabilityMaps/r_hippocampus_ProbabilityMap.nii.gz",
                   "probabilityMaps/l_putamen_ProbabilityMap.nii.gz",
                   "probabilityMaps/r_putamen_ProbabilityMap.nii.gz",
                   "probabilityMaps/l_thalamus_ProbabilityMap.nii.gz",
                   "probabilityMaps/r_thalamus_ProbabilityMap.nii.gz",

                   "spatialImages/phi.nii.gz",
                   "spatialImages/rho.nii.gz",
                   "spatialImages/theta.nii.gz",

                   "modelFiles/RandomForestAllSubcorticalsBalancedModel.txtD0060NT0060.gz"

                      ]
## Remove filename extensions for images, but replace . with _ for other file types
atlas_file_keys=[os.path.basename(fn).replace('.nii.gz','').replace('.','_') for fn in atlas_file_names]
atlas_outputs_filename_match = dict(zip(atlas_file_keys,atlas_file_names))

def MakeAtlasNode(atlasDirectory,AtlasNodeName):
    BAtlas = pe.Node(interface=nio.DataGrabber(outfields=atlas_file_keys),
                               run_without_submitting=True,
                               name=AtlasNodeName)
    BAtlas.inputs.base_directory = atlasDirectory
    BAtlas.inputs.template = '*'
    ## Prefix every filename with atlasDirectory
    atlas_search_paths=['{0}'.format(fn) for fn in atlas_file_names]
    BAtlas.inputs.field_template = dict(zip(atlas_file_keys,atlas_search_paths))
    ## Give 'atlasDirectory' as the substitution argument
    atlas_template_args_match=[ [[]] for i in atlas_file_keys ] ##build a list of proper lenght with repeated entries
    BAtlas.inputs.template_args = dict(zip(atlas_file_keys,atlas_template_args_match))
    return BAtlas
