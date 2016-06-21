from __future__ import print_function
from builtins import zip
from builtins import range
def MakeAtlasNode(atlasDirectory, name, atlasParts):
    """ Make an atlas node that contains the elements requested in the atlasParts section
        This will allow more fine grained data grabbers to be used, thereby allowing enhanced
        compartmentalization of algorithmic components.

        (S_) Static files that are relevant for any atlas
        (W_) Files that require warping to subjecgt specific atlas

        KEY:
          [S|W]_BRAINSABCSupport
          [S|W]_BRAINSABCSupport
          [S|W]_BRAINSCutSupport
          [S|W]_BCDSupport
          [S|W]_LabelMapsSupport
          [S|W]_ExtraSupport
    """

    import nipype.interfaces.io as nio  # Data i/o
    import nipype.pipeline.engine as pe  # pypeline engine
    import os

    valid_choices = [
        'S_BRAINSABCSupport',
        'S_BRAINSABCSupport',
        'S_BRAINSCutSupport',
        'S_BCDSupport',
        'S_LabelMapsSupport',
        'S_ExtraSupport',
        'W_BRAINSABCSupport',
        'W_BRAINSABCSupport',
        'W_BRAINSCutSupport',
        'W_BCDSupport',
        'W_LabelMapsSupport',
        'W_ExtraSupport'
    ]
    for ap in atlasParts:
        assert ap in valid_choices, "ERROR: Invalid choice: {0} not in {1}".format(ap, valid_choices)

    # Generate by running a file system list "ls -1 $AtlasDir *.nii.gz *.xml *.fcsv *.wgts"
    # atlas_file_names=atlas_file_list.split(' ')
    atlas_file_names = list()
    if 'S_BRAINSABCSupport' in atlasParts:
        atlas_file_names.extend([
            "ExtendedAtlasDefinition.xml.in"
        ])
    if 'W_BRAINSABCSupport' in atlasParts:
        atlas_file_names.extend([
            "template_headregion.nii.gz",
            "ExtendedAtlasDefinition.xml"
        ])
    if 'S_BRAINSCutSupport' in atlasParts:
        atlas_file_names.extend([
            "modelFiles/trainModelFile.txtD0060NT0060.gz"
        ])
    if 'W_BRAINSCutSupport' in atlasParts:
        atlas_file_names.extend([
            "hncma-atlas.nii.gz",
            "template_t1_denoised_gaussian.nii.gz",
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
        ])
    if 'S_BCDSupport' in atlasParts:
        atlas_file_names.extend([
            "20141004_BCD/LLSModel_50Lmks.h5",
            "20141004_BCD/T1_50Lmks.mdl",
            "20141004_BCD/template_weights_50Lmks.wts"
        ])
    if 'W_BCDSupport' in atlasParts:
        atlas_file_names.extend([
            "template_t1_denoised_gaussian.nii.gz",
            "20141004_BCD/template_landmarks_50Lmks.fcsv",
        ])
    if 'W_LabelMapsSupport' in atlasParts:
        atlas_file_names.extend([
            "hncma-atlas.nii.gz",
            "hncma-atlas-lut-mod2.ctbl",
            "template_rightHemisphere.nii.gz",
            "template_leftHemisphere.nii.gz",
            "template_WMPM2_labels.nii.gz",
            "template_WMPM2_labels.txt",
            "template_nac_labels.nii.gz",
            "template_nac_labels.txt",
            "template_ventricles.nii.gz",
            "template_headregion.nii.gz"
        ])
    if 'W_ExtraSupport' in atlasParts:
        atlas_file_names.extend([
            "tempNOTVBBOX.nii.gz",
            "template_ABC_labels.nii.gz",
            "avg_t1.nii.gz",
            "avg_t2.nii.gz",
            "template_brain.nii.gz",
            "template_cerebellum.nii.gz",
            "template_class.nii.gz",
            "template_headregion.nii.gz",
            "template_t1_denoised_gaussian.nii.gz",
            "template_t2_denoised_gaussian.nii.gz",
            "template_t1_clipped.nii.gz",
            "template_t2_clipped.nii.gz"
        ])
    atlas_file_names = list(set(atlas_file_names))  # Make a unique listing
    # # Remove filename extensions for images, but replace . with _ for other file types
    atlas_file_keys = [os.path.basename(fn).replace('.nii.gz', '').replace('.', '_').replace('-', '_') for fn in
                       atlas_file_names]
    atlas_outputs_filename_match = dict(list(zip(atlas_file_keys, atlas_file_names)))

    node = pe.Node(interface=nio.DataGrabber(force_output=False, outfields=atlas_file_keys),
                   run_without_submitting=True,
                   name=name)
    node.inputs.base_directory = atlasDirectory
    node.inputs.sort_filelist = False
    # node.inputs.raise_on_empty = True
    node.inputs.template = '*'
    ## Prefix every filename with atlasDirectory
    atlas_search_paths = ['{0}'.format(fn) for fn in atlas_file_names]
    node.inputs.field_template = dict(list(zip(atlas_file_keys, atlas_search_paths)))
    ## Give 'atlasDirectory' as the substitution argument
    atlas_template_args_match = [[[]] for i in atlas_file_keys]  # build a list of proper length with repeated entries
    node.inputs.template_args = dict(list(zip(atlas_file_keys, atlas_template_args_match)))
    # print "+" * 100
    # print node.inputs
    # print "-" * 100
    return node


def CreateAtlasXMLAndCleanedDeformedAverages(t1_image, deformed_list, AtlasTemplate, outDefinition):
    import os
    import sys
    import SimpleITK as sitk

    patternDict = {
        'AVG_WM.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_WM.nii.gz',
        'AVG_SURFGM.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_SURFGM.nii.gz',
        'AVG_BASAL.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_BASAL.nii.gz',
        'AVG_GLOBUS.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_GLOBUS.nii.gz',
        'AVG_THALAMUS.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_THALAMUS.nii.gz',
        'AVG_HIPPOCAMPUS.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_HIPPOCAMPUS.nii.gz',
        'AVG_CRBLGM.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_CRBLGM.nii.gz',
        'AVG_CRBLWM.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_CRBLWM.nii.gz',
        'AVG_CSF.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_CSF.nii.gz',
        'AVG_VB.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_VB.nii.gz',
        'AVG_NOTCSF.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_NOTCSF.nii.gz',
        'AVG_NOTGM.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_NOTGM.nii.gz',
        'AVG_NOTWM.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_NOTWM.nii.gz',
        'AVG_NOTVB.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_NOTVB.nii.gz',
        'AVG_AIR.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/GENERATED_AIR.nii.gz',
        'AVG_BRAINMASK.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_brain.nii.gz',
        'T1_RESHAPED.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_t1_denoised_gaussian.nii.gz',
        'AVG_T2.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_t2_denoised_gaussian.nii.gz',
        'AVG_PD.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_t2_denoised_gaussian.nii.gz',
        'AVG_FL.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_t2_denoised_gaussian.nii.gz',
        'AVG_hncma_atlas.nii.gz': 'IGNORED',
        'AVG_r_caudate_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_r_putamen_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_r_accumben_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_template_nac_labels.nii.gz': 'IGNORED',
        'AVG_l_hippocampus_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_theta.nii.gz': 'IGNORED',
        'AVG_l_accumben_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_phi.nii.gz': 'IGNORED',
        'AVG_l_thalamus_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_l_globus_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_template_ventricles.nii.gz': 'IGNORED',
        'AVG_template_headregion.nii.gz': 'IGNORED',
        'AVG_r_thalamus_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_l_putamen_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_rho.nii.gz': 'IGNORED',
        'AVG_r_hippocampus_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_r_globus_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_l_caudate_ProbabilityMap.nii.gz': 'IGNORED',
        'AVG_template_leftHemisphere.nii.gz': 'IGNORED',
        'AVG_template_WMPM2_labels.nii.gz': 'IGNORED',
        'AVG_template_rightHemisphere.nii.gz': 'IGNORED'
    }
    templateFile = open(AtlasTemplate, 'r')
    xmlAtlasFileContents = templateFile.read()  # read entire file into memory
    templateFile.close()

    # # Now clean up the posteriors based on anatomical knowlege.
    ## sometimes the posteriors are not relevant for priors
    ## due to anomolies around the edges.
    #print("\n\n\nALL_FILES: {0}\n\n\n".format(deformed_list))
    load_images_list = dict()
    for full_pathname in deformed_list:
        full_pathname=str(full_pathname)
        base_name = os.path.basename(full_pathname)
        if base_name in list(patternDict.keys()):
            load_images_list[base_name] = sitk.ReadImage(full_pathname)
        else:
            print("MISSING FILE FROM patternDict: {0}".format(base_name))
    ## Make binary dilated mask
    binmask = sitk.BinaryThreshold(load_images_list['AVG_BRAINMASK.nii.gz'], 1, 1000000)
    brainmask_dilatedBy5 = sitk.DilateObjectMorphology(binmask, 5)
    brainmask_dilatedBy5 = sitk.Cast(brainmask_dilatedBy5, sitk.sitkFloat32)  # Convert to Float32 for multiply

    inv_brainmask_erodedBy5 = 1 - sitk.ErodeObjectMorphology(binmask, 5)
    inv_brainmask_erodedBy5 = sitk.Cast(inv_brainmask_erodedBy5, sitk.sitkFloat32)  # Convert to Float32 for multiply

    ## Now clip the interior brain mask with brainmask_dilatedBy5
    interiorPriors = [
        'AVG_WM.nii.gz',
        'AVG_SURFGM.nii.gz',
        'AVG_BASAL.nii.gz',
        'AVG_CRBLGM.nii.gz',
        'AVG_CRBLWM.nii.gz',
        'AVG_CSF.nii.gz',
        'AVG_VB.nii.gz',
        'AVG_GLOBUS.nii.gz',
        'AVG_THALAMUS.nii.gz',
        'AVG_HIPPOCAMPUS.nii.gz',
    ]
    exteriorPriors = [
        'AVG_NOTWM.nii.gz',
        'AVG_NOTGM.nii.gz',
        'AVG_NOTCSF.nii.gz',
        'AVG_NOTVB.nii.gz',
        'AVG_AIR.nii.gz'
    ]
    extraFiles = [
        'AVG_hncma_atlas.nii.gz',
        'AVG_r_caudate_ProbabilityMap.nii.gz',
        'AVG_r_putamen_ProbabilityMap.nii.gz',
        'AVG_r_accumben_ProbabilityMap.nii.gz',
        'AVG_template_nac_labels.nii.gz',
        'AVG_l_hippocampus_ProbabilityMap.nii.gz',
        'AVG_theta.nii.gz',
        'AVG_l_accumben_ProbabilityMap.nii.gz',
        'AVG_phi.nii.gz',
        'AVG_l_thalamus_ProbabilityMap.nii.gz',
        'AVG_l_globus_ProbabilityMap.nii.gz',
        'AVG_template_ventricles.nii.gz',
        'AVG_template_headregion.nii.gz',
        'AVG_r_thalamus_ProbabilityMap.nii.gz',
        'AVG_l_putamen_ProbabilityMap.nii.gz',
        'AVG_rho.nii.gz',
        'AVG_r_hippocampus_ProbabilityMap.nii.gz',
        'AVG_r_globus_ProbabilityMap.nii.gz',
        'AVG_l_caudate_ProbabilityMap.nii.gz',
        'AVG_template_leftHemisphere.nii.gz',
        'AVG_template_WMPM2_labels.nii.gz',
        'AVG_template_rightHemisphere.nii.gz',
    ]
    clean_deformed_list = deformed_list
    T2File = None
    PDFile = None
    for index in range(0, len(deformed_list)):
        full_pathname = str(deformed_list[index])
        base_name = os.path.basename(full_pathname)
        if base_name == 'AVG_BRAINMASK.nii.gz':
            ### Make Brain Mask Binary
            clipped_name = 'CLIPPED_' + base_name
            patternDict[clipped_name] = patternDict[base_name]
            sitk.WriteImage(binmask, clipped_name)
            clean_deformed_list[index] = os.path.realpath(clipped_name)
        elif base_name == 'AVG_T2.nii.gz':
            T2File = full_pathname
        elif base_name == 'AVG_PD.nii.gz':
            PDFile = full_pathname
        elif base_name in interiorPriors:
            ### Make clipped posteriors for brain regions
            curr = sitk.Cast(sitk.ReadImage(full_pathname), sitk.sitkFloat32)
            curr = curr * brainmask_dilatedBy5
            clipped_name = 'CLIPPED_' + base_name
            patternDict[clipped_name] = patternDict[base_name]
            sitk.WriteImage(curr, clipped_name)
            clean_deformed_list[index] = os.path.realpath(clipped_name)
            #print "HACK: ", clean_deformed_list[index]
            curr = None
        elif base_name in exteriorPriors:
            ### Make clipped posteriors for brain regions
            curr = sitk.Cast(sitk.ReadImage(full_pathname), sitk.sitkFloat32)
            curr = curr * inv_brainmask_erodedBy5
            clipped_name = 'CLIPPED_' + base_name
            patternDict[clipped_name] = patternDict[base_name]
            sitk.WriteImage(curr, clipped_name)
            clean_deformed_list[index] = os.path.realpath(clipped_name)
            #print "HACK: ", clean_deformed_list[index]
            curr = None
        elif base_name in extraFiles:
            pass
        else:
            import sys

            print("ERROR: basename {0} not in list!! \n{1}".format(base_name, ['AVG_BRAINMASK.nii.gz', 'AVG_T2.nii.gz',
                                                                               'AVG_PD.nii.gz', interiorPriors,
                                                                               exteriorPriors]))
            sys.exit(-1)

    binmask = None
    brainmask_dilatedBy5 = None
    inv_brainmask_erodedBy5 = None

    for full_pathname in clean_deformed_list:
        base_name = os.path.basename(full_pathname)
        if base_name in list(patternDict.keys()):
            xmlAtlasFileContents = xmlAtlasFileContents.replace(patternDict[base_name], base_name)
    ## If there is no T2, then use the PD image
    if T2File is not None:
        xmlAtlasFileContents = xmlAtlasFileContents.replace('@ATLAS_INSTALL_DIRECTORY@/template_t2_denoised_gaussian.nii.gz',
                                                            os.path.basename(T2File))
    elif PDFile is not None:
        xmlAtlasFileContents = xmlAtlasFileContents.replace('@ATLAS_INSTALL_DIRECTORY@/template_t2_denoised_gaussian.nii.gz',
                                                            os.path.basename(PDFile))
    xmlAtlasFileContents = xmlAtlasFileContents.replace('@ATLAS_INSTALL_DIRECTORY@/template_t1_denoised_gaussian.nii.gz', 'AVG_T1.nii.gz')
    ## NOTE:  HEAD REGION CAN JUST BE T1 image.
    xmlAtlasFileContents = xmlAtlasFileContents.replace('@ATLAS_INSTALL_DIRECTORY@/template_headregion.nii.gz',
                                                        os.path.basename(t1_image))
    ## NOTE:  BRAIN REGION CAN JUST BE the label images.
    outAtlasFullPath = os.path.realpath(outDefinition)
    newFile = open(outAtlasFullPath, 'w')
    newFile.write(xmlAtlasFileContents)  # write the file with the text substitution
    newFile.close()
    return outAtlasFullPath, clean_deformed_list
