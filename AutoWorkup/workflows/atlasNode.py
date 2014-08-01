def MakeAtlasNode(atlasDirectory, name):
    import nipype.interfaces.io as nio   # Data i/o
    import nipype.pipeline.engine as pe  # pypeline engine

    from utilities import atlas_file_names, atlas_file_keys, atlas_outputs_filename_match

    node = pe.Node(interface=nio.DataGrabber(force_output=False, outfields=atlas_file_keys),
                     run_without_submitting=True,
                     name=name)
    node.inputs.base_directory = atlasDirectory
    node.inputs.sort_filelist = False
    node.inputs.template = '*'
    ## Prefix every filename with atlasDirectory
    atlas_search_paths = ['{0}'.format(fn) for fn in atlas_file_names]
    node.inputs.field_template = dict(zip(atlas_file_keys, atlas_search_paths))
    ## Give 'atlasDirectory' as the substitution argument
    atlas_template_args_match = [ [[]] for i in atlas_file_keys]  # build a list of proper length with repeated entries
    node.inputs.template_args = dict(zip(atlas_file_keys, atlas_template_args_match))
    # print "+" * 100
    # print node.inputs
    # print "-" * 100
    return node


def GetAtlasNode(previousresult, name):
    """ Guarantee that template experiment Atlas matches baseline experiment Atlas """
    import os.path
    from atlasNode import MakeAtlasNode

    previousAtlasDir = os.path.abspath(os.path.join(previousresult, 'Atlas'))
    assert os.path.exists(previousAtlasDir), "Previous experiment's Atlas directory cannot be found! {0}".format(previousAtlasDir)
    return MakeAtlasNode(previousAtlasDir, name)


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
        'T1_RESHAPED.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_t1.nii.gz',
        'AVG_T2.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_t2.nii.gz',
        'AVG_PD.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_t2.nii.gz',
        'AVG_FL.nii.gz': '@ATLAS_INSTALL_DIRECTORY@/template_t2.nii.gz'
    }
    templateFile = open(AtlasTemplate, 'r')
    content = templateFile.read()              # read entire file into memory
    templateFile.close()

    ## Now clean up the posteriors based on anatomical knowlege.
    ## sometimes the posteriors are not relevant for priors
    ## due to anomolies around the edges.
    #print("\n\n\nALL_FILES: {0}\n\n\n".format(deformed_list))
    load_images_list = dict()
    for full_pathname in deformed_list:
        base_name = os.path.basename(full_pathname)
        if base_name in patternDict.keys():
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
    clean_deformed_list = deformed_list
    T2File = None
    PDFile = None
    for index in range(0, len(deformed_list)):
        full_pathname = deformed_list[index]
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
        else:
            import sys
            print "ERROR: basename {0} not in list!! \n{1}".format(base_name,['AVG_BRAINMASK.nii.gz','AVG_T2.nii.gz','AVG_PD.nii.gz',interiorPriors,exteriorPriors])
            sys.exit(-1)

    binmask = None
    brainmask_dilatedBy5 = None
    inv_brainmask_erodedBy5 = None

    for full_pathname in clean_deformed_list:
        base_name = os.path.basename(full_pathname)
        if base_name in patternDict.keys():
            content = content.replace(patternDict[base_name], full_pathname)
    ## If there is no T2, then use the PD image
    if T2File is not None:
        content = content.replace('@ATLAS_INSTALL_DIRECTORY@/template_t2.nii.gz', T2File)
    elif PDFile is not None:
        content = content.replace('@ATLAS_INSTALL_DIRECTORY@/template_t2.nii.gz', PDFile)
    content = content.replace('@ATLAS_INSTALL_DIRECTORY@/template_t1.nii.gz', t1_image)
    ## NOTE:  HEAD REGION CAN JUST BE T1 image.
    content = content.replace('@ATLAS_INSTALL_DIRECTORY@/template_headregion.nii.gz', t1_image)
    ## NOTE:  BRAIN REGION CAN JUST BE the label images.
    outAtlasFullPath = os.path.realpath(outDefinition)
    newFile = open(outAtlasFullPath, 'w')
    newFile.write(content)  # write the file with the text substitution
    newFile.close()
    return outAtlasFullPath, clean_deformed_list
