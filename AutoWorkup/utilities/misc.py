from __future__ import print_function
from builtins import str
from builtins import range
FS_VARS = ['FREESURFER_HOME',
           'FSFAST_HOME',
           'FSF_OUTPUT_FORMAT',
           'SUBJECTS_DIR',
           'MNI_DIR',
           'FSL_DIR']

def CommonANTsRegistrationSettings(antsRegistrationNode,
            registrationTypeDescription,
            output_transform_prefix,
            output_warped_image,
            output_inverse_warped_image,
            save_state,
            invert_initial_moving_transform = None):
    """ Ants registration settings are difficult
    to get correct all the time.  This utility function
    is designed to assist with getting the common settings
    correct across different registration units.

    By placing these settings in one location, we can
    review all antsRegistration settings for the
    entire pipeline in one spot.
    """
    ## TODO: Consider registration masking

    if ( registrationTypeDescription == "5StagesSingleModal" ) or ( "5StagesMultiModal" ):
        local_num_stages=5
        if registrationTypeDescription == "5StagesSingleModal":
            antsRegistrationNode.inputs.metric = ['MI','MI','CC','CC','CC']
            antsRegistrationNode.inputs.metric_weight = [1.0,1.0,1.0,1.0,1.0]
            antsRegistrationNode.inputs.sampling_strategy = ['Regular','Regular',None,None,None]
            antsRegistrationNode.inputs.sampling_percentage = [0.27,0.27,1.0,1.0,1.0]
            antsRegistrationNode.inputs.radius_or_number_of_bins = [32,32,4,4,4]
        else:
            antsRegistrationNode.inputs.metric = ['MI',['MI','MI'],'CC','CC',['CC','CC']]
            antsRegistrationNode.inputs.metric_weight = [1.0,[1.0,1.0],1.0,1.0,[1.0,1.0]]
            antsRegistrationNode.inputs.sampling_strategy = ['Regular',['Regular','Regular'],None,None,[None,None]]
            antsRegistrationNode.inputs.sampling_percentage = [0.27,[0.27,0.27],1.0,1.0,[1.0,1.0]]
            antsRegistrationNode.inputs.radius_or_number_of_bins = [32,[32,32],4,4,[4,4]]

        antsRegistrationNode.inputs.transforms = ["Rigid","Affine","SyN","SyN","SyN"]
        antsRegistrationNode.inputs.transform_parameters = [[0.1],[0.1],[0.1,3.0,0.0],[0.1,3.0,0.0],[0.1,3.0,0.0]]
        antsRegistrationNode.inputs.number_of_iterations = [[1000,1000,1000,1000],[1000,1000,1000,1000],[1000,250],[140],[25]]
        antsRegistrationNode.inputs.convergence_threshold = [5e-8,5e-8,5e-7,5e-6,5e-5]
        antsRegistrationNode.inputs.shrink_factors =   [[8,4,2,1],[8,4,2,1],[8,4],[2],[1]]
        antsRegistrationNode.inputs.smoothing_sigmas = [[3,2,1,0],[3,2,1,0],[3,2],[1],[0]]
        #
        antsRegistrationNode.inputs.use_estimate_learning_rate_once = [False]*local_num_stages
        antsRegistrationNode.inputs.winsorize_lower_quantile = 0.01
        antsRegistrationNode.inputs.winsorize_upper_quantile = 0.99

    elif ( registrationTypeDescription == "JointFusionT1Only" ) or ( registrationTypeDescription == "JointFusionMultiModal" ):
        local_num_stages=4
        if registrationTypeDescription == "JointFusionT1Only":
            antsRegistrationNode.inputs.metric = ['MI','MI','CC','CC']
            antsRegistrationNode.inputs.metric_weight = [1.0,1.0,1.0,1.0]
            antsRegistrationNode.inputs.sampling_strategy = ['Regular','Regular',None,None]
            antsRegistrationNode.inputs.sampling_percentage = [.5,.5,1.0,1.0]
            antsRegistrationNode.inputs.radius_or_number_of_bins = [32,32,4,4]
        else:
            antsRegistrationNode.inputs.metric = ['MI',['MI','MI'],'CC',['CC','CC']]
            antsRegistrationNode.inputs.metric_weight = [1.0,[1.0,1.0],1.0,[1.0,1.0]]
            antsRegistrationNode.inputs.sampling_strategy = ['Regular',['Regular','Regular'],None,[None,None]]
            antsRegistrationNode.inputs.sampling_percentage = [.5,[.5,0.5],1.0,[1.0,1.0]]
            antsRegistrationNode.inputs.radius_or_number_of_bins = [32,[32,32],4,[4,4]]

        antsRegistrationNode.inputs.transforms = ["Affine","Affine","SyN","SyN"]
        antsRegistrationNode.inputs.transform_parameters = [[0.1],[0.1],[0.1, 3, 0],[0.1, 3, 0]]
        antsRegistrationNode.inputs.number_of_iterations = [[1000,1000,500],[500,500],[500,500],[500,70]]
        antsRegistrationNode.inputs.convergence_threshold = [1e-8,1e-6,1e-8,1e-6]
        antsRegistrationNode.inputs.shrink_factors = [[8, 4, 2],[2, 1],[8, 4],[2, 1]]
        antsRegistrationNode.inputs.smoothing_sigmas = [[3, 2, 1],[1, 0],[3, 2],[1, 0]]

    #if registrationTypeDescription == "FromWorkupT1T2ANTS.py":
    #    local_num_stages=1
    #    ## This is outdated and does not work
    #    antsRegistrationNode.inputs.metric = 'CC'  # This is a family of interfaces, CC,MeanSquares,Demons,GC,MI,Mattes
    #    antsRegistrationNode.inputs.transform = 'SyN[0.25,3.0,0.0]'
    #    antsRegistrationNode.inputs.number_of_iterations = [250, 100, 20]
    #    antsRegistrationNode.inputs.convergence_threshold = 1e-7
    #    antsRegistrationNode.inputs.smoothing_sigmas = [0, 0, 0]
    #    antsRegistrationNode.inputs.shrink_factors = [3, 2, 1]

    elif registrationTypeDescription == 'AtlasToSubjectANTsPreABC_Rigid':
        local_num_stages=1
        antsRegistrationNode.inputs.transforms = ["Affine",]
        antsRegistrationNode.inputs.transform_parameters = [[0.1]]
        antsRegistrationNode.inputs.metric = ['MI']
        antsRegistrationNode.inputs.sampling_strategy = ['Regular']
        antsRegistrationNode.inputs.sampling_percentage = [0.5]
        antsRegistrationNode.inputs.metric_weight = [1.0]
        antsRegistrationNode.inputs.radius_or_number_of_bins = [32]
        antsRegistrationNode.inputs.number_of_iterations = [[1000,1000, 500, 100]]

        antsRegistrationNode.inputs.convergence_threshold = [1e-8]

        antsRegistrationNode.inputs.shrink_factors = [[8, 4, 2, 1]]
        antsRegistrationNode.inputs.smoothing_sigmas = [[3, 2, 1, 0]]

    elif registrationTypeDescription == 'AtlasToSubjectANTsPreABC_SyN':
        local_num_stages = 2
        antsRegistrationNode.inputs.transforms = ["SyN","SyN"]
        antsRegistrationNode.inputs.transform_parameters = [[0.1, 3, 0],[0.1, 3, 0]]
        antsRegistrationNode.inputs.metric = ['CC','CC']
        antsRegistrationNode.inputs.sampling_strategy = [None,None]
        antsRegistrationNode.inputs.sampling_percentage = [1.0,1.0]
        antsRegistrationNode.inputs.metric_weight = [1.0,1.0]
        antsRegistrationNode.inputs.radius_or_number_of_bins = [4,4]
        antsRegistrationNode.inputs.number_of_iterations = [[500, 500], [500, 70]]

        antsRegistrationNode.inputs.convergence_threshold = [1e-8,1e-6]

        antsRegistrationNode.inputs.shrink_factors = [[8, 4], [2, 1]]
        antsRegistrationNode.inputs.smoothing_sigmas = [[3, 2], [1, 0]]

    elif registrationTypeDescription == "antsRegistrationNode":
        local_num_stages = 1
        antsRegistrationNode.inputs.transforms = ["SyN"]
        antsRegistrationNode.inputs.transform_parameters = [[0.1, 3, 0]]
        antsRegistrationNode.inputs.metric = ['CC']
        antsRegistrationNode.inputs.sampling_strategy = [None]
        antsRegistrationNode.inputs.sampling_percentage = [1.0]
        antsRegistrationNode.inputs.metric_weight = [1.0]
        antsRegistrationNode.inputs.radius_or_number_of_bins = [4]
        antsRegistrationNode.inputs.number_of_iterations = [[70]]

        antsRegistrationNode.inputs.convergence_threshold = [1e-6]

        antsRegistrationNode.inputs.shrink_factors = [[1]]
        antsRegistrationNode.inputs.smoothing_sigmas = [[0]]
    else:
        print("!!"*160 + "ERROR invalid registration description")
        raise NameError(registrationTypeDescription)

    ## COMMON SETTINGS
    antsRegistrationNode.inputs.interpolation = "Linear"
    antsRegistrationNode.inputs.dimension = 3
    antsRegistrationNode.inputs.float = True
    antsRegistrationNode.inputs.num_threads = -1

    antsRegistrationNode.inputs.write_composite_transform = True # Required for initialize_transforms_per_stage
    antsRegistrationNode.inputs.collapse_output_transforms = False # Mutually Exclusive with initialize_transforms_per_stage
    antsRegistrationNode.inputs.initialize_transforms_per_stage = True

    antsRegistrationNode.inputs.convergence_window_size = [12]*local_num_stages
    antsRegistrationNode.inputs.sigma_units = ["vox"]*local_num_stages
    antsRegistrationNode.inputs.use_histogram_matching = [True]*local_num_stages
    #
    antsRegistrationNode.inputs.use_estimate_learning_rate_once = [False]*local_num_stages
    antsRegistrationNode.inputs.winsorize_lower_quantile = 0.01
    antsRegistrationNode.inputs.winsorize_upper_quantile = 0.99

    if invert_initial_moving_transform is not None:
        antsRegistrationNode.inputs.invert_initial_moving_transform = invert_initial_moving_transform
    if output_transform_prefix is not None:
        antsRegistrationNode.inputs.output_transform_prefix = output_transform_prefix
    if output_warped_image is not None:
        antsRegistrationNode.inputs.output_warped_image = output_warped_image
    if output_inverse_warped_image is not None:
        antsRegistrationNode.inputs.output_inverse_warped_image = output_inverse_warped_image

def MakeOutFileList(T1List, T2List, PDList, FLList, OTHERList, postfix, PrimaryT1, ListOutType=False):
    #
    #for BABC: "_corrected.nii.gz"
    #for UNM Denoise: "_UNM_denoised.nii.gz"

    #
    #make image file list
    def GetExtBaseName(filename):
        '''
        Get the filename without the extension.  Works for .ext and .ext.gz
        '''
        import os
        currBaseName = os.path.basename(filename)
        currExt = os.path.splitext(currBaseName)[1]
        currBaseName = os.path.splitext(currBaseName)[0]
        if currExt == ".gz":
            currBaseName = os.path.splitext(currBaseName)[0]
            currExt = os.path.splitext(currBaseName)[1]
        return currBaseName

    all_files = list()
    all_files.extend(T1List)
    all_files.extend(T2List)
    all_files.extend(PDList)
    all_files.extend(FLList)
    all_files.extend(OTHERList)
    outImageList = []
    for i in all_files:
        out_name = GetExtBaseName(i) + postfix
        if ListOutType:
            out_name = [str(out_name)]
        else:
            out_name = str(out_name)
        outImageList.append( out_name )
    #
    #make type list
    imageTypeList = ["T1"] * len(T1List)
    imageTypeList.extend(["T2"] * len(T2List))
    imageTypeList.extend(["PD"] * len(PDList))
    imageTypeList.extend(["FL"] * len(FLList))
    imageTypeList.extend(["OTHER"] * len(OTHERList))

    #make input raw images single list
    inImageList = list()
    inImageList.extend(T1List)
    inImageList.extend(T2List)
    inImageList.extend(PDList)
    inImageList.extend(FLList)
    inImageList.extend(OTHERList)

    """ This function uses PrimaryT1 for the first T1, and the append the rest of the T1's and T2's """
    if PrimaryT1 is not None:
        inImageList[0]=PrimaryT1
    print ("inImageList:::")
    print (inImageList)
    print ("outImageList:::")
    print (outImageList)
    print ("imageTypeList:::")
    print (imageTypeList)
    return inImageList, outImageList, imageTypeList

def GenerateSeparateImageTypeList(inFileList, inTypeList):
    allListDict = dict()
    allListDict["T1"]=list()
    allListDict["T2"]=list()
    allListDict["PD"]=list()
    allListDict["FL"]=list()
    allListDict["OTHER"]=list()
    T1List=list()
    for i in range(0,len(inFileList)):
        allListDict[ inTypeList[i] ].append( inFileList[i] )

    return allListDict["T1"], allListDict["T2"], allListDict["PD"], allListDict["FL"], allListDict["OTHER"]


def add_dict(d1, d2, force=False):
    from copy import deepcopy
    retval = deepcopy(d1)
    if d2:
        if not force:
            try:
                print("d1.keys():::")
                print(list(d1.keys()))
                print("d2.keys():::")
                print(list(d2.keys()))
                assert set(d1.keys()).isdisjoint(set(d2.keys()))
            except AssertionError:
                raise ValueError("Dictionaries have one or more duplicate keys")
        for key in list(d2.keys()):
            if key in list(retval.keys()) and force:
                try:
                    retval[key] += d2[key]
                except:
                    raise
            else:
                retval[key] = deepcopy(d2[key])
    return retval


def GenerateWFName(projectid, subjectid, sessionid, processing_phase):
    return 'WF_' + str(subjectid) + "_" + str(sessionid) + "_" + str(projectid) + "_" + processing_phase


def GenerateSubjectOutputPattern(subjectid):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard """
    import os.path

    patternList = []
    find_pat = "_subject_" + subjectid + "/"
    replace_pat = ""
    patternList.append((find_pat,replace_pat))

    find_pat = 'ReshapeAverageImageWithShapeUpdate.nii.gz'
    replace_pat = r'AVG_T1.nii.gz'
    patternList.append((find_pat, replace_pat))

    # find_pat = os.path.join('Atlas',
    #                         r'_ReshapeAveragePassiveImageWithShapeUpdate[0-9]*/AVG_[A-Z0-9]*WARP_(?P<structure>AVG_[A-Z0-9]*.nii.gz)')
    find_pat = r'_ReshapeAveragePassiveImageWithShapeUpdate[0-9]*/AVG_(?P<structure>.*.nii.gz)'
    replace_pat = r'AVG_\g<structure>'
    patternList.append((find_pat, replace_pat))

    #find_pat = r'CLIPPED_AVG_(?P<structure>.*.nii.gz)'
    #replace_pat = r'AVG_\g<structure>'
    #patternList.append((find_pat, replace_pat))

    #print "HACK: ", patternList
    return patternList


def GenerateOutputPattern(projectid, subjectid, sessionid, DefaultNodeName):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    patternList = []
    find_pat = os.path.join(DefaultNodeName)
    replace_pat = os.path.join(projectid, subjectid, sessionid, DefaultNodeName)
    patternList.append((find_pat, replace_pat))
    #print "HACK: ", patternList
    return patternList
