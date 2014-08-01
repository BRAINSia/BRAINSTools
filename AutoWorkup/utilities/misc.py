FS_VARS = ['FREESURFER_HOME',
           'FSFAST_HOME',
           'FSF_OUTPUT_FORMAT',
           'SUBJECTS_DIR',
           'MNI_DIR',
           'FSL_DIR']

def MakeOutFileList(T1List, T2List, PDList, FLList, OtherList, postfix, PrimaryT1):
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
    all_files.extend(OtherList)
    outImageList = []
    for i in all_files:
        out_name = GetExtBaseName(i) + postfix
        outImageList.append(out_name)
    #
    #make type list
    imageTypeList = ["T1"] * len(T1List)
    imageTypeList.extend(["T2"] * len(T2List))
    imageTypeList.extend(["PD"] * len(PDList))
    imageTypeList.extend(["FL"] * len(FLList))
    imageTypeList.extend(["OTHER"] * len(OtherList))

    #make input raw images single list
    inImageList = list()
    inImageList.extend(T1List)
    inImageList.extend(T2List)
    inImageList.extend(PDList)
    inImageList.extend(FLList)
    inImageList.extend(OtherList)

    """ This funciton uses PrimaryT1 for the first T1, and the append the rest of the T1's and T2's """
    if PrimaryT1 is not None:
        inImageList[0]=PrimaryT1

    return inImageList, outImageList, imageTypeList

def GenerateSeparateImageTypeList(inFileList, inTypeList):
    allListDict = dict()
    allListDict["T1"]=list()
    allListDict["T2"]=list()
    allListDict["PD"]=list()
    allListDict["FL"]=list()
    allListDict["Other"]=list()
    T1List=list()
    for i in range(0,len(inFileList)):
        allListDict[ inTypeList[i] ].append( inFileList[i] )

    return allListDict["T1"], allListDict["T2"], allListDict["PD"], allListDict["FL"], allListDict["Other"]


def add_dict(d1, d2, force=False):
    from copy import deepcopy
    retval = deepcopy(d1)
    if not force:
        try:
            assert set(d1.keys()).isdisjoint(set(d2.keys()))
        except AssertionError:
            raise ValueError("Dictionaries have one or more duplicate keys")
    for key in d2.keys():
        if key in retval.keys() and force:
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
