FS_VARS = ['FREESURFER_HOME',
           'FSFAST_HOME',
           'FSF_OUTPUT_FORMAT',
           'SUBJECTS_DIR',
           'MNI_DIR',
           'FSL_DIR']

WFRUN = ['helium_all.q',
         'helium_all.q_graph',
         'ipl_OSX',
         'local_4',
         'local_12',
         'local',
         'ds_runner']


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

class DS_runner(object):
    def run(self, graph, **kwargs):
        for node in graph.nodes():
            if '_ds' in node.name.lower():
                node.run()


def get_cpus(option):
    assert option in WFRUN, "Unknown wfrun option"
    from multiprocessing import cpu_count
    import os
    total_cpus = cpu_count()
    suffix = option.rsplit('local', 1)[1]
    if suffix == '':
        assert option in ['local', 'ds_runner'], "wfrun parse error!  Current option: {0}".format(option)
        threads = 1
        if option == 'local':
            print "RUNNING WITHOUT POOL BUILDING"
    else:
        threads = int(suffix.strip('_'))
    return int(total_cpus / threads)


def nipype_plugin(wfrun):
    assert wfrun in WFRUN, "Unknown workflow run environment: {0}".format(wfrun)
    sge_flavor = 'SGEGraph'  # NOT 'SGE' per Satra 2014/07/18
    if wfrun in ['helium_all.q', 'helium_all.q_graph', 'ipl_OSX']:
        return sge_flavor
    elif wfrun in ['local_4', 'local_12']:
        return 'MultiProc'
    elif wfrun == 'ds_runner':
        return DS_runner()
    else:
        assert wfrun == 'local', "You must specify a valid run environment type.  Invalid: {0}".format(wfrun)
        return 'Linear'
    return None


def nipype_plugin_args(wfrun, cluster, template=''):
    assert wfrun in WFRUN, "Unknown workflow run environment: {0}".format(wfrun)
    qsub_args = "-S /bin/bash -cwd -pe smp 1- -l h_vmem=19G,mem_free=9G -o /dev/null -e /dev/null {0}"
    if wfrun in ['helium_all.q', 'helium_all.q_graph', 'ipl_OSX']:
        return {'template': template,
                'qsub_args': qsub_args.format(cluster['queue']),
                'qstatProgramPath': cluster['qstat'],
                'qstatCachedProgramPath': cluster['qstat_cached']}
    elif wfrun in ['local_4', 'local_12']:
        proc_count = int(wfrun.split('local_')[1])
        print "Running with {0} parallel processes on local machine".format(proc_count)
        return {'n_procs': proc_count}
    else:
        assert wfrun in ['local', 'ds_runner'], "You must specify a valid run environment type.  Invalid: {0}".format(wfrun)
    return None


def nipype_execution(plugin='Linear', stop_on_first_crash=False, stop_on_first_rerun=False):
    stop_crash = 'false'
    stop_rerun = 'false'
    if stop_on_first_crash:
        stop_crash = 'true'
    if stop_on_first_rerun:
        # This stops at first attempt to rerun, before running, and before deleting previous results
        stop_rerun = 'true'
    return {'plugin': nipype_plugin(plugin),
            'stop_on_first_crash': stop_crash,
            'stop_on_first_rerun': stop_rerun,
            'hash_method': 'timestamp',          # default
            'single_thread_matlab': 'true',      # default # Multi-core 2011a  multi-core for matrix multiplication.
            'use_relative_paths': 'false',       # default # relative paths should be on, require hash update when changed.
            'remove_node_directories': 'false',  # default
            'remove_unnecessary_outputs': 'false',
            'local_hash_check': 'true',          # default
            'job_finished_timeout': 25}


def nipype_logging(cachedir):
    return {'workflow_level': 'INFO', # possible options:
            'filemanip_level': 'INFO',#   INFO (default) | DEBUG
            'interface_level': 'INFO',
            'log_directory': cachedir}


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
