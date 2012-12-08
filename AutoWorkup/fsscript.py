import argparse
import os
import subprocess
import errno
import SimpleITK as sitk

def normalizeWM(t1, wm_prob):
    """This function will compute the mean value of wm and rescale the image so that the mean WM=100"""
    WM_MEAN_FINAL = 110.0
    WM_THRESHOLD = 0.66
    wm_mask = sitk.BinaryThreshold(wm_prob, WM_THRESHOLD)
    ls = sitk.LabelStatisticsImageFilter()
    ls.Execute(t1, wm_mask)
    wm_value = 1
    myMeasurementMap = ls.GetMeasurementMap(wm_value)
    MeanWM = myMeasurementMap['Mean']
    t1_new = sitk.Cast(sitk.Cast(t1, sitk.sitkFloat32) * WM_MEAN_FINAL / MeanWM, sitk.sitkUInt8)
    return t1_new

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def make_dummy_file(fn):
    """This function just makes a file with the correct name and time stamp"""
    import time
    mkdir_p(os.path.dirname(fn))
    ff=open(fn,'w')
    ff.write("DummyFile with Proper time stamp")
    time.sleep(1) # 1 second
    ff.close()

def IsFirstNewerThanSecond(firstFile, secondFile):
    if not os.path.exists(firstFile):
        print "ERROR: image missing", firstFile
        return True
    image_time=os.path.getmtime(firstFile)
    if not os.path.exists(secondFile):
        return True
    reference_time=os.path.getmtime(secondFile)
    return image_time > reference_time


def run_mri_convert_script(niftiVol,mgzVol,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT):
    FS_SCRIPT_FN=os.path.join(FREESURFER_HOME,FS_SCRIPT)
    mri_convert_script="""#!/bin/bash
    export FREESURFER_HOME={FSHOME}
    export SUBJECTS_DIR={FSSUBJDIR}
    source {SOURCE_SCRIPT}
    {FSHOME}/bin/mri_convert --out_data_type uchar {invol} {outvol}
""".format(SOURCE_SCRIPT=FS_SCRIPT_FN,FSHOME=FREESURFER_HOME,FSSUBJDIR=FS_SUBJECTS_DIR,invol=niftiVol,outvol=mgzVol)
    script_name=mgzVol+'_convert.sh'
    script=open(script_name,'w')
    script.write(mri_convert_script)
    script.close()
    os.chmod(script_name,0777)
    script_name_stdout=mgzVol+'_convert.out'
    script_name_stdout_fid=open(script_name_stdout,'w')
    print "Starting mri_convert"
    subprocess.check_call([script_name], stdout=script_name_stdout_fid, stderr=subprocess.STDOUT, shell='/bin/bash')
    print "Ending mri_convert"
    script_name_stdout_fid.close()
    return

def baw_Recon1(t1_fn,wm_fn,brain_fn,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT,subjID):
    base_subj_dir=os.path.join(FS_SUBJECTS_DIR,subjID,'mri')
    output_brainmask_fn=os.path.join(base_subj_dir,'brainmask.nii.gz')
    output_nu_fn=os.path.join(base_subj_dir,'nu.nii.gz')
    output_brainmask_fn_mgz=os.path.join(base_subj_dir,'brainmask.mgz')
    output_nu_fn_mgz=os.path.join(base_subj_dir,'nu.mgz')
    if IsFirstNewerThanSecond(t1_fn,output_brainmask_fn_mgz):
        print "PREPARING ALTERNATE recon-auto1 stage"
        mkdir_p(base_subj_dir)
        make_dummy_file(os.path.join(base_subj_dir,'orig/001.mgz') )
        make_dummy_file(os.path.join(base_subj_dir,'rawavg.mgz') )
        make_dummy_file(os.path.join(base_subj_dir,'orig.mgz') )
        make_dummy_file(os.path.join(base_subj_dir,'transforms/talairach.auto.xfm') )
        make_dummy_file(os.path.join(base_subj_dir,'transforms/talairach.xfm') )

        t1 = sitk.ReadImage(t1_fn)
        wm = sitk.ReadImage(wm_fn)
        t1_new = sitk.Cast(normalizeWM(t1, wm),sitk.sitkUInt8)
        sitk.WriteImage(t1_new,output_nu_fn)
        run_mri_convert_script(output_nu_fn,output_nu_fn_mgz,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT)

        make_dummy_file(os.path.join(base_subj_dir,'T1.mgz') )
        make_dummy_file(os.path.join(base_subj_dir,'transforms/talairach_with_skull.lta') )
        make_dummy_file(os.path.join(base_subj_dir,'brainmask.auto.mgz') )
        brain = sitk.ReadImage(brain_fn)
        blood = sitk.BinaryThreshold(brain, 5, 5)
        not_blood = 1 - blood
        clipping = sitk.BinaryThreshold(brain, 1, 1000000) - blood
        fill_size=2
        ## HACK: Unfortunately we need to hole fill because of a WM bug in BABC where
        ## some white matter is being classified as background when it it being avoid due
        ## to too strict of multi-modal thresholding.
        hole_filled=sitk.ErodeObjectMorphology(sitk.DilateObjectMorphology(clipping,fill_size),fill_size)
        clipped = sitk.Cast(t1_new * hole_filled * not_blood,sitk.sitkUInt8)
        sitk.WriteImage(clipped,output_brainmask_fn) ## brain_matter image with values normalized 0-110, no skull or surface blood
        run_mri_convert_script(output_brainmask_fn,output_brainmask_fn_mgz,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT)
    else:
        print "NOTHING TO BE DONE, SO SKIPPING."
        return # Nothing to be done, files are already up-to-date.

def runAutoReconStage(subjID,StageToRun,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT):
    FS_SCRIPT_FN=os.path.join(FREESURFER_HOME,FS_SCRIPT)
    auto_recon_script="""#!/bin/bash
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}
source {SOURCE_SCRIPT}
{FSHOME}/bin/recon-all -debug -subjid {SUBJID} -make autorecon{AUTORECONSTAGE}
""".format(SOURCE_SCRIPT=FS_SCRIPT_FN,FSHOME=FREESURFER_HOME,FSSUBJDIR=FS_SUBJECTS_DIR,AUTORECONSTAGE=StageToRun,SUBJID=subjID)
    script_name=os.path.join(FS_SUBJECTS_DIR,subjID,'run_autorecon_stage'+str(StageToRun)+'.sh')
    script=open(script_name,'w')
    script.write(auto_recon_script)
    script.close()
    os.chmod(script_name,0777)
    script_name_stdout=script_name+'_out'
    script_name_stdout_fid=open(script_name_stdout,'w')
    print "Starting auto_recon Stage: {0} for SubjectSession {1}".format(StageToRun,subjID)
    subprocess.check_call([script_name], stdout=script_name_stdout_fid, stderr=subprocess.STDOUT, shell='/bin/bash')
    print "Ending auto_recon Stage: {0} for SubjectSession {1}".format(StageToRun,subjID)
    script_name_stdout_fid.close()
    return

def runAutoRecon(t1_fn,wm_fn,brain_fn,subjID,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT):
    """Run all stages of AutoRecon For Freesurfer, including the custom BAW initialization."""
    baw_Recon1(t1_fn,wm_fn,brain_fn,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT,subjID)
    runAutoReconStage(subjID,2,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT)
    runAutoReconStage(subjID,3,FREESURFER_HOME,FS_SUBJECTS_DIR,FS_SCRIPT)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    DELETE LATER: This is an just example of the commands required to run FreeSurfer recon-all:
    /bin/bash
    setenv FREESURFER_HOME /ipldev/sharedopt/20110601/MacOSX_10.6/freesurfer
    source ${FREESURFER_HOME}/FreeSurferEnv.sh
    setenv SUBJECTS_DIR /IPLlinux/raid0/homes/jforbes/freesurfer/recon-all/autorecon1_copy
    recon-all -make all -subjid 0074_24832

    Link to recon-all i/o table:
    http://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllDevTable

    """)
    # TODO: Make parser group "Inputs"
    parser.add_argument('--T1image', action='store', dest='t1Image', help='Original T1 image')
    #parser.add_argument('--T2image', action='store', dest='t2Image', help='Original T2 image')
    parser.add_argument('--SubjID', action='store', dest='subjID', help='Subject_Session')
    parser.add_argument('--BrainLabel', action='store', dest='brainLabel',
                        help='The normalized T1 image with the skull removed. Normalized 0-110 where white matter=110.')
    parser.add_argument('--WMProb', action='store', dest='wmProbImage', help='')
    # TODO: Make parser group "Environment"
    parser.add_argument('--FSHomeDir', action='store', dest='FREESURFER_HOME',
                        default='/ipldev/sharedopt/20110601/MacOSX_10.6/freesurfer',
                        help='Location of FreeSurfer (differs for Mac and Linux environments')
    parser.add_argument('--FSSubjDir', action='store', dest='FS_SUBJECTS_DIR', help='FreeSurfer subjects directory')
    parser.add_argument('--FSSource', action='store', dest='FS_SCRIPT',
                        default='FreeSurferEnv.sh', help='')
    all_args=parser.parse_args()

    runAutoRecon(all_args.t1Image,
                 all_args.wmProbImage,
                 all_args.brainLabel,
                 all_args.subjID,
                 all_args.FREESURFER_HOME,
                 all_args.FS_SUBJECTS_DIR,
                 all_args.FS_SCRIPT)
