#!/usr/bin/env python

import argparse
import os
import shutil
import SimpleITK as sitk
import subprocess

from PipeLineFunctionHelpers import mkdir_p, make_dummy_file


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


def IsFirstNewerThanSecond(firstFile, secondFile):
    if not os.path.exists(firstFile):
        print "ERROR: image missing", firstFile
        return True
    if not os.path.exists(secondFile):
        print "Returning True because file is missing:  {0}".format(secondFile)
        return True
    image_time = os.path.getmtime(firstFile)
    reference_time = os.path.getmtime(secondFile)
    if image_time > reference_time:
        print "Returning True because {0} is newer than {1}".format(firstFile, secondFile)
        return True
    return False

def run_mri_convert_script(niftiVol, mgzVol, subjects_dir, subj_session_id, FREESURFER_HOME, FS_SCRIPT):
    FS_SCRIPT_FN = os.path.join(FREESURFER_HOME, FS_SCRIPT)
    mri_convert_script = """#!/bin/bash
#$ -o {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/mri_convert_qsub.out
#$ -e {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/mri_convert_qsub.err
#$ -cwd
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}
source {SOURCE_SCRIPT}
{FSHOME}/bin/mri_convert --conform --out_data_type uchar {invol} {outvol}
status=$?
exit $status
""".format(SOURCE_SCRIPT=FS_SCRIPT_FN,
                                                                                         FSHOME=FREESURFER_HOME,
                                                                                         FSSUBJDIR=subjects_dir,
                                                                                         SUBJ_SESSION_ID=subj_session_id,
                                                                                         invol=niftiVol,
                                                                                         outvol=mgzVol)
    script_name = mgzVol + '_convert.sh'
    script = open(script_name, 'w')
    script.write(mri_convert_script)
    script.close()
    os.chmod(script_name, 0777)
    script_name_stdout = mgzVol + '_convert.out'
    script_name_stdout_fid = open(script_name_stdout, 'w')
    print "Starting mri_convert"
    scriptStatus = subprocess.check_call([script_name], stdout=script_name_stdout_fid, stderr=subprocess.STDOUT, shell='/bin/bash')
    if scriptStatus != 0:
        sys.exit(scriptStatus)
    print "Ending mri_convert"
    script_name_stdout_fid.close()
    return


def run_mri_mask_script(output_brainmask_fn_mgz, output_custom_brainmask_fn_mgz, output_nu_fn_mgz, subjects_dir, subj_session_id, FREESURFER_HOME, FS_SCRIPT):
    FS_SCRIPT_FN = os.path.join(FREESURFER_HOME, FS_SCRIPT)
    mri_mask_script = """#!/bin/bash
#$ -o {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/mri_mask_scipt_qsub.out
#$ -e {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/mri_mask_scipt_qsub.err
#$ -cwd
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}
source {SOURCE_SCRIPT}
{FSHOME}/bin/mri_add_xform_to_header -c {CURRENT}/transforms/talairach.xfm {maskvol} {maskvol}
{FSHOME}/bin/mri_mask {invol} {maskvol} {outvol}
status=$?
exit $status
""".format(SOURCE_SCRIPT=FS_SCRIPT_FN,
           FSHOME=FREESURFER_HOME,
           FSSUBJDIR=subjects_dir,
           SUBJ_SESSION_ID=subj_session_id,
           CURRENT=os.path.dirname(output_custom_brainmask_fn_mgz),
           invol=output_nu_fn_mgz,
           maskvol=output_custom_brainmask_fn_mgz,
           # regmaskvol=output_custom_brainmask_fn_mgz.replace(".mgz", "_reg.mgz"),
           outvol=output_brainmask_fn_mgz)
    script_name = output_brainmask_fn_mgz + '_mask.sh'
    script = open(script_name, 'w')
    script.write(mri_mask_script)
    script.close()
    os.chmod(script_name, 0777)
    script_name_stdout = output_brainmask_fn_mgz + '_convert.out'
    script_name_stdout_fid = open(script_name_stdout, 'w')
    print "Starting mri_mask"
    scriptStatus = subprocess.check_call([script_name], stdout=script_name_stdout_fid, stderr=subprocess.STDOUT, shell='/bin/bash')
    if scriptStatus != 0:
        sys.exit(scriptStatus)
    print "Ending mri_mask"
    script_name_stdout_fid.close()
    return

def baw_FixBrainMask(brainmask, subjects_dir, FREESURFER_HOME, FS_SCRIPT, subj_session_id):
    base_subj_dir = os.path.join(subjects_dir, subj_session_id, 'mri')
    mkdir_p(base_subj_dir)
    output_brainmask_fn_mgz = os.path.join(base_subj_dir, 'brainmask.mgz')
    output_custom_brainmask_fn = os.path.join(base_subj_dir, 'custom_brain_mask.nii.gz')
    output_custom_brainmask_fn_mgz = os.path.join(base_subj_dir, 'custom_brain_mask.mgz')
    output_nu_fn_mgz = os.path.join(base_subj_dir, 'nu.mgz')
    if IsFirstNewerThanSecond(brainmask, output_brainmask_fn_mgz) \
        or IsFirstNewerThanSecond(brainmask, output_nu_fn_mgz) \
        or IsFirstNewerThanSecond(output_nu_fn_mgz, output_custom_brainmask_fn_mgz):
        print "Fixing BrainMask recon-auto1 stage"
        brain = sitk.ReadImage(brainmask)
        blood = sitk.BinaryThreshold(brain, 5, 5)
        not_blood = 1 - blood
        clipping = sitk.BinaryThreshold(brain, 1, 1000000) - blood
        fill_size = 2
        ## HACK: Unfortunately we need to hole fill because of a WM bug in BABC where
        ## some white matter is being classified as background when it it being avoid due
        ## to too strict of multi-modal thresholding.
        hole_filled = sitk.ErodeObjectMorphology(sitk.DilateObjectMorphology(clipping, fill_size), fill_size)
        final_mask = sitk.Cast(hole_filled * not_blood, sitk.sitkUInt8)
        ## Now make an mgz version of the binary custom brain mask
        sitk.WriteImage(final_mask, output_custom_brainmask_fn) ## brain_matter mask with blood zero'ed out, and no "NOT" regions.
        run_mri_convert_script(output_custom_brainmask_fn, output_custom_brainmask_fn_mgz, subjects_dir,subj_session_id, FREESURFER_HOME, FS_SCRIPT)
        os.rename(output_brainmask_fn_mgz, output_brainmask_fn_mgz.replace('.mgz','_orig_backup.mgz'))
        ## Multipy output_brainmask_fn_mgz =  output_custom_brainmask_fn_mgz * output_nu_fn_mgz
        run_mri_mask_script(output_brainmask_fn_mgz, output_custom_brainmask_fn_mgz, output_nu_fn_mgz, subjects_dir, subj_session_id, FREESURFER_HOME, FS_SCRIPT)
    else:
        print "NOTHING TO BE DONE, SO SKIPPING."
        return  # Nothing to be done, files are already up-to-date.


def removeDir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)

def runAutoReconStage(subj_session_id, StageToRun, t1_fn, subjects_dir, FREESURFER_HOME, FS_SCRIPT):
    FS_SCRIPT_FN = os.path.join(FREESURFER_HOME, FS_SCRIPT)
    base_subj_dir = os.path.join(subjects_dir, subj_session_id)
    orig_001_mgz_fn = os.path.join(base_subj_dir,'mri','orig','001.mgz')
    if IsFirstNewerThanSecond(t1_fn, orig_001_mgz_fn):
        if os.path.exists(base_subj_dir):
            removeDir(base_subj_dir)
        mkdir_p(os.path.dirname(orig_001_mgz_fn))
        run_mri_convert_script(t1_fn, orig_001_mgz_fn, subjects_dir, subj_session_id, FREESURFER_HOME, FS_SCRIPT)
    auto_recon_script="""#!/bin/bash
#$ -o {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/autorecon{AUTORECONSTAGE}_qsub.out
#$ -e {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/autorecon{AUTORECONSTAGE}_qsub.err
#$ -cwd
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}
source {SOURCE_SCRIPT}
## Need to delete the "IsRunning" flags, if nipype pipeline is running this, then nipype prevents duplicates
rm -f {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/Is*

{FSHOME}/bin/recon-all -debug -subjid {SUBJ_SESSION_ID} -make autorecon{AUTORECONSTAGE}
status=$?
exit $status
""".format(SOURCE_SCRIPT=FS_SCRIPT_FN,
               FSHOME=FREESURFER_HOME,
               FSSUBJDIR=subjects_dir,
               AUTORECONSTAGE=StageToRun,
               SUBJ_SESSION_ID=subj_session_id)
    base_run_dir = os.path.join(subjects_dir, subj_session_id, 'scripts')
    mkdir_p(base_run_dir)
    script_name = os.path.join(base_run_dir,'run_autorecon_stage'+str(StageToRun)+'.sh')
    script = open(script_name, 'w')
    script.write(auto_recon_script)
    script.close()
    os.chmod(script_name, 0777)
    script_name_stdout = script_name + '_out'
    script_name_stdout_fid = open(script_name_stdout, 'w')
    print "Starting auto_recon Stage: {0} for SubjectSession {1}".format(StageToRun, subj_session_id)
    scriptStatus = subprocess.check_call([script_name], stdout=script_name_stdout_fid, stderr=subprocess.STDOUT, shell='/bin/bash')
    if scriptStatus != 0:
        sys.exit(scriptStatus)
    print "Ending auto_recon Stage: {0} for SubjectSession {1}".format(StageToRun, subj_session_id)
    script_name_stdout_fid.close()
    return


def runSubjectTemplate(args, FREESURFER_HOME, FS_SCRIPT):
    """ Create the within-subject template """
    base_template_id = args.base_template_id
    list_all_subj_session_ids = args.list_all_subj_session_ids
    subjects_dir = args.subjects_dir
    print "X"*80
    print "base_template_id :{0}:".format(base_template_id)
    print "Input a list of list_all_subj_session_ids :{0}:".format(list_all_subj_session_ids)
    print "subjects_dir :{0}:".format(subjects_dir)
    print "X"*80
    assert isinstance(list_all_subj_session_ids, list), "Must input a list of list_all_subj_session_ids :{0}:".format(list_all_subj_session_ids)
    StageToRun = "Within-SubjectTemplate"
    FS_SCRIPT_FN = os.path.join(FREESURFER_HOME, FS_SCRIPT)
    allTimePointFlags = ""
    for subj_session_id in list_all_subj_session_ids:
        allTimePointFlags += " -tp {timepoint}".format(timepoint=subj_session_id)
    allTimePointFlags += " -all"
    auto_recon_script="""#!/bin/bash
#$ -o {FSSUBJDIR}/{TEMPLATEID}/scripts/base_{TEMPLATEID}_qsub.out
#$ -e {FSSUBJDIR}/{TEMPLATEID}/scripts/base_{TEMPLATEID}_qsub.err
#$ -cwd
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}
source {SOURCE_SCRIPT}
## Need to delete the "IsRunning" flags, if nipype pipeline is running this, then nipype prevents duplicates
rm -f {FSSUBJDIR}/{TEMPLATEID}/scripts/Is*

if [ -f {FSSUBJDIR}/{TEMPLATEID}/stats/rh.entorhinal_exvivo.stats ]; then
   echo "--- SKIPPING: {TEMPLATEID} sentinal file already exits: {FSSUBJDIR}/{TEMPLATEID}/stats/rh.entorhinal_exvivo.stats"
   status=$?
else
   {FSHOME}/bin/recon-all -debug -base {TEMPLATEID} {ALL_TIME_POINTS}
   status=$?
fi
exit $status
""".format(SOURCE_SCRIPT=FS_SCRIPT_FN,
               FSHOME=FREESURFER_HOME,
               FSSUBJDIR=subjects_dir,
               TEMPLATEID=base_template_id,
               ALL_TIME_POINTS=allTimePointFlags)
    base_run_dir = os.path.join(subjects_dir, base_template_id, 'scripts')
    mkdir_p(base_run_dir)
    script_name = os.path.join(base_run_dir,'run_autorecon_stage_'+str(StageToRun)+'.sh')
    script = open(script_name, 'w')
    script.write(auto_recon_script)
    script.close()
    os.chmod(script_name, 0777)
    script_name_stdout = script_name + '_out'
    script_name_stdout_fid = open(script_name_stdout, 'w')
    print "Starting auto_recon Stage: {0} for SubjectSession {1}".format(StageToRun, base_template_id)
    scriptStatus = subprocess.check_call([script_name], stdout=script_name_stdout_fid, stderr=subprocess.STDOUT, shell='/bin/bash')
    if scriptStatus != 0:
        sys.exit(scriptStatus)
    print "Ending auto_recon Stage: {0} for SubjectSession {1}".format(StageToRun, base_template_id)
    script_name_stdout_fid.close()
    return


def runLongitudinal(args, FREESURFER_HOME, FS_SCRIPT):
    """ Create the longitudinal analysis """
    subj_session_id = args.subj_session_id
    subjects_dir = args.subjects_dir
    base_template_id = args.base_template_id
    assert isinstance(subj_session_id, str), "Must input a singel subj_session_id as string :{0}:".format(subj_session_id)
    StageToRun = "Longitudinal"
    FS_SCRIPT_FN = os.path.join(FREESURFER_HOME, FS_SCRIPT)
    auto_recon_script = """#!/bin/bash
#$ -o {FSSUBJDIR}/{LONGSESSIONID}/scripts/long_{TEMPLATEID}_{LONGSESSIONID}_qsub.out
#$ -e {FSSUBJDIR}/{LONGSESSIONID}/scripts/long_{TEMPLATEID}_{LONGSESSIONID}_qsub.err
#$ -cwd
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}
source {SOURCE_SCRIPT}
## Need to delete the "IsRunning" flags, if nipype pipeline is running this, then nipype prevents duplicates
rm -f {FSSUBJDIR}/{LONGSESSIONID}/scripts/Is*

if [ -f {FSSUBJDIR}/{LONGSESSIONID}.long.{TEMPLATEID}/stats/rh.entorhinal_exvivo.stats ]; then
   echo "--- SKIPPING: {LONGSESSIONID}.long.{TEMPLATEID} file already exits: {FSSUBJDIR}/{LONGSESSIONID}.long.{TEMPLATEID}/stats/rh.entorhinal_exvivo.stats"
   status=$?
else
   echo "--- RUNNING: {LONGSESSIONID}.long.{TEMPLATEID} file already exits: {FSSUBJDIR}/{LONGSESSIONID}.long.{TEMPLATEID}/stats/rh.entorhinal_exvivo.stats"
   {FSHOME}/bin/recon-all -debug -long {LONGSESSIONID} {TEMPLATEID} -all
   status=$?
fi
exit $status
""".format(SOURCE_SCRIPT=FS_SCRIPT_FN,
               FSHOME=FREESURFER_HOME,
               FSSUBJDIR=subjects_dir,
               TEMPLATEID=base_template_id,
               LONGSESSIONID=subj_session_id)
    base_run_dir = os.path.join(subjects_dir, subj_session_id, 'scripts')
    mkdir_p(base_run_dir)
    script_name = os.path.join(base_run_dir,'run_autorecon_stage_'+str(StageToRun)+'.sh')
    script = open(script_name, 'w')
    script.write(auto_recon_script)
    script.close()
    os.chmod(script_name, 0777)
    script_name_stdout = script_name + '_out'
    script_name_stdout_fid = open(script_name_stdout, 'w')
    print "Starting auto_recon Stage: {0} for SubjectSession {1}".format(StageToRun, base_template_id)
    scriptStatus = subprocess.check_call([script_name], stdout=script_name_stdout_fid, stderr=subprocess.STDOUT, shell='/bin/bash')
    if scriptStatus != 0:
        sys.exit(scriptStatus)
    print "Ending auto_recon Stage: {0} for SubjectSession {1}".format(StageToRun, base_template_id)
    script_name_stdout_fid.close()
    return


def runAutoRecon(args, FREESURFER_HOME, FS_SCRIPT):
    """Run all stages of AutoRecon For FreeSurfer, including the custom BAW initialization."""
    runAutoReconStage(args.subj_session_id, 1, args.T1_files, args.subjects_dir, FREESURFER_HOME, FS_SCRIPT)
    baw_FixBrainMask(args.brainmask, args.subjects_dir, FREESURFER_HOME, FS_SCRIPT, args.subj_session_id)
    runAutoReconStage(args.subj_session_id, 2, args.T1_files, args.subjects_dir, FREESURFER_HOME, FS_SCRIPT)
    runAutoReconStage(args.subj_session_id, 3, args.T1_files, args.subjects_dir, FREESURFER_HOME, FS_SCRIPT)
    dirname = os.path.join(args.subjects_dir, args.subj_session_id)
    t1 = os.path.join(dirname, 'mri', 'brain.mgz')
    label1 = os.path.join(dirname, 'mri_nifti', 'aparc+aseg.nii.gz')
    label2 = os.path.join(dirname, 'mri_nifti', 'aparc.a2009+aseg.nii.gz')
    return t1, label1, label2, args.subj_session_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
    Run various FreeSurfer's recon-all methods
                                     """, epilog="""
    For more information:
    ---------------------
    * http://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllDevTable
    * http://surfer.nmr.mgh.harvard.edu/fswiki/OtherUsefulFlags
    """)
    ## For the current nipype processing, the environments are set prior to running this script, so this code is not
    ## needed for using this script within the baw running envirionment.
    # TODO: Make parser group "Environment"
    # TODO: parser.add_argument('--FSHomeDir', action='store', dest='FREESURFER_HOME',
    # TODO:                  default='/ipldev/sharedopt/20110601/MacOSX_10.6/freesurfer',
    # TODO:                  help='Location of FreeSurfer (differs for Mac and Linux environments')
    # TODO: parser.add_argument('--FSSource', action='store', dest='FS_SCRIPT',
    # TODO:                  default='FreeSurferEnv.sh', help='')
    # TODO: local_FREESURFER_HOME = all_args.FREESURFER_HOME
    # TODO: local_FS_SCRIPT = all_args.FS_SCRIPT
    ### HACK:
    try:
        local_FREESURFER_HOME = os.environ['FREESURFER_HOME']
        if not os.path.exists(local_FREESURFER_HOME):
            raise Exception("INVALID PATH FOR FREESURFER HOME :{0}:".format(local_FREESURFER_HOME))
        # local_FS_SCRIPT = os.path.join(local_FREESURFER_HOME,'FreeSurferEnv.sh')
        local_FS_SCRIPT = 'FreeSurferEnv.sh'
    except KeyError, err:
        print KeyError
        print err
        raise KeyError
    ### END HACK
    subparsers = parser.add_subparsers(help='Currently supported subprocesses: "autorecon", "template", "longitudinal"')
    # Create -make subparser
    autorecon = subparsers.add_parser('autorecon', help='Link to recon-all i/o table: http://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllDevTable')
    autorecon.add_argument('--subjects_dir', action='store', dest='subjects_dir', help='FreeSurfer subjects directory')
    autorecon.add_argument('--subj_session_id', action='store', dest='subj_session_id', required=True, help='Subject_Session')
    autorecon.add_argument('--T1_files', action='store', dest='T1_files', required=True, help='Original T1 image')
    autorecon.add_argument('--brainmask', action='store', dest='brainmask', required=True,
                           help='The normalized T1 image with the skull removed. Normalized 0-110 where white matter=110.')
    autorecon.set_defaults(func=runAutoRecon)
    # Create -base subparser
    template = subparsers.add_parser('template', help='Link to recon-all longitudinal processing: http://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalProcessing')
    template.add_argument('--subjects_dir', action='store', dest='subjects_dir',required=True, help='FreeSurfer subjects directory')
    template.add_argument('--base_template_id', action='store', dest='base_template_id', required=True, help='Subject_template')
    template.add_argument('--list_all_subj_session_ids', action='store', dest='list_all_subj_session_ids', nargs='+', required=True, help='List of sessions for a subject template')
    template.set_defaults(func=runSubjectTemplate)
    # Create -long subparser
    longitudinal = subparsers.add_parser('longitudinal', help='Link to recon-all longitudinal processing: http://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalProcessing')
    longitudinal.add_argument('--subjects_dir', action='store', dest='subjects_dir', required=True, help='FreeSurfer subjects directory')
    longitudinal.add_argument('--subj_session_id', action='store', dest='subj_session_id', required=True, help='Subject_Session')
    longitudinal.add_argument('--base_template_id', action='store', dest='base_template_id', required=True, help='Template folder name (--base_template_id from "template" option)')
    longitudinal.set_defaults(func=runLongitudinal)
    # Parse inputs and run correct function
    all_args = parser.parse_args()
    all_args.func(all_args, local_FREESURFER_HOME, local_FS_SCRIPT)
