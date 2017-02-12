from __future__ import print_function
from __future__ import absolute_import
__author__ = 'johnsonhj'


import os
import SessionDB
import SimpleITK as sitk

#==================================
## GLOBALS ARE EVIL, but so is searching for all the items necessary to complete an analysis

FREESURFER_HOME='/Shared/pinc/sharedopt/20131115/RHEL6/freesurfer'
#subjects_dir='/Shared/paulsen/Experiments/20150617_PREDICTHD_FS'
subjects_dir='/Shared/johnsonhj/TrackOn/Experiments/TRACKON_FS20150701_NoT2'
scripts_dir = os.path.join(subjects_dir,'scripts')
if not os.path.exists(scripts_dir):
    os.makedirs(scripts_dir)

#subject_data_file='/Shared/johnsonhj/HDNI/20131116_TrackOn/scripts/edited_without_T2_PD_15s_track_autoworkup_TRACKON.csv'
#subject_data_file='/Shared/sinapse/CACHE/20131228_PREDICTHD_FS/20131124_FS.csv'
subject_data_file=os.path.join(subjects_dir,'scripts','2015-05-27_trackon_autoworkup.csv')
#subject_data_file=os.path.join(subjects_dir,'scripts/short_test.csv')
#subject_data_file=os.path.join(subjects_dir,'scripts/edited_without_T2_PD_15s_predict_autoworkup_PREDICT.csv')
USE_T2_FOR_FREESURFER=False ## For TRACKHD, don't use T2's
USE_15T_SESSIONS=False ## For TRACKHD, don't use 1.5T sessions
subjectDatabaseFile=os.path.join(scripts_dir,'subject_inputs.db')
#==================================

required_longfs_files=[
'surf/rh.volume',
'surf/lh.volume',
'stats/aseg.stats',
'stats/lh.BA.stats',
'stats/lh.BA.thresh.stats',
'stats/lh.aparc.DKTatlas40.stats',
'stats/lh.aparc.a2009s.stats',
'stats/lh.aparc.stats',
'stats/lh.curv.stats',
'stats/lh.entorhinal_exvivo.stats',
'stats/lh.w-g.pct.stats',
'stats/rh.BA.stats',
'stats/rh.BA.thresh.stats',
'stats/rh.aparc.DKTatlas40.stats',
'stats/rh.aparc.a2009s.stats',
'stats/rh.aparc.stats',
'stats/rh.curv.stats',
'stats/rh.entorhinal_exvivo.stats',
'stats/rh.w-g.pct.stats',
'stats/wmparc.stats',
'stats/aseg.stats',
]

required_qcache_files=[
'surf/lh.w-g.pct.mgh',
'surf/rh.w-g.pct.mgh',
'surf/lh.thickness.fsaverage.mgh',
'surf/lh.thickness.fwhm0.fsaverage.mgh',
'surf/lh.thickness.fwhm5.fsaverage.mgh',
'surf/lh.thickness.fwhm10.fsaverage.mgh',
'surf/lh.thickness.fwhm15.fsaverage.mgh',
'surf/lh.thickness.fwhm20.fsaverage.mgh',
'surf/lh.thickness.fwhm25.fsaverage.mgh',
'surf/lh.area.fsaverage.mgh',
'surf/lh.area.fwhm0.fsaverage.mgh',
'surf/lh.area.fwhm5.fsaverage.mgh',
'surf/lh.area.fwhm10.fsaverage.mgh',
'surf/lh.area.fwhm15.fsaverage.mgh',
'surf/lh.area.fwhm20.fsaverage.mgh',
'surf/lh.area.fwhm25.fsaverage.mgh',
'surf/lh.area.pial.fsaverage.mgh',
'surf/lh.area.pial.fwhm0.fsaverage.mgh',
'surf/lh.area.pial.fwhm5.fsaverage.mgh',
'surf/lh.area.pial.fwhm10.fsaverage.mgh',
'surf/lh.area.pial.fwhm15.fsaverage.mgh',
'surf/lh.area.pial.fwhm20.fsaverage.mgh',
'surf/lh.area.pial.fwhm25.fsaverage.mgh',
'surf/lh.volume.fsaverage.mgh',
'surf/lh.volume.fwhm0.fsaverage.mgh',
'surf/lh.volume.fwhm5.fsaverage.mgh',
'surf/lh.volume.fwhm10.fsaverage.mgh',
'surf/lh.volume.fwhm15.fsaverage.mgh',
'surf/lh.volume.fwhm20.fsaverage.mgh',
'surf/lh.volume.fwhm25.fsaverage.mgh',
'surf/lh.curv.fsaverage.mgh',
'surf/lh.curv.fwhm0.fsaverage.mgh',
'surf/lh.curv.fwhm5.fsaverage.mgh',
'surf/lh.curv.fwhm10.fsaverage.mgh',
'surf/lh.curv.fwhm15.fsaverage.mgh',
'surf/lh.curv.fwhm20.fsaverage.mgh',
'surf/lh.curv.fwhm25.fsaverage.mgh',
'surf/lh.sulc.fsaverage.mgh',
'surf/lh.sulc.fwhm0.fsaverage.mgh',
'surf/lh.sulc.fwhm5.fsaverage.mgh',
'surf/lh.sulc.fwhm10.fsaverage.mgh',
'surf/lh.sulc.fwhm15.fsaverage.mgh',
'surf/lh.sulc.fwhm20.fsaverage.mgh',
'surf/lh.sulc.fwhm25.fsaverage.mgh',
'surf/lh.jacobian_white.fsaverage.mgh',
'surf/lh.jacobian_white.fwhm0.fsaverage.mgh',
'surf/lh.jacobian_white.fwhm5.fsaverage.mgh',
'surf/lh.jacobian_white.fwhm10.fsaverage.mgh',
'surf/lh.jacobian_white.fwhm15.fsaverage.mgh',
'surf/lh.jacobian_white.fwhm20.fsaverage.mgh',
'surf/lh.jacobian_white.fwhm25.fsaverage.mgh',
'surf/lh.w-g.pct.mgh.fsaverage.mgh',
'surf/lh.w-g.pct.mgh.fwhm0.fsaverage.mgh',
'surf/lh.w-g.pct.mgh.fwhm5.fsaverage.mgh',
'surf/lh.w-g.pct.mgh.fwhm10.fsaverage.mgh',
'surf/lh.w-g.pct.mgh.fwhm15.fsaverage.mgh',
'surf/lh.w-g.pct.mgh.fwhm20.fsaverage.mgh',
'surf/lh.w-g.pct.mgh.fwhm25.fsaverage.mgh',
'surf/rh.thickness.fsaverage.mgh',
'surf/rh.thickness.fwhm0.fsaverage.mgh',
'surf/rh.thickness.fwhm5.fsaverage.mgh',
'surf/rh.thickness.fwhm10.fsaverage.mgh',
'surf/rh.thickness.fwhm15.fsaverage.mgh',
'surf/rh.thickness.fwhm20.fsaverage.mgh',
'surf/rh.thickness.fwhm25.fsaverage.mgh',
'surf/rh.area.fsaverage.mgh',
'surf/rh.area.fwhm0.fsaverage.mgh',
'surf/rh.area.fwhm5.fsaverage.mgh',
'surf/rh.area.fwhm10.fsaverage.mgh',
'surf/rh.area.fwhm15.fsaverage.mgh',
'surf/rh.area.fwhm20.fsaverage.mgh',
'surf/rh.area.fwhm25.fsaverage.mgh',
'surf/rh.area.pial.fsaverage.mgh',
'surf/rh.area.pial.fwhm0.fsaverage.mgh',
'surf/rh.area.pial.fwhm5.fsaverage.mgh',
'surf/rh.area.pial.fwhm10.fsaverage.mgh',
'surf/rh.area.pial.fwhm15.fsaverage.mgh',
'surf/rh.area.pial.fwhm20.fsaverage.mgh',
'surf/rh.area.pial.fwhm25.fsaverage.mgh',
'surf/rh.volume.fsaverage.mgh',
'surf/rh.volume.fwhm0.fsaverage.mgh',
'surf/rh.volume.fwhm5.fsaverage.mgh',
'surf/rh.volume.fwhm10.fsaverage.mgh',
'surf/rh.volume.fwhm15.fsaverage.mgh',
'surf/rh.volume.fwhm20.fsaverage.mgh',
'surf/rh.volume.fwhm25.fsaverage.mgh',
'surf/rh.curv.fsaverage.mgh',
'surf/rh.curv.fwhm0.fsaverage.mgh',
'surf/rh.curv.fwhm5.fsaverage.mgh',
'surf/rh.curv.fwhm10.fsaverage.mgh',
'surf/rh.curv.fwhm15.fsaverage.mgh',
'surf/rh.curv.fwhm20.fsaverage.mgh',
'surf/rh.curv.fwhm25.fsaverage.mgh',
'surf/rh.sulc.fsaverage.mgh',
'surf/rh.sulc.fwhm0.fsaverage.mgh',
'surf/rh.sulc.fwhm5.fsaverage.mgh',
'surf/rh.sulc.fwhm10.fsaverage.mgh',
'surf/rh.sulc.fwhm15.fsaverage.mgh',
'surf/rh.sulc.fwhm20.fsaverage.mgh',
'surf/rh.sulc.fwhm25.fsaverage.mgh',
'surf/rh.jacobian_white.fsaverage.mgh',
'surf/rh.jacobian_white.fwhm0.fsaverage.mgh',
'surf/rh.jacobian_white.fwhm5.fsaverage.mgh',
'surf/rh.jacobian_white.fwhm10.fsaverage.mgh',
'surf/rh.jacobian_white.fwhm15.fsaverage.mgh',
'surf/rh.jacobian_white.fwhm20.fsaverage.mgh',
'surf/rh.jacobian_white.fwhm25.fsaverage.mgh',
'surf/rh.w-g.pct.mgh.fsaverage.mgh',
'surf/rh.w-g.pct.mgh.fwhm0.fsaverage.mgh',
'surf/rh.w-g.pct.mgh.fwhm5.fsaverage.mgh',
'surf/rh.w-g.pct.mgh.fwhm10.fsaverage.mgh',
'surf/rh.w-g.pct.mgh.fwhm15.fsaverage.mgh',
'surf/rh.w-g.pct.mgh.fwhm20.fsaverage.mgh',
'surf/rh.w-g.pct.mgh.fwhm25.fsaverage.mgh'
]

def mkfsscript(session,outscript,t1list,t2list,is3T,useT2):
    ## TODO:  Work on staging files in/out

    T1_FLAGS=''
    for t1 in t1list:
        T1_FLAGS+= " -i " + t1

    T2_FLAGS=''
    if len(t2list) > 0 and useT2 == True:
        ## HACK:  Rachael requested that this NOT be used
        T2_FLAGS=' -T2 '+ t2list[0]
        pass

    ThreeTFlag = ''
    if is3T:
        ThreeTFlag=' -3T '
        job_name='F3B{SUBJ_SESSION_ID}'.format(SUBJ_SESSION_ID=session)
    else:
        job_name='F1B{SUBJ_SESSION_ID}'.format(SUBJ_SESSION_ID=session)

    auto_recon_script = """#!/bin/bash
#$ -o {FSSUBJDIR}/scripts/base_{SUBJ_SESSION_ID}_qsub.out
#$ -e {FSSUBJDIR}/scripts/base_{SUBJ_SESSION_ID}_qsub.err
#$ -cwd
#$ -N {JOB_NAME}

mkdir -p {FSSUBJDIR}/scripts/
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}

source {FSHOME}/FreeSurferEnv.sh

MAX_SLOTS=$NSLOTS
if [ -z "$MAX_SLOTS" ]; then
  MAX_SLOTS=2
fi

if [ ! -f {FSSUBJDIR}/{SUBJ_SESSION_ID}/mri/orig/001.mgz ]; then
  rm -rf {FSSUBJDIR}/{SUBJ_SESSION_ID}
  #Some data sets have FOV greater than 256mm, so force clipping to 256 with the -cw flag
  {FSHOME}/bin/recon-all {THREETFLAG} -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T1_FLAGS} {T2_FLAGS} -cw256
fi
rm -rf {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/IsRunning.*

{FSHOME}/bin/recon-all {THREETFLAG} -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T2_FLAGS} -make autorecon1
recon1_stats=$?
if [ $recon1_stats -eq 0 ];then
  {FSHOME}/bin/recon-all {THREETFLAG} -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T2_FLAGS} -make autorecon2
  recon2_stats=$?
  if [ $recon2_stats -eq 0 ];then
    {FSHOME}/bin/recon-all {THREETFLAG} -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T2_FLAGS} -make autorecon3
    recon3_stats=$?
      if [ $recon3_stats -eq 0 ];then
        {FSHOME}/bin/recon-all {THREETFLAG} -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T2_FLAGS} -make qcache
        qcache_stats=$?
      else
        echo "SKIPPING STAGE3 because STAGE2 failed"
      fi
  else
    echo "SKIPPING STAGE3 because STAGE2 failed"
  fi
else
  echo "SKIPPING STAGE2 because STAGE1 failed"
fi

status=$?
exit $status
""".format(
           FSHOME=FREESURFER_HOME,
           FSSUBJDIR=subjects_dir,
           SUBJ_SESSION_ID=session,
           JOB_NAME=job_name,
           THREETFLAG=ThreeTFlag,
           T1_FLAGS=T1_FLAGS,T2_FLAGS=T2_FLAGS)
    outfd = open(outscript,'w')
    outfd.write(auto_recon_script)
    outfd.close()
    return job_name


def ValidateBaseTPS(base_tps_file,found_sessions,subject,templateID):
    return_status=True
    previous_list = list()
    if os.path.exists(base_tps_file):
        tpsFF=open(base_tps_file,'r')
        for tps_session in tpsFF.readlines():
            previous_list.append(tps_session.lstrip().rstrip())
        tpsFF.close()
    else:
        return_status=False

    if len(found_sessions) != len(previous_list):
        return_status=False
    else:
        found_sessions.sort()
        previous_list.sort()
        for session in found_sessions:
            if not session in previous_list:
                return_status=False

    if return_status == False:
        print("WARNING:  DON'T KNOW WHAT TO DO: {0}".format(base_tps_file))
        print("current   {0}".format( found_sessions ) )
        print("base-tps   {0}".format( previous_list ) )
        import shutil
        templ_dir=os.path.join(subjects_dir,templateID)
        if os.path.exists(templ_dir):
            print("REMOVE TEMPLATE: {0}".format(templ_dir) )
            try:
              shutil.rmtree(templ_dir)
            except:
              print("MANUALLY REMOVE: {0}".format( templ_dir ))
            pass
        else:
            print("NO NEED TO REMOVE TEMPLATE: ".format( templ_dir))

        for session in found_sessions:
            long_sess_dir=os.path.join(subjects_dir,session+".long."+templateID)
            if os.path.exists(long_sess_dir):
                print("REMOVE LONG: ", long_sess_dir)
                try:
                    shutil.rmtree(long_sess_dir)
                except:
                    print("MANUALLY REMOVE: ", long_sess_dir)
                pass
            else:
                print("LONG already deleted: ", long_sess_dir)

    return return_status


def mktemplatescript(templateID, sessionList, outscript, dependantJobNames):
    timePoints=""
    for session in sessionList:
        timePoints += " -tp " + session

    if len(dependantJobNames) > 0:
        hold_jid='#$ -hold_jid '+",".join(dependantJobNames)
    else:
        hold_jid=''
    job_name = 'TPL{TEMP_ID}'.format(TEMP_ID=templateID)
    auto_recon_script = """#!/bin/bash
#$ -o {FSSUBJDIR}/scripts/temp_{SUBJ_SESSION_ID}_qsub.out
#$ -e {FSSUBJDIR}/scripts/temp_{SUBJ_SESSION_ID}_qsub.err
#$ -cwd
#$ -N {JOB_NAME}
{HOLD_FOR_JOBS}

mkdir -p {FSSUBJDIR}/scripts/
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}

source {FSHOME}/FreeSurferEnv.sh

MAX_SLOTS=$NSLOTS
if [ -z "$MAX_SLOTS" ]; then
  MAX_SLOTS=2
fi

rm -rf {FSSUBJDIR}/{TEMP_ID}/scripts/IsRunning.*
{FSHOME}/bin/recon-all -debug -openmp $MAX_SLOTS -base {TEMP_ID} {SESSIONS_FLAGS}  -all
recon_long_stat=$?

exit $recon_long_stat
""".format(
           FSHOME=FREESURFER_HOME,
           FSSUBJDIR=subjects_dir,
           SESSIONS_FLAGS=timePoints,
           TEMP_ID=templateID,
           JOB_NAME=job_name,
           SUBJ_SESSION_ID=thisSubject,
           HOLD_FOR_JOBS=hold_jid
           )
    outfd = open(outscript,'w')
    outfd.write(auto_recon_script)
    return job_name

def mklongscript(templateID, session, outscript,dependantJobNames,mode,is3T,useT2):
    if len(dependantJobNames) > 0:
        hold_jid='#$ -hold_jid '+",".join(dependantJobNames)
    else:
        hold_jid=''
    if is3T:
        job_SUFF='30'
    else:
        job_SUFF='15'

    if mode == "all":
        job_name='L{SUBJ_SESSION_ID}_{job_SUFF}'.format(SUBJ_SESSION_ID=session,job_SUFF=job_SUFF)
        mode_flag="all"
    else:
        job_name='Q{SUBJ_SESSION_ID}_{job_SUFF}'.format(SUBJ_SESSION_ID=session,job_SUFF=job_SUFF)
        mode_flag="qcache"

    ThreeTFlag = ''
    if is3T:
        ThreeTFlag=' -3T '


    auto_recon_script = """#!/bin/bash
#$ -o {FSSUBJDIR}/scripts/{JOB_NAME}_qsub.out
#$ -e {FSSUBJDIR}/scripts/{JOB_NAME}_qsub.err
#$ -cwd
#$ -N {JOB_NAME}
{HOLD_FOR_JOBS}

mkdir -p {FSSUBJDIR}/scripts/
export FREESURFER_HOME={FSHOME}
export SUBJECTS_DIR={FSSUBJDIR}

source {FSHOME}/FreeSurferEnv.sh

MAX_SLOTS=$NSLOTS
if [ -z "$MAX_SLOTS" ]; then
  MAX_SLOTS=2
fi

rm -rf {FSSUBJDIR}/{THIS_SESSION}.long.{TEMP_ID}/scripts/IsRunning.*
recon-all  {THREETFLAG} -debug -openmp $MAX_SLOTS -long {THIS_SESSION} {TEMP_ID} -{MODE_FLAG}

exit $recon_long_stat
""".format(
           FSHOME=FREESURFER_HOME,
           FSSUBJDIR=subjects_dir,
           TEMP_ID=templateID,
           SUBJ_SESSION_ID=session,
           JOB_NAME=job_name,
           THIS_SESSION=session,
           HOLD_FOR_JOBS=hold_jid,
           THREETFLAG=ThreeTFlag,
           MODE_FLAG=mode_flag
           )
    outfd = open(outscript,'w')
    outfd.write(auto_recon_script)
    outfd.close()
    return job_name


##################
## MAIN WORK #####
##################
single_subject=['all']
#single_subject=['000517510']
mountPrefix=''

if (not os.path.exists(subjectDatabaseFile)) or (
        os.path.getmtime(subjectDatabaseFile) < os.path.getmtime(subject_data_file)):
    ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
    ExperimentDatabase.MakeNewDB(subject_data_file, mountPrefix)
    ExperimentDatabase = None
    ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
else:
    print("Single_subject {0}: Using cached database, {1}".format(single_subject,subjectDatabaseFile))
    ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)

import pickle
pickled_good_list_fn = 'good_list.obj'
if os.path.exists(pickled_good_list_fn):
  good_list = pickle.load(open(pickled_good_list_fn,'r'))
else:
  good_list = dict()


#  pickle.dump(obj, file[, protocol])

def GetBaseSize(filename):
    filename = getInputFileName(filename)
    tmp = good_list
    base_size = good_list.get(filename)
    #print "#"*40+filename
    if not base_size:
      base_size = sitk.ReadImage(filename.encode('ascii','replace')).GetSpacing()
      good_list[filename] = base_size
      pf=open(pickled_good_list_fn,'w')
      pickle.dump(good_list,pf)
      pf.close()
    return base_size

def getInputFileName(filename):
    """ Does conversion to nifti if necessary, because NRRD files are not supported by FreeSurfer"""
    outfn=filename
    if filename.find('.nrrd'):
       outfn = filename.replace('.nrrd','.nii.gz')
       if not os.path.exists(outfn):
           tempIm = sitk.ReadImage(filename.encode('ascii','replace'))
           sitk.WriteImage( tempIm, outfn.encode('ascii','replace'))
    return outfn

def find_mgz(inlist_withNrrd):
    inlist = list()
    for ff in inlist_withNrrd:
        inlist.append(getInputFileName(ff))
    outlist = list()
    if len(inlist) == 0:
        return outlist

    same_size = list()
    base_size = GetBaseSize(inlist[0])
    for ff in inlist:
        new_size = GetBaseSize(ff)
        if new_size != base_size:
            continue
        same_size.append(ff)

    for ff in same_size:
        outfn=getInputFileName(ff)
        testmgz=outfn.replace('.nii.gz','.mgz')
        if os.path.exists(testmgz):
            outlist.append(testmgz)
        else:
            ## TODO: fix to handle this better
            #print("WARNING: Missing MGZ version so using nii.gz version: {0}".format(ff))
            outlist.append(outfn)
            pass
    return outlist

def GetMissingFilesList(subjects_dir,session_name,required_files_list):
    fs_full_paths = [ os.path.join(subjects_dir,session_name,fs_file) for fs_file in required_files_list ]
    missing_files = [ fs_file_fpath for fs_file_fpath in fs_full_paths if not os.path.exists( fs_file_fpath ) ]
    return fs_full_paths,missing_files


#==================================

type_report = ""
base_done = 0
temp_done = 0
long_done = 0

all_subjects=ExperimentDatabase.getAllSubjects()
all_subjects.sort()
all_subjects.reverse()
for thisSubject in all_subjects:
    #if thisSubject != '0152':
    #    continue
    thisSubject_sessions=ExperimentDatabase.getSessionsFromSubject(thisSubject)

    ## -------------------------
    ## Do initial timepoints
    ThreeT_sessions=list()
    OneT_sessions=list()
    base3T_job_names=list()
    base1T_job_names=list()
    for session in thisSubject_sessions:
        # NEVER DO THIS! if session != '77574':
        # NEVER DO THIS!   continue
        T1_files_30=ExperimentDatabase.getFilenamesByScantype(session,['T1-30'])
        T2_files=ExperimentDatabase.getFilenamesByScantype(session,['T2-30'])

        if len(T1_files_30) > 0:
            is3T = True
        else:
            if USE_15T_SESSIONS == False:
                continue  #Skip this session
            T1_files_15=ExperimentDatabase.getFilenamesByScantype(session,['T1-15'])
            is3T = False

        if is3T:
            fsscript = os.path.join(scripts_dir,'f3_'+session+'.sh')
        else:
            fsscript = os.path.join(scripts_dir,'f1_'+session+'.sh')
        T1_files = list()
        T2_files = list()
        this_session_base_done=False
        if is3T:#  len(T1_files_30) > 0:
            T1_files=find_mgz(T1_files_30)
            T2_files_30=ExperimentDatabase.getFilenamesByScantype(session,['T2-30'])
            T2_files=find_mgz(T2_files_30)
            if len(T2_files)  != len(T2_files_30):
                print("\n\nWARNING:\n\n: some mgz files missing: {0}".format( T2_files_30 ))
                print("\n {0}\n".format(T2_files))
                #continue
        else:
            T1_files=find_mgz(T1_files_15)
            T2_files=list()

        if len(T1_files) > 0:
            if is3T:
                ThreeT_sessions.append(session)
            else:
                OneT_sessions.append(session)
                pass
        else:
            print("SKIPPING: No T1 files missing")
            continue

        sentinal_file=os.path.join(subjects_dir,session,'surf/rh.jacobian_white.fwhm15.fsaverage.mgh')
        if os.path.exists(sentinal_file):
            print("1DONE:",  session, ":", sentinal_file,":")
            base_done += 1
            if os.path.exists(fsscript):
                os.unlink(fsscript)
            this_session_base_done=True
        else:
          print("1TODO:",  session, ":", sentinal_file,":")
          job_name=mkfsscript(session,fsscript,T1_files,T2_files,is3T,USE_T2_FOR_FREESURFER)
          if is3T:
              base3T_job_names.append(job_name)
          else:
              base1T_job_names.append(job_name)
        type_report += "{session},{is3T},{numT1s},{numT2s}\n".format(session=session,is3T=is3T,numT1s=len(T1_files),numT2s=len(T2_files))
        #print T1_files
        #print T2_files
        #print type_report

    if len(ThreeT_sessions) > 0:
        ## -------------------------
        ## Do 3T template building
        template_job_names=list()
        fsscript = os.path.join(scripts_dir,'tpl30_'+thisSubject+'.sh')
        templateID = thisSubject+'.template'
        list_of_sentinal_files=list()
        master_sentinal_file=os.path.join(subjects_dir,templateID,'surf/rh.volume')
        list_of_sentinal_files.append(sentinal_file)

        for session in ThreeT_sessions:
           sentinal_file=os.path.join(subjects_dir,templateID,'mri','transforms',
                                      thisSubject+".template_to_"+session+".lta")
           list_of_sentinal_files.append(sentinal_file)
        for session in ThreeT_sessions:
           sentinal_file=os.path.join(subjects_dir,templateID,'mri','transforms',
                                      session+"_to_"+thisSubject+".template.lta")
           list_of_sentinal_files.append(sentinal_file)
        file_missing=False
        for candiate_file in list_of_sentinal_files:
            if not os.path.exists(candiate_file):  ##If any files missing, break
                print("MISSING: {0}".format(candiate_file))
                file_missing=True
            else:
                print("OK: {0}".format(candiate_file))

        #/Shared/paulsen/Experiments/20150617_PREDICTHD_FS/2739.template/mri/transforms/2739.template_to_96241.lta
        #/Shared/paulsen/Experiments/20150617_PREDICTHD_FS/2739.template/mri/transforms/96241_to_2739.template.lta

        base_tps_file=os.path.join(subjects_dir,templateID,'base-tps')
        validatedsamelength = ValidateBaseTPS(base_tps_file,ThreeT_sessions,thisSubject,templateID)

        if not validatedsamelength or file_missing:
            if os.path.exists(master_sentinal_file):
                os.unlink(master_sentinal_file)
                pass;
        else:
            print("TPS OK for: ", templateID)

        if validatedsamelength and os.path.exists(master_sentinal_file):
            print("2DONE:", templateID,":")
            temp_done += 1
            if os.path.exists(fsscript):
                os.unlink(fsscript)
        else:
            print("2TODO:", templateID,":")
            template_job_name=mktemplatescript(templateID, ThreeT_sessions, fsscript, base3T_job_names)
            template_job_names.append(template_job_name)
            ## Do secondary run
            #make_template
            #template_reference_flag ="-t " + thisSubject
            #mkfsscript(session,fsscript,T1_files,T2_files,"baseline",is3T,USE_T2_FOR_FREESURFER)

        ## -------------------------
        ## Do 3T Longitudinal timepoints
        for session in ThreeT_sessions:
            T1_files=ExperimentDatabase.getFilenamesByScantype(session,['T1-30'])
            T2_files=ExperimentDatabase.getFilenamesByScantype(session,['T2-30'])

            T1_files=find_mgz(T1_files)
            T2_files=find_mgz(T2_files)

            long_job_names = list()
            fsscript = os.path.join(scripts_dir,'lg30_'+session+'.sh')
            sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.volume')

            fs_full_paths,missing_files = GetMissingFilesList(subjects_dir,session+'.long.'+templateID,required_longfs_files)
            #print("XXXXXXXX {0}\nYYYYYYYY {1}\nZZZZZZZ {2}\n".format(fs_full_paths,missing_files, len(missing_files) ))

            long_completed=False
            if len(missing_files) == 0:
                print("3DONE:",  session, "found :", len(fs_full_paths),": files")
                long_done += 1
                long_completed=True
                if os.path.exists(fsscript):
                    os.unlink(fsscript)
            else:
                print("3TODO:",  session, "missing :", len(missing_files), ": files", missing_files)
                long_job_name=mklongscript(templateID,session,fsscript,template_job_names,"all",is3T,USE_T2_FOR_FREESURFER)
                long_job_names.append(long_job_name)

            ## -------------------------
            ## Do 3T Longitudinal qcache timepoints
            qcache_job_names = list()
            fsscript = os.path.join(scripts_dir,'qc30_'+session+'.sh')
            sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.w-g.pct.mgh.fwhm25.fsaverage.mgh')
            fs_full_paths,missing_files = GetMissingFilesList(subjects_dir,session+'.long.'+templateID, required_qcache_files)
            if len(missing_files) == 0 and long_completed:
                print("3DONE:",  session, " found :", len(fs_full_paths) ,": files")
                if os.path.exists(fsscript):
                    os.unlink(fsscript)
            else:
                print("3TODO:",  session, " missing :", len(missing_files), ": files")
                qcache_job_name=mklongscript(templateID,session,fsscript,long_job_names,"qcache",is3T,USE_T2_FOR_FREESURFER)
                qcache_job_names.append(qcache_job_name)

    if len(OneT_sessions) > 0:
        ## -------------------------
        ## Do 1T template building
        template_job_names=list()
        fsscript = os.path.join(scripts_dir,'tpl15_'+thisSubject+'.sh')
        templateID = thisSubject+'.template15'
        sentinal_file=os.path.join(subjects_dir,templateID,'surf/rh.volume')

        base_tps_file=os.path.join(subjects_dir,templateID,'base-tps')
        validatedsamelength = ValidateBaseTPS(base_tps_file,OneT_sessions,thisSubject,templateID)

        if not validatedsamelength:
            if os.path.exists(sentinal_file):
                os.unlink(sentinal_file)
            pass;
        else:
            print("TPS OK for: ", templateID)

        if os.path.exists(sentinal_file):
            print("15TDONE:", templateID,":",sentinal_file)
            temp_done += 1
            if os.path.exists(fsscript):
                os.unlink(fsscript)
        else:
            print("15TTODO:", templateID,":",OneT_sessions)
            template_job_name=mktemplatescript(templateID, OneT_sessions, fsscript, base1T_job_names)
            template_job_names.append(template_job_name)
            ## Do secondary run
            #make_template
            #template_reference_flag ="-t " + thisSubject
            #mkfsscript(session,fsscript,T1_files,T2_files,"baseline")

        ## -------------------------
        ## Do 1T Longitudinal timepoints
        for session in OneT_sessions:
            T1_files=ExperimentDatabase.getFilenamesByScantype(session,['T1-15'])

            T1_files=find_mgz(T1_files)

            long_job_names = list()
            fsscript = os.path.join(scripts_dir,'lg15_'+session+'.sh')
            sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.volume')
            fs_full_paths,missing_files = GetMissingFilesList(subjects_dir,
                                                              session+'.long.'+templateID,
                                                              required_longfs_files)

            long_completed=False
            if len(missing_files) == 0:
                print("15LDONE:",  session, ":", sentinal_file,":")
                long_done += 1
                long_completed=True
                if os.path.exists(fsscript):
                    os.unlink(fsscript)
            else:
                print("15LTODO:",  session, "missing :", len(missing_files), ": files", missing_files)
                long_job_name=mklongscript(templateID,session,fsscript,template_job_names,"all",False,USE_T2_FOR_FREESURFER)
                long_job_names.append(long_job_name)


            ## -------------------------
            ## Do 1T Longitudinal qcache timepoints
            qcache_job_names = list()
            fsscript = os.path.join(scripts_dir,'qc15_'+session+'.sh')
            sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.w-g.pct.mgh.fwhm25.fsaverage.mgh')
            fs_full_paths,missing_files = GetMissingFilesList(subjects_dir,session+'.long.'+templateID, required_qcache_files)
            if len(missing_files) == 0 and long_completed:
                print("15QDONE:",  session, ":", sentinal_file,":")
                if os.path.exists(fsscript):
                    os.unlink(fsscript)
            else:
                print("15QTODO:",  session, ":", sentinal_file,":")
                qcache_job_name=mklongscript(templateID,session,fsscript,long_job_names,"qcache",is3T,USE_T2_FOR_FREESURFER)
                qcache_job_names.append(qcache_job_name)


ff = open('type_report.csv','w')
ff.write(type_report)
ff.close()

print("BASE COMPLETED: {0}".format(base_done))
print("TEMPLATE COMPLETED: {0}".format(temp_done))
print("LONGITUDINAL COMPLETED: {0}".format(long_done))
