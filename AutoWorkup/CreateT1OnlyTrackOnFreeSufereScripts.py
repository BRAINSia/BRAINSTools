__author__ = 'johnsonhj'


import os
import SessionDB

#==================================
## GLOBALS ARE EVIL, but so is searching for all the items necessary to complete an analysis

FREESURFER_HOME='/Shared/sinapse/sharedopt/20131115/RHEL6/freesurfer'
#subjects_dir='/Shared/sinapse/CACHE/TRACKON_FS20131220'
subjects_dir='/Shared/sinapse/CACHE/TRACKON_FS20131220_NoT2'
scripts_dir = os.path.join(subjects_dir,'scripts')
if not os.path.exists(scripts_dir):
    os.makedirs(scripts_dir)

#subject_data_file='/Shared/johnsonhj/HDNI/20131116_TrackOn/scripts/edited_without_T2_PD_15s_track_autoworkup_TRACKON.csv'
subject_data_file='/Shared/sinapse/CACHE/TRACKON_FS20131220_NoT2/scripts/FS_Alterations.csv'
subjectDatabaseFile=os.path.join(scripts_dir,'subject_inputs.db')
#==================================

def mkfsscript(session,outscript,t1list,t2list):

    T1_FLAGS=''
    for t1 in t1list:
        T1_FLAGS+= " -i " + t1

    T2_FLAGS=''
    if len(t2list) > 0 and subjects_dir == '/Shared/sinapse/CACHE/TRACKON_FS20131220' :
        ## HACK:  Rachael requested that this NOT be used
        T2_FLAGS=' -T2 '+ t2list[0]
        pass

    job_name='FSB{SUBJ_SESSION_ID}'.format(SUBJ_SESSION_ID=session)
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

rm -f {FSSUBJDIR}/{SUBJ_SESSION_ID}/scripts/IsRunning.lh+rh
if [ ! -f {FSSUBJDIR}/{SUBJ_SESSION_ID}/mri/orig/001.mgz ]; then
  #Some data sets have FOV greater than 256mm, so force clipping to 256 with the -cw flag
  {FSHOME}/bin/recon-all -3T -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T1_FLAGS} {T2_FLAGS} -cw256
fi

{FSHOME}/bin/recon-all -3T -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T2_FLAGS} -make autorecon1
recon1_stats=$?
if [ $recon1_stats -eq 0 ];then
  {FSHOME}/bin/recon-all -3T -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T2_FLAGS} -make autorecon2
  recon2_stats=$?
  if [ $recon2_stats -eq 0 ];then
    {FSHOME}/bin/recon-all -3T -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T2_FLAGS} -make autorecon3
    recon3_stats=$?
      if [ $recon3_stats -eq 0 ];then
        {FSHOME}/bin/recon-all -3T -debug -openmp $MAX_SLOTS  -subjid {SUBJ_SESSION_ID} {T2_FLAGS} -make qcache
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
           T1_FLAGS=T1_FLAGS,T2_FLAGS=T2_FLAGS)
    outfd = open(outscript,'w')
    outfd.write(auto_recon_script)
    return job_name



def mktemplatescript(templateID, sessionList, outscript, dependantJobNames):
    timePoints=""
    for session in sessionList:
        timePoints+= " -tp " + session

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

rm -f {FSSUBJDIR}/{TEMP_ID}/scripts/IsRunning.lh+rh
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

def mklongscript(templateID, session, outscript,dependantJobNames,mode):
    if len(dependantJobNames) > 0:
        hold_jid='#$ -hold_jid '+",".join(dependantJobNames)
    else:
        hold_jid=''

    if mode == "all":
        job_name='LNG{SUBJ_SESSION_ID}'.format(SUBJ_SESSION_ID=session)
        mode_flag="all"
    else:
        job_name='Q{SUBJ_SESSION_ID}'.format(SUBJ_SESSION_ID=session)
        mode_flag="qcache"

    auto_recon_script = """#!/bin/bash
#$ -o {FSSUBJDIR}/scripts/{MODE_FLAG}_{SUBJ_SESSION_ID}_qsub.out
#$ -e {FSSUBJDIR}/scripts/{MODE_FLAG}_{SUBJ_SESSION_ID}_qsub.err
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

rm -f {FSSUBJDIR}/{THIS_SESSION}.long.{TEMP_ID}/scripts/IsRunning.lh+rh
recon-all  -3T -debug -openmp $MAX_SLOTS -long {THIS_SESSION} {TEMP_ID} -{MODE_FLAG}

exit $recon_long_stat
""".format(
           FSHOME=FREESURFER_HOME,
           FSSUBJDIR=subjects_dir,
           TEMP_ID=templateID,
           SUBJ_SESSION_ID=session,
           JOB_NAME=job_name,
           THIS_SESSION=session,
           HOLD_FOR_JOBS=hold_jid,
           MODE_FLAG=mode_flag
           )
    outfd = open(outscript,'w')
    outfd.write(auto_recon_script)
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


def find_mgz(inlist):
    outlist = list()
    for ff in inlist:
        testmgz=ff.replace('.nii.gz','.mgz')
        ## TODO:  Search for .nrrd files and figure out what to do.
        if os.path.exists(testmgz):
            outlist.append(testmgz)
        else:
            outlist.append(ff)
    return outlist


#==================================

all_subjects=ExperimentDatabase.getAllSubjects()
for thisSubject in all_subjects:
    thisSubject_sessions=ExperimentDatabase.getSessionsFromSubject(thisSubject)

    ## -------------------------
    ## Do initial timepoints
    ThreeT_sessions=list()
    base_job_names=list()
    for session in thisSubject_sessions:
        T1_files=ExperimentDatabase.getFilenamesByScantype(session,['T1-30'])
        T2_files=ExperimentDatabase.getFilenamesByScantype(session,['T2-30'])

        T1_files=find_mgz(T1_files)
        T2_files=find_mgz(T2_files)
        #print T1_files
        #print T2_files

        if len(T1_files) > 0:
            ThreeT_sessions.append(session)
        else:
            continue

        fsscript = os.path.join(scripts_dir,'fs'+session+'.sh')
        sentinal_file=os.path.join(subjects_dir,session,'surf/rh.jacobian_white.fwhm15.fsaverage.mgh')
        if os.path.exists(sentinal_file):
            print "1DONE:",  session, ":", sentinal_file,":"
            if os.path.exists(fsscript):
                os.unlink(fsscript)
        else:
            print "1TODO:",  session, ":", sentinal_file,":"
            job_name=mkfsscript(session,fsscript,T1_files,T2_files)
            base_job_names.append(job_name)

    ## -------------------------
    ## Do template building
    template_job_names=list()
    fsscript = os.path.join(scripts_dir,'tpl'+thisSubject+'.sh')
    templateID = thisSubject+'.template'
    sentinal_file=os.path.join(subjects_dir,templateID,'surf/rh.volume')
    if os.path.exists(sentinal_file):
        print "2DONE:", templateID,":"
        if os.path.exists(fsscript):
            os.unlink(fsscript)
    else:
        print "2TODO:", templateID,":"
        template_job_name=mktemplatescript(templateID, ThreeT_sessions, fsscript, base_job_names)
        template_job_names.append(template_job_name)
        ## Do secondary run
        #make_template
        #template_reference_flag ="-t " + thisSubject
        #mkfsscript(session,fsscript,T1_files,T2_files,"baseline")

    ## -------------------------
    ## Do Longitudinal timepoints
    for session in ThreeT_sessions:
        T1_files=ExperimentDatabase.getFilenamesByScantype(session,['T1-30'])
        T2_files=ExperimentDatabase.getFilenamesByScantype(session,['T2-30'])

        T1_files=find_mgz(T1_files)
        T2_files=find_mgz(T2_files)

        long_job_names = list()
        fsscript = os.path.join(scripts_dir,'lg'+session+'.sh')
        sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.volume')
        if os.path.exists(sentinal_file):
            print "3DONE:",  session, ":", sentinal_file,":"
            if os.path.exists(fsscript):
                os.unlink(fsscript)
        else:
            print "3TODO:",  session, ":", sentinal_file,":"
            long_job_name=mklongscript(templateID,session,fsscript,template_job_names,"all")
            long_job_names.append(long_job_name)

        ## -------------------------
        ## Do Longitudinal qcache timepoints
        qcache_job_names = list()
        fsscript = os.path.join(scripts_dir,'qc'+session+'.sh')
        sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.w-g.pct.mgh.fwhm25.fsaverage.mgh')
        if os.path.exists(sentinal_file):
            print "3DONE:",  session, ":", sentinal_file,":"
            if os.path.exists(fsscript):
                os.unlink(fsscript)
        else:
            print "3TODO:",  session, ":", sentinal_file,":"
            qcache_job_name=mklongscript(templateID,session,fsscript,long_job_names,"qcache")
            qcache_job_names.append(qcache_job_name)
