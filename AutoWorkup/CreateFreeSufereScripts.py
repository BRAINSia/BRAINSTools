__author__ = 'johnsonhj'


import os
import SessionDB

## GLOBALS ARE EVIL, but so is searching for all the items necessary to complete an analysis

FREESURFER_HOME='/Shared/sinapse/sharedopt/20131115/RHEL6/freesurfer'
#subjects_dir='/Shared/sinapse/CACHE/TRACKON_FS20131220'
subjects_dir='/Shared/sinapse/CACHE/TRACKON_FS20131220_NoT2'
scripts_dir = os.path.join(subjects_dir,'scripts')
if not os.path.exists(scripts_dir):
    os.makedirs(scripts_dir)

subject_data_file='/Shared/johnsonhj/HDNI/20131116_TrackOn/scripts/edited_without_T2_PD_15s_track_autoworkup_TRACKON.csv'
subjectDatabaseFile=os.path.join(scripts_dir,'subject_inputs.db')

single_subject=['all']
#single_subject=['012366774']
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

all_sessions=ExperimentDatabase.getAllSessions()

print "sessions: ",all_sessions

def find_mgz(inlist):
    outlist = list()
    for ff in inlist:
        testmgz=ff.replace('.nii.gz','.mgz')
        if os.path.exists(testmgz):
            outlist.append(testmgz)
        else:
            outlist.append(ff)
    return outlist

def mkfsscript(session,outscript,t1list,t2list):

    T1_FLAGS=''
    for t1 in t1list:
        T1_FLAGS+= " -i " + t1

    T2_FLAGS=''
    if len(t2list) > 0 and subjects_dir == '/Shared/sinapse/CACHE/TRACKON_FS20131220' :
        ## HACK:  Rachael requested that this NOT be used
        T2_FLAGS=' -T2 '+ t2list[0]
        pass

    auto_recon_script = """#!/bin/bash
#$ -o {FSSUBJDIR}/scripts/base_{SUBJ_SESSION_ID}_qsub.out
#$ -e {FSSUBJDIR}/scripts/base_{SUBJ_SESSION_ID}_qsub.err
#$ -cwd
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
           T1_FLAGS=T1_FLAGS,T2_FLAGS=T2_FLAGS)
    outfd = open(outscript,'w')

    outfd.write(auto_recon_script)

for session in all_sessions:
    T1_files=ExperimentDatabase.getFilenamesByScantype(session,['T1-30'])
    T2_files=ExperimentDatabase.getFilenamesByScantype(session,['T2-30'])

    T1_files=find_mgz(T1_files)
    T2_files=find_mgz(T2_files)
    print T1_files
    print T2_files
    fsscript = os.path.join(scripts_dir,'fs'+session+'.sh')
    mkfsscript(session,fsscript,T1_files,T2_files)
