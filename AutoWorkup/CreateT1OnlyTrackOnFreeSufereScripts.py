from __future__ import print_function
from __future__ import absolute_import
__author__ = 'johnsonhj'


import os
from . import SessionDB

#==================================
## GLOBALS ARE EVIL, but so is searching for all the items necessary to complete an analysis

FREESURFER_HOME='/Shared/pinc/sharedopt/20131115/RHEL6/freesurfer'
#subjects_dir='/Shared/sinapse/CACHE/TRACKON_FS20131220'
subjects_dir='/Shared/johnsonhj/TrackOn/Experiments/TRACKON_FS20131220_NoT2'
scripts_dir = os.path.join(subjects_dir,'scripts')
if not os.path.exists(scripts_dir):
    os.makedirs(scripts_dir)

subject_data_file=os.path.join(subjects_dir,'scripts/FS_Alterations2015.csv')
subjectDatabaseFile=os.path.join(scripts_dir,'subject_inputs.db')
#==================================

def mkfsscript(session,outscript,t1list,t2list,is3T):

    T1_FLAGS=''
    for t1 in t1list:
        T1_FLAGS+= " -i " + t1

    T2_FLAGS=''
    if len(t2list) > 0 and subjects_dir == '/Shared/sinapse/CACHE/TRACKON_FS20131220' :
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
        print("WARNING:  DON'T KNOW WHAT TO DO: ", base_tps_file)
        print("current   ", found_sessions)
        print("base-tps   ", previous_list)
        import shutil
        templ_dir=os.path.join(subjects_dir,templateID)
        if os.path.exists(templ_dir):
            print("REMOVE TEMPLATE: ", templ_dir)
            try:
              shutil.rmtree(templ_dir)
            except:
              print("MANUALLY REMOVE: ", templ_dir)
            pass
        else:
            print("NO NEED TO REMOVE TEMPLATE: ", templ_dir)

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

def mklongscript(templateID, session, outscript,dependantJobNames,mode,is3T):
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

    ThreeTFlag = ''
    if is3T:
        ThreeTFlag=' -3T '


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
    tmp = good_list
    base_size = good_list.get(filename)
    if not base_size:
      base_size = sitk.ReadImage(filename.encode('ascii','replace')).GetSpacing()
      good_list[filename] = base_size
      pf=open(pickled_good_list_fn,'w')
      pickle.dump(good_list,pf)
      pf.close()
    return base_size

def find_mgz(inlist):
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
        if ff.find('.nrrd'):
            outfn = ff.replace('.nrrd','.nii.gz')
            if not os.path.exists(outfn):
                tempIm = sitk.ReadImage(ff.encode('ascii','replace'))
                sitk.WriteImage( tempIm, outfn.encode('ascii','replace'))
                ff = outfn
        testmgz=ff.replace('.nii.gz','.mgz')
        ## TODO:  Search for .nrrd files and figure out what to do.
        if os.path.exists(testmgz):
            outlist.append(testmgz)
        else:
            ## TODO: fix to handle this better
            #print("WARNING: Missing MGZ version so using nii.gz version: {0}".format(ff))
            outlist.append(ff)
            pass
    return outlist


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
                print("SKIPPING: mgz files missing: {0}".format( T2_files_30 ))
                continue
        else:
            #print "HACK: ",T1_files_15
            T1_files=find_mgz(T1_files_15)
            #print "HACK: ",T1_files
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
          job_name=mkfsscript(session,fsscript,T1_files,T2_files,is3T)
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
        sentinal_file=os.path.join(subjects_dir,templateID,'surf/rh.volume')

        base_tps_file=os.path.join(subjects_dir,templateID,'base-tps')
        validatedsamelength = ValidateBaseTPS(base_tps_file,ThreeT_sessions,thisSubject,templateID)

        if not validatedsamelength:
            if os.path.exists(sentinal_file):
                os.unlink(sentinal_file)
            pass;
        else:
            print("TPS OK for: ", templateID)

        if validatedsamelength and os.path.exists(sentinal_file):
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
            #mkfsscript(session,fsscript,T1_files,T2_files,"baseline")

        ## -------------------------
        ## Do 3T Longitudinal timepoints
        for session in ThreeT_sessions:
            T1_files=ExperimentDatabase.getFilenamesByScantype(session,['T1-30'])
            T2_files=ExperimentDatabase.getFilenamesByScantype(session,['T2-30'])

            T1_files=find_mgz(T1_files)
            T2_files=find_mgz(T2_files)

            long_job_names = list()
            fsscript = os.path.join(scripts_dir,'lg'+session+'.sh')
            sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.volume')
            if os.path.exists(sentinal_file):
                print("3DONE:",  session, ":", sentinal_file,":")
                long_done += 1
                if os.path.exists(fsscript):
                    os.unlink(fsscript)
            else:
                print("3TODO:",  session, ":", sentinal_file,":")
                long_job_name=mklongscript(templateID,session,fsscript,template_job_names,"all",is3T)
                long_job_names.append(long_job_name)

            ## -------------------------
            ## Do 3T Longitudinal qcache timepoints
            qcache_job_names = list()
            fsscript = os.path.join(scripts_dir,'qc'+session+'.sh')
            sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.w-g.pct.mgh.fwhm25.fsaverage.mgh')
            if os.path.exists(sentinal_file):
                print("3DONE:",  session, ":", sentinal_file,":")
                if os.path.exists(fsscript):
                    os.unlink(fsscript)
            else:
                print("3TODO:",  session, ":", sentinal_file,":")
                qcache_job_name=mklongscript(templateID,session,fsscript,long_job_names,"qcache",is3T)
                qcache_job_names.append(qcache_job_name)

    if len(OneT_sessions) > 0:
        ## -------------------------
        ## Do 1T template building
        template_job_names=list()
        fsscript = os.path.join(scripts_dir,'tpl15_'+thisSubject+'.sh')
        templateID = thisSubject+'.template15'
        sentinal_file=os.path.join(subjects_dir,templateID,'surf/rh.volume')
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
            fsscript = os.path.join(scripts_dir,'lg15'+session+'.sh')
            sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.volume')
            if os.path.exists(sentinal_file):
                print("15LDONE:",  session, ":", sentinal_file,":")
                long_done += 1
                if os.path.exists(fsscript):
                    os.unlink(fsscript)
            else:
                print("15LTODO:",  session, ":", sentinal_file,":")
                long_job_name=mklongscript(templateID,session,fsscript,template_job_names,"all",False)
                long_job_names.append(long_job_name)

            ## -------------------------
            ## Do 1T Longitudinal qcache timepoints
            qcache_job_names = list()
            fsscript = os.path.join(scripts_dir,'qc15'+session+'.sh')
            sentinal_file=os.path.join(subjects_dir,session+'.long.'+templateID,'surf/rh.w-g.pct.mgh.fwhm25.fsaverage.mgh')
            if os.path.exists(sentinal_file):
                print("15QDONE:",  session, ":", sentinal_file,":")
                if os.path.exists(fsscript):
                    os.unlink(fsscript)
            else:
                print("15QTODO:",  session, ":", sentinal_file,":")
                qcache_job_name=mklongscript(templateID,session,fsscript,long_job_names,"qcache",is3T)
                qcache_job_names.append(qcache_job_name)


ff = open('type_report.csv','w')
ff.write(type_report)
ff.close()

print("BASE COMPLETED: {0}".format(base_done))
print("TEMPLATE COMPLETED: {0}".format(temp_done))
print("LONGITUDINAL COMPLETED: {0}".format(long_done))
