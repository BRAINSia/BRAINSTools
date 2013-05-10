#
# script to download DICOM files, convert them and compare them
# with the already converted files in /paulsen/MRx

#imports
import glob
import os
import re
import subprocess
import dicom
import shutil
from pyxnat import Interface

# Directory root for finding current converted files
PREDICT_BASE = '/paulsen/MRx'
# Directory to hold test data
DEST_BASE = '/scratch/kent/DWI_test'
# User ID for contacting xnat instance
HAWKEYEID = 'williamsnk'
# URL for XNAT server
XNATURL = 'https://www.predict-hd.net/xnat'
# XNat cache directory
CACHEDIR = '/scratch/kent/DWI_test/cache'
# DWIconvert program path
DWICONVERT = '/scratch/kent/BRAINSTools/build/bin/DWIConvert'
# Comparison program
DWICOMPARE = '/scratch/kent/BRAINSTools/build/bin/DWISimpleCompare'

#
# caches the list of all NRRD files in the PREDICT base dir
def write_nrrds_to_file(file_name='/scratch/kent/DWI_test/kent_nrrds.txt'):
    nrrds = glob.glob(PREDICT_BASE + '/*/*/*/ANONRAW/*.nrrd')
    with open(file_name, 'w') as fn:
        for nr in nrrds:
            fn.write(nr + '\n')

#
# read in cached file list
def get_all_nrrds_from_file(file_name='/scratch/kent/DWI_test/kent_nrrds.txt'):
    nrrds = []
    with open(file_name, 'rU') as fn:
        for f in fn:
            nrrds.append(f.rstrip())
    return nrrds

#
# extract series # from filename
def get_series_number(nrrd_file):
    """
    (str) -> str
    Return series number for the given nrrd_file.
    >>>nrrd_file = '/paulsen/MRx/PHD_032/0454/31774/ANONRAW/0454_31774_DWI-31_7.nrrd'
    >>>get_scan_type(nrrd_file)
    7
    """
    return nrrd_file.split('_')[-1].replace('.nrrd', '')


#
# use pyxnat to download the DICOM dataset and unpack it
def download_scan(scans, dest_dir, series_number):
    """
    (pyxnat scan object,str) -> images in the given dest direcotry

    Download images in the given scan to the given directory.

    >>>destDir = '/IPLlinux/raid0/homes/cyang3/Desktop/yangc/DICOMTEST/FMRI_Conn/'
    >>>project = xnat.select.project('FMRI_HD_120')
    >>>subject = project.subject('2530')
    >>>experiment = subject.experiment('23953')
    >>>scans = experiment.scans()
    >>>downloadScan(scans, destDir)

    """

    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
#    else:
#        print 'The path exists!'
    scans.download(dest_dir, type=str(series_number), extract=True)


def get_project_subject_session(nrrd_file):
    """
    (str) -> tuple(str,str,str)
    Return a tuple which contains (project,subject,session) from the given nrrd fiel.
    >>>nrrd_file = '/paulsen/MRx/PHD_032/0454/31774/ANONRAW/0454_31774_DWI-31_7.nrrd'
    >>>get_scan_type(nrrd_file)
    (PHD_032,0450,31774)
    """
    nfs = nrrd_file.split('/')
    return (nfs[3], nfs[4], nfs[5])

def dicom_to_nrrds(dicom, dwi_convert, dest):
    pass

def compare_nrrds(nrrd1, nrrd2):
    pass

def get_scans(XNAT, project_subject_session):
    """
    (XNAT interface,tuple) -> scans object.
    Return scans for the given (project,subject,session) tuple.

    """
    pl = project_subject_session[0]
    sul = project_subject_session[1]
    sl = project_subject_session[2]
    scans = XNAT.select.project(pl).subject(sul).experiment(sl).scans()
    return scans

def get_scan_type(nrrd_name):
    rval = nrrd_name.split('_')[-2]
    rval = rval.split('-')
    rval = '_'.join(rval)
    return rval

def remove_zip(dir,number):
    pat = dir
    pat += '/*'
    pat += number
    pat += '.zip'
    for f in glob.glob(pat):
        os.remove(f)

def main(XNAT, nrrds):
    errFName = DEST_BASE
    errFName += '/ErrorScans.txt'
    errFile = open(errFName,'w')

    for nr in nrrds:
        pss = get_project_subject_session(nr)
        # print pss
        sn = get_series_number(nr)
        scans = get_scans(XNAT, pss)
        print 'Downloading DICOM for ', nr
        download_scan(scans, DEST_BASE, str(sn))
        remove_zip(DEST_BASE,sn)

        dicomdir = DEST_BASE
        dicomdir += '/'
        dicomdir += pss[-1]
        dicomdir += '/scans/'
        dicomdir += sn
        dicomdir += '-'
        dicomdir += get_scan_type(nr)
        dicomdir += '/resources/DICOM/files'

        outvol = os.path.basename(nr)
        convertcmd = [ DWICONVERT, '--inputDicomDirectory', \
                           dicomdir, '--outputVolume', outvol ]

        #
        # find out vendor
        firstDicom = dicomdir
        firstDicom += '/'
        firstDicom += os.listdir(dicomdir)[0]
        ds = dicom.read_file(firstDicom,stop_before_pixels=True)
        vendor = ds[0x0008,0x0070].value

        if re.search("[Ss][Ii][Ee][Mm][Ee][Nn][Ss]", vendor) is None:
            pass
        else:
            convertcmd.append('--useBMatrixGradientDirections')
            print 'Using B-Matrix for gradients'

        print 'Converting ', dicomdir, ' to ', outvol

        try:
            subprocess.check_output(convertcmd, stderr=subprocess.STDOUT, env = os.environ)
        except subprocess.CalledProcessError:
            print "can't convert ", dicomdir
            continue

        print 'Comparing ', nr, ' and ', outvol
        comparecmd = [ DWICOMPARE, '--inputVolume1', nr, \
                           '--inputVolume2', outvol, '--checkDWIData' ]
        convertresult = 0
        try:
            subprocess.check_call(comparecmd)
        except subprocess.CalledProcessError,e:
            convertresult = e.returncode

        if convertresult == 0:
            print 'conversion of ', nr, ' matches'
            #
            # clean up
            projdir = DEST_BASE
            projdir += '/'
            projdir += pss[-1]
            shutil.rmtree(projdir)
            os.remove(outvol)
        else:
            print 'conversion of ',nr, ' does not match'
            message = nr + " doesn't match " + outvol + '\n'
            errFile.write(message)
    errFile.close()



if __name__ == '__main__':
    xnatUrl = 'https://www.predict-hd.net/xnat'
    XNAT = Interface(server=xnatUrl, user='williamsnk', cachedir='/scratch/kent/DWI_test/cache')

    if not os.path.isfile('kent_nrrds.txt'):
        write_nrrds_to_file()

    nrrds = get_all_nrrds_from_file()
    main(XNAT, nrrds)
    print len(nrrds)
    pass
