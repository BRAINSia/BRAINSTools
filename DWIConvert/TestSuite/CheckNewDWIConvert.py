
from builtins import str
#
# script to download DICOM files, convert them and compare them
# with the already converted files in /paulsen/MRx
#
# NOTE: this will only work on the PINC computer network
#
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
DEST_BASE = '/scratch/kent/DWI_test/TestResults'

# User ID for contacting xnat instance
HAWKEYEID = 'williamsnk'

# URL for XNAT server
XNATURL = 'https://www.predict-hd.net/xnat'

# XNat cache directory
CACHE_DIR = DEST_BASE + '/cache'

# DWIconvert program path
DWICONVERT = '/scratch/kent/BRAINSTools/release/bin/DWIConvert'

# Comparison program
DWICOMPARE = '/scratch/kent/BRAINSTools/release/bin/DWISimpleCompare'

# cached list of all already-converted scans
ALL_SCANS = '/scratch/kent/DWI_test/TestResults/AllScans.txt'

#
# caches the list of all NRRD files in the PREDICT base dir
def write_nrrds_to_file(file_name):
    nrrds = glob.glob(PREDICT_BASE + '/*/*/*/ANONRAW/*.nrrd')
    with open(file_name, 'w') as fn:
        for nr in nrrds:
            fn.write(nr + '\n')

#
# read in cached file list
def get_all_nrrds_from_file(file_name):
    nrrds = []
    with open(file_name, 'rU') as fn:
        for f in fn:
            nrrds.append(f.rstrip())
    print((file_name, " size ", len(nrrds)))
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

def remove_zip(dir, number):
    pat = dir
    pat += '/*'
    pat += number
    pat += '.zip'
    for f in glob.glob(pat):
        os.remove(f)

def remove_already_processed(nrrdList, processedFname):
    processed = get_all_nrrds_from_file(processedFname)
    # if any files have been processed
    if processed:
        origsize = len(nrrdList)
        processedsize = len(processed)
        rval = [x for x in nrrdList if x not in processed]
        print(processedFname)
        print(("Original # of files ", origsize, \
            " # processed ", processedsize, \
            "# new list len ", len(rval)))

    else:
        rval = nrrdList
    return rval

def main(XNAT, nrrds):
    # keep track of files that did not pass
    errFName = DEST_BASE + '/ErrorScans.txt'
    if os.path.isfile(errFName):
        nrrds = remove_already_processed(nrrds, errFName)
        errFile = open(errFName, 'a', 0)
    else:
        errFile = open(errFName, 'w', 0)

    # also track all files that did pass
    passFName = DEST_BASE + "/PassedScans.txt"
    if os.path.isfile(passFName):
        nrrds = remove_already_processed(nrrds, passFName)
        passFile = open(passFName, 'a', 0)
    else:
        passFile = open(passFName, 'w', 0)

    failedConversionsFName = DEST_BASE + "/FailedConversions.txt"
    if os.path.isfile(failedConversionsFName):
        nrrds = remove_already_processed(nrrds, failedConversionsFName)
        failedConversions = open(failedConversionsFName, "a", 0)
    else:
        failedConversions = open(failedConversionsFName, "w", 0)


    try:
        for nr in nrrds:
            pss = get_project_subject_session(nr)
            # print pss
            sn = get_series_number(nr)
            scans = get_scans(XNAT, pss)
            print(('Downloading DICOM for ', nr))
            try:
                download_scan(scans, DEST_BASE, str(sn))
                remove_zip(DEST_BASE, sn)
            except KeyboardInterrupt:
                print('Keyboard Interrupt')
                break
            except:
                print(("Error downloading files for ", nr))
                continue

            dicomdir = DEST_BASE + '/' \
                + pss[-1] \
                + '/scans/' \
                + sn \
                + '-' \
                + get_scan_type(nr) \
                + '/resources/DICOM/files'

            outvol = DEST_BASE + '/' + os.path.basename(nr)
            convertcmd = [ DWICONVERT, '--inputDicomDirectory', \
                               dicomdir, '--outputVolume', outvol ]

            #
            # find out vendor
            firstDicom = dicomdir + '/' + os.listdir(dicomdir)[0]

            ds = dicom.read_file(firstDicom, stop_before_pixels=True)
            vendor = ds[0x0008, 0x0070].value
            model = ds[0x0008, 0x1090].value
            print(("Scanner vendor ", vendor, " model ", model))

            if re.search("SIEMENS", vendor.upper()) is None:
                pass
            elif re.search("ALLEGRA", model.upper()) is None and \
                    re.search("TRIOTIM", mode.upper()) is None:
                convertcmd.append('--useBMatrixGradientDirections')
                print('Using B-Matrix for gradients')

            print(('Converting ', dicomdir, ' to ', outvol))
            print(("Command line", convertcmd))
            try:
                subprocess.check_output(convertcmd, stderr=subprocess.STDOUT, env = os.environ)
#                subprocess.call(convertcmd, stderr=subprocess.STDOUT, env = os.environ)
            except subprocess.CalledProcessError as error:
                print(("can't convert ", dicomdir))
                print(('\n', error.output))
                failedConversions.write(nr + '\n')
                continue

            print(('Comparing ', nr, ' and ', outvol))
            comparecmd = [ DWICOMPARE, '--inputVolume1', nr, \
                               '--inputVolume2', outvol, '--checkDWIData' ]
            convertresult = 0
            try:
                subprocess.check_call(comparecmd)
            except subprocess.CalledProcessError as e:
                convertresult = e.returncode

            if convertresult == 0:
                print(('conversion of ', nr, ' matches'))
                #
                # clean up
                projdir = DEST_BASE
                projdir += '/'
                projdir += pss[-1]
                shutil.rmtree(projdir)
                os.remove(outvol)
                passFile.write(nr + '\n')
            else:
                print(('conversion of ', nr, ' does not match'))
                errFile.write(nr + '\n')
    except KeyboardInterrupt:
        print('Keyboard Interrupt')
        pass
    except:
        print('')
        print('Error during processing')
    passFile.close()
    errFile.close()
    failedConversions.close()


if __name__ == '__main__':
    xnatUrl = 'https://www.predict-hd.net/xnat'
    XNAT = Interface(server=xnatUrl, user='williamsnk', cachedir=CACHE_DIR)

    if not os.path.isfile(ALL_SCANS):
        write_nrrds_to_file(ALL_SCANS)

    nrrds = get_all_nrrds_from_file(ALL_SCANS)
    main(XNAT, nrrds)
    pass
