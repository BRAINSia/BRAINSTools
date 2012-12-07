import argparse
from nipype.interfaces.freesurfer.base import Info

class ReconContainer(object):
    def __init__(self, *args, **kwds):
        self.required = kwds['required']
        self._checkFiles()

    def _checkFiles(self):
        import os
        for dirname, flist in required.iteritems():
            for fname in flist:
                fullname = os.path.join(os.getcwd(), dirname, fname)
                if not os.path.exists(fullname):
                    self._touchFile()

    def _touchFile(self, fullname):
        import subprocess as sp
        touch = sp.Popen(['touch' fullname], stdout=sp.PIPE, stderr=sp.STDOUT)
        touch.stout.close()

def normalizeWMTo100(t1, wm_prob):
    """This function will compute the mean value of wm and rescale the image so that the mean WM=100"""
    WM_THRESHOLD = 0.66
    wm_mask = sitk.BinaryThreshold(wm_prob, WM_THRESHOLD)
    ls = sitk.LabelStatisticsImageFilter()
    ls.Execute(t1, wm_mask)
    wm_value = 1
    myMeasurementMap = ls.GetMeasurementMap(wm_value)
    MeanWM = myMeasurementMap['Mean']
    t1_new = sitk.Cast(sitk.Cast(t1, sitk.sitkFloat32) * 100.0 / MeanWM, sitk.sitkUInt8)
    return t1_new

def autoRecon1(*arg, **kwds):
    """
    TODO: Traditional Freesurfer stage (for comparison)
    """
    pass

def baw_Recon1(*arg, **kwds):
    required = {'orig':['001.mgz','002.mgz'],
                'transforms':['talairach.auto.xfm', 'talairach.xfm', 'talairach_avi.log', 'talairach_with_skull.lta'],
                '.':'brainmask.auto.mgz', 'brainmask.mgz', 'nu.mgz', 'T1.mgz', 'rawavg.mgz', 'orig.mgz']}
    container1 = ReconContainer(required=required)
    t1 = sitk.ReadImage(kwds['t1Image'])
    t2 = sitk.ReadImage(kwds['t2Image'])
    wm = sitk.ReadImage(kwds['wmProbImage'])
    brain = sitk.ReadImage(kwds['brainLabel'])
    t1_new = normalizeWMTo100(t1, wm)
    clipping = sitk.BinaryThreshold(brain, 1, 1000000) - sitk.BinaryThreshold(brain, 5, 5)
    clipped = t1_new * clipping

def autoRecon2(*arg, **kwds):
    pass

def autoRecon3(*arg, **kwds):
    pass

def run(*arg, **kwds):
    import os
    os.environ['SUBJECTS_DIR'] = kwds['FS_SUBJECTS_DIR']
    os.environ['FREESURFER_HOME'] = kwds['FREESURFER_HOME']

    # autoRecon1(*arg, **kwds)
    baw_Recon1(*arg, **kwds)
    autoRecon2(*arg, **kwds)
    autoRecon3(*arg, **kwds)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
                                     """

    """DELETE LATER: This is an just example of the commands required to run FreeSurfer recon-all:
    /bin/tcsh
    setenv FREESURFER_HOME /ipldev/sharedopt/20110601/MacOSX_10.6/freesurfer
    source ${FREESURFER_HOME}/FreeSurferEnv.csh
    setenv SUBJECTS_DIR /IPLlinux/raid0/homes/jforbes/freesurfer/recon-all/autorecon1_copy
    recon-all -make all -subjid 0074_24832

    Link to recon-all i/o table:
    http://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllDevTable

    """
                                     )
    parser.add_argument('-t1', '--T1image', action='store', dest='t1Image', help='Original T1 image')
    parser.add_argument('-t2', '--T2image', action='store', dest='t2Image', help='Original T2 image')
    parser.add_argument('-h', '--FSHomeDir', action='store', dest='FREESURFER_HOME',
                        default='/ipldev/sharedopt/20110601/MacOSX_10.6/freesurfer',
                        help='Location of FreeSurfer (differs for Mac and Linux environments')
    parser.add_argument('-d', '--FSSubjDir', action='store', dest='FS_SUBJECTS_DIR', help='FreeSurfer subjects directory')
    parser.add_argument('-s', '--FSSource', action='store', dest='FS_SOURCE',
                        default='FreeSurferEnv.csh', help='')
    parser.add_argument('-i', '--SubjID', action='store', dest='SubjID', help='Subject_Session')
    parser.add_argument('-b', '--BrainLabel', action='store', dest='brainLabel',
                        help='The normalized T1 image with the skull removed. Normalized 0-110 where white matter=110.')
    parser.add_argument('-w', '--WMProb', action='store', dest='wmProbImage', help='')
    kwds = vars(parser.parse_args())
    return runAutoRecon(**kwds)
