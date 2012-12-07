import argparse

def doIt():
    pass


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
    parser.add_argument('-t', '--T1image', action='store', dest='T1image', help='Original T1 image')
    parser.add_argument('-h', '--FSHomeDir', action='store', dest='FSHomeDir',
                        default='/ipldev/sharedopt/20110601/MacOSX_10.6/freesurfer',
                        help='Location of FreeSurfer (differs for Mac and Linux environments')
    parser.add_argument('-d', '--FSSubjDir', action='store', dest='FSSubjDir', help='FreeSurfer subjects directory')
    parser.add_argument('-s', '--FSSource', action='store', dest='FSSource',
                        default='source ${FREESURFER_HOME}/FreeSurferEnv.csh', help='')
    parser.add_argument('-i', '--SubjID', action='store', dest='SubjID', help='Subject_Session')
    parser.add_argument('-b', '--Brainmask', action='store', dest='Brainmask',
                        help='The normalized T1 image with the skull removed. Normalized 0-110 where white matter=110.')
    parser.add_argument('-w', '--WMProb', action='store', dest='WMProb', help='')
    inputArguments = parser.parse_args()