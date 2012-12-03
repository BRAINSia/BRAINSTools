##############################################################################
import os
import re
import sys
##############################################################################
##############################################################################
def iter_flatten(iterable):
  it = iter(iterable)
  for e in it:
    if isinstance(e, (list, tuple)):
      for f in iter_flatten(e):
        yield f
    else:
      yield e

##############################################################################
def get_global_sge_script(pythonPathsList,binPathsList,customEnvironment={}):
    """This is a wrapper script for running commands on an SGE cluster
so that all the python modules and commands are pathed properly"""

    custEnvString=""
    for key,value in customEnvironment.items():
        custEnvString+="export "+key+"="+value+"\n"

    PYTHONPATH=":".join(pythonPathsList)
    BASE_BUILDS=":".join(binPathsList)
    GLOBAL_SGE_SCRIPT="""#!/bin/bash
echo "STARTED at: $(date +'%F-%T')"
echo "Ran on: $(hostname)"
export PATH={BINPATH}
export PYTHONPATH={PYTHONPATH}

echo "========= CUSTOM ENVIORNMENT SETTINGS =========="
echo "export PYTHONPATH={PYTHONPATH}"
echo "export PATH={BINPATH}"
echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

echo "With custom environment:"
echo {CUSTENV}
{CUSTENV}
## NOTE:  nipype inserts the actual commands that need running below this section.
""".format(PYTHONPATH=PYTHONPATH,BINPATH=BASE_BUILDS,CUSTENV=custEnvString)
    return GLOBAL_SGE_SCRIPT

#{ --------------------------------------------------------------------------------------- #
def findInputImagesForSubject( inputT1, BRAINSAtlasDir, outputDir ):
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print inputT1
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    import os
    import sys

    TissueClassifyDir = os.path.dirname( inputT1 )
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print TissueClassifyDir
    outputT2 = TissueClassifyDir + "/t2_average_BRAINSABC.nii.gz"
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print outputT2



    ScanDir = os.path.dirname( os.path.dirname( TissueClassifyDir ))
    outputInitialTransform = ScanDir + "/ACPCAlign/landmarkInitializer_atlas_to_subject_transform.mat"

    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print ScanDir
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print outputInitialTransform

    result = {'t1':inputT1, 't2':outputT2, 'transform':outputInitialTransform}
    # --------------------------------------------------------------------------------------- #
    # subject workflow
    #
    scan = os.path.basename( ScanDir )
    subjectID = os.path.basename(  os.path.dirname( ScanDir ) )
    siteID    = os.path.basename(  os.path.dirname( os.path.dirname( ScanDir ) ))
    ## find site in case of TRACK
    #subjectID_date_postFix  = os.path.basename( ScanDir )   # ex) 427997062_20090826_30

    #siteID    = os.path.basename( os.path.dirname( os.path.dirname( ScanDir )))
    #subjectID = subjectID_date_postFix.split( '_' )[0]
    #scan      = subjectID_date_postFix.split( '_' )[1]

    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print "site::" + siteID +", subjectID::" + subjectID + ",scan::"+ scan
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

    import WFPerSubject   
    
    #inputTemplateDir = "/ipldev/scratch/eunyokim/src/BRAINSStandAlone/build20121015/ReferenceAtlas-build/Atlas/Atlas_20120830/"

    WFPerSubject.WFPerSubjectDef( result,
                                  BRAINSAtlasDir,
                                  outputDir+"_Cache/"+siteID+"/"+subjectID+"/"+scan, 
                                  outputDir+"_Result/"+siteID+"/"+subjectID+"/"+scan)
# --------------------------------------------------------------------------------------- }#

#{--------------------------------------------------------------------------------------- #
def main(argv=None):

    import nipype.pipeline.engine as pe
    from nipype.interfaces.utility import Function
    import nipype.interfaces.io as nio
    import argparse
    
    if argv == None:
      argv = sys.argv
    # PARSE Command Line Arguments
    #{ ------------------------------------------------------------------------------------- #
    argParser = argparse.ArgumentParser( description = 'Run Statistical Analysis')
    group = argParser.add_argument_group('Required')

    group.add_argument( '--runOption', action="store", 
                        dest='runOption', required=True,
                            help='Basic will run on your local machine while cluster uses OsX')
    group.add_argument( '--inputBAWDir', action="store", 
                            dest='inputBAWDir', required=True,
                            help='The input BAW directory to look for t1/t2 average images')
    group.add_argument( '--outputBaseDir', action="store", 
                            dest='outputBaseDir', required=True,
                            help='The output base directory to store results/caches from nipype')
    group.add_argument( '--BRAINSStandAloneBuildDir', action="store",
                            dest='BRAINSStandAloneBuildDir', required=True, 
                            help='BRAINS stand alone BUILD directory to be used')
    group.add_argument( '--BRAINSStandAloneSrcDir', action="store",
                            dest='BRAINSStandAloneSrcDir', required=True, 
                            help='BRAINS stand alone SOURCE directory to be used')
    group.add_argument( '--PythonBinDir', action="store",
                            dest='PythonBinDir', required=True, 
                            help='Python Binary directory to be used')
    group.add_argument( '--BRAINSAtlasDir', action="store",
                            dest='BRAINSAtlasDir', required=True, 
                            help='BRAINS Atlas Directory')
    group.add_argument( '--NipypeBinDir', action="store",
                            dest='NipypeBinDir', required=False, 
                            help='Nipype Binary directory to be used')
    group.add_argument( '--NipypeLibDir', action="store",
                            dest='NipypeLibDir', required=False, 
                            help='Nipype Libary directory to be used')
    group.add_argument( '--inputSubjectListFile', action="store",
                        dest='inputSubjectListFile', required=True,
                        help='Input sujbect list file')

    inputArg = argParser.parse_args()
    print( inputArg.PythonBinDir)

    # ------------------------------------------------------------------------------------- }#
    # --------------------------------------------------------------------------------------- #
    # workflow
    #
    # outputBaseDir = "/scratch/eunyokim/LabelStatistics/RobustStats/TrackOn_Analysis/"
    # outputBaseDir = "/nfsscratch/PREDICT/regina/LabelStatistics/RobustStats/TrackOn_Analysis/"
    
    myWF = pe.Workflow(name="Analysis")
    myWF.base_dir = inputArg.outputBaseDir+"/_Cache"
    
    # --------------------------------------------------------------------------------------- #
    # input subject list file 
    #
    subjectList = []

    import csv
    subjectListFile = open( inputArg.inputSubjectListFile, 'rb')
    reader = csv.reader( subjectListFile )
    nrow = 0
    for row in reader:
        if nrow == 0:
            header=row
        else:
            ncol = 0
            for col in row:
                print '%-8s: %s' % ( header[ncol],col)
                if header[ncol] == "subject":
                  subjectList.append( col );
                ncol += 1
        nrow += 1
    
    print( subjectList )

    # --------------------------------------------------------------------------------------- #
    # data src : looking for t1_average images
    #
    
    dataSrc =  nio.DataGrabber( infields = ['subjectid'], 
                                outfields = ['t1'])
    
    #dataSrc.inputs.base_directory = "/hjohnson/TrackOn/Experiments/TrackOn_2012_Results/"
    dataSrc.inputs.base_directory =  inputArg.inputBAWDir
    dataSrc.inputs.template = '*'
    dataSrc.inputs.field_template = dict( t1 = '*/%s/*/TissueClassify/BABC//t1_average_BRAINSABC.nii.gz' )
    dataSrc.inputs.subjectid = subjectList
    results=dataSrc.run()
    
    t1List = [i for i in iter_flatten( results.outputs.t1)]
    print t1List
    # --------------------------------------------------------------------------------------- #
    
    # --------------------------------------------------------------------------------------- #
    #connect discovered dataSrc to subject workflow
    #
    
    findRestInputsFromT1 = pe.MapNode( name = "FindInputsForSubject",
                                    interface= Function ( input_names = ['inputT1','BRAINSAtlasDir','outputDir'],
                                                          output_names = [ 'outputListOfSubjectVolumes' ],
                                                          function = findInputImagesForSubject ),
                                    iterfield=['inputT1']
                                  )
    findRestInputsFromT1.inputs.outputDir = inputArg.outputBaseDir 
    #findRestInputsFromT1.iterables = ('inputT1', results.outputs.t1) # run the pipeline for each subject
    findRestInputsFromT1.inputs.inputT1 = t1List 
    findRestInputsFromT1.inputs.BRAINSAtlasDir = inputArg.BRAINSAtlasDir
    
    # if cluster
    #BAWSrcDir="/nfsscratch/PREDICT/regina/src/BRAINSStandAloneSrc//BRAINSStandAlone/"
    #BAWBuildDir="/nfsscratch/PREDICT/regina/src/BRAINSStandAloneSrc/build_20121016/"
    
    pythonPath = inputArg.BRAINSStandAloneSrcDir + "/BRAINSCut/BRAINSFeatureCreators/RobustStatisticComputations:" + inputArg.BRAINSStandAloneSrcDir + "/AutoWorkup/:" + inputArg.BRAINSStandAloneSrcDir + "/AutoWorkup/BRAINSTools/:" + inputArg.BRAINSStandAloneBuildDir + "/SimpleITK-build/bin/" + inputArg.BRAINSStandAloneBuildDir + "/SimpleITK-build/lib:" + inputArg.PythonBinDir
    #+ ":" + inputArg.NipypeBinDir + ":" + inputArg.NipypeLibDir  
    binPath = inputArg.BRAINSStandAloneBuildDir + "/bin:" + inputArg.BRAINSStandAloneBuildDir+ "/lib"

    myWF.add_nodes( [findRestInputsFromT1] )
    ############################################

    # Platform specific information
    #     Prepend the python search paths
    PYTHON_AUX_PATHS= pythonPath
    PYTHON_AUX_PATHS=PYTHON_AUX_PATHS.split(':')                                                                                  
    PYTHON_AUX_PATHS.extend(sys.path)                                                                                             
    sys.path=PYTHON_AUX_PATHS                                                                                                     
    #print sys.path
    import SimpleITK as sitk
    #     Prepend the shell environment search paths
    PROGRAM_PATHS= binPath
    PROGRAM_PATHS=PROGRAM_PATHS.split(':')
    PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))                                                                           
    os.environ['PATH']=':'.join(PROGRAM_PATHS)

    ###########################################
    
    if inputArg.runOption == "cluster":
        Cluster_Script = get_global_sge_script( PYTHON_AUX_PATHS, 
                                                PROGRAM_PATHS,
                                                {}
                                              )
        myWF.run( plugin='SGE',
                  plugin_args = dict( template = Cluster_Script, 
                                      qsub_args = "-S /bin/bash -pe smp1 4-8 -o /dev/null -q OSX "))
    else:
        myWF.run()
#--------------------------------------------------------------------------------------- }#

if __name__ == "__main__":
    sys.exit(main())  



# Python Path ex. from Hans

#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages/nibabel-1.2.0-py2.7.egg:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages/pydicom-0.9.7-py2.7.egg:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages/pudb-2012.2-py2.7.egg:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages/urwid-1.0.1-py2.7-linux-x86_64.egg:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages/SimpleITK-0.5.0-py2.7.egg:
#/opt/ortheus/src/python:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python27.zip:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/plat-linux2:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/lib-tk:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/lib-old:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/lib-dynload:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages/PIL:
#/nfsscratch/PREDICT/opt/epd-7.2-1-rh5-x86_64/lib/python2.7/site-packages/IPython/extensions
