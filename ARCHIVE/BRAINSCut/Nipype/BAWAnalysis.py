from builtins import str

## Analysis interface for BAW version
## input: segmentations from BCut
##        manual traces
## output: a set of similarity measures
##         mean/median of them
##         ICC measures
##############################################################################
##############################################################################

import computeSimilarityWF

#########################################################################################
def get_global_sge_script(pythonPathsList, binPathsList, customEnvironment={}):
    print("""get_global_sge_script""")
    """This is a wrapper script for running commands on an SGE cluster
so that all the python modules and commands are pathed properly"""

    custEnvString = ""
    for key, value in list(customEnvironment.items()):
        custEnvString += "export " + key + "=" + value + "\n"

    PYTHONPATH = ":".join(pythonPathsList)
    BASE_BUILDS = ":".join(binPathsList)
    GLOBAL_SGE_SCRIPT = """#!/bin/bash
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
""".format(
        PYTHONPATH=PYTHONPATH, BINPATH=BASE_BUILDS, CUSTENV=custEnvString
    )
    return GLOBAL_SGE_SCRIPT


#########################################################################################


def writeCSV(dataDictList, outputCSVFilename):
    import csv

    f = open(outputCSVFilename, "wb")
    w = csv.DictWriter(f, list(dataDictList[0].keys()))
    w.writeheader()
    for row in dataDictList:
        w.writerow(row)
    f.close()

    return outputCSVFilename


#########################################################################################


def getData(
    ResultDir, normalization, methodParameter, sessionID, optionalString=""
):  # ex = ANNLabel_seg.nii.gz
    import nipype.interfaces.io as nio
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    #
    # get label map
    #
    DG = nio.DataGrabber(
        infields=["normalization", "method", "sessionID"], outfields=["outputList"]
    )
    DG.inputs.base_directory = ResultDir
    DG.inputs.template = "Test*/%s/RF/%s/%s*" + optionalString
    DG.inputs.template_args = OrderedDict(
        outputList=[["normalization", "method", "sessionID"]]
    )
    DG.inputs.normalization = normalization
    DG.inputs.method = methodParameter
    DG.inputs.sessionID = sessionID
    dt = DG.run()
    print(
        (
            """Data grabber with {opt}:::
           {str}
           """.format(
                opt=optionalString, str=dt.outputs.outputList
            )
        )
    )

    return dt.outputs.outputList


#########################################################################################
def experimentAnalysis(
    resultDir,
    outputCSVFilename,
    normalization,
    methodParameter,
    manualDict,
    roiList,
    sessionIDList,
    doRawComparison,
):
    import nipype.interfaces.io as nio
    import ast
    import analysis as this

    roiLabel = {}
    label = 1
    for roi in sorted(set(roiList)):
        roiLabel[roi] = label
        print(("""{r} ===> {l}""".format(r=roi, l=label)))
        label = label + 1

    autoFileList = []
    defFileList = []
    refFileList = []
    autoLabelList = []
    identityDictList = []
    roiSeqList = []
    sessionList = []
    for sessionID in sessionIDList:
        #
        # get label map
        #

        roiManualDict = ast.literal_eval(manualDict[sessionID]["roiList"])
        identityDict = {}
        for roi in roiList:
            identityDict["sessionID"] = sessionID
            identityDict["roi"] = roi
            identityDictList.append(identityDict)
            roiSeqList.append(roi)
            sessionList.append(sessionID)

            if doRawComparison:
                labelDT = this.getData(
                    resultDir,
                    normalization,
                    methodParameter,
                    sessionID,
                    roi + "_ANNContinuous.nii.gz",
                )
                autoLabelList.append(1)

            else:
                labelDT = this.getData(
                    resultDir,
                    normalization,
                    methodParameter,
                    sessionID,
                    "ANNLabel_seg.nii.gz",
                )
                autoLabelList.append(roiLabel[roi])

            autoFileList.append(labelDT)

            print(
                (
                    """compute roi of ::
                   {s}
                   """.format(
                        s=roi
                    )
                )
            )
            #
            # get get deformation
            #
            defFileList.append(
                this.getData(
                    resultDir,
                    normalization,
                    methodParameter,
                    sessionID,
                    roi + "*.nii.gzdef.nii.gz",
                )
            )
            # read manual image
            refFileList.append(roiManualDict[roi])

    print(
        (
            """LENGTH:::
        length( roiSeqList) = {roiL}
        length( sessionList ) = {sessionL}
        length( autoFileList ) = {autuL}
        length( refFilename ) = {refL}
    """.format(
                roiL=len(roiSeqList),
                sessionL=len(sessionList),
                autuL=len(autoFileList),
                refL=len(refFileList),
            )
        )
    )

    import nipype.pipeline.engine as pe
    import os

    workFlowName = "experimentAnalysis"

    exp = pe.Workflow(name=workFlowName)
    outputCSVFilename = os.path.abspath(outputCSVFilename)
    exp.base_dir = os.path.dirname(outputCSVFilename)

    from nipype.interfaces.utility import Function

    computeSimilarityND = pe.MapNode(
        name="computeSimilarityND",
        interface=Function(
            input_names=[
                "autoFilename",
                "refFilename",
                "autoLabel",
                "roi",
                "session",
                "defFilename",
            ],
            output_names=["outDict"],
            function=computeSimilarityWF.computeSimilarity,
        ),
        iterfield=[
            "autoFilename",
            "defFilename",
            "refFilename",
            "autoLabel",
            "roi",
            "session",
        ],
    )
    computeSimilarityND.inputs.autoFilename = autoFileList
    computeSimilarityND.inputs.defFilename = defFileList
    computeSimilarityND.inputs.refFilename = refFileList
    computeSimilarityND.inputs.autoLabel = autoLabelList
    computeSimilarityND.inputs.roi = roiSeqList
    computeSimilarityND.inputs.session = sessionList

    exp.add_nodes([computeSimilarityND])

    writeCSVFileND = pe.Node(
        name="writeCSVFileND",
        interface=Function(
            input_names=["dataDictList", "outputCSVFilename"],
            output_names=["outputCSVFilename"],
            function=this.writeCSV,
        ),
    )
    writeCSVFileND.inputs.outputCSVFilename = outputCSVFilename
    exp.connect(computeSimilarityND, "outDict", writeCSVFileND, "dataDictList")

    exp.run()
    return outputCSVFilename


#########################################################################################


def similarityComputeWorkflow(
    ResultDir,
    OutputDir,
    ExperimentalConfigurationFile,
    runOption,
    PythonBinDir,
    BRAINSToolsSrcDir,
    BRAINSToolsBuildDir,
):

    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering
    import sys
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    #
    # get normalization option from experimental configuration file
    #
    import ConfigurationParser as configParser
    import analysis as this

    configMap = configParser.ConfigurationSectionMap(ExperimentalConfigurationFile)
    normalizationOptions = configMap["Options"]["normalization"]
    print(
        (
            """ Normalization Option:::
          {str}
          """.format(
                str=normalizationOptions
            )
        )
    )

    #
    # get methods
    #
    import ast

    methodOptionsDictList = configMap["Options"]["modelParameter".lower()]
    methodOptions = []
    print(methodOptionsDictList)
    for option in methodOptionsDictList:
        methodStr = (
            "TreeDepth"
            + str(option["--randomTreeDepth"])
            + "_TreeNumber"
            + str(option["--numberOfTrees"])
        )
        methodOptions.append(methodStr)
    print(
        (
            """ Method Option:::
          {str}
          """.format(
                str=methodStr
            )
        )
    )

    #
    # get roiList
    #
    roiList = list(configMap["Options"]["roiBooleanCreator".lower()].keys())
    print(
        (
            """ ROIList:::
          {str}
          """.format(
                str=roiList
            )
        )
    )

    #
    # get sessionList and manualDict
    #
    import XMLConfigurationGenerator

    subjectListFilename = configMap["ListFiles"]["subjectListFilename".lower()]
    manualDict = XMLConfigurationGenerator.combineCSVs(subjectListFilename, {})
    sessionList = list(manualDict.keys())

    #
    # workflow
    #
    workFlowName = "outputDataCollector"
    workflow = pe.Workflow(name=workFlowName)
    workflow.base_dir = OutputDir

    from nipype.interfaces.utility import Function

    experimentalND = pe.Node(
        name="experimentalND",
        interface=Function(
            input_names=[
                "resultDir",
                "outputCSVFilename",
                "normalization",
                "methodParameter",
                "manualDict",
                "roiList",
                "sessionIDList",
                "doRawComparison",
            ],
            output_names="outputCSVFilename",
            function=this.experimentAnalysis,
        ),
    )
    experimentalND.inputs.resultDir = ResultDir
    experimentalND.inputs.outputCSVFilename = "experimentalResult.csv"
    experimentalND.inputs.roiList = roiList
    experimentalND.inputs.manualDict = manualDict
    experimentalND.inputs.sessionIDList = sessionList
    # experimentalND.inputs.doRawComparison = doRawComparison
    experimentalND.iterables = [
        ("normalization", normalizationOptions),
        ("methodParameter", methodOptions),
        ("doRawComparison", [True, False]),
    ]
    workflow.add_nodes([experimentalND])

    summaryND = pe.Node(
        name="summaryND",
        interface=Function(
            input_names=["inputCSVFilename", "outputCSVPrefix"],
            output_names=["outputCSVList"],
            function=computeSimilarityWF.computeSummaryFromCSV,
        ),
    )

    summaryND.inputs.outputCSVPrefix = "summaryOutput"
    workflow.connect(experimentalND, "outputCSVFilename", summaryND, "inputCSVFilename")

    if runOption == "cluster":
        ############################################
        # Platform specific information
        #     Prepend the python search paths
        pythonPath = (
            BRAINSToolsSrcDir
            + "/BRAINSCut/BRAINSFeatureCreators/RobustStatisticComputations:"
            + BRAINSToolsSrcDir
            + "/AutoWorkup/:"
            + BRAINSToolsSrcDir
            + "/AutoWorkup/BRAINSTools/:"
            + BRAINSToolsBuildDir
            + "/SimpleITK-build/bin/"
            + BRAINSToolsBuildDir
            + "/SimpleITK-build/lib:"
            + PythonBinDir
        )
        binPath = BRAINSToolsBuildDir + "/bin:" + BRAINSToolsBuildDir + "/lib"

        PYTHON_AUX_PATHS = pythonPath
        PYTHON_AUX_PATHS = PYTHON_AUX_PATHS.split(":")
        PYTHON_AUX_PATHS.extend(sys.path)
        sys.path = PYTHON_AUX_PATHS
        # print sys.path
        import SimpleITK as sitk

        #     Prepend the shell environment search paths
        PROGRAM_PATHS = binPath
        PROGRAM_PATHS = PROGRAM_PATHS.split(":")
        import os

        PROGRAM_PATHS.extend(os.environ["PATH"].split(":"))
        os.environ["PATH"] = ":".join(PROGRAM_PATHS)

        Cluster_Script = get_global_sge_script(PYTHON_AUX_PATHS, PROGRAM_PATHS, {})
        workflow.run(
            plugin="SGE",
            plugin_args=OrderedDict(
                template=Cluster_Script,
                qsub_args="-S /bin/bash -pe smp 4-8 -o /dev/null ",
            ),
        )
    else:
        workflow.run()


def main(argv=None):
    import os
    import sys

    from nipype import config

    config.enable_debug_mode()

    # -------------------------------- argument parser
    import argparse

    argParser = argparse.ArgumentParser(
        description="""****************************
        10-cross validation analysis
        """
    )
    # workup arguments
    argWfGrp = argParser.add_argument_group(
        "argWfGrp",
        """****************************
        auto workflow arguments for cross validation
        """,
    )
    argWfGrp.add_argument(
        "--experimentalConfigurationFile",
        help="""experimentalConfigurationFile
        Configuration file name with FULL PATH""",
        dest="experimentalConfigurationFile",
        required=True,
    )
    argWfGrp.add_argument(
        "--expDir",
        help="""expDir
        """,
        dest="expDir",
        required=False,
        default=".",
    )
    argWfGrp.add_argument(
        "--baseDir",
        help="""baseDir
        """,
        dest="baseDir",
        required=False,
        default=".",
    )
    argWfGrp.add_argument(
        "--runOption",
        help="""runOption [local/cluster]
        """,
        dest="runOption",
        required=False,
        default="local",
    )
    argWfGrp.add_argument(
        "--PythonBinDir",
        help="""PythonBinDir [local/cluster]
        """,
        dest="PythonBinDir",
        required=False,
        default="NA",
    )
    argWfGrp.add_argument(
        "--BRAINSToolsSrcDir",
        help="""BRAINSToolsSrcDir [local/cluster]
        """,
        dest="BRAINSToolsSrcDir",
        required=False,
        default="NA",
    )
    argWfGrp.add_argument(
        "--BRAINSToolsBuildDir",
        help="""BRAINSToolsBuildDir [local/cluster]
        """,
        dest="BRAINSToolsBuildDir",
        required=False,
        default="NA",
    )

    args = argParser.parse_args()
    similarityComputeWorkflow(
        args.expDir,
        args.baseDir,
        args.experimentalConfigurationFile,
        args.runOption,
        args.PythonBinDir,
        args.BRAINSToolsSrcDir,
        args.BRAINSToolsBuildDir,
    )


import sys

if __name__ == "__main__":
    sys.exit(main())
#########################################################################################
# unit test
#
# ResultDir = "/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec/Experiment_20121222/Result/Labels/"
# outputDir = '/ipldev/scratch/eunyokim/src/BRAINSTools/BRAINSTools/BRAINSCut/Nipype/'
#
# outputCSVFilename = '/ipldev/scratch/eunyokim/src/BRAINSTools/BRAINSTools/BRAINSCut/Nipype/output.csv'
# normalization = 'Linear'
# methodParameter = 'TreeDepth50_TreeNumber50'
# configFilename = '/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec/Dec22/model.config'
#
# similarityComputeWorkflow( ResultDir, outputDir, configFilename )
#
# import ConfigurationParser as configParser
# import XMLConfigurationGenerator
# configMap = configParser.ConfigurationSectionMap( configFilename )
# subjectListFilename = configMap[ 'ListFiles' ]['subjectListFilename'.lower() ]
# manualDict  = XMLConfigurationGenerator.combineCSVs( subjectListFilename, {} )
# sessionIDList = manualDict.keys()
# print( sessionIDList )
#
# testDictList = [ {'roi':'l_accumben','vol':1000,'NF':3},
#                 {'roi':'r_accumben','vol':100,'NF':3} ]
# writeCSV(  testDictList, 'out.csv' )
# roiList =  [
#            'l_thalamus'
#           ]
# experimentAnalysis( ResultDir,
#                    outputCSVFilename,
#                    normalization,
#                    methodParameter,
#                    manualDict,
#                    roiList,
#                    sessionIDList )
#
