#!/usr/bin/env python
#
# TODO
# :: copy model file into md5 repository
# :: connect input/output in the BAW
"""
BRAINSCutCMD.py
=================
Description:

Author:

Usage:

"""


import argparse
import subprocess
import sys

if False:
    print(("^" * 100))
    print(("PYTHON EXEC: ", sys.executable))
    print(("VERSION: ", sys.version))
    # HACK
    dummy = False
    if sys.version[:3] != "2.7":
        print(__file__)
        dummy = True
    # END HACK
    print("ARGV: [\n")
    print((sys.argv[0]))
    for __argv in sys.argv[1:]:
        if __argv.startswith("-"):
            print(("    {0}= ".format(__argv)))
        else:
            print(__argv)
    print("]")
    del __argv
    print(("^" * 100))

    if dummy:  # debugging
        import os

        print("PATH: [")
        __PATH = os.environ["PATH"]
        for p in __PATH.split(":"):
            print((p, ","))
        print("]")
        __PPATH = ""
        try:
            print("PYTHONPATH: [")
            __PPATH = os.environ["PYTHONPATH"]
            for p in __PPATH.split(":"):
                print((p, ","))
            print("]")
        except KeyError:
            pass
        del __PATH, __PPATH
        print(("^" * 100))
        raise NotImplementedError
        sys.exit(-1)
        # END HACK


def add_probability_map_element(probabilityMap, maskName, outputStream):
    """
    This function...

    :param probabilityMap:
    :param maskName:
    :param outputStream:
    :return:
    """
    outputStream.write('  <ProbabilityMap StructureID    = "' + maskName + '"\n')
    outputStream.write('      Gaussian       = "1.0"\n')
    outputStream.write('      GenerateVector = "true"\n')
    outputStream.write('      Filename       = "' + probabilityMap + '"\n')
    outputStream.write("   />\n")


def xml_generator(args, roi=""):
    """
    This function...

    :param args:
    :param roi:
    :return:
    """
    xmlFilename = args.xmlFilename + roi + ".xml"
    outputStream = open(xmlFilename, "w")
    registrationID = "BSpline_ROI"

    outputStream.write("<AutoSegProcessDescription>\n")

    #
    # template
    #
    outputStream.write('  <DataSet Name="template" Type="Atlas" >\n')
    outputStream.write(
        '      <Image Type="T1" Filename="{fn}" />\n'.format(fn=args.inputTemplateT1)
    )
    if args.inputSubjectT2Filename is not None:
        outputStream.write(
            '      <Image Type="T2" Filename="{fn}" />\n'.format(fn="na")
        )
        outputStream.write(
            '      <Image Type="GadSG" Filename="{fn}" />\n'.format(fn="na")
        )
    # outputStream.write( "      <Image Type=\"TotalGM\" Filename=\"{fn}\" />\n".format(fn="na"))
    # outputStream.write( "      <Mask  Type=\"RegistrationROI\" Filename=\"{fn}\" />\n".format(fn=args.inputTemplateRegistrationROIFilename))

    outputStream.write(
        '      <SpatialLocation Type="rho" Filename="'
        + args.inputTemplateRhoFilename
        + '" />\n'
    )
    outputStream.write(
        '      <SpatialLocation Type="phi" Filename="'
        + args.inputTemplatePhiFilename
        + '" />\n'
    )
    outputStream.write(
        '      <SpatialLocation Type="theta" Filename="'
        + args.inputTemplateThetaFilename
        + '" />\n'
    )
    outputStream.write("  </DataSet>\n")

    #
    # Registration
    #
    outputStream.write("  <RegistrationConfiguration \n")
    outputStream.write('          ImageTypeToUse  = "T1"\n')
    outputStream.write('          ID              = "' + registrationID + '"\n')
    outputStream.write('          BRAINSROIAutoDilateSize= "1"\n')
    outputStream.write("   />\n")

    #
    # training vector configuration  (feature vector)
    #

    outputStream.write('   <NeuralNetParams MaskSmoothingValue     = "0.0"\n')
    outputStream.write('          GradientProfileSize    = "1"\n')
    outputStream.write(
        '          TrainingVectorFilename = "' + args.trainingVectorFilename + '"\n'
    )
    outputStream.write('          TrainingModelFilename  = "na"\n')
    outputStream.write('          TestVectorFilename     = "na"\n')
    if roi == "caudate":
        outputStream.write('          Normalization          = "' + "Linear" + '"\n')
    else:
        outputStream.write(
            '          Normalization          = "' + args.vectorNormalization + '"\n'
        )

    outputStream.write("   />\n")

    #
    # random forest parameters
    #
    outputStream.write("   <RandomForestParameters \n")
    outputStream.write('      MaxDepth= "1"\n')  # dummy
    outputStream.write('      MaxTreeCount= "1"\n')  # dummy
    outputStream.write('      MinSampleCount= "5"\n')
    outputStream.write('      UseSurrogates= "false"\n')
    outputStream.write('      CalcVarImportance= "false"\n')
    outputStream.write("      />\n")

    #
    # ANN Parameters
    #
    outputStream.write('   <ANNParameters Iterations             = "5"\n')
    outputStream.write('                     MaximumVectorsPerEpoch = "70000"\n')
    outputStream.write('                     EpochIterations        = "100"\n')
    outputStream.write('                     ErrorInterval          = "1"\n')
    outputStream.write('                     DesiredError           = "0.000001"\n')
    outputStream.write('                     NumberOfHiddenNodes    = "100"\n')
    outputStream.write('                     ActivationSlope        = "1.0"\n')
    outputStream.write('                     ActivationMinMax       = "1.0"\n')
    outputStream.write("    />\n")

    #
    # apply conditions
    #
    outputStream.write('<ApplyModel         CutOutThresh           = "0.05"\n')
    outputStream.write('                    MaskThresh             = "0.5"\n')
    outputStream.write('                    GaussianSmoothingSigma = "0.0"\n')
    outputStream.write("   />\n")

    #
    # add probability maps (ROIs)
    #
    if roi == "caudate":
        add_probability_map_element(
            args.probabilityMapsLeftCaudate, "l_caudate", outputStream
        )
        add_probability_map_element(
            args.probabilityMapsRightCaudate, "r_caudate", outputStream
        )
    elif roi == "putamen":
        add_probability_map_element(
            args.probabilityMapsLeftPutamen, "l_putamen", outputStream
        )
        add_probability_map_element(
            args.probabilityMapsRightPutamen, "r_putamen", outputStream
        )
    elif roi == "thalamus":
        add_probability_map_element(
            args.probabilityMapsLeftThalamus, "l_thalamus", outputStream
        )
        add_probability_map_element(
            args.probabilityMapsRightThalamus, "r_thalamus", outputStream
        )
    elif roi == "hippocampus":
        add_probability_map_element(
            args.probabilityMapsLeftHippocampus, "l_hippocampus", outputStream
        )
        add_probability_map_element(
            args.probabilityMapsRightHippocampus, "r_hippocampus", outputStream
        )
    elif roi == "accumben":
        add_probability_map_element(
            args.probabilityMapsLeftAccumben, "l_accumben", outputStream
        )
        add_probability_map_element(
            args.probabilityMapsRightAccumben, "r_accumben", outputStream
        )
    elif roi == "globus":
        add_probability_map_element(
            args.probabilityMapsLeftGlobus, "l_globus", outputStream
        )
        add_probability_map_element(
            args.probabilityMapsRightGlobus, "r_globus", outputStream
        )

    #
    # subject
    #
    outputStream.write('  <DataSet Name="subject" Type="Apply"')
    outputStream.write('      OutputDir="./" >\n')
    outputStream.write(
        '    <Image Type="T1" Filename="' + args.inputSubjectT1Filename + '" />\n'
    )
    if args.inputSubjectT2Filename is not None:
        outputStream.write(
            '    <Image Type="T2" Filename="' + args.inputSubjectT2Filename + '" />\n'
        )
        outputStream.write(
            '    <Image Type="GadSG" Filename="'
            + args.inputSubjectGadSGFilename
            + '" />\n'
        )
    if roi == "caudate":
        outputStream.write(
            '    <Image Type="candiateRegion" Filename="'
            + args.candidateRegion
            + '" />\n'
        )
    # outputStream.write( "    <Image Type=\"TotalGM\" Filename=\"{fn}\" />\n".format(fn=args.inputSubjectTotalGMFilename))
    # outputStream.write( "    <Mask  Type=\"RegistrationROI\" Filename=\"{fn}\" />\n".format(fn=args.inputSubjectRegistrationROIFilename))

    # outputStream.write( "    <Mask Type=\"l_caudate\" Filename=\""+args.outputBinaryLeftCaudate+"\" />\n")
    # outputStream.write( "    <Mask Type=\"r_caudate\" Filename=\""+args.outputBinaryRightCaudate+"\" />\n")
    # outputStream.write( "    <Mask Type=\"l_putamen\" Filename=\""+args.outputBinaryLeftPutamen+"\" />\n")
    # outputStream.write( "    <Mask Type=\"r_putamen\" Filename=\""+args.outputBinaryRightPutamen+"\" />\n")
    # outputStream.write( "    <Mask Type=\"l_thalamus\" Filename=\""+args.outputBinaryLeftThalamus+"\" />\n")
    # outputStream.write( "    <Mask Type=\"r_thalamus\" Filename=\""+args.outputBinaryRightThalamus+"\" />\n")
    # outputStream.write( "    <Mask Type=\"l_hippocampus\" Filename=\""+args.outputBinaryLeftHippocampus+"\" />\n")
    # outputStream.write( "    <Mask Type=\"r_hippocampus\" Filename=\""+args.outputBinaryRightHippocampus+"\" />\n")
    # outputStream.write( "    <Mask Type=\"l_accumben\" Filename=\""+args.outputBinaryLeftAccumben+"\" />\n")
    # outputStream.write( "    <Mask Type=\"r_accumben\" Filename=\""+args.outputBinaryRightAccumben+"\" />\n")
    # outputStream.write( "    <Mask Type=\"l_globus\" Filename=\""+args.outputBinaryLeftGlobus+"\" />\n")
    # outputStream.write( "    <Mask Type=\"r_globus\" Filename=\""+args.outputBinaryRightGlobus+"\" />\n")

    # if args.inputSubjectBrainMaskFilename != "NA":
    #    outputStream.write( "    <Mask Type=\"RegistrationROIi\"  Filename=\""+args.inputSubjectBrainMaskFilename+"\" />\n")

    if args.deformationFromSubjectToTemplate is not None:
        outputStream.write(
            '    <Registration SubjToAtlasRegistrationFilename="'
            + args.deformationFromSubjectToTemplate
            + '"\n'
        )
    else:
        outputStream.write('    <Registration SubjToAtlasRegistrationFilename="" \n')
    outputStream.write(
        '       AtlasToSubjRegistrationFilename="'
        + args.deformationFromTemplateToSubject
        + '"\n'
    )
    outputStream.write('       ID="' + registrationID + '" /> \n')
    outputStream.write("  </DataSet>\n")

    outputStream.write("</AutoSegProcessDescription>\n")
    outputStream.close()

    return xmlFilename


if __name__ == "__main__":
    brainscutParser = argparse.ArgumentParser(
        description="BRAINSCut command line argument parser"
    )

    # HACK:  This is to allow special treatment of caudates with masking
    brainscutParser.add_argument(
        "--candidateRegion",
        help="Specify the valid candidate region for caudate",
        required=True,
    )

    #
    # input arguments
    #
    brainscutParser.add_argument(
        "--inputSubjectT1Filename", help="T1 subject filename", required=True
    )
    brainscutParser.add_argument(
        "--inputSubjectT2Filename", help="T2 subject filename", required=False
    )
    # brainscutParser.add_argument('--inputSubjectTotalGMFilename', help='TotalGM filename', required=True )
    brainscutParser.add_argument(
        "--inputSubjectGadSGFilename", help="GadSG subject filename", required=False
    )
    # brainscutParser.add_argument('--inputSubjectBrainMaskFilename', help='BrainMask subject filename' )
    # brainscutParser.add_argument('--inputSubjectRegistrationROIFilename', help='template brain mask filename' )

    brainscutParser.add_argument(
        "--inputTemplateT1", help="template T1-weighted filename", required=True
    )
    # brainscutParser.add_argument('--inputTemplateRegistrationROIFilename', help='template brain region filename', required=True )

    brainscutParser.add_argument(
        "--inputTemplateRhoFilename", help="template rho filename", required=True
    )
    brainscutParser.add_argument(
        "--inputTemplatePhiFilename", help="template phi filename", required=True
    )
    brainscutParser.add_argument(
        "--inputTemplateThetaFilename", help="template theta filename", required=True
    )

    brainscutParser.add_argument(
        "--trainingVectorFilename", help="training vector filename", default="NA"
    )
    # brainscutParser.add_argument('--modelFileBasename', help='model filei base name for net configuration file (xml).', default="NA" )
    brainscutParser.add_argument(
        "--modelFilename", help="model filename", default="NA", required=True
    )
    brainscutParser.add_argument(
        "--vectorNormalization",
        help="feature vector normalization (IQR,Linear,Sigmoid_Q01,Sigmoid_Q05,ZScore,NONE)",
        required=True,
    )

    # probability maps
    brainscutParser.add_argument(
        "--probabilityMapsLeftCaudate",
        help="model probability maps for left caudate",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsRightCaudate",
        help="model probability maps for right caudate",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsLeftPutamen",
        help="model probability maps for left putamen",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsRightPutamen",
        help="model probability maps for right putamen",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsLeftThalamus",
        help="model probability maps for left thalamus",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsRightThalamus",
        help="model probability maps for right thalamus",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsLeftHippocampus",
        help="model probability maps for left hippocampus",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsRightHippocampus",
        help="model probability maps for right hippocampus",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsLeftAccumben",
        help="model probability maps for left accumben",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsRightAccumben",
        help="model probability maps for right accumben",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsLeftGlobus",
        help="model probability maps for left globus",
        required=True,
    )
    brainscutParser.add_argument(
        "--probabilityMapsRightGlobus",
        help="model probability maps for right globus",
        required=True,
    )

    brainscutParser.add_argument(
        "--deformationFromTemplateToSubject", help="deformationFromTemplateToSubject"
    )
    brainscutParser.add_argument(
        "--deformationFromSubjectToTemplate", help="deformationFromSubjectToTemplate"
    )

    #
    # output arguments
    #
    brainscutParser.add_argument(
        "--outputBinaryLeftCaudate", help="output binary file name for left caudate"
    )
    brainscutParser.add_argument(
        "--outputBinaryRightCaudate", help="output binary file name for right caudate"
    )
    brainscutParser.add_argument(
        "--outputBinaryLeftPutamen", help="output binary file name for left putamen"
    )
    brainscutParser.add_argument(
        "--outputBinaryRightPutamen", help="output binary file name for right putamen"
    )
    brainscutParser.add_argument(
        "--outputBinaryLeftThalamus", help="output binary file name for left thalamus"
    )
    brainscutParser.add_argument(
        "--outputBinaryRightThalamus", help="output binary file name for right thalamus"
    )
    brainscutParser.add_argument(
        "--outputBinaryLeftHippocampus",
        help="output binary file name for left hippocampus",
    )
    brainscutParser.add_argument(
        "--outputBinaryRightHippocampus",
        help="output binary file name for right hippocampus",
    )
    brainscutParser.add_argument(
        "--outputBinaryLeftAccumben", help="output binary file name for left accumben"
    )
    brainscutParser.add_argument(
        "--outputBinaryRightAccumben", help="output binary file name for right accumben"
    )
    brainscutParser.add_argument(
        "--outputBinaryLeftGlobus", help="output binary file name for left globus"
    )
    brainscutParser.add_argument(
        "--outputBinaryRightGlobus", help="output binary file name for right globus"
    )

    brainscutParser.add_argument(
        "--xmlFilename",
        help="BRAINSCut xml configuration filename",
        default="output.xml",
    )

    args = brainscutParser.parse_args()
    ## HACK:  DOUBLE CHECK THAT IQR IS USED
    if args.vectorNormalization != "IQR":
        print("ERROR:   ONLY IQR SUPPORTED AT THE MOMENT")
        exit - 1

    print(args)
    roiList = ["accumben", "caudate", "putamen", "globus", "thalamus", "hippocampus"]

    for roi in roiList:
        currentXmlFilename = xml_generator(args, roi)
        if roi == "caudate":
            currentModelFilename = (
                args.modelFilename[:-3] + "_" + roi + "_LinearWithMask.gz"
            )  # trainModelFile.txtD0060NT0060_caudate_LinearWithMask.gz
        else:
            currentModelFilename = (
                args.modelFilename[:-3] + "_" + roi + ".gz"
            )  # trainModelFile.txtD0060NT0060_accumben.gz

        BRAINSCutCommand = [
            "BRAINSCut"
            + " --applyModel "
            + " --netConfiguration "
            + currentXmlFilename
            + " --modelFilename "
            + currentModelFilename
            + " --method RandomForest"
            + " --numberOfTrees 60  --randomTreeDepth 60"
        ]
        print(("HACK:  BRAINCUT COMMAND: {0}".format(BRAINSCutCommand)))
        subprocess.check_call(BRAINSCutCommand, shell=True)
    """
    script to be run
      BRAINSCut  --applyModel --netConfiguration BRAINSTools-build/BRAINSCut/TestSuite/TestSuite/NetConfigurations/output.xml --modelFilename TrainedModels/20110919ANNModel_allSubcorticals.txtD0050NT0050   --method RandomForest
    """
