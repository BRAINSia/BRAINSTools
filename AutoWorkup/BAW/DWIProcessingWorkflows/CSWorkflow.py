## \author Ali Ghayoor
##
## This workflow runs "compressed Sensing" in Matlab on a DWI scan,
## that is already corrected and aligned to t2 image space.
##
"""
CSWorkflow.py
============================
Description:
    The purpose of this is to...
    
Author:

Usage:

"""
import os

import nipype
import nipype.interfaces.io as nio  # Data i/oS
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces import ants
from nipype.interfaces.base import (
    CommandLine,
    CommandLineInputSpec,
    TraitedSpec,
    File,
    Directory,
)
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface


def create_cs_from_workflow(WFname, PYTHON_AUX_PATHS):
    """
    This Function runs "compressed Sensing" in Matlab on a DWI scan, that is already correct and aligned to t2 image spaces

    :param WFname:
    :param PYTHON_AUX_PATHS:
    :return:
    """
    if PYTHON_AUX_PATHS is not None:
        Path_to_Matlab_Func = os.path.join(
            PYTHON_AUX_PATHS[0], "DWIProcessingWorkflows"
        )
        assert os.path.exists(Path_to_Matlab_Func), (
            "Path to CS matlab function is not found: %s" % Path_to_Matlab_Func
        )

    #### Utility function ####
    def run_cs_by_matlab(inputScan, inputMask, CSScanFileName, Path_to_Matlab_Func):
        """
        This Function takes in...

        :param inputScan:
        :param inputMask:
        :param CSScanFileName:
        :param Path_to_Matlab_Func:
        :return: outputCSFilename
        """
        import os
        import nipype.interfaces.matlab as matlab

        script = (
            "runCS('" + inputScan + "','" + inputMask + "','" + CSScanFileName + "')"
        )
        mlab = matlab.MatlabCommand()
        mlab.set_default_matlab_cmd("matlab")
        mlab.inputs.single_comp_thread = False
        mlab.inputs.nodesktop = True
        mlab.inputs.nosplash = True
        mlab.inputs.paths = Path_to_Matlab_Func
        mlab.inputs.script = script
        mlab.run()
        outputCSFilename = os.path.join(
            os.getcwd(), CSScanFileName
        )  # return output CS filename
        assert os.path.isfile(outputCSFilename), (
            "CS file is not found: %s" % outputCSFilename
        )
        return outputCSFilename

    #########################

    CSWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(
        interface=IdentityInterface(fields=["DWI_Corrected_Aligned", "DWIBrainMask"]),
        name="inputsSpec",
    )

    outputsSpec = pe.Node(
        interface=IdentityInterface(fields=["DWI_Corrected_Aligned_CS"]),
        name="outputsSpec",
    )

    # Running a matlab node directly could be a nicer way, but I doesn't have any outputspec,
    # so I chose to run the matlab code using a python function.
    """
    runCS=pe.Node(interface=matlab.MatlabCommand(),name="run_cs_by_matlab")
    matlab.MatlabCommand().set_default_matlab_cmd("matlab")
    runCS.inputs.single_comp_thread = False
    runCS.inputs.nodesktop = True
    runCS.inputs.nosplash = True
    runCS.inputs.paths = Path_to_Matlab_Func
    CSWF.connect(createMatlabScriptNode,'matlabScript',runCS,'script')
    CSWF.connect(runCS,'out',outputsSpec,'DWI_Corrected_Aligned_CS') #This line cause problem
    """
    runCS = pe.Node(
        interface=Function(
            function=run_cs_by_matlab,
            input_names=[
                "inputScan",
                "inputMask",
                "CSScanFileName",
                "Path_to_Matlab_Func",
            ],
            output_names=["outputCSFilename"],
        ),
        name="runCS",
    )
    runCS.inputs.Path_to_Matlab_Func = Path_to_Matlab_Func
    runCS.inputs.CSScanFileName = "DWI_Corrected_Aligned_CS.nrrd"
    CSWF.connect(
        [
            (
                inputsSpec,
                runCS,
                [("DWI_Corrected_Aligned", "inputScan"), ("DWIBrainMask", "inputMask")],
            )
        ]
    )
    CSWF.connect(runCS, "outputCSFilename", outputsSpec, "DWI_Corrected_Aligned_CS")

    return CSWF
