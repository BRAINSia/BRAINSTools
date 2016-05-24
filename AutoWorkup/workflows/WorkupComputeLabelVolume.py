#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from nipype.interfaces.base import CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined
from nipype.interfaces.utility import Merge, Function, Rename, IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine

from utilities.measureVolumes import *
from utilities.misc import *

def CreateVolumeMeasureWorkflow( WFname, master_config):
    volumeMeasureWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['subj_t1_image', #Input T1 image
                                                             'subj_label_image' #Input Label image
                                                            ]),
                         run_without_submitting=True,
                         name='inputspec')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['csvFilename','jsonFilename']),
                          run_without_submitting=True,
                          name='outputspec')

    """
    Measure volumes according to
    1) label image
    2) label look up table (following format for 3D Slicer color lookup table)
    3) Reference image

    and produce measured volumes in
    1) CSV format
    2) JSON Format
    """

    volumeMeasureWF = pe.Workflow(name=WFname)
    makeDictND = pe.Node(Function(function=MakeLabelDictionary,
                                  input_names=['inputColorLookUpTableFilename'],
                                  output_names = ['labelDictionary']),
                         run_without_submitting = True,
                         name = 'makeLabelDict')
    makeDictND.inputs.inputColorLookUpTableFilename = master_config['labelmap_colorlookup_table']

    getVolumesND = pe.Node(Function(function=GetLabelVolumes,
                                    input_names = ['labelVolume','RefVolume','labelDictionary'],
                                    output_names = ['outputLabelVolumes']),
                           run_without_submitting = False,
                           name='getVolumes')
    volumeMeasureWF.connect( makeDictND, 'labelDictionary',
                             getVolumesND, 'labelDictionary')
    volumeMeasureWF.connect( inputsSpec, 'subj_t1_image',
                             getVolumesND, 'RefVolume')
    volumeMeasureWF.connect( inputsSpec, 'subj_label_image',
                             getVolumesND, 'labelVolume')

    writeCSVND = pe.Node(Function(function=WriteDictionaryToCSV,
                                  input_names = ['inputList', 'outputFilename'],
                                   output_names = ['outputFilename']),
                         run_without_submitting = False,
                         name='writeCSV')
    volumeMeasureWF.connect( getVolumesND, 'outputLabelVolumes',
                             writeCSVND, 'inputList')
    writeCSVND.inputs.outputFilename = 'labelVolume.csv'
    volumeMeasureWF.connect( writeCSVND, 'outputFilename',
                             outputsSpec, 'csvFilename')




    writeJSONND = pe.Node(Function(function=WriteDictionaryToJson,
                                  input_names = ['inputList', 'outputFilename'],
                                   output_names = ['outputFilename']),
                         run_without_submitting = False,
                         name='writeJSON')
    volumeMeasureWF.connect( getVolumesND, 'outputLabelVolumes',
                             writeJSONND, 'inputList')
    writeJSONND.inputs.outputFilename = 'labelVolume.json'
    volumeMeasureWF.connect( writeJSONND, 'outputFilename',
                             outputsSpec, 'jsonFilename')

    return volumeMeasureWF
