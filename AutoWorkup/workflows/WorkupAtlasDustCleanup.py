#!/usr/bin/env python

from nipype.interfaces.utility import Function, IdentityInterface
import nipype.pipeline.engine as pe  # pypeline engine

def MakeVector(inFN1, inAtlas, outAtlas, maxIslandCount,
               useFullyConnected, forceLabelChange,
               inFN2=None, includeList=None, excludeList=None):
    arguments = {'--inputT1Path': inFN1,
                 '--inputT2Path': inFN2,
                 '--inputAtlasPath': inAtlas,
                 '--outputAtlasPath': outAtlas,
                 '--maximumIslandVoxelCount': maxIslandCount,
                 '--useFullyConnectedInConnectedComponentFilter': useFullyConnected,
                 '--forceSuspiciousLabelChange': forceLabelChange,
                 '--includeLabelsList': includeList,
                 '--excludeLabelsList': excludeList
               }
    return arguments

def runAutomaticCleanupScript(arguments):
  from atlasSmallIslandCleanup import DustCleanup
  print arguments
  localDustCleanupObject = DustCleanup(arguments=arguments)
  localDustCleanupObject.main()
  return arguments['--outputAtlasPath']

def CreateDustCleanupWorkflow(WFname, onlyT1, master_config,
                              runFixFusionLabelMap=True):

    if onlyT1:
      n_modality = 1
    else:
      n_modality = 2
    CLUSTER_QUEUE = master_config['queue']
    CLUSTER_QUEUE_LONG = master_config['long_q']

    dustCleanupWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['subj_t1_image', #Input T1 image
                                                             'subj_t2_image', #Input T2 image
                                                             'subj_label_atlas' #Input label atlas image
                                                            ]),
                         run_without_submitting=True,
                         name='inputspec')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['dustCleanedOutputAtlas_label']),
                          run_without_submitting=True,
                          name='outputspec')

    """
    multimodal atlas dust cleanup if t2 exists
    """
    sessionMakeMultimodalInput = pe.Node(Function(function=MakeVector,
                                                  input_names=['inFN1', 'inFN2', 'inAtlas', 'outAtlas',
                                                               'maxIslandCount', 'useFullyConnected',
                                                               'forceLabelChange', 'noDilation',
                                                               'includeList', 'excludeList'],
                                                  output_names=['arguments']),
                                run_without_submitting=True, name="sessionMakeMultimodalInput")
    dustCleanupWF.connect(inputsSpec, 'subj_t1_image', sessionMakeMultimodalInput, 'inFN1')
    if not onlyT1:
        dustCleanupWF.connect(inputsSpec, 'subj_t2_image', sessionMakeMultimodalInput, 'inFN2')
    else:
        pass

    dustCleanupWF.connect(inputsSpec, 'subj_label_atlas', sessionMakeMultimodalInput, 'inAtlas')
    sessionMakeMultimodalInput.inputs.outAtlas = 'dustCleanedLabel_maxIsland5.nii.gz'
    sessionMakeMultimodalInput.inputs.maxIslandCount = 5
    sessionMakeMultimodalInput.inputs.useFullyConnected = True
    sessionMakeMultimodalInput.inputs.forceLabelChange = True
    sessionMakeMultimodalInput.inputs.noDilation = False
    sessionMakeMultimodalInput.inputs.includeList = None
    sessionMakeMultimodalInput.inputs.excludeList = '4,5,14,15,21,24,31,43,44,63,72,85,98,128,219,15000'


    sessionRunDustCleanupInput = pe.Node(Function(function=runAutomaticCleanupScript,
                                                  input_names=['arguments'],
                                                  output_names=['cleanedLabelImage']),
                                run_without_submitting=True, name="sessionRunDustCleanupInput")

    dustCleanupWF.connect(sessionMakeMultimodalInput, 'arguments', sessionRunDustCleanupInput, 'arguments')

    dustCleanupWF.connect(sessionRunDustCleanupInput, 'cleanedLabelImage', outputsSpec, 'dustCleanedOutputAtlas_label')

    return dustCleanupWF
