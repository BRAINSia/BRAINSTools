#!/usr/bin/env python
"""
crossValidate.py
================
This is used to cross validate Multi-Atlas Label Fusion using ANTs JointFusion tool

Usage:
  crossValidate.py -h | --help
  crossValidate.py -t | --test
  crossValidate.py --testSamplesize SIZE [--header] --sampleFile listFILE --workphase WORKPHASE [--wfrun PLGIN] --pe ENV  --ExperimentConfig FILE [--rewrite-datasinks]

Options:
  -h, --help    Show this help and exit
  -t, --test    Run doctests
  --testSamplesize=SIZE    a size of test sample. Should be less the the total size list in sample File
  --sampleFile=listFILE        A data list including all the information.
  --header                 Give this flag if CSV file has a header line
  --pe=ENV                 The processing environment to use from configuration file
  --wfrun=PLUGIN            The name of the workflow plugin option (default: 'local')
  --ExperimentConfig=FILE   The configuration file
  --workphase WORKPHASE The type of processing to be done [cross-validation]
  --rewrite-datasinks   Turn on the Nipype option to overwrite all files in the 'results' directory


PIPELINE
        CSVreader() <--- FILE
          -> createTests(iterables=[(sampleLists, [...]), (...)]  <--- SIZE
          -> MapNode(workflow)
             -> SelectTest/SampleT1s, SelectTest/SampleT2s, SelectTest/SampleLabels
             -> JointFusion(), iterfield=[sample_list, test_list])
             -> DataSink()
"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
import os.path

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.base import BaseInterface, traits, TraitedSpec, DynamicTraitedSpec, File, BaseInterfaceInputSpec, isdefined
from nipype.interfaces.io import add_traits, SelectFiles, DataSink
from nipype.interfaces.utility import Select, Merge, Split, Function, Rename, IdentityInterface
from nipype.interfaces.ants.segmentation import JointFusion


def subsample_crossValidationSet(length, test_size):
    """
    >>> print zip(*subsample_crossValidationSet(10, 2))  #doctest: +NORMALIZE_WHITESPACE
    [([0, 1], [2, 3, 4, 5, 6, 7, 8, 9]),
     ([2, 3], [0, 1, 4, 5, 6, 7, 8, 9]),
     ([4, 5], [0, 1, 2, 3, 6, 7, 8, 9]),
     ([6, 7], [0, 1, 2, 3, 4, 5, 8, 9]),
     ([8, 9], [0, 1, 2, 3, 4, 5, 6, 7])]

    >>> print zip(*subsample_crossValidationSet(9, 3))  #doctest: +NORMALIZE_WHITESPACE
    [([0, 1, 2], [3, 4, 5, 6, 7, 8]),
     ([3, 4, 5], [0, 1, 2, 6, 7, 8]),
     ([6, 7, 8], [0, 1, 2, 3, 4, 5])]
    """
    test_size = int(test_size)
    subsample_data_index = []
    base_train = list(range(test_size))
    for x in range(0, length, test_size):
        test = [y + x for y in base_train]
        train = list(range(length))
        for y in test:
            try:
                train.remove(y)
            except ValueError:
                raise ValueError("List test size is not evenly divisible by N({0})".format(test_size))
        subsample_data_index.append( {'train':train, 'test':test } )
    print("="*80)
    print(subsample_data_index)
    return subsample_data_index

def writeCVSubsetFile( environment, experiment, pipeline, cluster, csv_file, test_size, hasHeader):

    from utilities.misc import add_dict
    master_config = {}
    for configDict in [environment, experiment, pipeline, cluster]:
        master_config = add_dict(master_config, configDict)

    """
    read in csv file
    """
    import csv
    csv_data=[]
    with open(csv_file, mode='r') as infile:
          reader = csv.DictReader(infile, skipinitialspace=True)
          for row in reader:
            csv_data.append(row)
    print(csv_data)

    totalSampleSize = len(csv_data)
    print(totalSampleSize)
    cv_subsets = subsample_crossValidationSet( totalSampleSize, test_size)

    """
    global variable
    """
    ## HACK FOR NOW SHOULD BE MORE ELEGANT FROM THE .config file
    BASE_DATA_GRABBER_DIR='/Shared/johnsonhj/HDNI/ReferenceData/Neuromorphometrics/2012Subscription'
    #master_config = {'queue':'HJ',
    #    'long_q':'HJ'}

    """
    workflow
    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from .WorkupJointFusion import CreateJointFusionWorkflow
    CV_JointFusion_WF = pe.Workflow(name="CV_JointFusion")
    CV_JointFusion_WF.base_dir = master_config['cachedir']


    subset_no = 1
    for subset in cv_subsets:
        print("-"*80)
        print(" Creat a subset workflow Set " + str(subset_no))
        print("-"*80)
        trainData = [ csv_data[i] for i in subset['train'] ]
        testData = [ csv_data[i] for i in subset['test'] ]

        print([ (trainData[i])['id'] for i in range( len(trainData))])

        for testSession in testData:
            JointFusionWFName = "JointFusion_Set{0}_{1}".format(subset_no, testSession['id'])
            myJointFusion = CreateJointFusionWorkflow( JointFusionWFName,
                                         master_config,
                                         [ (trainData[i])['id'] for i in range( len(trainData))],
                                         BASE_DATA_GRABBER_DIR,
                                         runFixFusionLabelMap=False)

            testSessionName= "testSessionSpec_Set{0}_{1}".format(subset_no, testSession['id'])
            testSessionSpec = pe.Node( interface=IdentityInterface( fields=['t1_average',
                                                                         'tissueLabel',
                                                                         'template_leftHemisphere',
                                                                         'landmarkInACPCAlignedSpace',
                                                                         'template_weights_50Lmks_wts',
                                                                         'labelFilename']),
                                    run_without_submitting = True,
                                    name=testSessionName)

            CV_JointFusion_WF.connect(testSessionSpec,'t1_average', myJointFusion,'inputspec.subj_t1_image')
            CV_JointFusion_WF.connect(testSessionSpec,'tissueLabel',myJointFusion,'inputspec.subj_fixed_head_labels')

            CV_JointFusion_WF.connect(testSessionSpec,'template_leftHemisphere', myJointFusion,'inputspec.subj_left_hemisphere')
            CV_JointFusion_WF.connect(testSessionSpec,'landmarkInACPCAlignedSpace', myJointFusion,'inputspec.subj_lmks')
            CV_JointFusion_WF.connect(testSessionSpec,'template_weights_50Lmks_wts', myJointFusion,'inputspec.atlasWeightFilename')
            CV_JointFusion_WF.connect(testSessionSpec, 'labelFilename', myJointFusion, 'inputspec.labelBaseFilename')

            """ set test image information
            """
            print(testSession)
            testSessionSpec.inputs.t1_average = testSession['t1']
            testSessionSpec.inputs.tissueLabel = testSession['fixed_head_label']
            testSessionSpec.inputs.template_leftHemisphere = testSession['warpedAtlasLeftHemisphere']
            testSessionSpec.inputs.landmarkInACPCAlignedSpace = testSession['lmk']
            testSessionSpec.inputs.template_weights_50Lmks_wts = "/Shared/sinapse/scratch/eunyokim/src/NamicExternal/build_Mac_201501/bin/Atlas/Atlas_20131115/20141004_BCD/template_landmarks_50Lmks.fcsv"
            testSessionSpec.inputs.labelFilename='FS_wmparc.nii.gz'

            """
            DataSink
            """
            dsName = "DataSink_DS_Set{0}_{1}".format(subset_no,testSession['id'])
            DataSink = pe.Node(name=dsName, interface=nio.DataSink())
            DataSink.overwrite = master_config['ds_overwrite']
            DataSink.inputs.container = 'CV_Set{0}/{1}'.format(subset_no, testSession['id'])
            DataSink.inputs.base_directory = master_config['resultdir']

            CV_JointFusion_WF.connect(myJointFusion, 'outputspec.JointFusion_neuro2012_labelmap',
                               DataSink, 'Segmentation.@JointFusion_neuro2012_labelmap')

            subset_no=subset_no+1

    #CV_JointFusion_WF.write_graph()
    CV_JointFusion_WF.run( plugin=master_config['plugin_name'],
                    plugin_args=master_config['plugin_args'])



class CrossValidationJointFusionWorkflow(Workflow):
    """ Nipype workflow for Multi-Label Atlas Fusion cross-validation experiment """
    csv_file = None
    hasHeader = None
    sample_size = None

    def __init__(self, csv_file=None, size=0, hasHeader=False, name='CrossValidationJointFusionWorkflow', **kwargs):
        super(CrossValidationJointFusionWorkflow, self).__init__(name=name, **kwargs)
        self.csv_file = File(value=os.path.abspath(csv_file), exists=True)
        self.hasHeader = traits.Bool(hasHeader)
        self.sample_size = traits.Int(size)
        self.config['execution'] = {'remove_unnecessary_outputs': 'true'}

    def create(self):  #, **kwargs):
        """ Create the nodes and connections for the workflow """
        # Preamble
        csvReader = CSVReader()
        csvReader.inputs.in_file = self.csv_file.default_value
        csvReader.inputs.header = self.hasHeader.default_value
        csvOut = csvReader.run()

        print("="*80)
        print(csvOut.outputs.__dict__)
        print("="*80)

        iters = {}
        label = list(csvOut.outputs.__dict__.keys())[0]
        result = eval("csvOut.outputs.{0}".format(label))
        iters['tests'], iters['trains'] = subsample_crossValidationSet(result, self.sample_size.default_value)
        # Main event
        out_fields = ['T1', 'T2', 'Label', 'trainindex', 'testindex']
        inputsND = Node(interface=IdentityInterface(fields=out_fields),
                       run_without_submitting=True, name='inputs')
        inputsND.iterables = [('trainindex', iters['trains']),
                             ('testindex', iters['tests'])]
        if not self.hasHeader.default_value:
            inputsND.inputs.T1 = csvOut.outputs.column_0
            inputsND.inputs.Label = csvOut.outputs.column_1
            inputsND.inputs.T2 = csvOut.outputs.column_2
        else:
            inputsND.inputs.T1 = csvOut.outputs.__dict__['t1']
            inputsND.inputs.Label = csvOut.outputs.__dict__['label']
            inputsND.inputs.T2 = csvOut.outputs.__dict__['t2']
            pass #TODO
        metaflow = Workflow(name='metaflow')
        metaflow.config['execution'] = {
            'plugin': 'Linear',
            'stop_on_first_crash': 'false',
            'stop_on_first_rerun': 'false',  # This stops at first attempt to rerun, before running, and before deleting previous results.
            'hash_method': 'timestamp',
            'single_thread_matlab': 'true',  # Multi-core 2011a  multi-core for matrix multiplication.
            'remove_unnecessary_outputs': 'true',
            'use_relative_paths': 'false',  # relative paths should be on, require hash update when changed.
            'remove_node_directories': 'false',  # Experimental
            'local_hash_check': 'false'
        }

        metaflow.add_nodes([inputsND])
        """import pdb; pdb.set_trace()"""
        fusionflow = FusionLabelWorkflow()
        self.connect([(metaflow, fusionflow, [('inputs.trainindex', 'trainT1s.index'), ('inputs.T1',    'trainT1s.inlist')]),
                      (metaflow, fusionflow, [('inputs.trainindex', 'trainLabels.index'), ('inputs.Label', 'trainLabels.inlist')]),
                      (metaflow, fusionflow, [('inputs.testindex',  'testT1s.index'), ('inputs.T1',    'testT1s.inlist')])
                      ])

    # def _connect_subworkflow(self):
    #     labelFusion = MapNode(FusionLabelWorkflow(),
    #                           iterfield=['trainindex', 'testindex'],
    #                           name='FusionLabelWorkflow')
    #     self.connect([(self.get_node('csvReader'), labelFusion.inputspec, [('column_0', 'T1'),
    #                                                                        ('column_2', 'T2'),
    #                                                                        ('column_1', 'Label')]),
    #                   (self.get_node('createTests'), labelFusion.inputspec, [('trains', 'trainindex'),
    #                                                                          ('tests', 'testindex')])
    #                  ])
    def _connect_subworkflow(self, node):
        self.connect(createTests, 'trains', node, 'trainindex')
        self.connect(createTests, 'tests', node, 'testindex')
        self.connect(csvReader, 'T1', node.inputspec, 'T1')
        self.connect(csvReader, 'T2', node.inputspec, 'T2')
        self.connect(csvReader, 'LabelMaps', node.inputspec, 'Label')


class FusionLabelWorkflow(Workflow):
    """ Subworkflow to use with MapNode """
    def __init__(self, name='FusionLabelWorkflow', **kwargs):
        super(FusionLabelWorkflow, self).__init__(name=name, **kwargs)
        self.create()
        # self.connect = None  # Don't allow instances to add to workflow

    def connect(self, *args, **kwargs):
        try:
            super(FusionLabelWorkflow, self).connect(*args, **kwargs)
        except:
            from pprint import pprint
            pprint(args)
            pprint(kwargs)
            raise


    def create(self):
        trainT1s    = Node(interface=Select(), name='trainT1s')
        trainT2s    = Node(interface=Select(), name='trainT2s')
        trainLabels = Node(interface=Select(), name='trainLabels')
        testT1s      = Node(interface=Select(), name='testT1s')

        #intensityImages = Node(interface=Merge(2), name='intensityImages')

        jointFusion = Node(interface=JointFusion(), name='jointFusion')
        jointFusion.inputs.num_threads = -1
        jointFusion.inputs.dimension = 3
        jointFusion.inputs.modalities = 1  #TODO: verify 2 for T1/T2
        jointFusion.inputs.method =  "Joint[0.1,2]" # this does not work
        jointFusion.inputs.output_label_image = 'fusion_neuro2012_20.nii.gz'

        outputs = Node(interface=IdentityInterface(fields=['output_label_image']),
                       run_without_submitting=True, name='outputspec')

        self.connect([# Don't worry about T2s now per Regina
                      # (trainT1s, intensityImages, [('out', 'in1')]),
                      # (trainT2s, intensityImages, [('out', 'in2')]),
                      (testT1s, jointFusion, [('out', 'target_image')]),
                      (trainT1s, jointFusion, [('out', 'warped_intensity_images')]),
                      (trainLabels, jointFusion, [('out', 'warped_label_images')]),
                      (jointFusion, outputs, [('output_label_image', 'output_label_image')]),
                      ])
        ## output => jointFusion.outputs.output_label_image


def main(environment, experiment, pipeline, cluster, **kwargs):
    from utilities.configFileParser import nipype_options

    print("Copying Atlas directory and determining appropriate Nipype options...")
    pipeline = nipype_options(kwargs, pipeline, cluster, experiment, environment)  # Generate Nipype options
    print("Getting session(s) from database...")

    writeCVSubsetFile( environment,
                       experiment,
                       pipeline,
                       cluster,
                       csv_file=kwargs['--sampleFile'],
                       test_size=kwargs['--testSamplesize'],
                       hasHeader=kwargs['--header'] )


def _test():
    import doctest
    doctest.testmod(verbose=True)
    return 0


if __name__ == "__main__":
    import sys

    from docopt import docopt

    argv = docopt(__doc__, version='1.1')
    print(argv)
    print('=' * 100)
    from AutoWorkup import setup_environment
    if argv['--test']:
        sys.exit(_test())
    environment, experiment, pipeline, cluster = setup_environment(argv)
    # from AutoWorkup import setup_environment
    # environment, experiment, pipeline, cluster = setup_envir""
    sys.exit(main(environment, experiment, pipeline, cluster, **argv))
