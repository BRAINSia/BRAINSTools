#!/usr/bin/env python
"""
crossValidate.py
================
This is used to cross validate Multi-Atlas Label Fusion using ANTs JointFusion tool

Usage:
  crossValidate.py -h | --help
  crossValidate.py -t | --test
  crossValidate.py SIZE [--header] FILE

Arguments:
  SIZE          Sample size
  FILE          comma-seperated file in the format:
                     labelmap, t1_average, t2_average
Options:
  -h, --help    Show this help and exit
  -t, --test    Run doctests
  --header      Give this flag if CSV file has a header line


PIPELINE
        CSVreader() <--- FILE
          -> createTests(iterables=[(sampleLists, [...]), (...)]  <--- SIZE
          -> MapNode(workflow)
             -> SelectTest/SampleT1s, SelectTest/SampleT2s, SelectTest/SampleLabels
             -> JointFusion(), iterfield=[sample_list, test_list])
             -> DataSink()
"""
import os.path

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.base import BaseInterface, traits, TraitedSpec, DynamicTraitedSpec, File, BaseInterfaceInputSpec, isdefined
from nipype.interfaces.io import add_traits, SelectFiles, DataSink
from nipype.interfaces.utility import Select, Merge, Split, Function, Rename, IdentityInterface
from nipype.interfaces.ants.segmentation import JointFusion


class CSVReaderInputSpec(DynamicTraitedSpec, TraitedSpec):
    in_file = File(exists=True, mandatory=True, desc='Input comma-seperated value (CSV) file')
    header = traits.Bool(False, usedefault=True, desc='True if the first line is a column header')


class CSVReader(BaseInterface):
    """
    Example
    -------

    >>> def testNoHeader():
    ...   import os.path
    ...   import crossValidate as cv
    ...   reader = cv.CSVReader()
    ...   reader.inputs.in_file = os.path.abspath('test_noHeader.csv')
    ...   # reader.inputs.header = False
    ...   out = reader.run()
    ...   assert len(out.outputs.column_0) == 4
    ...   for output in out.outputs.column_0:
    ...       assert os.path.basename(output) == 't1_average_BRAINSABC.nii.gz'
    ...   assert len(out.outputs.column_1) == 4
    ...   for output in out.outputs.column_1:
    ...       assert os.path.basename(output) == 'neuro2012_20fusion_merge_seg.nii.gz'
    ...   assert len(out.outputs.column_2) == 4
    ...   for output in out.outputs.column_2:
    ...       assert os.path.basename(output) == 't2_average_BRAINSABC.nii.gz'

    >>> testNoHeader()

    >>> def testHeader():
    ...   import os.path
    ...   import crossValidate as cv
    ...   reader = cv.CSVReader()
    ...   reader.inputs.in_file = os.path.abspath('test_Header.csv')
    ...   reader.inputs.header = True
    ...   out = reader.run()
    ...   assert len(out.outputs.t1) == 4
    ...   for output in out.outputs.t1:
    ...       assert os.path.basename(output) == 't1_average_BRAINSABC.nii.gz'
    ...   assert len(out.outputs.t2) == 4
    ...   for output in out.outputs.t2:
    ...       assert os.path.basename(output) == 't2_average_BRAINSABC.nii.gz'
    ...   assert len(out.outputs.label) == 4
    ...   for output in out.outputs.label:
    ...       assert os.path.basename(output) == 'neuro2012_20fusion_merge_seg.nii.gz'

    >>> testHeader()


    """
    input_spec = CSVReaderInputSpec
    output_spec = DynamicTraitedSpec
    _always_run = True

    def _append_entry(self, outputs, entry):
        for key, value in zip(self._outfields, entry):
            outputs[key].append(value)
        return outputs

    def _parse_line(self, line):
        line = line.replace('\n', '')
        entry = [x.strip() for x in line.split(',')]
        return entry

    def _get_outfields(self):
        with open(self.inputs.in_file, 'r') as fid:
            entry = self._parse_line(fid.readline())
            if self.inputs.header:
                self._outfields = tuple(entry)
            else:
                self._outfields = tuple(['column_' + str(x) for x in range(len(entry))])
        return self._outfields

    def _run_interface(self, runtime):
        self._get_outfields()
        return runtime

    def _outputs(self):
        return self._add_output_traits(super(CSVReader, self)._outputs())

    def _add_output_traits(self, base):
        return add_traits(base, self._get_outfields())

    def _list_outputs(self):
        outputs = self.output_spec().get()
        isHeader = True
        for key in self._outfields:
            outputs[key] = []  # initialize outfields
        with open(self.inputs.in_file, 'r') as fid:
            for line in fid.readlines():
                if self.inputs.header and isHeader:  # skip header line
                    isHeader = False
                    continue
                entry = self._parse_line(line)
                outputs = self._append_entry(outputs, entry)
        return outputs




def sample_test_lists(in_list, size):
    """
    >>> print zip(*sample_test_lists(range(10), 2))  #doctest: +NORMALIZE_WHITESPACE
    [([0, 1], [2, 3, 4, 5, 6, 7, 8, 9]),
     ([2, 3], [0, 1, 4, 5, 6, 7, 8, 9]),
     ([4, 5], [0, 1, 2, 3, 6, 7, 8, 9]),
     ([6, 7], [0, 1, 2, 3, 4, 5, 8, 9]),
     ([8, 9], [0, 1, 2, 3, 4, 5, 6, 7])]

    >>> print zip(*sample_test_lists(range(9), 3))  #doctest: +NORMALIZE_WHITESPACE
    [([0, 1, 2], [3, 4, 5, 6, 7, 8]),
     ([3, 4, 5], [0, 1, 2, 6, 7, 8]),
     ([6, 7, 8], [0, 1, 2, 3, 4, 5])]
    """
    size = int(size)
    samples = list()
    tests = list()
    length = len(in_list)
    base_sample = range(size)
    for x in range(0, length, size):
        test = [y + x for y in base_sample]
        sample = range(length)
        for y in test:
            try:
                sample.remove(y)
            except ValueError:
                raise ValueError("List size is not evenly divisible by N({0})".format(size))
        samples.append(sample)
        tests.append(test)
    return tests, samples


class CrossValidationWorkflow(Workflow):
    """ Nipype workflow for Multi-Label Atlas Fusion cross-validation experiment """
    csv_file = None
    hasHeader = None
    sample_size = None

    def __init__(self, csv_file=None, size=0, hasHeader=False, name='CrossValidationWorkflow', **kwargs):
        super(CrossValidationWorkflow, self).__init__(name=name, **kwargs)
        self.csv_file = File(value=os.path.abspath(csv_file), exists=True)
        self.hasHeader = traits.Bool(hasHeader)
        self.sample_size = traits.Int(size)
        self.config['execution'] = {'remove_unnecessary_outputs': 'false'}

    def create(self):  #, **kwargs):
        """ Create the nodes and connections for the workflow """
        # Preamble
        csvReader = CSVReader()
        csvReader.inputs.in_file = self.csv_file.default_value
        csvReader.inputs.header = self.hasHeader.default_value
        csvOut = csvReader.run()

        iters = {}
        label = csvOut.outputs.__dict__.keys()[0]
        result = eval("csvOut.outputs.{0}".format(label))
        iters['tests'], iters['samples'] = sample_test_lists(result, self.sample_size.default_value)
        # Main event
        out_fields = ['T1', 'T2', 'Label', 'sampleindex', 'testindex']
        inputs = Node(interface=IdentityInterface(fields=out_fields),
                       run_without_submitting=True, name='inputs')
        inputs.iterables = [('sampleindex', iters['samples']),
                             ('testindex', iters['tests'])]
        if not self.hasHeader.default_value:
            inputs.inputs.T1 = csvOut.outputs.column_0
            inputs.inputs.Label = csvOut.outputs.column_1
            inputs.inputs.T2 = csvOut.outputs.column_2
        else:
            pass #TODO
        metaflow = Workflow(name='metaflow')
        metaflow.add_nodes([inputs])
        import pdb; pdb.set_trace()
        fusionflow = FusionLabelWorkflow()
        self.connect([(metaflow, fusionflow, [('inputs.sampleindex', 'sampleT1s.index'), ('inputs.T1',    'sampleT1s.inlist')]),
                      (metaflow, fusionflow, [('inputs.sampleindex', 'sampleT2s.index'), ('inputs.T2',    'sampleT2s.inlist')]),
                      (metaflow, fusionflow, [('inputs.sampleindex', 'sampleLabels.index'), ('inputs.Label', 'sampleLabels.inlist')]),
                      (metaflow, fusionflow, [('inputs.testindex',   'testT1s.index'), ('inputs.T1',    'testT1s.inlist')]),
                      (metaflow, fusionflow, [('inputs.testindex',   'testT2s.index'), ('inputs.T2',    'testT2s.inlist')]),
                      (metaflow, fusionflow, [('inputs.testindex',   'testLabels.index'), ('inputs.Label', 'testLabels.inlist')])
                      ])

    # def _connect_subworkflow(self):
    #     labelFusion = MapNode(FusionLabelWorkflow(),
    #                           iterfield=['sampleindex', 'testindex'],
    #                           name='FusionLabelWorkflow')
    #     self.connect([(self.get_node('csvReader'), labelFusion.inputspec, [('column_0', 'T1'),
    #                                                                        ('column_2', 'T2'),
    #                                                                        ('column_1', 'Label')]),
    #                   (self.get_node('createTests'), labelFusion.inputspec, [('samples', 'sampleindex'),
    #                                                                          ('tests', 'testindex')])
    #                  ])
    def _connect_subworkflow(self, node):
        self.connect(createTests, 'samples', node, 'sampleindex')
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
        sampleT1s    = Node(interface=Select(), name='sampleT1s')
        sampleT2s    = Node(interface=Select(), name='sampleT2s')
        sampleLabels = Node(interface=Select(), name='sampleLabels')

        testT1s      = Node(interface=Select(), name='testT1s')
        testT2s      = Node(interface=Select(), name='testT2s')
        testLabels   = Node(interface=Select(), name='testLabels')

        intensityImages = Node(interface=Merge(2), name='intensityImages')

        jointFusion = Node(interface=JointFusion(), name='jointFusion')
        jointFusion.inputs.dimension = 3
        jointFusion.inputs.modalities = 1  #TODO: verify 2 for T1/T2
        jointFusion.inputs.method = 'Joint[0.1, 2]'
        jointFusion.inputs.output_label_image = 'fusion_neuro2012_20.nii.gz'

        outputs = Node(interface=IdentityInterface(fields=['output_label_image']),
                       run_without_submitting=True, name='outputspec')

        self.connect([# Don't worry about T2s now per Regina
                      # (sampleT1s, intensityImages, [('out', 'in1')]),
                      # (sampleT2s, intensityImages, [('out', 'in2')]),
                      # (intensityImages, jointFusion, [('out', 'warped_intensity_images')]),
                      (sampleT1s, jointFusion, [('out', 'warped_intensity_images')]),
                      #END: per Regina
                      (sampleLabels, jointFusion, [('out', 'warped_label_images')]),
                      (jointFusion, outputs, [('output_label_image', 'output_label_image')]),
                      ])
        ## output => jointFusion.outputs.output_label_image


def main(**kwargs):
    workflow = CrossValidationWorkflow(csv_file=kwargs['FILE'], size=kwargs['SIZE'], hasHeader=kwargs['--header'])
    workflow.create()
    return workflow.write_graph()


def _test():
    import doctest
    doctest.testmod(verbose=True)
    return 0


if __name__ == "__main__":
    import sys

    from docopt import docopt

    argv = docopt(__doc__, version='1.1')
    print argv
    print '=' * 100
    if argv['--test']:
        sys.exit(_test())
    # from AutoWorkup import setup_environment
    # environment, experiment, pipeline, cluster = setup_environment(argv)
    sys.exit(main(**argv))
