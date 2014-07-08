#! /usr/bin/env python
"""
longitudinal.py
=========
This program is used to generate the subject- and session-specific workflows for BRAINSTool processing

Usage:
  longitudinal.py [--rewrite-datasinks] [--wfrun PLUGIN] --subject ID... --pe ENV --ExperimentConfig FILE
  longitudinal.py -v | --version
  longitudinal.py -h | --help

Arguments:


Options:
  -h, --help            Show this help and exit
  -v, --version         Print the version and exit
  --rewrite-datasinks   Turn on the Nipype option to overwrite all files in the 'results' directory
  --pe ENV              The processing environment to use from configuration file
  --subject ID          The subject ID to process
  --wfrun PLUGIN        The name of the workflow plugin option (default: 'local')
  --ExperimentConfig FILE   The configuration file

Examples:
  $ longitudinal.py --subject 1058 --pe OSX --ExperimentConfig my_baw.config
  $ longitudinal.py --wfrun helium_all.q --subject 1058 --pe OSX --ExperimentConfig my_baw.config
  $ longitudinal.py --rewrite-datasinks --subject 1058 --pe OSX --ExperimentConfig my_baw.config

"""

def create_longitudinal(project, subject, session, master_config, interpMode='Linear', pipeline_name=''):
    """
    create longitudinal workflow on a single session

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """
    from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, Directory
    from nipype.interfaces.base import traits, isdefined, BaseInterface
    from nipype.interfaces.utility import Split, Rename, IdentityInterface, Function
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio

    from baseline import baseline_workflow as create_baseline

    baw201 = create_baseline(project, subject, session, master_config,
                               phase='longitudinal',
                               interpMode=interpMode,
                               pipeline_name=pipeline_name)
    inputsSpec = baw201.get_node('inputspec')

    if 'segmentation' in master_config['components']:
        from workflows.segmentation import segmentation
        sname = 'segmentation'
        # sname = GenerateWFName(project, subject, session, 'segmentation')
        onlyT1 = not(len(inputsSpec.inputs.T2s) > 0)
        segWF = segmentation(project, subject, session, master_config, onlyT1, pipeline_name=sname)
        outputSpec = baw201.get_node('outputspec')
        baw201.connect([(outputSpec, segWF, [('t1_average', 'inputspec.t1_average'),
                                             ('LMIatlasToSubject_tx', 'inputspec.LMIatlasToSubject_tx'),
                                             ('outputLabels', 'inputspec.inputLabels'),
                                             ('posteriorImages', 'inputspec.posteriorImages'),
                                             ('tc_atlas2sessionInverse_tx',
                                              'inputspec.TissueClassifyatlasToSubjectInverseTransform'),
                                             ('UpdatedPosteriorsList', 'inputspec.UpdatedPosteriorsList'),
                                             ('outputHeadLabels', 'inputspec.inputHeadLabels')])
                                ])
        if not onlyT1:
            baw201.connect([(outputSpec, segWF, [('t1_average', 'inputspec.t2_average')])])

    return baw201

def main():
    # from atlasNode import MakeAtlasNode
    # atlasNode = MakeAtlasNode(master_config['atlascache'], 'BAtlas')
    # workflow = create_longitudial(...)
    # workflow.connect([(atlas -> baseline), (atlas -> segmentation)]) # See singleSubject for more details...
    pass


if __name__ == '__main__':
    from docopt import docopt
    from AutoWorkup import setup, run

    argv = docopt(__doc__, version='1.1')
    # print argv
    # print '=' * 100
    configs = setup(argv)
    exit = run(argv, *configs)
    sys.exit(exit)
