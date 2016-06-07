#! /usr/bin/env python

"""
configFileParser.py
===================

Usage:
  configFileParser.py [--debug] ENV PHASE FILE
  configFileParser.py (-h | --help | -v | --version)

Arguments:
  ENV            Environment name in configuration file
  PHASE          Workflow phase to use: atlas-based-reference, subject-template-generation or subject-based-reference
  FILE           Configuration file

Options:
  -h, --help     Print this and exit
  -v, --version  Print file version and exit
  --debug        Run doctests for file  # TODO
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from past.utils import old_div
from builtins import object
from configparser import ConfigParser
import os
import sys
import copy

from .pathHandling import *
from .distributed import modify_qsub_args
from . import misc

# http://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python
def str2bool(v):
  if str(v).lower() in ("yes", "true", "t", "1"):
      return True
  elif str(v).lower() in ("no", "false", "f", "0"):
      return False
  raise ValueError( "ERROR: INVALID String to bool conversion for '{0}'".format(v) )

def parseEnvironment(parser, environment):
    """ Parse the environment environment given by 'section' and return a dictionary
        Values are shell-centric, i.e. PYTHONPATH is a colon-seperated string

    """
    retval = dict()
    if parser.has_option(environment, 'ENVAR_DICT'):
        retval['env'] = eval(parser.get(environment, 'ENVAR_DICT'))
    else:
        retval['env'] = dict()
    if 'PYTHONPATH' in list(retval['env'].keys()):
        pythonpath = appendPathList(parser.get(environment, 'APPEND_PYTHONPATH'), retval['env']['PYTHONPATH'])
        retval['env']['PYTHONPATH'] = pythonpath  # Create append to PYTHONPATH
    else:
        retval['env']['PYTHONPATH'] = parser.get(environment, 'APPEND_PYTHONPATH')
    if 'PATH' in list(retval['env'].keys()):
        envpath = appendPathList(parser.get(environment, 'APPEND_PATH'), retval['env']['PATH'])
        retval['env']['PATH'] = envpath  # Create append to PATH
    else:
        retval['env']['PATH'] = parser.get(environment, 'APPEND_PATH')

    retval['prefix'] = validatePath(parser.get(environment, 'MOUNT_PREFIX'), True, True)
    if retval['prefix'] is None:
        retval['prefix'] = ''
    if parser.has_option(environment, 'VIRTUALENV_DIR'):
        retval['virtualenv_dir'] = validatePath(parser.get(environment, 'VIRTUALENV_DIR'), False, True)
    else:
        retval['virtualenv_dir'] = None
    retval_cluster = dict()
    retval_cluster['modules'] = eval(parser.get(environment, 'MODULES'))
    retval_cluster['queue'] = parser.get(environment, 'QUEUE')
    retval_cluster['long_q'] = parser.get(environment, 'QUEUE_LONG')
    retval_cluster['qstat'] = parser.get(environment, 'QSTAT_IMMEDIATE')
    retval_cluster['qstat_cached'] = parser.get(environment, 'QSTAT_CACHED')

    return retval, retval_cluster


def create_experiment_dir(dirname, name, suffix, verify=False):
    """ Construct directories given the base directory, the experiment name, and the suffix ['CACHE', 'Results'] """
    basename = name + '_' + suffix
    fullpath = os.path.join(dirname, basename)
    if verify:
        return validatePath(fullpath, False, True)
    else:
        if os.path.isdir(fullpath):
            print("WARNING: Experiment directory already exists.  Continuing will overwrite the previous results...")
            print("   Path: {0}".format(fullpath))
            return fullpath
        try:
            os.makedirs(fullpath)
        except OSError:
            raise
    return fullpath


def parseExperiment(parser, workflow_phase):
    """ Parse the experiment section and return a dictionary """
    retval = dict()
    dirname = validatePath(parser.get('EXPERIMENT', 'BASE_OUTPUT_DIR'), False, True)
    if workflow_phase == 'atlas-based-reference':
        current_suffix = '_BASE'
    elif workflow_phase == 'subject-template-generation':
        current_suffix = '_TEMP'
    elif workflow_phase == 'subject-based-reference':
        current_suffix = '_LONG'
    elif workflow_phase == 'cross-validation':
        current_suffix = '_CV'
    else:
        assert 0 == 1, "ERROR INVALID workflow_phase"
    current = parser.get('EXPERIMENT', 'EXPERIMENT' + current_suffix)

    """ output directory """
    retval['cachedir'] = create_experiment_dir(dirname, current, 'CACHE')
    retval['resultdir'] = create_experiment_dir(dirname, current, 'Results')

    """ any previous run HACK: DO WE EVER USE THIS?"""
    if parser.has_option('EXPERIMENT', 'EXPERIMENT' + current_suffix + '_INPUT'):
        # If this is the initial run, there will be no previous experiment
        previous = parser.get('EXPERIMENT', 'EXPERIMENT' + current_suffix + '_INPUT')
        retval['previousresult'] = create_experiment_dir(dirname, previous, 'Results', verify=True)

    useRegistrationMasking = True
    try:
        regMasking = parser.get('EXPERIMENT', 'USE_REGISTRATION_MASKING')
        useRegistrationMasking = str2bool(regMasking)
    except:
        pass
    retval['use_registration_masking'] = useRegistrationMasking

    atlas = validatePath(parser.get('EXPERIMENT', 'ATLAS_PATH'), False, True)
    retval['atlascache'] = clone_atlas_dir(retval['cachedir'], atlas)

    if workflow_phase == 'cross-validation':
        retval['components'] = ['']
    else:
        retval['dbfile'] = validatePath(parser.get('EXPERIMENT', 'SESSION_DB' + current_suffix), False, False)
        retval['components'] = [x.lower() for x in eval(parser.get('EXPERIMENT', 'WORKFLOW_COMPONENTS' + current_suffix))]
        if 'jointfusion_2015_wholebrain' in retval['components']:
            print("'jointFusion_2015_wholebrain' will be run with a specified 'jointfusion_atlas_db_base'.")
            """ HACK: warp_atlas_to_subject is coupled with jointFusion????"""
            retval['jointfusion_atlas_db_base'] = validatePath(parser.get('EXPERIMENT', 'JointFusion_ATLAS_DB_BASE'),
                                                       allow_empty=False,
                                                       isDirectory=False)
            retval['labelmap_colorlookup_table'] = validatePath(parser.get('EXPERIMENT', 'LABELMAP_COLORLOOKUP_TABLE'),
                                                       allow_empty=False,
                                                       isDirectory=False)
            retval['relabel2lobes_filename'] = validatePath(parser.get('EXPERIMENT', 'RELABEL2LOBES_FILENAME'),
                                                       allow_empty=True,
                                                       isDirectory=False)
        retval['workflow_phase'] = workflow_phase
    return retval


def parseNIPYPE(parser):
    """ Parse the nipype section and return a dictionary """
    retval = dict()
    retval['ds_overwrite'] = parser.getboolean('NIPYPE', 'GLOBAL_DATA_SINK_REWRITE')

    if parser.has_option('NIPYPE', 'CRASHDUMP_DIR'):
        retval['CRASHDUMP_DIR'] = parser.get('NIPYPE', 'CRASHDUMP_DIR')
    else:
        retval['CRASHDUMP_DIR'] = None

    return retval


#def parseCluster(parser, env):
#    """ Parse the cluster section and return a dictionary """
#    retval = dict()
#    retval['modules'] = eval(parser.get(env, 'MODULES'))
#    retval['queue'] = parser.get(env, 'QUEUE')
#    retval['long_q'] = parser.get(env, 'QUEUE_LONG')
#    retval['qstat'] = parser.get(env, 'QSTAT_IMMEDIATE')
#    retval['qstat_cached'] = parser.get(env, 'QSTAT_CACHED')
#    return retval


def parseFile(configFile, env, workphase):

    configFile = os.path.realpath(configFile)
    assert os.path.exists(configFile), "Configuration file could not be found: {0}".format(configFile)
    parser = ConfigParser(allow_no_value=True)  # Parse configuration file parser = ConfigParser()
    parser.read(configFile)
    assert (parser.has_option(env, '_BUILD_DIR') or parser.has_option('DEFAULT', '_BUILD_DIR')
            ), "BUILD_DIR option not in {0}".format(env)
    environment, cluster = parseEnvironment(parser, env)
    experiment = parseExperiment(parser, workphase)
    pipeline = parseNIPYPE(parser)
    return environment, experiment, pipeline, cluster


def resolveDataSinkOption(args, pipeline):
    if args["--rewrite-datasinks"] or pipeline['ds_overwrite']:  # GLOBAL_DATA_SINK_REWRITE
        return True
    return False


class _create_DS_runner(object):

    def run(self, graph, **kwargs):
        for node in graph.nodes():
            if '_ds' in node.name.lower():
                node.run()


_WFRUN_VALID_TYPES = ['SGE',
                      'SGEGraph',
                      'local_4',
                      'local_12',
                      'local',
                      'ds_runner']


def get_cpus(option):
    assert option in _WFRUN_VALID_TYPES, "Unknown wfrun option"
    from multiprocessing import cpu_count
    import os
    total_cpus = cpu_count()
    suffix = option.rsplit('local', 1)[1]
    if suffix == '':
        assert option in ['local', 'ds_runner'], "wfrun parse error!  Current option: {0}".format(option)
        threads = 1
        if option == 'local':
            print("RUNNING WITHOUT POOL BUILDING")
    else:
        threads = int(suffix.strip('_'))
    return int(old_div(total_cpus, threads))


def _nipype_plugin_config(wfrun, cluster, template=''):
    assert wfrun in _WFRUN_VALID_TYPES, "Unknown workflow run environment: {0}".format(wfrun)
    if wfrun in ['SGEGraph', 'SGE']:
        plugin_name = wfrun
        plugin_args = {'template': template,
                       'qsub_args': modify_qsub_args(cluster['queue'],2,1,1),
                       'qstatProgramPath': cluster['qstat'],
                       'qstatCachedProgramPath': cluster['qstat_cached']}
    elif wfrun in ['local_4', 'local_12']:
        plugin_name = 'MultiProc'
        proc_count = int(wfrun.split('local_')[1])
        print("Running with {0} parallel processes on local machine".format(proc_count))
        plugin_args = {'n_procs': proc_count}
    elif wfrun == 'ds_runner':
        plugin_name = _create_DS_runner()
        plugin_args = {}
    else:
        assert wfrun in ['local',
                         'ds_runner'], "You must specify a valid run environment type.  Invalid: {0}".format(wfrun)
        plugin_name = 'Linear'
        plugin_args = {}

    return plugin_name, plugin_args


def _nipype_execution_config(stop_on_first_crash=False, stop_on_first_rerun=False, crashdumpTempDirName=None):
    stop_crash = 'false'
    stop_rerun = 'false'
    if stop_on_first_crash:
        stop_crash = 'true'
    if stop_on_first_rerun:
        # This stops at first attempt to rerun, before running, and before deleting previous results
        stop_rerun = 'true'

    if crashdumpTempDirName is None:
        import tempfile
        crashdumpTempDirName=tempfile.gettempdir()
    print( "*** Note")
    print( "    Crash file will be written to '{0}'".format(crashdumpTempDirName))
    return {
        'stop_on_first_crash': stop_crash,
        'stop_on_first_rerun': stop_rerun,
        'hash_method': 'timestamp',          # default
        'single_thread_matlab': 'true',      # default # Multi-core 2011a  multi-core for matrix multiplication.
        # default # relative paths should be on, require hash update when changed.
        'use_relative_paths': 'false',
        'remove_node_directories': 'false',  # default
        'remove_unnecessary_outputs': 'true', #remove any interface outputs not needed by the workflow
        'local_hash_check': 'true',          # default
        'job_finished_timeout': 25,
        'crashdump_dir':crashdumpTempDirName}


def _nipype_logging_config(cachedir):
    return {'workflow_level': 'INFO',  # possible options:
            'filemanip_level': 'INFO',  # INFO (default) | DEBUG
            'interface_level': 'INFO',
            'log_directory': cachedir}


def nipype_options(args, pipeline, cluster, experiment, environment):
    """
    Chicken-egg problem: cannot create pipeline dictionary with Nipype defaults until Nipype is found by environment
    # from nipype.utils.config import homedir, default_cfg

    # retval = eval(default_cfg)
    # for key, value in kwds.items():
    #     retval['execution'][key] = value
    """
    retval = copy.deepcopy(pipeline)
    from .distributed import create_global_sge_script
    template = create_global_sge_script(cluster, environment)
    #else:
    #    template = None
    plugin_name, plugin_args = _nipype_plugin_config(args['--wfrun'], cluster, template)
    retval['plugin_name'] = plugin_name
    retval['plugin_args'] = plugin_args
    retval['execution'] = _nipype_execution_config(stop_on_first_crash=True, stop_on_first_rerun=False, crashdumpTempDirName=pipeline['CRASHDUMP_DIR'])
    retval['logging'] = _nipype_logging_config(experiment['cachedir'])
    return retval


if __name__ == "__main__":
    from docopt import docopt
    SILENT = True
    args = docopt(__doc__, version='0.1')  # Get argv as dictionary
    assert args["PHASE"] in ['atlas-based-reference',
                             'subject-template-generation', 'subject-based-reference'], "Unknown phase!"
    if args['--debug']:
        # TODO: Add and run doctests!
        pass
    output = parseFile(args["FILE"], args["ENV"], args["PHASE"])
    from pprint import pprint
    print("")
    print("**** OUTPUT ****")
    for d in output:
        pprint(d)
        print("")
