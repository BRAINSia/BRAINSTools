#! /usr/bin/env python

"""
configFileParser.py
===================

Usage:
  configFileParser.py [--debug] ENV FILE
  configFileParser.py (-h | --help | -v | --version)

Arguments:
  ENV            Environment name in configuration file
  FILE           Configuration file

Options:
  -h, --help     Print this and exit
  -v, --version  Print file version and exit
  --debug        Run doctests for file
"""
from ConfigParser import ConfigParser
import os
import sys

from pathHandling import *
import misc


def parseEnvironment(parser, section):
    """ Parse the environment section given by 'section' and return a dictionary
        Values are shell-centric, i.e. PYTHONPATH is a colon-seperated string

    """
    retval = dict()
    if parser.has_option(section, 'ENVAR_DICT'):
        retval['env'] = eval(parser.get(section, 'ENVAR_DICT'))
    else:
        retval['env'] = dict()
    if 'PYTHONPATH' in retval['env'].keys():
        pythonpath = appendPathList(parser.get(section, 'APPEND_PYTHONPATH'), retval['env']['PYTHONPATH'])
        retval['env']['PYTHONPATH'] = pythonpath  # Create append to PYTHONPATH
    else:
        retval['env']['PYTHONPATH'] = parser.get(section, 'APPEND_PYTHONPATH')
    if 'PATH' in retval['env'].keys():
        envpath = appendPathList(parser.get(section, 'APPEND_PATH'), retval['env']['PATH'])
        retval['env']['PATH'] = envpath  # Create append to PATH
    else:
        retval['env']['PATH'] = parser.get(section, 'APPEND_PATH')
    retval['cluster'] = parser.getboolean(section, 'CLUSTER')
    retval['prefix'] = validatePath(parser.get(section, 'MOUNT_PREFIX'), True, True)
    if retval['prefix'] is None:
        retval['prefix'] = ''
    retval['virtualenv_dir'] = validatePath(parser.get(section, 'VIRTUALENV_DIR'), False, True)
    return retval


def create_experiment_dir(dirname, name, suffix, verify=False):
    """ Construct directories given the base directory, the experiment name, and the suffix ['CACHE', 'Results'] """
    basename = name + '_' + suffix
    fullpath = os.path.join(dirname, basename)
    if verify:
        return validatePath(fullpath, False, True)
    else:
        if os.path.isdir(fullpath):
            print "WARNING: Experiment directory already exists.  Continuing will overwrite the previous results..."
            print "   Path: {0}".format(fullpath)
            return fullpath
        try:
            os.makedirs(fullpath)
        except OSError:
            raise
    return fullpath


def parseExperiment(parser):
    """ Parse the experiment section and return a dictionary """
    retval = dict()
    dirname = validatePath(parser.get('EXPERIMENT', 'BASE_OUTPUT_DIR'), False, True)
    current = parser.get('EXPERIMENT', 'EXPERIMENT')
    retval['cachedir'] = create_experiment_dir(dirname, current, 'CACHE')
    retval['resultdir'] = create_experiment_dir(dirname, current, 'Results')
    if parser.has_option('EXPERIMENT', 'PREVIOUS_EXPERIMENT'):
        # If this is the initial run, there will be no previous experiment
        previous = parser.get('EXPERIMENT', 'PREVIOUS_EXPERIMENT')
        retval['previouscache'] = create_experiment_dir(dirname, previous, 'CACHE', verify=True)
        retval['previousresult'] = create_experiment_dir(dirname, previous, 'Results', verify=True)
    atlas = validatePath(parser.get('EXPERIMENT', 'ATLAS_PATH'), False, True)
    retval['atlascache'] = clone_atlas_dir(retval['cachedir'], atlas)
    retval['dbfile'] = validatePath(parser.get('EXPERIMENT', 'SESSION_DB'), False, False)
    retval['components'] = [x.lower() for x in eval(parser.get('EXPERIMENT', 'WORKFLOW_COMPONENTS'))]
    return retval


def parsePipeline(parser):
    """ Parse the pipeline section and return a dictionary """
    retval = dict()
    retval['ds_overwrite'] = parser.getboolean('NIPYPE', 'GLOBAL_DATA_SINK_REWRITE')
    return retval


def parseCluster(parser):
    """ Parse the cluster section and return a dictionary """
    retval = dict()
    retval['modules'] = eval(parser.get('CLUSTER', 'MODULES'))
    retval['queue'] = parser.get('CLUSTER', 'QUEUE')
    retval['long_q'] = parser.get('CLUSTER', 'QUEUE_LONG')
    retval['qstat'] = parser.get('CLUSTER', 'QSTAT_IMMEDIATE')
    retval['qstat_cached'] = parser.get('CLUSTER', 'QSTAT_CACHED')
    return retval


def parseFile(configFile, env):
    configFile = os.path.realpath(configFile)
    assert os.path.exists(configFile), "Configuration file could not be found: {0}".format(configFile)
    parser = ConfigParser(allow_no_value=True)  # Parse configuration file parser = ConfigParser()
    parser.read(configFile)
    assert (parser.has_option(env, '_BUILD_DIR') or parser.has_option('DEFAULT', '_BUILD_DIR')), "BUILD_DIR option not in {0}".format(env)
    environment = parseEnvironment(parser, env)
    experiment = parseExperiment(parser)
    pipeline = parsePipeline(parser)
    if environment['cluster']:
        cluster = parseCluster(parser)
        return environment, experiment, pipeline, cluster
    return environment, experiment, pipeline, None


def resolveDataSinkOption(args, pipeline):
    if args["--rewrite-datasinks"] or pipeline['ds_overwrite']:  # GLOBAL_DATA_SINK_REWRITE
        return True
    return False


def nipype_options(args, pipeline, cluster, template, experiment):
    retval = {}
    ## Chicken-egg problem: cannot create pipeline dictionary with Nipype defaults until Nipype is found by environment
    # from nipype.utils.config import homedir, default_cfg

    # retval = eval(default_cfg)
    # for key, value in kwds.items():
    #     retval['execution'][key] = value
    print pipeline
    retval['ds_overwrite'] = pipeline['ds_overwrite']  # resolveDataSinkOption(args, pipeline)
    retval['execution'] = misc.nipype_execution(plugin=args['--wfrun'])
    retval['plugin_args'] = misc.nipype_plugin_args(args['--wfrun'], cluster, template)
    retval['logging'] = misc.nipype_logging(experiment['cachedir'])
    return retval


if __name__ == "__main__":
    from docopt import docopt

    args = docopt(__doc__, version='0.1')  # Get argv as dictionary
    # TODO: Add and run doctests!
    print parseFile(args["FILE"], args["ENV"])
