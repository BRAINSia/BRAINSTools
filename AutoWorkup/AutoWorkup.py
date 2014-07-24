#! /usr/bin/env python
from baw_exp import OpenSubjectDatabase

def load_modules(modules):
    """ The command 'module' is actually a script call in bash:

    module=() {eval `/opt/modules/Modules/$MODULE_VERSION/bin/modulecmd bash $* }`

    So _template_runnerning os.execvp() on it doesn't work without the correct file path to the module executable """
    import os
    for module in modules:
        os.system(" ".join(['module load', module]))  # os.execv(module_exe, 'bash', 'load', module])


def setup(argv):
    print "Configuring environment..."
    import os, os.path
    from utilities.configFileParser import resolveDataSinkOption, parseFile
    from utilities.pathHandling import validatePath
    from utilities import misc
    environment, experiment, pipeline, cluster = parseFile(argv["--ExperimentConfig"], argv["--pe"])
    pipeline['ds_overwrite'] = resolveDataSinkOption(argv, pipeline)

    if cluster is None:
        assert argv["--wf_template_runner"] in misc.WFRUN, \
          "wf_template_runner options for clusters can only be given when the configuration file's CLUSTER option == True"
        os.environ['NSLOTS'] = str(misc.get_cpus(argv["--wf_template_runner"]))
    else:
        load_modules(cluster['modules'])  # Load modules if not already done  ## MODS PATH
        # print os.environ['LOADEDMODULES']
    if environment['virtualenv_dir']:  ## MODS PATH
        activate_this = validatePath(os.path.join(environment['virtualenv_dir'], 'bin', 'activate_this.py'), False, False)
        execfile(activate_this, dict(__file__=activate_this))
    utilities_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utilities')
    configure_env = validatePath(os.path.join(utilities_path, 'configure_env.py'), False, False)
    # Add the AutoWorkup directory to the PYTHONPATH every time - REQUIRED FOR CLUSTER DISPATCHING
    environment['env']['PYTHONPATH'] =  environment['env']['PYTHONPATH'] + ":" + os.path.dirname(__file__)
    execfile(configure_env, dict(__file__=__file__,
                                 append_os_path=environment['env']['PATH'],
                                 append_sys_path=environment['env']['PYTHONPATH'])
        )  ## MODS PATH
    from nipype import config
    config.enable_debug_mode()
    from utilities.package_check import verify_packages
    verify_packages()
    if 'FREESURFER' in experiment['components']:  ## FREESURFER MODS
        configure_FS = validatePath(os.path.join(utilities_path, 'utilities', 'configure_FS.py'), False, False)
        execfile(configure_FS, dict(FS_VARS=misc.FS_VARS, env=environment['env']))
        print "FREESURFER needs to check for sane environment here!"  # TODO: raise warning, write method, what???
    for key, value in environment['env'].items():
        if key in ['PATH', 'PYTHONPATH'] + misc.FS_VARS:
            pass
        else:
            os.environ[key] = value  # Do not use os.putenv (see Python documentation)
    return environment, experiment, pipeline, cluster


def get_subjects(argv, cache, prefix, dbfile, shuffle=True):
    import random
    _temp = OpenSubjectDatabase(cache, ['all'], prefix, dbfile)
    subjects = argv["--subjects"].split(',')
    if "all" in subjects:
        subjects = _temp.getAllSubjects()
    if shuffle:
        random.shuffle(subjects)  # randomly shuffle to get max
    return subjects
