#! /usr/bin/env python
"""
AutoWorkup.py
===========================
Description:

Author:

Usage:

"""


def load_modules(modules):
    """ The command 'module' is actually a script call in bash:

    module=() {eval `/opt/modules/Modules/$MODULE_VERSION/bin/modulecmd bash $* }`

    So _template_runnerning os.execvp() on it doesn't work without the correct file path to the module executable

    :param modules:
    :return: None
    """
    import os

    for module in modules:
        os.system(
            " ".join(["module load", module])
        )  # os.execv(module_exe, 'bash', 'load', module])


def setup_environment(argv):
    """
    This function...

    :param argv:
    :return: environment, experiment, pipeline, cluster
    """
    print("Configuring environment...")
    import os
    import os.path
    from BAW.utilities.configFileParser import resolve_data_sink_option, parse_file
    from BAW.utilities.pathHandling import validate_path
    from BAW.utilities import misc
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    environment, experiment, pipeline, cluster = parse_file(
        argv["--ExperimentConfig"], argv["--pe"], argv["--workphase"]
    )
    pipeline["ds_overwrite"] = resolve_data_sink_option(argv, pipeline)
    if cluster is None:
        print("Running on local")
        # raise NotImplementedError("Running local has old code and has not been tested!")
        # assert argv["--wfrun"] in argvWFRUN, \
        #    "wfrun  options for clusters can only be given when the configuration file's CLUSTER option == True"
        # os.environ['NSLOTS'] = str(misc.get_cpus(argv["--wf_template_runner"]))
    else:
        load_modules(
            cluster["modules"]
        )  # Load modules if not already done  ## MODS PATH
        # print os.environ['LOADEDMODULES']
        # if environment['virtualenv_dir'] is not None:  # MODS PATH
        # activate_this = validate_path(
        #    os.path.join(environment['virtualenv_dir'], 'bin', 'activate_this.py'), False, False)
        # if os.path.exists( activate_this ) :
        #    exec(open(activate_this).read(), OrderedDict(__file__=activate_this))
    utilities_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "utilities"
    )
    configure_env = validate_path(
        os.path.join(utilities_path, "configure_env.py"), False, False
    )
    # Add the AutoWorkup directory to the PYTHONPATH every time - REQUIRED FOR CLUSTER DISPATCHING
    environment["env"]["PYTHONPATH"] = (
        environment["env"]["PYTHONPATH"] + ":" + os.path.dirname(__file__)
    )

    exec(
        open(configure_env).read(),
        OrderedDict(
            __file__=__file__,
            append_os_path=environment["env"]["PATH"],
            append_sys_path=environment["env"]["PYTHONPATH"],
        ),
    )  # MODS PATH

    print(("@" * 80))
    print((environment["env"]["PYTHONPATH"]))
    print(("@" * 80))
    print((environment["env"]["PATH"]))
    print(("@" * 80))

    from nipype import config

    config.enable_debug_mode()
    # config.enable_provenance()

    from BAW.utilities.package_check import verify_packages

    verify_packages()
    if "FREESURFER" in experiment["components"]:  # FREESURFER MODS
        configure_FS = validate_path(
            os.path.join(utilities_path, "utilities", "configure_FS.py"), False, False
        )
        exec(
            open(configure_FS).read(),
            OrderedDict(FS_VARS=misc.FS_VARS, env=environment["env"]),
        )
        print(
            "FREESURFER needs to check for sane environment here!"
        )  # TODO: raise warning, write method, what???
    for key, value in list(environment["env"].items()):
        if key in ["PATH", "PYTHONPATH"] + misc.FS_VARS:
            pass
        else:
            os.environ[key] = value  # Do not use os.putenv (see Python documentation)
    return environment, experiment, pipeline, cluster
