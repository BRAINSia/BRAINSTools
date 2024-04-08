"""
distributed.py
============================
Description:
    The purpose of this is to...

Author:

Usage:

"""

import math
from past.utils import old_div
from collections import (
    OrderedDict,
)  # Need OrderedDict internally to ensure consistent ordering


def load_cluster(modules=[]):
    """
    This Function takes in...

    :param modules:
    :return:
    """
    if len(modules) > 0:
        module_list = []
        for module in modules:
            module_list.append(f"module load {module}")
        assert len(modules) == len(module_list)
        return "\n".join(module_list)
    return ""


def source_virtualenv(virtualenv_dir=""):
    """
    This Function takes in...

    :param virtualenv_dir:
    :return:
    """
    if virtualenv_dir is None:
        return ""
    assert virtualenv_dir != ""
    return f"source {virtualenv_dir}"


def prepend_env(environment=OrderedDict()):
    """
    This Function takes in...

    :param environment:
    :return:
    """
    import os

    export_list = []
    for key, value in list(environment.items()):
        export_list.append(
            "export {key}={value}{sep}${key}".format(
                key=key, value=value, sep=os.pathsep
            )
        )  # Append to variable
    return "\n".join(export_list)


def create_global_sge_script(cluster, environment):
    """
    This is a wrapper script for running commands on an SGE cluster
    so that all the python modules and commands are pathed properly

    :param cluster:
    :param environment:
    :return:
    """
    #    >>> import os
    #    >>> nomodules = open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'TestSuite', 'node.sh.template.nomodules'), 'r')
    #    >>> create_global_sge_script({'modules':[]}, {'virtualenv_dir':'/path/to/virtualenv_dir', 'env': os.environ}).split('\n')[0]
    #    True
    #    >>> create_global_sge_script({'modules':[]}, {'virtualenv_dir':'/path/to/virtualenv_dir', 'env': os.environ}).split('\n')[0] == '#!/bin/bash FAIL'

    import os
    from string import Template
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    sub_dict = OrderedDict(
        LOAD_MODULES=load_cluster(cluster["modules"]),
        VIRTUALENV_DIR=source_virtualenv(environment["virtualenv_dir"]),
        EXPORT_ENV=prepend_env(environment["env"]),
    )
    with open(os.path.join(os.path.dirname(__file__), "node.sh.template")) as fid:
        tpl = fid.read()
    retval = Template(tpl).substitute(sub_dict)
    return retval


def modify_qsub_args(
    queue, memoryGB, minThreads, maxThreads, stdout="/dev/null", stderr="/dev/null"
):
    """
    Outputs qsub_args string for Nipype nodes
    queue is the string to specify the queue "-q all.q | -q HJ,ICTS,UI"
    memoryGB is a numeric in gigabytes to be given (ie 2.1 will result in "-l mem_free=2.1G")
    if memoryGB = 0, then it is automatically computed.
    minThreads The fewest number of threads to use (if an algorithm has benifits from more than 1 thread)
    maxThreads The max number of threads to use (if an algorithm is not multi-threaded, then just use 1)
    stdout Where to put stdout logs
    stderr Where to put stderr logs

    >>> modify_qsub_args('test', 2, 5, None)
    -S /bin/bash -cwd -pe smp 5 -l mem_free=2G -o /dev/null -e /dev/null test FAIL
    >>> modify_qsub_args('test', 2, 5, -1 )
    -S /bin/bash -cwd -pe smp 5- -l mem_free=2G -o /dev/null -e /dev/null test FAIL
    >>> modify_qsub_args('test', 8, 5, 7)
    -S /bin/bash -cwd -pe smp 5-7 -l mem_free=8G -o /dev/null -e /dev/null test FAIL
    >>> modify_qsub_args('test', 8, 5, 7, -1)
    -S /bin/bash -cwd -pe smp 5-7 -l mem_free=8G -o /dev/null -e /dev/null test FAIL
    >>> modify_qsub_args('test', 1, 5, 7, stdout='/my/path', stderr='/my/error')
    -S /bin/bash -cwd -pe smp 5-7 -l mem_free=1G -o /my/path -e /my/error test FAIL

    :param queue:
    :param memoryGB:
    :param minThreads:
    :param maxThreads:
    :param stdout: '/dev/null'
    :param stderr: '/dev/null'
    :return:

    """
    assert (
        memoryGB <= 48
    ), "Memory must be supplied in GB, so anything more than 24 seems not-useful now."

    ## NOTE: At least 1 thread needs to be requested per 2GB needed
    memoryThreads = int(
        math.ceil(old_div(math.ceil(memoryGB), 2))
    )  # Ensure that threads are integers
    minThreads = max(minThreads, memoryThreads)
    maxThreads = max(maxThreads, memoryThreads)
    maxThreads = int(maxThreads)  # Ensure that threads are integers
    minThreads = int(minThreads)  # Ensure that threads are integers

    if maxThreads is None or minThreads == maxThreads:
        threadsRangeString = f"{minThreads}"
        maxThreads = minThreads
    elif maxThreads == -1:
        threadsRangeString = f"{minThreasds}-"
        maxThreads = 12345  # HUGE NUMBER!
    else:
        threadsRangeString = f"{minThreads}-{maxThreads}"

    if maxThreads < minThreads:
        assert (
            maxThreads > minThreads
        ), f"Must specify maxThreads({minThreads}) > minThreads({maxThreads})"

    ## INFO:  May need to figure out how to set memory and threads for cluster.
    ## for now just let the number of threads requested take care of this because
    ## the job manager on helium is really slow with lots of constraints
    ##  -l mem_free={mem}

    ## format_str = '-S /bin/bash -cwd -pe smp {mint}{maxt} -o {stdout} -e {stderr} {queue}'
    format_str = "-S /bin/bash -cwd -pe smp {totalThreads} -o {stdout} -e {stderr} {queue}".format(
        totalThreads=threadsRangeString,
        stdout=stdout,
        stderr=stderr,
        queue=queue,
    )
    return format_str
