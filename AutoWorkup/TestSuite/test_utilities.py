"""
test_utilities.py
===================
Description:

Author:

Usage:

"""

from AutoWorkup import utilities


def configure_env_test():
    """
    This function...
    """
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    config_env = os.path.join(os.path.dirname(utilities.__file__), "configure_env.py")
    for p in range(10):
        file_template = f"/my/test/path/{p}"

    exec(
        compile(open(config_env).read(), config_env, "exec"),
        OrderedDict(
            __file__=__file__,
            append_os_path=[
                "/my/test/path/1:/my/test/path/2",
                "/my/test/path/3:/my/test/path/4",
            ],
        ),
    )
    assert os.environ["PATH"].split(os.pathsep)[0]
