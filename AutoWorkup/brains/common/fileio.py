"""
fileio.py
=================
Description:

Author:

Usage:

"""


def check_file(path):
    """
    This function..

    :param path:
    :return:
    """
    import os.path

    full = os.path.abspath(path)
    if os.path.exists(full):
        return full
    return None


def parse_labels_file():
    """
    This function...

    :return:
    """
    import os.path
    from ..config import _config
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    build_tree = _config.get("Resources", "build_directory")
    filename = _config.get("Resources", "label_template")
    fullname = check_file(os.path.join(build_tree, filename))
    labelDict = OrderedDict()
    with open(fullname, "r") as fid:
        for line in fid.readlines():
            if line[0] != "#":
                parts = line.split(" ")
                labelDict[int(parts[0])] = parts[1]
    return labelDict
