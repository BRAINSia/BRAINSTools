"""
pathHandling.py
============================
Description:
    The purpose of this is to..

Usage:

"""

import os.path
import re


def validate_path(path, allow_empty, isDirectory):
    """ Check if a path exists and return the path with all variables expanded or raise AssertionError

    >>> validate_path('/usr/bin')
    /usr/bin
    >>> import os
    >>> os.environ['HOME'] == validate_path('~')
    True
    >>> validate_path('/dev/null')
    Traceback (most recent call last):
        ...
    AssertionError: Path could not be found! /dev/null

    :param path:
    :param allow_empty:
    :param isDirectory:
    :return:
    """
    msg = "Path could not be found! {0}"
    if (path is None or path == "") and allow_empty:
        return None
    full = os.path.realpath(os.path.abspath(path))
    assert os.path.exists(full), msg.format(full)
    if isDirectory:
        assert os.path.isdir(full), msg.format(full)
    return full


def validate_paths(pathString):
    """ Run validate_path() on all paths in a ':' seperated string

    >>> validate_paths('/:/usr/bin')
    /:/usr/bin
    >>> validate_paths('/:/usr/bin:/dev/null')
    Traceback (most recent call last):
        ...
    AssertionError: Path could not be found! /dev/null

    :param pathString:
    :return:
    """
    return ":".join([validate_path(path, False, True) for path in ":".split(pathString)])


def append_path_list(new, old=None):
    """ Join the new and old ":"-seperated path strings and return the result

    >>> append_path_list('/usr', '/bin:/usr/bin')
    /usr:/bin:/usr/bin
    >>> append_path_list('/usr:/usr/bin')
    /usr:/usr/bin
    >>> append_path_list('')
    Traceback (most recent call last):
        ...
    AssertionError: Path could not be found!
    >>> append_path_list('/dev/null')
    Traceback (most recent call last):
        ...
    AssertionError: Path could not be found! /dev/null

    :param new:
    :param old:
    :return:
    """
    new = validate_paths(new, False, True)
    if old is None or old == "":
        return new
    old = validate_paths(old, False, True)
    return ":".join([new, old])


def file_replace(in_file, out_file, pattern, repl):
    """
    This function...

    :param in_file:
    :param out_file:
    :param patterm:
    :param repl:
    :return:
    """
    # From http://stackoverflow.com/questions/1597649/replace-strings-in-files-by-python
    from platform import system

    if system().lower() == "linux":
        pass  # assert not os.path.samefile(in_file, out_file), "Input and output files refer to the same file!"
    assert in_file != out_file, "Input and output files are the same!"
    with open(in_file) as f:
        assert any(
            re.search(pattern, line) for line in f
        ), "Pattern not found in input file!"
    with open(in_file) as f, open(out_file, "w") as out:
        for line in f:
            out.write(re.sub(pattern, repl, line))
    return True


def clone_atlas_dir(cachedir, atlasdir):
    """
    This function...

    :param cachedir:
    :param atlasdir:
    :return:
    """
    from distutils.dir_util import copy_tree, remove_tree
    import stat

    new_dir = os.path.join(cachedir, "Atlas")
    print("Searching for atlas directory in cache...")
    if not os.path.exists(new_dir):
        old_dir = validate_path(atlasdir, False, True)
        print(("Copying new atlas {0} to cache directory...".format(old_dir)))
        newfiles = copy_tree(
            old_dir, new_dir, preserve_mode=1, preserve_times=1, verbose=True
        )
        xml_file = "ExtendedAtlasDefinition.xml"
        old_xml = os.path.join(old_dir, xml_file + ".in")
        new_xml = os.path.join(new_dir, xml_file)
        assert file_replace(old_xml, new_xml, "@ATLAS_INSTALL_DIRECTORY@", new_dir)
        newfiles.append(new_xml)
        for fname in newfiles:
            os.chmod(fname, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
            assert os.path.isfile(fname)
    else:
        print(("Atlas directory found at {0}".format(new_dir)))
    return new_dir
