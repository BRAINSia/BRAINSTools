from __future__ import print_function

import os.path
import re


def validatePath(path, allow_empty, isDirectory):
    """ Check if a path exists and return the path with all variables expanded or raise AssertionError

    >>> validatePath('/usr/bin')
    /usr/bin
    >>> import os
    >>> os.environ['HOME'] == validatePath('~')
    True
    >>> validatePath('/dev/null')
    Traceback (most recent call last):
        ...
    AssertionError: Path could not be found! /dev/null
    """
    msg = "Path could not be found! {0}"
    if (path is None or path == '') and allow_empty:
        return None
    full = os.path.realpath(os.path.abspath(path))
    assert os.path.exists(full), msg.format(full)
    if isDirectory:
        assert os.path.isdir(full), msg.format(full)
    return full


def validatePaths(pathString):
    """ Run validatePath() on all paths in a ':' seperated string

    >>> validatePaths('/:/usr/bin')
    /:/usr/bin
    >>> validatePaths('/:/usr/bin:/dev/null')
    Traceback (most recent call last):
        ...
    AssertionError: Path could not be found! /dev/null
    """
    return ':'.join([validatePath(path, False, True) for path in ':'.split(pathString)])


def appendPathList(new, old=None):
    """ Join the new and old ":"-seperated path strings and return the result

    >>> appendPathList('/usr', '/bin:/usr/bin')
    /usr:/bin:/usr/bin
    >>> appendPathList('/usr:/usr/bin')
    /usr:/usr/bin
    >>> appendPathList('')
    Traceback (most recent call last):
        ...
    AssertionError: Path could not be found!
    >>> appendPathList('/dev/null')
    Traceback (most recent call last):
        ...
    AssertionError: Path could not be found! /dev/null
    """
    new = validatePaths(new, False, True)
    if old is None or old == '':
        return new
    old = validatePaths(old, False, True)
    return ':'.join([new, old])


def file_replace(in_file, out_file, pattern, repl):
    # From http://stackoverflow.com/questions/1597649/replace-strings-in-files-by-python
    from platform import system
    if system().lower() == 'linux':
        pass  # assert not os.path.samefile(in_file, out_file), "Input and output files refer to the same file!"
    assert in_file != out_file, "Input and output files are the same!"
    with open(in_file) as f:
        assert any(re.search(pattern, line) for line in f), "Pattern not found in input file!"
    with open(in_file) as f, open(out_file, "w") as out:
        for line in f:
            out.write(re.sub(pattern, repl, line))
    return True


def clone_atlas_dir(cachedir, atlasdir):
    from distutils.dir_util import copy_tree, remove_tree
    import stat

    new_dir = os.path.join(cachedir, 'Atlas')
    print("Searching for atlas directory in cache...")
    if not os.path.exists(new_dir):
        old_dir = validatePath(atlasdir, False, True)
        print("Copying new atlas {0} to cache directory...".format(old_dir))
        newfiles = copy_tree(old_dir, new_dir, preserve_mode=1, preserve_times=1, verbose=True)
        xml_file = 'ExtendedAtlasDefinition.xml'
        old_xml = os.path.join(old_dir, xml_file + '.in')
        new_xml = os.path.join(new_dir, xml_file)
        assert file_replace(old_xml, new_xml, "@ATLAS_INSTALL_DIRECTORY@", new_dir)
        newfiles.append(new_xml)
        for fname in newfiles:
            os.chmod(fname, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
            assert os.path.isfile(fname)
    else:
        print("Atlas directory found at {0}".format(new_dir))
    return new_dir
