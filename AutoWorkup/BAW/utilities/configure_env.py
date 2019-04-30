"""
configure_env.py
============================
Description:
    The purpose of this is to...
Author:

Usage:


By using
         execfile(this_file, OrderedDict(__file__=__file__,
                                  append_os_path=['list', 'of', 'PATH'],
                                  append_sys_path=['list', 'of', 'PYTHONPATH']))
    you will activate the Autoworkup environment.
    __file__ is the current Python file that exists in the BRAINSTools/AutoWorkup/ directory.
"""


from builtins import str

try:
    __file__
    append_os_path
    append_sys_path
except NameError:
    raise AssertionError(
        "You must run this like execfile('path/to/utilities/configure_env.py',\n\
        \tOrderedDict(__file__=__file__,\n\
        \tappend_os_path=['list', 'of', 'PATH'],\n\
        \tappend_sys_path=['list', 'of', 'PYTHONPATH']))")
import sys
import os

# PATH
append_os_path = [os.path.abspath(os.path.expanduser(os.path.expandvars(item))) for item in
                  append_os_path.split(os.pathsep)]
old_os_path = list(os.environ['PATH'].split(os.pathsep))
package_path = os.path.dirname(os.path.abspath(__file__))
append_os_path.insert(0, package_path)  # AutoWorkup/, if called from AutoWorkup.py
append_os_path.insert(1, (os.path.join(package_path, 'bin')))  # AutoWorkup/bin, if ...
# Move the added items to the front of the path:
new_os_path = []
for item in append_os_path:
    if item not in old_os_path:
        new_os_path.append(item)
os.environ['PATH'] = os.pathsep.join(new_os_path + old_os_path)

# PYTHONPATH
append_sys_path = [os.path.abspath(os.path.expanduser(os.path.expandvars(item))) for item in
                   append_sys_path.split(os.pathsep)]
# Move the added items to the front of the path:
old_sys_path = list(sys.path)
# old_sys_path.remove('')
for item in list(old_sys_path):
    if item in append_sys_path:
        old_sys_path.remove(item)
sys.path = [''] + append_sys_path + old_sys_path

print(("=" * 100))
print(("=" * 100))
print(("=" * 100))
print(("=" * 100))
print(("=" * 100))
print(("NEW PATH env " + os.environ['PATH']))
print(("=" * 100))
print(("NEW PYTHONPATH env " + str(sys.path)))
print(("=" * 100))
