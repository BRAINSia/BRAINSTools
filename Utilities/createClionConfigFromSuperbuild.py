#
# \author Hans J. Johnson
#
# This script will assist with generating the cmake command line options used by CLion
# in the dialog box "CLion"->"Preferences"->"Build,Excecution,Deployment"->"CMake"->"[CMake options]"
#
# The effect should be a CLion environment that mirrors the configuration of the SuperBuild environment.
#
import sys
import os
import shlex

def usage():
    print("ERROR: First argument must be a cmake generated cache file")
    print("       i.e. /scratch/johnsonhj/src/BT-11/BRAINSTools-prefix/tmp/BRAINSTools-cache-Release.cmake")
    sys.exit(-1)

if len(sys.argv) > 1:
  input_cache_file_from_superbuild=sys.argv[1]
  if not os.path.exists(input_cache_file_from_superbuild):
      usage()
else:
  usage()

with open(input_cache_file_from_superbuild,'r') as fid:
    lines=fid.readlines()


## By placing all values into a dictionary, only the last processed value is stored
options_dict =dict()
for line in lines:
    ll = line.replace("set(","").replace(")","")
    elems=shlex.split(ll)
    if len(elems) > 5:
        options_dict[elems[0]] = (elems[3],elems[1])

## A few items should be under the control of CLion
skip_cmake_directives_list=["CMAKE_BUILD_TYPE","CMAKE_GENERATOR","MAKECOMMAND","EXTERNAL_PROJECT_BUILD_TYPE"];
for k,v in options_dict.items():
  if k in skip_cmake_directives_list:
    continue
  if len(v[1]) > 0 and v[1] != '""':
     print("-D{0}:{1}=\"{2}\"".format(k,v[0],v[1]))
