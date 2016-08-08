MEX FILE HINTS
==============
The compiler configuration for building mex files
often needs fine tuning (especially on mac)

Compiler flags/ settings for mex building are set in "${MATLAB\_ROOT}/bin/mexopts.sh"
Often a google search will provide hints for how to set values for your
compiler/os/matlabversion configuration.

OS X -- with VTK there will be undefined externals because of references to an OS X framework in GDCM.  Mex doesn't allow the -framework command line option needed to link to the CoreFoundation framework.

STEPS TO FIX THIS

1. Run "mex -setup"
answer '1' to the prompt, this will ask to create or overwrite the mexopts.sh in ~/.matlab/R2013b/mexopts.sh
Hit enter to agree with overwriting (if you haven't already made local modifications).

2. ~/.matlab/R2013b/mexopts.sh is created read-only, so run "chmod u+w ~/.matlab/R2013b/mexopts.sh"

3. Edit ~/.matlab/R2013b/mexopts.sh and add "-framework CoreFoundation" to the LDFLAGS variable.

4. HINTS: utils/\*_mexopts.sh includes some files that may help identify useful flags and settings.

LDFLAGS is in the file twice, once under the glnxa64 architecture case, and once under the maci64 case.  You can change both, but the default is the "maci64" architecture.

TODO
===============
Eric Pahl Working on updates
