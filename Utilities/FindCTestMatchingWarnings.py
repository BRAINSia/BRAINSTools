#!/bin/python


import termcolor  # <-- used to make pretty printing in color
import regex  # <-- support Perl style regular expressions that closely match what is used in CTest
import sys

USAGE_COMMENT = r"""
# \author Hans J. Johnson
# \date 2013-04-04
# This script it to help identify regular epsression
# patterns # for the CTEST_CUSTOM_WARNING_EXCEPTION
# variable in a file like ${CMAKE_SOURCE_DIR}/CMake/CTestCustom.cmake.in
# for the purpose of supressing irrelavant warnings.
#
# It can be frustrating to find the correct regular
# expressions that are needed to avoid a polluted
# CDash submission
#
# The intended usage paradigm is to do the following:
# 1. In a completely clean build directory, run cmake
#    to your preferred settings
# 2. Run "make -j1 > CDASH_WARNING_LOGGER 2>&1" to make a
#    static log of all the warning messages that
#    you want to suppress
# 3. Modify the suppressions file used by this build system
#    that is usually in ${CMAKE_SOURCE_DIR}/CMake/CTestCustom.cmake.in
# 4. python FindCTestMatchingWarnings.py CDASH_WARNING_LOGGER ${CMAKE_SOURCE_DIR}/CMake/CTestCustom.cmake.in
# 5. Review output from (4), if suppressable warnings
#    identified, return to (3), else goto (6)
# 6. CELEBRATE AN UNSULLIED DASHBOARD!
#
"""

################
################

if len(sys.argv) != 3:
    print(f"{USAGE_COMMENT}\n\n")
    print(
        "USAGE {} <CDASH_WARNING_LOGGER> <$CMAKE_SOURCE_DIR/CMake/CTestCustom.cmake.in>\n".format(
            sys.argv[0]
        )
    )
    sys.exit(0)


def GetListOfCustomExludePatterns(patterns_file_name):
    patterns_file = open(patterns_file_name)
    ## NOTE: This is not very robust, but if the first line is by itself,
    ## each regex is on its own line, andthe final closing ) is on it's own line,
    ## then this works.
    start_excepts_pat = regex.compile(r"^ *set *\( *CTEST_CUSTOM_WARNING_EXCEPTION *")
    comment_pat = regex.compile(r"^ *#")
    string_pattern_pat = regex.compile(r"^ *\".*\" *$")
    end_excepts_pat = regex.compile(r"^ *\) *$")

    in_excepts_block = False

    custom_except_pattern_list = list()
    linenum = 0
    for line in patterns_file.readlines():
        linenum += 1
        iscomment_matchresult = comment_pat.search(line)
        if iscomment_matchresult:
            continue
        if in_excepts_block is False:
            start_matchresult = start_excepts_pat.search(line)
            if start_matchresult:
                in_excepts_block = True
                print(("Starting Block on line ", linenum))
        else:
            end_matchresult = end_excepts_pat.search(line)
            if end_matchresult:
                in_excepts_block = False
                print(("Ending Block on line ", linenum))
            else:
                ispat_matchresult = string_pattern_pat.search(line)
                if ispat_matchresult:
                    custom_except_pattern_list.append(
                        line.rstrip().lstrip().rstrip('"').lstrip('"')
                    )
    patterns_file.close()
    return custom_except_pattern_list


#######
##
## You must modify the following regular expressions
## so that it matches your list of CTEST_CUSTOM_WARNING_EXCEPTION
## that often are stored in a file like CMake/CTestCustom.cmake.in
##
######
print(
    termcolor.colored(f"Generating Custom exceptions from file: {sys.argv[2]}", "green")
)
CTEST_CUSTOM_WARNING_EXCEPTION = GetListOfCustomExludePatterns(sys.argv[2])
print(
    termcolor.colored(
        "count= {} {}".format(
            len(CTEST_CUSTOM_WARNING_EXCEPTION), CTEST_CUSTOM_WARNING_EXCEPTION
        ),
        "green",
    )
)

CTEST_CUSTOM_WARNING_EXCEPTION_compiled = list()
for cc in CTEST_CUSTOM_WARNING_EXCEPTION:
    CTEST_CUSTOM_WARNING_EXCEPTION_compiled.append(regex.compile(cc))


####  These are patterns that replicate the find regex from
####  115 static const char* cmCTestWarningMatches[] = {}
####  from Cmake/Source/CTest/cmCTestBuildHandler.cxx
warn_regex_patterns = [
    r"([^ :]+):([0-9]+): warning:.*",
    r"([^ :]+):([0-9]+): note:",
    r"^cc[^C]*CC: WARNING File = ([^,]+), Line = ([0-9]+)",
    r"^ld([^:])*:([ \\t])*WARNING([^:])*:",
    r"([^:]+): warning ([0-9]+):",
    r"^\"[^\"]+\", line [0-9]+: [Ww](arning|arnung)",
    r"([^:]+): warning[ \\t]*[0-9]+[ \\t]*:",
    r"^(Warning|Warnung) ([0-9]+):",
    r"^(Warning|Warnung)[ :]",
    r"WARNING: ",
    r"([^ :]+) : warning",
    r"([^:]+): warning",
    r"\", line [0-9]+\\.[0-9]+: [0-9]+-[0-9]+ \\([WI]\\)",
    r"^cxx: Warning:",
    r".*file: .* has no symbols",
    r"([^ :]+):([0-9]+): (Warning|Warnung)",
    r"\\([0-9]*\\): remark #[0-9]*",
    r"\".*\", line [0-9]+: remark\\([0-9]*\\):",
    r"cc-[0-9]* CC: REMARK File = .*, Line = [0-9]*",
    r"^CMake Warning.*:",
    r"^\\[WARNING\\]",
]

warn_regex_compiled = list()
for cc in warn_regex_patterns:
    warn_regex_compiled.append(regex.compile(cc))

print(termcolor.colored(f"Processing file: {sys.argv[1]}", "green"))
ff = open(sys.argv[1])
all_lines = ff.readlines()
ff.close
for idx in range(0, len(warn_regex_patterns)):
    print("")
    print(termcolor.colored(warn_regex_patterns[idx], "blue"))
    print("=" * 80)

    line_count = 0
    for line in all_lines:
        line_count += 1
        mresult = warn_regex_compiled[idx].search(line)
        if mresult:
            ignore_found = False
            for ignore_pat_idx in range(
                0, len(CTEST_CUSTOM_WARNING_EXCEPTION_compiled)
            ):
                ignore_result = CTEST_CUSTOM_WARNING_EXCEPTION_compiled[
                    ignore_pat_idx
                ].search(line)
                if ignore_result and ignore_found is False:
                    # print("{0}:
                    # {1}".format(CTEST_CUSTOM_WARNING_EXCEPTION[ignore_pat_idx],line))
                    ignore_found = True
                    continue
            if ignore_found is False:
                print(termcolor.colored(f"No Match :{line_count}: {line}    ", "red"))
