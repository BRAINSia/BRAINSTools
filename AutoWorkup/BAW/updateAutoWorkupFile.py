"""
updateAutoWorkupFile.py
============================
Description:

Author:

Usage:

"""

import argparse
import csv
import os
import textwrap
from builtins import object


class UpdateAutoWorkup(object):
    """This class represents a..."""

    def _get_black_list(self):
        """
        This function...
        :return: blackListDict, list(blackListDict.keys())
        """
        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        handle = csv.reader(
            open(inputArguments.blackList, "rb"), delimiter=",", quotechar='"'
        )
        blackListDict = OrderedDict()
        for row in handle:
            if len(row) == 3:
                blackListDict[row[0]] = row[1]
            else:
                print(
                    (
                        "WARNING: WRONG # of columns in csv (should be 3): {0}".format(
                            row
                        )
                    )
                )
        return blackListDict, list(blackListDict.keys())

    def _generate_new_path_name(self):
        """
        This function...

        :return:
        """
        dirname = os.path.dirname(inputArguments.autoWorkupFile)
        basename = os.path.basename(inputArguments.autoWorkupFile)
        newPath = os.path.join(dirname, "edited_{0}".format(basename))
        return newPath

    def update_auto_workup(self):
        """
        This function...
        """
        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        newPath = self._generate_new_path_name()
        newFile = csv.writer(open(newPath, "wb"), quoting=csv.QUOTE_ALL)
        col_name_list = ["project", "subject", "session", "imagefiles"]
        newFile.writerow(col_name_list)
        oldFile = csv.reader(
            open(inputArguments.autoWorkupFile, "rb"), delimiter=",", quotechar='"'
        )
        blackListDict, blackListKeys = self._get_black_list()
        print((blackListDict, blackListKeys))
        for row in oldFile:
            ## skip header
            if oldFile.line_num > 1:
                scanDict = eval(row[3])
                newScanDict = OrderedDict()
                for scan in list(scanDict.keys()):
                    filepaths = scanDict[scan]
                    for path in filepaths:
                        if path in blackListKeys:
                            newPath = blackListDict[path]
                            if newPath == "":
                                continue
                            if scan not in list(newScanDict.keys()):
                                newScanDict[scan] = [newPath]
                            else:
                                newScanDict[scan].append(newPath)
                        else:
                            if scan not in list(newScanDict.keys()):
                                newScanDict[scan] = [path]
                            else:
                                newScanDict[scan].append(path)
                if newScanDict != OrderedDict():
                    project = row[0]
                    subject = row[1]
                    session = row[2]
                    line = (project, subject, session, newScanDict)
                    newFile.writerow(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
This program is used to update the Auto Workup input csv file with the black list. \n
The blacklist needs to be a comma separated file.
    Example line:
    file/path/you/desire/to/change, new/file/path \n
If a file path is to deleted from the auto workup file, leave the
new/file/path blank in the blacklist file, but keep the comma:
    Example line:
    filepath/you/desire/to/change,

Common Usage: \n
$ python updateAutoWorkupFile.py -a autoworkupfile -b blacklistfile

Example:
$ python updateAutoWorkupFile.py -a example_autoworkup.csv -b example_blacklist.csv"""
        ),
    )
    parser.add_argument(
        "-a", "--autoWorkupFile", action="store", dest="autoWorkupFile", help=""
    )
    parser.add_argument("-b", "--blackList", action="store", dest="blackList", help="")
    inputArguments = parser.parse_args()
    Object = UpdateAutoWorkup()
    Object.update_auto_workup()
