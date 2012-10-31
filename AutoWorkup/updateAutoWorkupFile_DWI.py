import argparse
import csv
import os
import textwrap

class UpdateAutoWorkup():

    def _generateNewPathName(self):
        dirname = os.path.dirname(inputArguments.autoWorkupFile)
        basename = os.path.basename(inputArguments.autoWorkupFile)
        newPath = os.path.join(dirname, "edited_{0}".format(basename))
        return newPath

    def updateAutoWorkup(self):
        newPath = self._generateNewPathName()
        newFile = csv.writer(open(newPath, 'wb'), quoting=csv.QUOTE_ALL)
        col_name_list = ["project", "subject", "session", "imagefiles"]
        newFile.writerow(col_name_list)
        oldFile = csv.reader(open(inputArguments.autoWorkupFile, 'rb'), delimiter=',', quotechar='\"')
        NewImageDict = MakeNewImageDict()
        for row in oldFile:
            ## skip header
            if oldFile.line_num > 1:
                project = row[0]
                subject = row[1]
                session = row[2]
                scanDict = eval(row[3])
                newScanDict = dict()
                if project =
                for scan in scanDict.keys():
                    filepaths = scanDict[scan]
                    for path in filepaths:
                        if path in blackListKeys:
                            newPath = blackListDict[path]
                            if newPath == '':
                                continue
                            if scan not in newScanDict.keys():
                                newScanDict[scan] = [newPath]
                            else:
                                newScanDict[scan].append(newPath)
                        else:
                            if scan not in newScanDict.keys():
                                newScanDict[scan] = [path]
                            else:
                                newScanDict[scan].append(path)
                if newScanDict != {}:
                    project = row[0]
                    subject = row[1]
                    session = row[2]
                    line = (project, subject, session, newScanDict)
                    newFile.writerow(line)

class MakeNewImageDict():

    def __init__(self):
        self.newImageDict = dict()
        self._populateNewImageDict()

    def _populateNewImageDict(self):
        handle = csv.reader(open(inputArguments.blackList, 'rb'), delimiter=',', quotechar='\"')
        for row in handle:
            if len(row) == 5:
                #imageInfo = {'modality': row[0],'project': row[1], 'subject': row[2],
                #             'session': row[3], 'imagefile': row[4]}
                modality = row[0]
                project = row[1]
                subject = row[2]
                session = row[3]
                imagefile = row[4]
                self._addModality(modality)
                self._addProject(modality, project)
                self._addSubject(modality, project, subject)
                self._addSession(modality, project, subject, session)
                self._addImageFile(modality, project, subject, session, imagefile)
            else:
                print("WARNING: WRONG # of columns in csv (should be 5): {0}".format(row))

    def _addModality(self, modality):
        if modality not in self.newImageDict.keys():
            self.newImageDict[modality] = dict()

    def _addProject(self, modality, project):
        if project not in self.newImageDict[modality].keys():
            self.newImageDict[modality][project] = dict()

    def _addSubject(self, modality, project, subject):
        if subject not in self.newImageDict[modality][project].keys():
            self.newImageDict[modality][project][subject] = dict()

    def _addSession(self, modality, project, subject, session):
        if session not in self.newImageDict[modality][project][subject].keys():
            self.newImageDict[modality][project][subject][session] = list()

    def _addImageFile(self, modality, project, subject, session, imagefile):
        if imagefile not in self.newImageDict[modality][project][subject][session]:
            self.newImageDict[modality][project][subject][session].append(imagefile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""
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
$ python updateAutoWorkupFile.py -a example_autoworkup.csv -b example_blacklist.csv"""))
    parser.add_argument('-a', '--autoWorkupFile', action='store', dest='autoWorkupFile', help='')
    parser.add_argument('-b', '--blackList', action='store', dest='blackList', help='')
    parser.add_argument('-m', '--modality', action='store', dest='modality', help='T1-30, T2-30, or DWI')
    inputArguments = parser.parse_args()
    Object = UpdateAutoWorkup()
    Object.updateAutoWorkup()
