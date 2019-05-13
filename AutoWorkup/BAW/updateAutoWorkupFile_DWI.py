"""
updateAutoWorkupFile_DWI.py
===============================
Description:

Author:

Usage:

"""


import argparse
import csv
import os
import sqlite3 as lite
import textwrap
from builtins import object
from builtins import str


class UpdateAutoWorkup(object):
    """This class represents a..."""

    def update_auto_workup(self):
        """
        This function...
        """
        newPath = self._generate_new_path_name()
        newFile = csv.writer(open(newPath, "wb"), quoting=csv.QUOTE_ALL)
        col_name_list = ["project", "subject", "session", "imagefiles"]
        newFile.writerow(col_name_list)
        oldFile = csv.reader(
            open(inputArguments.autoWorkupFile, "rb"), delimiter=",", quotechar='"'
        )
        NewImageDict = MakeNewImageDict()
        for row in oldFile:
            ## skip header
            if oldFile.line_num > 1:
                project = row[0]
                subject = row[1]
                session = row[2]
                newImagesList = NewImageDict.get_new_image_list(project, subject, session)
                scanDict = eval(row[3])
                if newImagesList != []:
                    if inputArguments.modality not in list(scanDict.keys()):
                        scanDict[inputArguments.modality] = newImagesList
                    else:
                        for filepath in newImagesList:
                            if filepath not in scanDict[inputArguments.modality]:
                                scanDict[inputArguments.modality].append(filepath)
                line = (project, subject, session, scanDict)
                newFile.writerow(line)

    def _generate_new_path_name(self):
        """
        This function...

        :return:
        """
        dirname = os.path.dirname(inputArguments.autoWorkupFile)
        basename = os.path.basename(inputArguments.autoWorkupFile)
        newPath = os.path.join(
            dirname, "{}_{}".format(inputArguments.modality, basename)
        )
        return newPath


class MakeNewImageDict(object):
    """This class represents a..."""

    def __init__(self):
        """This function..."""
        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        self.newImageDict = OrderedDict()
        self.commandList = list()
        self.dbName = "NewImages.db"
        self.dbTableName = "NewImages"
        self.newImagesFilepath = "{}_Images.list".format(inputArguments.modality)
        self._make_new_images_file()
        self._make_db()
        self._create_command_list()
        self._fill_db()

    def _make_new_images_file(self):
        """
        This function...
        """
        command = (
            'find %s -name "*_DWI_CONCAT_QCed.nrrd" |awk -F/ \'{print "%s," $5 "," $6 "," $7 "," $0}\' |tee %s'
            % (inputArguments.inputDir, inputArguments.modality, self.newImagesFilepath)
        )
        os.system(command)

    def _create_command_list(self):
        """
        This function...
        """
        handle = csv.reader(
            open(self.newImagesFilepath, "rb"), delimiter=",", quotechar='"'
        )
        for row in handle:
            if handle.line_num > 1:
                if len(row) == 5:
                    imageInfo = {
                        "modality": row[0],
                        "project": row[1],
                        "subject": row[2],
                        "session": row[3],
                        "filepath": row[4],
                    }
                    sqlCommand = self._make_sqlite_command(imageInfo)
                    self._append_command(sqlCommand)
                else:
                    print(
                        (
                            "WARNING: Wrong number of columns in csv file (should be 5): {0}".format(
                                row
                            )
                        )
                    )

    def _make_db(self):
        """This function..."""
        if os.path.exists(self.dbName):
            os.remove(self.dbName)
        dbColTypes = (
            "modality TEXT, project TEXT, subject TEXT, session TEXT, filepath TEXT"
        )
        con = lite.connect(self.dbName)
        dbCur = con.cursor()
        dbCur.execute(
            "CREATE TABLE {dbTableName}({dbColTypes});".format(
                dbColTypes=dbColTypes, dbTableName=self.dbTableName
            )
        )
        dbCur.close()

    def _fill_db(self):
        """This function..."""
        con = lite.connect(self.dbName)
        dbCur = con.cursor()
        for command in self.commandList:
            dbCur.execute(command)
            con.commit()
        dbCur.close()

    def _make_sqlite_command(self, imageDict):
        """
        This function..

        :param imageDict:
        """
        keys = list(imageDict.keys())
        vals = list(imageDict.values())
        col_names = ",".join(keys)
        values = ", ".join(["'" + x + "'" for x in vals])
        sqlCommand = "INSERT INTO {dbTableName} ({col_names}) VALUES ({values});".format(
            dbTableName=self.dbTableName, col_names=col_names, values=values
        )
        return sqlCommand

    def _append_command(self, val):
        """
        This function...
        :param val:
        """
        self.commandList.append(val)

    def get_new_image_list(self, project, subject, session):
        """
        This function..

        :param project:
        :param subject:
        :param session:
        :return:
        """
        sqlQuery = self._make_db_query(project, subject, session)
        dbInfo = self._get_info_from_db(sqlQuery)
        newImages = list()
        for item in dbInfo:
            newImages.append(str(item[0]))
        return newImages

    def _get_info_from_db(self, sqlQuery):
        """
        This function...

        :param sqlQuery:
        :return:
        """
        con = lite.connect(self.dbName)
        dbCur = con.cursor()
        dbCur.execute(sqlQuery)
        dbInfo = dbCur.fetchall()
        dbCur.close()
        return dbInfo

    def _make_db_query(self, project, subject, session):
        """
        This function...

        :param project:
        :param subject:
        :param session:
        :return:
        """
        return "SELECT filepath FROM {dbTableName} WHERE modality='{modality}' AND project='{project}' AND subject='{subject}' AND session='{session}';".format(
            dbTableName=self.dbTableName,
            modality=inputArguments.modality,
            project=project,
            subject=subject,
            session=session,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
This program is used to update the Auto Workup input csv file with new filepaths. \n

Example Usage: \n
$ python updateAutoWorkupFile_DWI.py -a example_autoworkup.csv -m DWI -d /paulsen/Experiments/20120814_DTIPREP
"""
        ),
    )
    parser.add_argument(
        "-a", "--autoWorkupFile", action="store", dest="autoWorkupFile", help=""
    )
    parser.add_argument(
        "-m",
        "--modality",
        action="store",
        dest="modality",
        help="T1-15, T1-30, T2-15, T2-30, PD-15, or DWI",
    )
    parser.add_argument(
        "-d",
        "--inputDir",
        action="store",
        dest="inputDir",
        help="Directory to search through for the new filepaths",
    )
    inputArguments = parser.parse_args()
    Object = UpdateAutoWorkup()
    Object.update_auto_workup()
