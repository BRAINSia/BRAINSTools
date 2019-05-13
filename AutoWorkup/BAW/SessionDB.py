"""
SessionDB.py
=================
Description:

Author:

Usage:

"""


import csv
import os
import sqlite3 as lite
import sys
from builtins import object
from builtins import range
from builtins import str
from collections import OrderedDict


class SessionDB(object):
    """This class represents a..."""

    def __init__(self, defaultDBName="TempFileForDB.db", subject_list=[]):
        """This function represents a"""
        self.MasterTableName = "MasterDB"
        self.dbName = defaultDBName
        self.subjectList = subject_list
        self.cursor = None
        self.connection = None
        subject_filter = "( "
        for curr_subject in subject_list:
            subject_filter += "'" + curr_subject + "',"
        subject_filter = subject_filter.rstrip(",")  # Remove last ,
        subject_filter += " )"
        if subject_list[0] == "all":
            self.MasterQueryFilter = "SELECT * FROM {_tablename}".format(
                _tablename=self.MasterTableName
            )
        else:
            self.MasterQueryFilter = "SELECT * FROM {_tablename} WHERE subj IN {_subjid}".format(
                _tablename=self.MasterTableName, _subjid=subject_filter
            )

    def open_connection(self):
        """This function represents a """
        self.connection = lite.connect(self.dbName)
        self.cursor = self.connection.cursor()

    def close_connection(self):
        """This function represents a """
        if not self.cursor is None:
            self.cursor.close()
        if not self.connection is None:
            self.connection.close()

    def _local_fill_db_and_close(self, sqlCommandList):
        """This function represents a

        :param sqlCommandList:
        """
        print("Filling SQLite database SessionDB.py")
        for sqlCommand in sqlCommandList:
            self.cursor.execute(sqlCommand)
        self.connection.commit()
        print("Finished filling SQLite database SessionDB.py")

    def make_new_db(self, subject_data_file, mountPrefix):
        """
        This function...

        :param subject_data_file:
        :param mountPrefix:
        :return:
        """
        ## First close so that we can delete.
        self.close_connection()
        if os.path.exists(self.dbName):
            os.remove(self.dbName)
        self.open_connection()
        dbColTypes = (
            "project TEXT, subj TEXT, session TEXT, type TEXT, Qpos INT, filename TEXT"
        )
        self.cursor.execute(
            "CREATE TABLE {tablename}({coltypes});".format(
                tablename=self.MasterTableName, coltypes=dbColTypes
            )
        )
        self.connection.commit()
        sqlCommandList = list()
        missingFilesLog = self.dbName + "_MissingFiles.log"
        missingCount = 0
        print(("MISSING FILES RECORED IN {0}".format(missingFilesLog)))
        missingFiles = open(missingFilesLog, "w")
        print(("Building Subject returnList: " + subject_data_file))
        subjData = csv.reader(
            open(subject_data_file, "rt"), delimiter=",", quotechar='"'
        )
        allEntriesOK = True
        for row in subjData:
            if len(row) < 1:
                # contine of it is an empty row
                continue
            if row[0][0] == "#":
                # if the first character is a #, then it is commented out
                continue
            if row[0] == "project":
                # continue if header line
                continue
            currDict = OrderedDict()
            validEntry = True
            if len(row) == 4:
                currDict = {"project": row[0], "subj": row[1], "session": row[2]}
                rawDict = OrderedDict(eval(row[3]))
                dictionary_keys = list(rawDict.keys())
                if not (("T1-15" in dictionary_keys) or ("T1-30" in dictionary_keys)):
                    print(
                        (
                            "ERROR: Skipping session {0} due to missing T1's: {1}".format(
                                currDict, dictionary_keys
                            )
                        )
                    )
                    print("REMOVE OR FIX BEFORE CONTINUING")
                    allEntriesOK = False
                for imageType in dictionary_keys:
                    currDict["type"] = imageType
                    fullPaths = [mountPrefix + i for i in rawDict[imageType]]
                    if len(fullPaths) < 1:
                        print(("Invalid Entry!  {0}".format(currDict)))
                        validEntry = False
                    for i in range(len(fullPaths)):
                        imagePath = fullPaths[i]
                        if not os.path.exists(imagePath):
                            print(
                                (
                                    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Missing File: {0}\n".format(
                                        imagePath
                                    )
                                )
                            )
                            missingFiles.write("Missing File: {0}\n".format(imagePath))
                            validEntry = False
                            missingCount += 1
                        else:
                            print(("Found file {0}".format(imagePath)))
                        if validEntry == True:
                            currDict["Qpos"] = str(i)
                            currDict["filename"] = imagePath
                            sqlCommand = self.make_sqlite_command(currDict)
                            sqlCommandList.append(sqlCommand)
            else:
                print("ERROR:  Invalid number of elements in row")
                print(row)
        sqlCommandList
        self._local_fill_db_and_close(sqlCommandList)
        if (missingCount > 0) or (allEntriesOK == False):
            if os.path.exists(self.dbName):
                os.remove(self.dbName)
            missingFiles.close()
            print(("ABORTING: At least 1 missing file\n" * 20))
            sys.exit(-1)
        else:
            missingFiles.write("NO_MISSING_FILES")
        missingFiles.close()
        self.close_connection()

    def get_subject_filter(self):
        """
        This function...

        :return:
        """
        return self.MasterQueryFilter

    def make_sqlite_command(self, imageDict):
        """
        This function...

        :param imageDict
        :return:
        """
        keys = list(imageDict.keys())
        vals = list(imageDict.values())
        col_names = ",".join(keys)
        values = ", ".join(["'" + x + "'" for x in vals])
        sqlCommand = "INSERT INTO {_tablename} ({_col_names}) VALUES ({_values});".format(
            _tablename=self.MasterTableName, _col_names=col_names, _values=values
        )
        return sqlCommand

    def get_info_from_db(self, sqlCommand):
        """
        This function...

        :param sqlCommand:
        :return:
        """
        # print("get_info_from_db({0})".format(sqlCommand))
        self.open_connection()
        self.cursor.execute(sqlCommand)
        dbInfo = self.cursor.fetchall()
        self.close_connection()
        return dbInfo

    def get_first_scan(self, sessionid, scantype):
        """
        This function...

        :param sessionid:
        :param scantype:
        :return:
        """
        sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' AND type='{_scantype}' AND Qpos='0';".format(
            _master_query=self.MasterQueryFilter,
            _sessionid=sessionid,
            _scantype=scantype,
        )
        val = self.get_info_from_db(sqlCommand)
        filename = str(val[0][0])
        return filename

    def get_first_t1(self, sessionid):
        """
        This function...

        :param sessionid:
        :return:
        """
        scantype = "T1-30"
        sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' AND type='{_scantype}' AND Qpos='0';".format(
            _master_query=self.MasterQueryFilter,
            _sessionid=sessionid,
            _scantype=scantype,
        )
        val = self.get_info_from_db(sqlCommand)
        # print "HACK: ",sqlCommand
        # print "HACK: ", val
        filename = str(val[0][0])
        return filename

    def get_filenames_by_scan_type(self, sessionid, scantypelist):
        """
        This function...

        :param sessionid:
        :param scantypelist:
        :return:
        """
        returnList = list()
        for currScanType in scantypelist:
            sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' AND type='{_scantype}' ORDER BY Qpos ASC;".format(
                _master_query=self.MasterQueryFilter,
                _sessionid=sessionid,
                _scantype=currScanType,
            )
            val = self.get_info_from_db(sqlCommand)
            for i in val:
                returnList.append(str(i[0]))
        return returnList

    def find_scan_type_length(self, sessionid, scantypelist):
        """
        This function...

        :param sessionid:
        :param scantypelist:
        :return:
        """
        countList = self.get_filenames_by_scan_type(sessionid, scantypelist)
        return len(countlist)

    def get_t1s_t2s(self, sessionid):
        """
        This function...
        :param sessionid:
        :return: returnList
        """
        sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' ORDER BY type ASC, Qpos ASC;".format(
            _master_query=self.MasterQueryFilter, _sessionid=sessionid
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def get_all_projects(self):
        """
        This function..

        :return:
        """
        sqlCommand = "SELECT DISTINCT project FROM ({_master_query});".format(
            _master_query=self.MasterQueryFilter
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def get_all_subjects(self):
        """
        This function..

        :return:
        """
        sqlCommand = "SELECT DISTINCT subj FROM ({_master_query});".format(
            _master_query=self.MasterQueryFilter
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def get_all_sessions(self):
        """
        This function..

        :return:
        """
        # print("HACK:  This is a temporary until complete re-write")
        sqlCommand = "SELECT DISTINCT session FROM ({_master_query});".format(
            _master_query=self.MasterQueryFilter
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def get_sessions_from_subject(self, subj):
        """
        This function..

        :param subj:
        :return:
        """
        sqlCommand = "SELECT DISTINCT session FROM ({_master_query}) WHERE subj='{_subjid}';".format(
            _master_query=self.MasterQueryFilter, _subjid=subj
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def get_everything(self):
        """
        This function..

        :return:
        """
        sqlCommand = "SELECT * FROM ({_master_query});".format(
            _master_query=self.MasterQueryFilter
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(i)
        return returnList

    def get_subjects_from_project(self, project):
        """
        This function..

        :param project:
        :return:
        """
        sqlCommand = "SELECT DISTINCT subj FROM ({_master_query}) WHERE project='{_projectid}';".format(
            _master_query=self.MasterQueryFilter, _projectid=project
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def get_subj_from_session(self, session):
        """
        This function..

        :param session:
        :return:
        """
        sqlCommand = "SELECT DISTINCT subj FROM ({_master_query}) WHERE session='{_sessionid}';".format(
            _master_query=self.MasterQueryFilter, _sessionid=session
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        if len(returnList) != 1:
            print("ERROR: More than one subject found")
            sys.exit(-1)
        return returnList[0]

    def get_proj_from_session(self, session):
        """
        This function..

        :param session:
        :return:
        """
        sqlCommand = "SELECT DISTINCT project FROM ({_master_query}) WHERE session='{_sessionid}';".format(
            _master_query=self.MasterQueryFilter, _sessionid=session
        )
        val = self.get_info_from_db(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        if len(returnList) != 1:
            print("ERROR: More than one project found")
            sys.exit(-1)
        return returnList[0]


#
# import SessionDB
# a=SessionDB.SessionDB()
# a=SessionDB.SessionDB('predict_autoworkup.csv',''))
# a.get_first_scan('42245','T1-30')
# a.get_info_from_db("SELECT filename FROM SessionDB WHERE session=42245 ORDER BY type ASC, Qpos ASC;")
# a.get_info_from_db("SELECT DISTINCT subj FROM SessionDB;")
