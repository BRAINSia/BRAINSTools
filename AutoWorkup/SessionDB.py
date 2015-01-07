import os
import sys
import sqlite3 as lite
import csv

DEBUG = False

class SessionDB():
    def __init__(self, defaultDBName='TempFileForDB.db', subject_list=[]):
        self.MasterTableName = "MasterDB"
        self.dbName = defaultDBName
        self.subjectList = subject_list
        self.cursor = None
        self.connection = None
        subject_filter = "( "
        for curr_subject in subject_list:
            subject_filter += "'" + curr_subject + "',"
        subject_filter = subject_filter.rstrip(',')  # Remove last ,
        subject_filter += " )"
        if  subject_list[0] == 'all':
            self.MasterQueryFilter = "SELECT * FROM {_tablename}".format(
                _tablename=self.MasterTableName)
        else:
            self.MasterQueryFilter = "SELECT * FROM {_tablename} WHERE subj IN {_subjid}".format(
                _tablename=self.MasterTableName,
                _subjid=subject_filter)

    def open_connection(self):
        self.connection = lite.connect(self.dbName)
        self.cursor = self.connection.cursor()

    def close_connection(self):
        if not self.cursor is None:
            self.cursor.close()
        if not self.connection is None:
            self.connection.close()

    def _local_fillDB_AndClose(self, sqlCommandList):
        if DEBUG: print "Filling SQLite database SessionDB.py"
        for sqlCommand in sqlCommandList:
            self.cursor.execute(sqlCommand)
        self.connection.commit()
        if DEBUG: print "Finished filling SQLite database SessionDB.py"

    def MakeNewDB(self, subject_data_file, mountPrefix):
        ## First close so that we can delete.
        self.close_connection()
        if os.path.exists(self.dbName):
            os.remove(self.dbName)
        self.open_connection()
        dbColTypes = "project TEXT, subj TEXT, session TEXT, type TEXT, Qpos INT, filename TEXT"
        self.cursor.execute("CREATE TABLE {tablename}({coltypes});".format(tablename=self.MasterTableName, coltypes=dbColTypes))
        self.connection.commit()
        sqlCommandList = list()
        missingFilesLog = self.dbName + "_MissingFiles.log"
        missingCount = 0
        if DEBUG: print("MISSING FILES RECORED IN {0}".format(missingFilesLog))
        missingFiles = open(missingFilesLog, 'w')
        if DEBUG: print "Building Subject returnList: " + subject_data_file
        subjData = csv.reader(open(subject_data_file, 'rb'), delimiter=',', quotechar='"')
        for linenumber, row in enumerate(subjData):
            if len(row) < 1 or row[0][0] == '#' or row[0] == 'project':
                continue  # empty row, comment, or header line
            elif len(row) != 4:
                raise RuntimeError("Invalid number of elements in subject data file {0}, line {0}", subject_data_file, linenumber)
            currDict = dict()
            validEntry = True
            currDict = {'project': row[0],
                        'subj': row[1],
                        'session': row[2]}
            rawDict = eval(row[3])
            for imageType in rawDict.keys():
                currDict['type'] = imageType
                fullPaths = [mountPrefix + i for i in rawDict[imageType]]
                if len(fullPaths) < 1:
                    if DEBUG: print("Invalid Entry!  {0}".format(currDict))
                    validEntry = False
                for i in range(len(fullPaths)):
                    imagePath = fullPaths[i]
                    if not os.path.exists(imagePath):
                        if DEBUG: print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Missing File: {0}\n".format(imagePath))
                        missingFiles.write("Missing File: {0}\n".format(imagePath))
                        validEntry = False
                        missingCount += 1
                    else:
                        if DEBUG: print("Found file {0}".format(imagePath))
                    if validEntry == True:
                        currDict['Qpos'] = str(i)
                        currDict['filename'] = imagePath
                        sqlCommand = self.makeSQLiteCommand(currDict)
                        sqlCommandList.append(sqlCommand)
        # sqlCommandList
        self._local_fillDB_AndClose(sqlCommandList)
        if missingCount > 0:
            if os.path.exists(self.dbName):
                os.remove(self.dbName)
            missingFiles.close()
            if DEBUG: print("ABORTING: At least 1 missing file\n"*20)
            sys.exit(-1)
        else:
            missingFiles.write("NO_MISSING_FILES")
        missingFiles.close()
        self.close_connection()

    def getSubjectFilter(self):
        return self.MasterQueryFilter

    def makeSQLiteCommand(self, imageDict):
        keys = imageDict.keys()
        vals = imageDict.values()
        col_names = ",".join(keys)
        values = ', '.join(map(lambda x: "'" + x + "'", vals))
        sqlCommand = "INSERT INTO {_tablename} ({_col_names}) VALUES ({_values});".format(
            _tablename=self.MasterTableName,
            _col_names=col_names, _values=values)
        return sqlCommand

    def getInfoFromDB(self, sqlCommand):
        #print("getInfoFromDB({0})".format(sqlCommand))
        self.open_connection()
        self.cursor.execute(sqlCommand)
        dbInfo = self.cursor.fetchall()
        self.close_connection()
        return dbInfo

    def getFirstScan(self, sessionid, scantype):
        sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' AND type='{_scantype}' AND Qpos='0';".format(
            _master_query=self.MasterQueryFilter,
            _sessionid=sessionid, _scantype=scantype)
        val = self.getInfoFromDB(sqlCommand)
        filename = str(val[0][0])
        return filename

    def getFirstT1(self, sessionid):
        scantype = 'T1-30'
        sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' AND type='{_scantype}' AND Qpos='0';".format(
            _master_query=self.MasterQueryFilter,
            _sessionid=sessionid, _scantype=scantype)
        val = self.getInfoFromDB(sqlCommand)
        # print "HACK: ",sqlCommand
        # print "HACK: ", val
        filename = str(val[0][0])
        return filename

    def getFilenamesByScantype(self, sessionid, scantypelist):
        returnList = list()
        for currScanType in scantypelist:
            sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' AND type='{_scantype}' ORDER BY Qpos ASC;".format(
                _master_query=self.MasterQueryFilter,
                _sessionid=sessionid, _scantype=currScanType)
            val = self.getInfoFromDB(sqlCommand)
            for i in val:
                returnList.append(str(i[0]))
        return returnList

    def findScanTypeLength(self, sessionid, scantypelist):
        countList = self.getFilenamesByScantype(sessionid, scantypelist)
        return len(countlist)

    def getT1sT2s(self, sessionid):
        sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' ORDER BY type ASC, Qpos ASC;".format(
            _master_query=self.MasterQueryFilter,
            _sessionid=sessionid)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getAllProjects(self):
        sqlCommand = "SELECT DISTINCT project FROM ({_master_query});".format(_master_query=self.MasterQueryFilter)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getAllSubjects(self):
        sqlCommand = "SELECT DISTINCT subj FROM ({_master_query});".format(_master_query=self.MasterQueryFilter)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getAllSessions(self):
        # print("HACK:  This is a temporary until complete re-write")
        sqlCommand = "SELECT DISTINCT session FROM ({_master_query});".format(_master_query=self.MasterQueryFilter)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getSessionsFromSubject(self, subj):
        sqlCommand = "SELECT DISTINCT session FROM ({_master_query}) WHERE subj='{_subjid}';".format(
            _master_query=self.MasterQueryFilter,
            _subjid=subj)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getSessionsFromGroup(self, groupfile):
        """ Create a dictionary of {group: {subject: [sessions]}, ...}

        Keyword Arguments:
        self   -- SessionDB object
        groupfile -- path to group dictionary file

        NOTA BENE: groupfile MUST have a Python dictionary named 'groups' defined within it

        >>> def test():
        ...   import os, os.path
        ...   import SessionDB
        ...   db = SessionDB.SessionDB(subject_list=['all'])
        ...   dbfile = os.path.join(os.getcwd(), 'test_SessionDB.csv')
        ...   db.MakeNewDB(dbfile, '')
        ...   groupfile = os.path.abspath('test_groupfile.pydict')
        ...   return db.getSessionsFromGroup(groupfile)
        >>> test()
        {'groupA': ['63819', '59911'], 'groupC': [], 'groupB': ['29876']}
        """
        execfile(groupfile, globals())
        assert 'groups' in globals(), "'groups' not found in {0}".format(groupfile)
        allSessions = set(self.getAllSessions())
        for k, v in groups.items():
            valid = list(set(v).intersection(allSessions))
            if len(valid) == 0:
                if DEBUG: print "WARNING: Group {0} is empty!".format(k)
            groups[k] = valid
        return groups


    def getEverything(self):
        sqlCommand = "SELECT * FROM ({_master_query});".format(_master_query=self.MasterQueryFilter)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(i)
        return returnList

    def getSubjectsFromProject(self, project):
        sqlCommand = "SELECT DISTINCT subj FROM ({_master_query}) WHERE project='{_projectid}';".format(
            _master_query=self.MasterQueryFilter,
            _projectid=project)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getSubjFromSession(self, session):
        sqlCommand = "SELECT DISTINCT subj FROM ({_master_query}) WHERE session='{_sessionid}';".format(
            _master_query=self.MasterQueryFilter,
            _sessionid=session)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        if len(returnList) != 1:
            print("ERROR: More than one subject found")
            sys.exit(-1)
        return returnList[0]

    def getProjFromSession(self, session):
        sqlCommand = "SELECT DISTINCT project FROM ({_master_query}) WHERE session='{_sessionid}';".format(
            _master_query=self.MasterQueryFilter,
            _sessionid=session)
        val = self.getInfoFromDB(sqlCommand)
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
# a.getFirstScan('42245','T1-30')
# a.getInfoFromDB("SELECT filename FROM SessionDB WHERE session=42245 ORDER BY type ASC, Qpos ASC;")
# a.getInfoFromDB("SELECT DISTINCT subj FROM SessionDB;")

if __name__ == '__main__':
    import os.path
    import doctest

    doctest.testmod(verbose=True, raise_on_error=True)
    # dirname = os.path.dirname(__file__)
    # doctest.testfile(os.path.join(dirname, 'tests/SessionDB.test'))
