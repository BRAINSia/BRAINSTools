import os
import sqlite3 as lite
import csv

class SessionDB():

    def __init__(self, defaultDBName='TempFileForDB.db',subject_filter='00000000000'):
        self.MasterTableName = "MasterDB"
        self.dbName = defaultDBName
        self._local_openDB()
        self.MasterQueryFilter = "SELECT * FROM {_tablename} WHERE subj='{_subjid}'".format(
          _tablename=self.MasterTableName,
          _subjid=subject_filter)

    def _local_openDB(self):
        self.connection = lite.connect(self.dbName)
        self.cursor = self.connection.cursor()

    def _local_fillDB(self, sqlCommandList):
        print "Filling SQLite database SessionDB.py"
        for sqlCommand in sqlCommandList:
            self.cursor.execute(sqlCommand)
        self.connection.commit()
        print "Finished filling SQLite database SessionDB.py"

    def MakeNewDB(self, subject_data_file, mountPrefix):
        ## First close so that we can delete.
        self.cursor.close()
        self.connection.close()
        if os.path.exists(self.dbName):
            os.remove(self.dbName)
        self._local_openDB()

        dbColTypes =  "project TEXT, subj TEXT, session TEXT, type TEXT, Qpos INT, filename TEXT"
        self.cursor.execute("CREATE TABLE {tablename}({coltypes});".format(tablename=self.MasterTableName,coltypes=dbColTypes))
        self.connection.commit()
        sqlCommandList = list()
        print "Building Subject returnList: " + subject_data_file
        subjData=csv.reader(open(subject_data_file,'rb'), delimiter=',', quotechar='"')
        for row in subjData:
            if len(row) < 1:
                # contine of it is an empty row
                continue
            if row[0][0] == '#':
                # if the first character is a #, then it is commented out
                continue
            if row[0] == 'project':
                # continue if header line
                continue
            currDict=dict()
            validEntry=True
            if len(row) == 4:
                currDict = {'project': row[0],
                            'subj': row[1],
                            'session': row[2]}
                rawDict=eval(row[3])
                for imageType in rawDict.keys():
                    currDict['type'] = imageType
                    fullPaths=[ mountPrefix+i for i in rawDict[imageType] ]
                    if len(fullPaths) < 1:
                        print("Invalid Entry!  {0}".format(currDict))
                        validEntry=False
                    for i in range(len(fullPaths)):
                        imagePath = fullPaths[i]
                        if not os.path.exists(imagePath):
                            print("Missing File: {0}".format(imagePath))
                            validEntry=False
                        if validEntry == True:
                            currDict['Qpos'] = str(i)
                            currDict['filename'] = imagePath
                            sqlCommand = self.makeSQLiteCommand(currDict)
                            sqlCommandList.append(sqlCommand)
            else:
                print "ERROR:  Invalid number of elements in row"
                print row
        sqlCommandList
        self._local_fillDB(sqlCommandList)

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
        self.cursor.execute(sqlCommand)
        dbInfo = self.cursor.fetchall()
        return dbInfo

    def getFirstScan(self, sessionid, scantype):
        sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' AND type='{_scantype}' AND Qpos='0';".format(
            _master_query=self.MasterQueryFilter,
            _sessionid=sessionid, _scantype=scantype)
        val = self.getInfoFromDB(sqlCommand)
        filename = str(val[0][0])
        return filename

    def getFirstT1(self, sessionid):
        scantype='T1-30'
        sqlCommand = "SELECT filename FROM ({_master_query}) WHERE session='{_sessionid}' AND type='{_scantype}' AND Qpos='0';".format(
          _master_query=self.MasterQueryFilter,
          _sessionid=sessionid, _scantype=scantype)
        val = self.getInfoFromDB(sqlCommand)
        #print "HACK: ",sqlCommand
        #print "HACK: ", val
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
        countList=self.getFilenamesByScantype(sessionid,scantypelist)
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
        #print("HACK:  This is a temporary until complete re-write")
        sqlCommand = "SELECT DISTINCT session FROM ({_master_query});".format(_master_query=self.MasterQueryFilter)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getSessionsFromSubject(self,subj):
        sqlCommand = "SELECT DISTINCT session FROM ({_master_query}) WHERE subj='{_subjid}';".format(
          _master_query=self.MasterQueryFilter,
          _subjid=subj)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getEverything(self):
        sqlCommand = "SELECT * FROM ({_master_query});".format(_master_query=self.MasterQueryFilter)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(i)
        return returnList

    def getSubjectsFromProject(self,project):
        sqlCommand = "SELECT DISTINCT subj FROM ({_master_query}) WHERE project='{_projectid}';".format(
          _master_query=self.MasterQueryFilter,
          _projectid=project)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList


    def getSubjFromSession(self,session):
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

    def getProjFromSession(self,session):
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
