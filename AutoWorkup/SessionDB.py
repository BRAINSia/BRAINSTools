import os
import sqlite3 as lite
import csv

class SessionDB():

    def __init__(self, defaultDBName='TempFileForDB.db'):
        self.dbName = defaultDBName

    def _local_createDB(self):
        if os.path.exists(self.dbName):
            os.remove(self.dbName)
        local_connection = lite.connect(self.dbName)
        dbCur = local_connection.cursor()
        dbColTypes =  "project TEXT, subj TEXT, session TEXT, type TEXT, Qpos INT, filename TEXT"
        dbCur.execute("CREATE TABLE SessionDB({0});".format(dbColTypes))
        dbCur.close()
        local_connection.commit()
        local_connection.close()
        local_connection=None

    def _local_fillDB(self, sqlCommandList):
        local_connection = lite.connect(self.dbName)
        dbCur = local_connection.cursor()
        print "Filling SQLite database SessionDB.py"
        for sqlCommand in sqlCommandList:
            dbCur.execute(sqlCommand)
            local_connection.commit()
        dbCur.close()
        local_connection.close()
        local_connection=None
        print "Finished filling SQLite database SessionDB.py"


    def MakeNewDB(self, subject_data_file, mountPrefix):
        self._local_createDB()
        sqlCommandList = list()
        print "Building Subject returnList: " + subject_data_file
        subjData=csv.reader(open(subject_data_file,'rb'), delimiter=',', quotechar='"')
        for row in subjData:
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
                            print("Missing File: {0}".format(i))
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

    def makeSQLiteCommand(self, imageDict):
        keys = imageDict.keys()
        vals = imageDict.values()
        col_names = ",".join(keys)
        values = ', '.join(map(lambda x: "'" + x + "'", vals))
        sqlCommand = "INSERT INTO SessionDB ({0}) VALUES ({1});".format(col_names, values)
        return sqlCommand

    def getInfoFromDB(self, sqlCommand):
        local_connection = lite.connect(self.dbName)
        dbCur = local_connection.cursor()
        #print("getInfoFromDB({0})".format(sqlCommand))
        dbCur.execute(sqlCommand)
        dbInfo = dbCur.fetchall()
        dbCur.close()
        local_connection.close()
        local_connection=None
        return dbInfo

    def getFirstScan(self, sessionid, scantype):
        sqlCommand = "SELECT filename FROM SessionDB WHERE session='{0}' AND type='{1}' AND Qpos='0';".format(sessionid, scantype)
        val = self.getInfoFromDB(sqlCommand)
        filename = str(val[0][0])
        return filename

    def getFirstT1(self, sessionid):
        scantype='T1-30'
        sqlCommand = "SELECT filename FROM SessionDB WHERE session='{0}' AND type='{1}' AND Qpos='0';".format(sessionid, scantype)
        val = self.getInfoFromDB(sqlCommand)
        print "HACK: ",sqlCommand
        print "HACK: ", val
        filename = str(val[0][0])
        return filename

    def getFilenamesByScantype(self, sessionid, scantype):
        sqlCommand = "SELECT filename FROM SessionDB WHERE session='{0}' AND type='{1}' ORDER BY Qpos ASC;".format(sessionid, scantype)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def findScanTypeLength(self, sessionid, scantype):
        sqlCommand = "SELECT COUNT(DISTINCT filename) FROM SessionDB WHERE session='{0}' AND type='{1}';".format(sessionid, scantype)
        val = self.getInfoFromDB(sqlCommand)
        count = val[0][0]
        return count

    def getT1sT2s(self, sessionid):
        sqlCommand = "SELECT filename FROM SessionDB WHERE session='{0}' ORDER BY type ASC, Qpos ASC;".format(sessionid)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getAllProjects(self):
        sqlCommand = "SELECT DISTINCT project FROM SessionDB;"
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getAllSubjects(self):
        sqlCommand = "SELECT DISTINCT subj FROM SessionDB;"
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getAllSessions(self):
        #print("HACK:  This is a temporary until complete re-write")
        sqlCommand = "SELECT DISTINCT session FROM SessionDB;"
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getSessionsFromProjectSubject(self,proj,subj):
        sqlCommand = "SELECT DISTINCT session FROM SessionDB WHERE project='{0}' AND subj='{1}';".format(proj,subj)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList

    def getEverything(self):
        sqlCommand = "SELECT * FROM SessionDB;"
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(i)
        return returnList

    def getSubjectsFromProject(self,project):
        sqlCommand = "SELECT DISTINCT subj FROM SessionDB WHERE project='{0}';".format(project)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        return returnList


    def getSubjFromSession(self,session):
        sqlCommand = "SELECT DISTINCT subj FROM SessionDB WHERE session='{0}';".format(session)
        val = self.getInfoFromDB(sqlCommand)
        returnList = list()
        for i in val:
            returnList.append(str(i[0]))
        if len(returnList) != 1:
            print("ERROR: More than one subject found")
            sys.exit(-1)
        return returnList[0]

    def getProjFromSession(self,session):
        sqlCommand = "SELECT DISTINCT project FROM SessionDB WHERE session='{0}';".format(session)
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
