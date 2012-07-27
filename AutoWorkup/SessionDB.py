import os
import sqlite3 as lite
import csv

class SessionDB():

    def __init__(self, autoworkupFile, mountPrefix):
        self.dbName = 'SessionDB.db'
        sqlCommandList = self.parseSubjDataFile(autoworkupFile, mountPrefix)
        self.makeDB()
        self.fillDB(sqlCommandList)

    def parseSubjDataFile(self, subject_data_file, mountPrefix):
        sqlCommandList = list()
        print "Building Subject List: " + subject_data_file
        subjData=csv.reader(open(subject_data_file,'rb'), delimiter=',', quotechar='"')
        for row in subjData:
            if row[0] == 'project':
                # continue if header line
                continue
            currDict=dict()
            validEntry=True
            if len(row) == 4:
                currDict = {'site': row[0],
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
        return sqlCommandList

    def makeDB(self):
        if os.path.exists(self.dbName):
            os.remove(self.dbName)
        con = lite.connect(self.dbName)
        dbCur = con.cursor()
        dbColTypes =  "site TEXT, subj TEXT, session TEXT, type TEXT, Qpos INT, filename TEXT"
        dbCur.execute("CREATE TABLE SessionDB({0});".format(dbColTypes))
        dbCur.close()

    def fillDB(self, sqlCommandList):
        con = lite.connect(self.dbName)
        dbCur = con.cursor()
        print "Filling SQLite database SessionDB.py"
        for sqlCommand in sqlCommandList:
            dbCur.execute(sqlCommand)
            con.commit()
        dbCur.close()
        print "Finished filling SQLite database SessionDB.py"

    def makeSQLiteCommand(self, imageDict):
        keys = imageDict.keys()
        vals = imageDict.values()
        col_names = ",".join(keys)
        values = ', '.join(map(lambda x: "'" + x + "'", vals))
        sqlCommand = "INSERT INTO SessionDB ({0}) VALUES ({1});".format(col_names, values)
        return sqlCommand

    def getInfoFromDB(self, sqlCommand):
        con = lite.connect(self.dbName)
        dbCur = con.cursor()
        dbCur.execute(sqlCommand)
        dbInfo = dbCur.fetchall()
        dbCur.close()
        return dbInfo

    def getFirstScan(self, uid, scantype):
        sqlCommand = "SELECT filename FROM SessionDB WHERE session={0} AND type='{1}' AND Qpos=0;".format(uid, scantype)
        val = self.getInfoFromDB(sqlCommand)
        filename = str(val[0][0])
        return filename

    def getFilenamesByScantype(self, uid, scantype):
        sqlCommand = "SELECT filename FROM SessionDB WHERE session={0} AND type='{1}' ORDER BY Qpos ASC;".format(uid, scantype)
        val = self.getInfoFromDB(sqlCommand)
        List = list()
        for i in val:
            List.append(str(i[0]))
        return List

    def findScanTypeLength(self, uid, scantype):
        sqlCommand = "SELECT COUNT(DISTINCT filename) FROM SessionDB WHERE session={0} AND type='{1}';".format(uid, scantype)
        val = self.getInfoFromDB(sqlCommand)
        count = val[0][0]
        return count

    def getT1sT2s(self, uid):
        sqlCommand = "SELECT filename FROM SessionDB WHERE session={0} ORDER BY type ASC, Qpos ASC;".format(uid)
        val = self.getInfoFromDB(sqlCommand)
        List = list()
        for i in val:
            List.append(str(i[0]))
        return List
