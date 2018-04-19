# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Import logger base classes
from response_logger_base import ResponseLogger

# Import additional libraries
import csv
import json
import time as timer
from copy import deepcopy

# ==============================================================================
class ResponseLoggerProjectedPosition( ResponseLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, communicator, optimizationSettings, timer ):
        self.communicator = communicator
        self.optimizationSettings = optimizationSettings
        self.timer = timer

        # columns to be printed on the left of the csv file
        self.__namePrintFirst = ["iterCounter","stepLength","sDx","sInit","lInequality","lEquality","inequality","equality","objective","objectivePercent"]


    def __FilterNamePrintFirst(self,p):
        #remove unused names
        if len(p.inequalityIds)==0:
            for name in self.__namePrintFirst:
                if "inequality" in name.lower():
                    self.__namePrintFirst.remove(name)
        if len(p.equalityIds)==0:
            for name in self.__namePrintFirst:
                if "equality" in name.lower():
                    self.__namePrintFirst.remove(name)

    # --------------------------------------------------------------------------
    def __GetPathCsvFile( self, optimizationSettings ):
        path = "../results/"+self.optimizationSettings.runId+"/"+self.optimizationSettings.runId+".csv"
        return path

    # --------------------------------------------------------------------------
    def __GetPathJsonFile( self, optimizationSettings):
        path = "../results/"+self.optimizationSettings.runId+"/"+self.optimizationSettings.runId+".json"
        return path

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        pass

    # return "-" if not defined
    def __GetValueForCsv(self,i,key,p):
        if key not in p.r:
            return ["-"]
        val = p.r[key][i]
        if not isinstance(val,list):
            val = [val]
        if len(val)==0:
            return ["-"]
        for i in range(len(val)):
            if val[i] is None or (isinstance(val[i],list) and len(val[i])==0):
                val[i] = "-"
        return val

    def __GetLength(self,name,p):
        if name not in p.r:
            return 1
        if len(p.r[name])>0 and isinstance(p.r[name][0],list):
            return len(p.r[name][0])
        return 1

    def __FillValues(self,p):
        r = p.r

        # initialize
        for name in ["time","inequalityFull","equalityFull","objectivePercent","lObjectiveActual","lInequalityActual","lEqualityActual"]:
            if name not in r:
                r[name] = []
        r["inequalityFull"].append([None for _ in r["inequality"][0]])
        r["lInequalityActual"].append([None for _ in r["inequality"][0]])
        r["equalityFull"].append([None for _ in r["equality"][0]])
        r["lEqualityActual"].append([None for _ in r["equality"][0]])

        r["time"].append(int(round(self.timer.GetTotalTime())))

        # add inequalityFull
        for i,(vconstr,vref) in enumerate(zip(r["inequality"][-1],p.inequalityReferences)):
            if p.inequalityTypes[i]=="<":
                r["inequalityFull"][-1][i] = vconstr + vref
            if p.inequalityTypes[i]==">":
                r["inequalityFull"][-1][i] = vref - vconstr

            if len(r["iterCounter"])>=2:
                r["lInequalityActual"][-2][i] = (r["inequality"][-1][i]-r["inequality"][-2][i])/r["gradInequalityL"][-2][i]

        # add equalityFull
        for i, (veq,vref) in enumerate( zip(r["equality"][-1],p.equalityReferences) ):
            r["equalityFull"][-1][i] = veq+vref

            if len(r["iterCounter"])>=2:
                r["lEqualityActual"][-2][i] = (r["equality"][-1][i]-r["equality"][-2][i])/r["gradEqualityL"][-2][i]

        objectiveInit = r["objective"][0]
        objectivePercent = (objectiveInit-r["objective"][-1])/objectiveInit*100
        r["objectivePercent"].append("{:.0f}%".format(objectivePercent))

        r["lObjectiveActual"].append(None)
        if len(r["iterCounter"])>=2:
            r["lObjectiveActual"][-2] = (r["objective"][-1]-r["objective"][-2])/r["gradObjectiveL"][-2]


    # write the database in an array with the right columns and no list
    def __DictToArray(self,p):
        ar = []

        names = self.__namePrintFirst + ["others-->"] +  list( set(p.r.keys())-set(self.__namePrintFirst) )

        # headers
        headerRow = []
        for name in names:
            length = self.__GetLength(name,p)
            headerRow += [name] + ["" for i in range(length-1)] # "" place holders for further columns
        ar.append(headerRow)

        # fill array with values, and put the columns in the right order
        for i in range(len(p.r["iterCounter"])):
            row = []
            for name in names:
                value = self.__GetValueForCsv(i,name,p)
                length = self.__GetLength(name,p)
                row += value + ["-" for i in range(length-len(value))] # fill with "-" if values are missing
            ar.append(row)
        return ar

    # take array of single values (no list) as an input, first row is the name of the columns
    # convert to string with the right format and write the csv file
    def __WriteArrayToCsv(self,ar):
        with open(self.__GetPathCsvFile(self.optimizationSettings), 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            headerRow = ar[0][:]
            #convert values to string (only num and string, no list)
            for row in ar:
                pr = row[:]
                for i,v in enumerate(row):
                    width = str(6) if headerRow[i]=="" else str(max(12,len(headerRow[i])+1))
                    if isinstance(v,(str,int)):
                        pr[i] =  ("{:>"+width+"s}").format(str(v))
                    elif isinstance(v,float):
                        if abs(v)<1000 and abs(v)>0.01:
                            pr[i] = ("{:>"+width+".2f}").format(v)
                        else:
                            pr[i] = ("{:>"+width+".2g}").format(v)
                historyWriter.writerow(pr)

    def __WriteParamResultsToJson(self,p):
        dataDict = {"results": p.r, "param": p.param}
        with open(self.pathJsonFile, 'w') as jsonfile:
            json.dump(dataDict, jsonfile, sort_keys=True, indent=4)

    def LogIteration(self,dicoIter,p):
        self.pathJsonFile = self.__GetPathJsonFile(self.optimizationSettings)

        self.__FilterNamePrintFirst(p)
        self.__FillValues(p)


        if len(p.r["iterCounter"])!=p.r["iterCounter"][-1]:
            print(p.r["iterCounter"])
            err

        self.__WriteArrayToCsv(self.__DictToArray(p))
        self.__WriteParamResultsToJson(p)

    def FinalizeLogging( self ):
        pass # No finalization necessary here

