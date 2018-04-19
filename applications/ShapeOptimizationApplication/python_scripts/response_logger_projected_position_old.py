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
import time as timer

# ==============================================================================
class ResponseLoggerProjectedPositionOld( ResponseLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, communicator, optimizationSettings, timer ):
        self.communicator = communicator
        self.optimizationSettings = optimizationSettings
        self.timer = timer

        ##################################
        # get response function params
        self.objectiveName = optimizationSettings["objectives"][0]["identifier"].GetString()
        self.objectiveType = optimizationSettings["objectives"][0]["type"].GetString()

        self.constraintNames = []
        self.constraintType = []
        self.constraintReference = []
        self.equalityNames = []
        self.equalityTarget = []
        self.equalityTolerance = []

        for i in range(self.optimizationSettings["constraints"].size()):
            cname = optimizationSettings["constraints"][i]["identifier"].GetString()
            ctype = optimizationSettings["constraints"][i]["type"].GetString()
            if ctype == "=":
                self.equalityNames.append(cname)
                if optimizationSettings["constraints"][i]["reference"].IsDouble():
                    self.equalityTarget.append(optimizationSettings["constraints"][i]["reference"].GetDouble())
                else:
                    self.equalityTarget.append(optimizationSettings["constraints"][i]["reference"].GetString()) # for "initial value"
                # self.equalityTolerance.append(optimizationSettings["constraints"][i]["tolerance"].GetDouble())
            else:
                self.constraintNames.append(cname)
                self.constraintType.append(ctype)
                if optimizationSettings["constraints"][i]["reference"].IsDouble():
                    self.constraintReference.append(optimizationSettings["constraints"][i]["reference"].GetDouble())
                else:
                    self.constraintReference.append(optimizationSettings["constraints"][i]["reference"].GetString())
        ####################################
        ## database and print
        self.__timeFirstLog = None

        self.completeResponseLogFileName = self.__CreateCompleteResponseLogFilename( optimizationSettings )

        self.__namePrintFirst = ["iterOuter","iterInner","stepLength","sDx","sInit","lConstraint","lEquality","constraint","equality","objective","objectivePercent","endIter"]
        #remove unused names
        if len(self.constraintNames)==0:
            for name in self.__namePrintFirst:
                if "constraint" in name.lower():
                    self.__namePrintFirst.remove(name)
        if len(self.equalityNames)==0:
            for name in self.__namePrintFirst:
                if "equality" in name.lower():
                    self.__namePrintFirst.remove(name)

        self.__db = []
        self.__nColumnByName = {}

    # --------------------------------------------------------------------------
    def __CreateCompleteResponseLogFilename( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["output"]["output_directory"].GetString()
        responseLogFilename = optimizationSettings["output"]["response_log_filename"].GetString()
        completeResponseLogFilename = resultsDirectory+"/"+responseLogFilename+".csv"
        return completeResponseLogFilename

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        pass

    # return "-" if not defined
    def __GetValue(self,i,key):
        dico = self.__db[i]
        if key not in dico:
            return ["-"]
        return dico[key]

    def __AddToDb(self,dicoIter,p):
        dicoIter["time"] = int(round(self.timer.GetTotalTime()))

        # complete the values of dicoIter
        # if "constraint" in dicoIter:
            # add constraintFull
            # dicoIter["constraintFull"] = dicoIter["constraint"][:]
            # for i,(vconstr,vref) in enumerate(zip(dicoIter["constraintFull"],p.constraintReference)):
            #     if vconstr is not None:
            #         if self.constraintType[i]=="<":
            #             dicoIter["constraintFull"][i] = vconstr + vref
            #         if self.constraintType[i]==">":
            #             dicoIter["constraintFull"][i] = vref - vconstr
            # add gradConstraintLActual
            # if "lConstraint" in dicoIter:
            #     dicoIter["gradConstraintLActual"] = dicoIter["constraint"][:]
            #     for i,(v,l) in enumerate(zip(dicoIter["constraint"],dicoIter["lConstraint"])):
            #         if v is not None:
            #             dicoIter["gradConstraintLActual"][i] = v/l if l!=0 else None

        # if "equality" in dicoIter:
            # add equalityFull
            # dicoIter["equalityFull"] = dicoIter["equality"][:]
            # for i, (veq,target,tolerance) in enumerate( zip(dicoIter["equalityFull"],p.equalityTarget,p.equalityTolerance) ):
            #     if veq is not None:
            #         dicoIter["equalityFull"][i] = veq * (target*tolerance) + target
            # add gradEqualityLActual
            # if "lEquality" in dicoIter and "equalityOuter" in dicoIter:
            #     dicoIter["gradEqualityLActual"] = dicoIter["equality"][:]
            #     for i,(v,vOuter,l) in enumerate(zip(dicoIter["equality"],dicoIter["equalityOuter"],dicoIter["lEquality"])):
            #         if v is not None:
            #             dicoIter["gradEqualityLActual"][i] = (v-vOuter)/l if l!=0 else None
            # if "iterOuter" in dicoIter:
            #     dicoIter["equalityTolerances"] = dicoIter["equality"][:]
            #     for i in range(len(dicoIter["equality"])):
            #         dicoIter["equalityTolerances"][i] = self.__GetEqualityTolerance(i,dicoIter["iterOuter"])

        if ("objective" in dicoIter) and (dicoIter["objective"] is not None) and len(self.__db)>0:
            objectiveInit = self.__db[0]["objective"][0]
            dicoIter["objectivePercent"] = "{:.0f}%".format((objectiveInit-dicoIter["objective"])/objectiveInit*100)
            if "objectiveOuter" in dicoIter and "lObjective" in dicoIter :
                dicoIter["gradobjectiveLActual"] = (dicoIter["objective"]-dicoIter["objectiveOuter"])/(dicoIter["lObjective"] if dicoIter["lObjective"]!=0 else 1)

        for k,v in dicoIter.items():
            # all values are list
            if not isinstance(v,list):
                dicoIter[k] = [v]
                v = [v]

            # replace None by "-"
            for i in range(len(v)):
                if v[i] is None:
                    dicoIter[k][i] = "-"

            if k not in self.__nColumnByName:
                self.__nColumnByName[k] = max(len(v),1)
            else:
                self.__nColumnByName[k] = max(len(v),self.__nColumnByName[k])

        self.__db.append(dicoIter)



    # write the database in an array with the right columns and no list
    def __DbToArray(self):
        ar = []
        names = self.__namePrintFirst + ["others-->"] +  list( set(self.__nColumnByName.keys())-set(self.__namePrintFirst) )

        # headers
        headerRow = []
        for name in names:
            length = 1 if name not in self.__nColumnByName else self.__nColumnByName[name]
            headerRow += [name] + ["" for i in range(length-1)] # "" place holders for further columns
        ar.append(headerRow)

        # fill array with values, and put the columns in the right order
        for i in range(len(self.__db)):
            row = []
            for name in names:
                value = self.__GetValue(i,name)
                length = 1 if name not in self.__nColumnByName else self.__nColumnByName[name]
                row += value + ["-" for i in range(length-len(value))] # fill with "-" if values are missing
            ar.append(row)
        return ar

    # take array of single values (no list) as an input, first row is the name of the columns
    # convert to string with the right format and write the csv file
    def __WriteArrayToCsv(self,ar):
        with open(self.completeResponseLogFileName, 'w') as csvfile:
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

    def LogIteration(self,dicoIter,p):
        self.__AddToDb(dicoIter,p)
        self.__WriteArrayToCsv(self.__DbToArray())
        pass

    def FinalizeLogging( self ):
        pass # No finalization necessary here

