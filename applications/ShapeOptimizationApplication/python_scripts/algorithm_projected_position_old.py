# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Rousseau Pierre, https://github.com/pr-pierre-rousseau
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Import algorithm base classes
from algorithm_base import OptimizationAlgorithm

# Additional imports
import timer_factory as timer_factory
import logging
import contextlib

from projected_position_modules.matrix import norminf3d
from projected_position_modules.projection import plusSphereOld

# ==============================================================================
class AlgorithmProjectedPositionOld( OptimizationAlgorithm ) :

    # --------------------------------------------------------------------------
    def __init__( self,
                  ModelPartController,
                  Analyzer,
                  Communicator,
                  Mapper,
                  DataLogger,
                  OptimizationSettings ):

        self.ModelPartController = ModelPartController
        self.Analyzer = Analyzer
        self.Communicator = Communicator
        self.Mapper = Mapper
        self.DataLogger = DataLogger
        self.OptimizationSettings = OptimizationSettings

        self.OptimizationModelPart = ModelPartController.GetOptimizationModelPart()
        self.DesignSurface = ModelPartController.GetDesignSurface()

        # self.projectionOnNormalsIsSpecified = OptimizationSettings["optimization_algorithm"]["project_gradients_on_surface_normals"].GetBool()
        self.dampingIsSpecified = OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()

        self.GeometryUtilities = GeometryUtilities( self.DesignSurface )
        self.OptimizationUtilities = OptimizationUtilities( self.DesignSurface, OptimizationSettings )
        if self.dampingIsSpecified:
            damping_regions = self.ModelPartController.GetDampingRegions()
            self.DampingUtilities = DampingUtilities( self.DesignSurface, damping_regions, self.OptimizationSettings )

        Communicator.p = self

        ######################################################################
        # Algorithm parameters
        osa = OptimizationSettings["optimization_algorithm"]
        self.iterOuterMax = osa["max_outer_iterations"].GetInt()
        self.iterInnerMax = osa["max_inner_iterations"].GetInt()
        self.maxStepLength = osa["max_step_length"].GetDouble()

        self.dlInequalityIncrRel = osa["delta_l_inequality"].GetDouble()
        self.dlEqualityIncrRel = osa["delta_l_equality"].GetDouble()
        self.dlUnfeasibleRel = osa["delta_l_infeasible"].GetDouble()

        self.enforceDecreasingObjectiveIfUnfeasible = osa["decreasing_objective_if_infeasible"].GetBool()
        self.optimizeInitOfLInner = osa["efficient_initialization_of_l"].GetBool()
        self.lEqualityInitPersistence = osa["l_equality_init_persistence"].GetDouble()
        self.lEqualityInitInterval = osa["l_equality_init_interval"].GetDouble()
        self.stepLengthDecreaseFactor = osa["step_length_decrease_factor"].GetDouble()
        self.stepLengthIncreaseFactor = osa["step_length_increase_factor"].GetDouble()
        self.stepLengthTermination = osa["step_length_termination"].GetDouble()

        defaultParam = {
            "sphere_norm": "geom_norminf3d",
            "sphere_sum": "addition"
        }

        if osa.Has("sphere_norm"):  # "geom_norminf3d" or "coord_norm1"
            self.sphereNorm = osa["sphere_norm"].GetString()
        else:
            self.sphereNorm = defaultParam["sphere_norm"]

        if osa.Has("sphere_sum"): # "addition" or "projection"
            self.sphereSum = osa["sphere_sum"].GetString()
        else:
            self.sphereSum = defaultParam["sphere_sum"]

        ##########################################################################
        # Response functions

        self.objectiveId = OptimizationSettings["objectives"][0]["identifier"].GetString()
        self.inequalityConstraintIds = []
        self.constraintTypes = []
        self.constraintReferences = []
        self.equalityConstraintIds = []
        # self.equalityTolerances = []
        self.equalityReferences = []

        # read the constraint and equality response function in the params
        for i in range(self.OptimizationSettings["constraints"].size()):
            cname = OptimizationSettings["constraints"][i]["identifier"].GetString()
            ctype = OptimizationSettings["constraints"][i]["type"].GetString()
            if ctype == "=":
                self.equalityConstraintIds.append(cname)
                # self.equalityTolerances.append(OptimizationSettings["constraints"][i]["tolerance"].GetDouble())
                if OptimizationSettings["constraints"][i]["reference"].IsDouble():
                    self.equalityReferences.append(OptimizationSettings["constraints"][i]["reference"].GetDouble())
                else:
                    self.equalityReferences.append(OptimizationSettings["constraints"][i]["reference"].GetString())
            else:
                self.inequalityConstraintIds.append(cname)
                self.constraintTypes.append(ctype)
                if OptimizationSettings["constraints"][i]["reference"].IsDouble():
                    self.constraintReferences.append(OptimizationSettings["constraints"][i]["reference"].GetDouble())
                else:
                    self.constraintReferences.append(OptimizationSettings["constraints"][i]["reference"].GetString())

        self.n = 3*self.DesignSurface.NumberOfNodes()
        self.nc = len(self.inequalityConstraintIds)
        self.ne = len(self.equalityConstraintIds)
        self.xInit = self.__ReadDesignSurfaceToList()
        self.objectiveInit = None

    # --------------------------------------------------------------------------
    def initializeOptimizationLoop( self ):
        self.Analyzer.initializeBeforeOptimizationLoop()
        self.ModelPartController.InitializeMeshController()
        self.DataLogger.StartTimer()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def runOptimizationLoop( self ):

        # init variables
        gradInequalityOuter = [0 for _ in range(self.nc)] # item will be replaced by lists
        gradEqualityOuter = [0 for _ in range(self.ne)] # item will be replaced by lists
        eConstraint = [0 for _ in range(self.nc)] # item will be replaced by lists
        eEquality = [0 for _ in range(self.ne)] # item will be replaced by lists
        gradInequalityL = [0 for _ in range(self.nc)]
        gradEqualityL = [0 for _ in range(self.ne)]
        inequalityOuter = [0 for _ in range(self.nc)]
        equalityOuter = [0 for _ in range(self.ne)]
        lConstraint = [0 for _ in range(self.nc)]
        lEquality = [0 for _ in range(self.ne)]
        lInequality0 = [0 for _ in range(self.nc)]
        lEquality0 = [0 for _ in range(self.ne)]
        timeConstraint = [0 for _ in range(self.nc)]
        timeEquality = [0 for _ in range(self.ne)]

        self.iterOuter = 0
        self.iterCounter = 0 # cumulative number of inner iterations
        xOuter = self.xInit[:]
        xInner = self.xInit[:]
        self.DataLogger.DesignLogger.LogCurrentDesign(0)

        # Compute response for initial design
        timer = timer_factory.CreateTimer()
        timer.StartNewLap()
        objectiveOuter = self.__GetObjective()
        timeObjective = timer.GetLapTime()
        for i in range(self.nc):
            timer.StartNewLap()
            inequalityOuter[i] = self.__GetInequalityConstraint(i)
            timeConstraint[i] = timer.GetLapTime()
        for i in range(self.ne):
            timer.StartNewLap()
            equalityOuter[i] = self.__GetEqualityConstraint(i)
            timeEquality[i] = timer.GetLapTime()
        self.callOrder = self.__GetCallOrder(timeObjective,timeConstraint,timeEquality) # build callOrder by increasing computation time
        print("Algo> callOrder = ",self.callOrder)
        self.DataLogger.ResponseLogger.LogIteration({"iterOuter": 0, "inequality":inequalityOuter,"equality": equalityOuter, "objective": objectiveOuter},self)

        # get the constraints and equality violated by xInit
        iInequalityInfeasible = [i for i,vconstr in enumerate(inequalityOuter) if vconstr>0]
        iEqualityInfeasible = [i for i,veq in enumerate(equalityOuter) if abs(veq)>1]
        print("Algo> iInequalityInfeasible=",iInequalityInfeasible,"   iEqualityInfeasible=",iEqualityInfeasible)

        stepLength = self.maxStepLength

        iterOuter = 0
        while True:
            iterOuter = iterOuter+1
            self.iterOuter = iterOuter
            print("\n\nAlgo> iterOuter {}".format(iterOuter))

            # compute gradients on xOuter
            gradObjectiveOuter = self.__GetObjectiveGradient()
            for i in range(self.nc):
                gradInequalityOuter[i] = self.__GetInequalityConstraintGradient(i)
            for i in range(self.ne):
                gradEqualityOuter[i] = self.__GetEqualityConstraintGradient(i)

            self.__RevertPossibleShapeModificationsDuringAnalysis()

            # convert from value-gradient format to length-direction format
            _,eObjective,gradObjectiveL = self.__ValueFormatToLengthFormat(1,gradObjectiveOuter)
            for i in range(self.nc):
                lInequality0[i],eConstraint[i],gradInequalityL[i] = self.__ValueFormatToLengthFormat(constraintOuter[i],gradInequalityOuter[i])
            for i in range(self.ne):
                lEquality0[i],eEquality[i],gradEqualityL[i] = self.__ValueFormatToLengthFormat(equalityOuter[i],gradEqualityOuter[i])

            ## Initialization for the first inner iteration
            stepLength = min(self.stepLengthIncreaseFactor*stepLength,self.maxStepLength)

            if self.optimizeInitOfLInner: # make a guess on lConstraint and lEquality to reduce the number of inner iterations (might be counterproductive in some case, for instance: strain equality)
                for i in range(self.nc):
                    lConstraint[i] = max( lInequality0[i]+self.dlInequalityIncrRel*stepLength , lConstraint[i]-self.dlInequalityIncrRel*stepLength )
                for i in range(self.ne):
                    lEquality[i] = (1-self.lEqualityInitPersistence)*self.__InInterval(lEquality0[i],self.lEqualityInitInterval*stepLength)+self.lEqualityInitPersistence*lEquality[i]
            else: # safe assumptions
                lConstraint = lInequality0
                lEquality = lEquality0
            # set l for the unfeasible constraints and equalities
            for i in iInequalityInfeasible:
                lConstraint[i] = min( lInequality0[i] , self.dlUnfeasibleRel*stepLength )
            for i in iEqualityInfeasible:
                lEquality[i] = self.dlUnfeasibleRel*stepLength if equalityOuter[i]>0 else -self.dlUnfeasibleRel*stepLength #self.__InInterval(lEquality[i],self.dlUnfeasibleRel*stepLength)

            #######################################################################
            ## Inner iterations
            iterInner = 0
            while True:
                iterInner += 1
                self.iterCounter += 1
                print("\nAlgo> iterOuter {}: iterInner {}".format(iterOuter,iterInner))
                print("Algo> stepLength: {:.3g}".format(stepLength))

                lConstraintRel = [l/stepLength for l in lConstraint]
                lEqualityRel = [l/stepLength for l in lEquality]

                dxRel,lObjectiveRel,sInitRel,sDxRel,lXConstraintRel = plusSphereOld(eObjective,lConstraint,eConstraint,lEquality,eEquality,self)

                dx = [dxRel[i]*stepLength for i in range(self.n)]
                lObjective = lObjectiveRel*stepLength
                sInit = sInitRel*stepLength
                sDx = sDxRel*stepLength
                lXConstraint = [lXConstraintRel[i]*stepLength for i in range(self.nc)]

                # update mesh
                dxLastInner = [xOuter[i]+dx[i]-xInner[i] for i in range(self.n)]
                xInner = [xOuter[i]+dx[i] for i in range(self.n)]
                self.__WriteListOfValuesToNodalVariable(dxLastInner,SHAPE_UPDATE)
                self.ModelPartController.UpdateMeshAccordingInputVariable( SHAPE_UPDATE )
                self.ModelPartController.SetReferenceMeshToMesh()

                # logging
                dxAbsolute = [xInner[i] - self.xInit[i] for i in range(self.n)]
                self.__WriteListOfValuesToNodalVariable(dxAbsolute,SHAPE_CHANGE)
                self.DataLogger.DesignLogger.LogCurrentDesign(iterOuter*100+iterInner)
                responseLoggerIterationData = {"iterOuter": iterOuter, "iterInner": iterInner, "stepLength": stepLength, "sDx": sDx, "sInit": sInit, "lConstraint": lConstraint[:], "lEquality": lEquality[:], "gradConstraintL": gradInequalityL[:], "gradEqualityL": gradEqualityL[:], "gradObjectiveL": gradObjectiveL, "lObjective": lObjective, "objectiveOuter": objectiveOuter }

                # Check response values and update (stepLength,lConstraint,lEquality) if xInner is not feasible
                isFeasible = True
                if sInit>stepLength: # check if is inside ball when projecting
                    isFeasible = False
                    endIter = "sInit" if sInit>stepLength else "sDx"
                    stepLength = stepLength/self.stepLengthDecreaseFactor
                    for i in range(self.nc):
                        lConstraint[i] = min(lInequality0[i],-self.dlInequalityIncrRel*stepLength)
                    for i in range(self.ne):
                        lEquality[i] = self.__InInterval(lEquality0[i],self.dlEqualityIncrRel*stepLength)

                objectiveInner = None
                constraintInner = [None for i in range(self.nc)]
                equalityInner = [None for i in range(self.ne)]
                if isFeasible:
                    #remove the initially unfeasible constraints and inequalities
                    callOrder = self.callOrder[:]
                    for i in iInequalityInfeasible:
                        callOrder.remove(i+1)
                    for i in iEqualityInfeasible:
                        callOrder.remove(-i-1)
                    # evaluate response function following the callOrder
                    for iCall in callOrder:
                        if iCall>0: # constraint
                            i = iCall-1 # get the right index to access the lists
                            constraintInner[i] = self.__GetInequalityConstraint(i)
                            if constraintInner[i]>0:
                                isFeasible = False
                                endIter = "inequality"
                                dlConstraint = -constraintInner[i]/gradInequalityL[i] + self.dlInequalityIncrRel*stepLength
                                lConstraint[i] = lXConstraint[i]+dlConstraint
                                break
                        elif iCall<0: # equality
                            i = -iCall-1 # get the right index to access the lists
                            equalityInner[i] = self.__GetEqualityConstraint(i)
                            if abs(equalityInner[i])>1:
                                isFeasible = False
                                endIter = "equality"
                                dlEquality = self.__InInterval(-equalityInner[i]/gradEqualityL[i], self.dlEqualityIncrRel*stepLength)
                                lEquality[i] = lEquality[i]+dlEquality
                                break
                        elif iCall==0: # objective
                            objectiveInner = self.__GetObjective()
                            # may not enforce decreasing objective if there are initially unfeasible constraints, and if said so in the parameters
                            enforceDecreasingObjective = self.enforceDecreasingObjectiveIfUnfeasible or (len(iInequalityInfeasible)+len(iEqualityInfeasible)==0)
                            if objectiveInner>=objectiveOuter and enforceDecreasingObjective:
                                isFeasible = False
                                endIter = "objective"
                                stepLength = stepLength/self.stepLengthDecreaseFactor
                                for k in range(self.nc):
                                    lConstraint[k] = min(lInequality0[k],-self.dlInequalityIncrRel*stepLength)
                                for k in range(self.ne):
                                    lEquality[k] = self.__InInterval(lEquality0[k],self.dlEqualityIncrRel*stepLength)
                                break

                # End of inner iteration
                if isFeasible or iterInner==self.iterInnerMax or stepLength==self.stepLengthTermination: # end of outer iteration
                    endIter = "feasible" if isFeasible else ("innerMax" if iterInner==self.iterInnerMax else "s<termination" )
                    print("Algo> endIter: {}".format(endIter))
                    break
                else: # next inner iteration
                    print("Algo> endIter: {}".format(endIter))
                    responseLoggerIterationData.update( {"constraint": constraintInner, "equality": equalityInner, "objective": objectiveInner, "endIter": endIter, "equalityOuter": equalityOuter} )
                    self.DataLogger.ResponseLogger.LogIteration(responseLoggerIterationData,self)

            ## End of outer iteration: prepare next outer iteration

            # check if initially unfeasible constraints and equalities are now satisfied. If so, remove it from i*Unfeasible
            for i in iInequalityInfeasible:
                constraintInner[i] = self.__GetInequalityConstraint(i)
                if constraintInner[i]<=0:
                    iInequalityInfeasible.remove(i)
            for i in iEqualityInfeasible:
                equalityInner[i] = self.__GetEqualityConstraint(i)
                if abs(equalityInner[i])<=1:
                    iEqualityInfeasible.remove(i)

            # add unfullfilled constraint to infeasible (reached maximum nb of inner iteration)
            for i in range(self.nc):
                if constraintInner[i]>0 and i not in iInequalityInfeasible:
                    iInequalityInfeasible.append(i)
            for i in range(self.ne):
                if abs(equalityInner[i])>1 and i not in iEqualityInfeasible:
                    iEqualityInfeasible.append(i)
            print("Algo> iInequalityInfeasible=",iInequalityInfeasible,"   iEqualityInfeasible=",iEqualityInfeasible)

            xOuter = xInner[:]
            objectiveOuter = objectiveInner
            constraintOuter = constraintInner
            equalityOuter = equalityInner

            responseLoggerIterationData.update( {"constraint": constraintInner, "equality": equalityInner, "objective": objectiveInner} )

            if iterOuter==self.iterOuterMax or stepLength<self.stepLengthTermination:# or iterInner==self.iterInnerMax: # terminate algorithm
                responseLoggerIterationData["endIter"] = "outerMax" if iterOuter==self.iterOuterMax else ("s<termination" if stepLength<self.stepLengthTermination else "innerMax")
                self.DataLogger.ResponseLogger.LogIteration(responseLoggerIterationData,self)
                break
            else: # next outer iteration
                responseLoggerIterationData.update({"endIter": endIter, "iInequalityInfeasible": iInequalityInfeasible, "iEqualityInfeasible": iEqualityInfeasible})
                self.DataLogger.ResponseLogger.LogIteration(responseLoggerIterationData,self)


    # --------------------------------------------------------------------------
    def finalizeOptimizationLoop( self ):
        self.DataLogger.FinalizeDataLogging()
        self.Analyzer.finalizeAfterOptimizationLoop()


#################################################################################################
# Utilities

    # --------------------------------------------------------------------------
    def __GetObjective( self ):
        value = self.__GetValue( self.objectiveId )
        if self.objectiveInit is None:
            self.objectiveInit = value
            percentInit = 0
        else:
            percentInit = (self.objectiveInit-value)/self.objectiveInit*100
        print("Algo> value objective {} = {:.2g} -> {:.0f}%".format(self.objectiveId,value,percentInit))
        return value

    # --------------------------------------------------------------------------
    def __GetInequalityConstraint( self,i ):
        value = self.__GetValue( self.inequalityConstraintIds[i] )
        print("Algo> value inequality {} = {:.2g}".format(self.inequalityConstraintIds[i],value))
        return value

    # --------------------------------------------------------------------------
    def __GetEqualityConstraint( self,i ):
        value = self.__GetValue( self.equalityConstraintIds[i] )
        tol = self.__GetEqualityTolerance(i)
        ref = self.Communicator.getReferenceValue(self.equalityConstraintIds[i])
        abs_tolerance = abs(tol*ref)
        value = value/abs_tolerance
        print("Algo> value equality {} = {:.2g}".format(self.equalityConstraintIds[i],value))
        return value
    # --------------------------------------------------------------------------
    def __GetObjectiveGradient( self ):
        return self.__GetGradient(self.objectiveId,OBJECTIVE_SENSITIVITY)

    # --------------------------------------------------------------------------
    def __GetInequalityConstraintGradient( self,i ):
        return self.__GetGradient(self.inequalityConstraintIds[i],CONSTRAINT_SENSITIVITY)

    # --------------------------------------------------------------------------
    def __GetEqualityConstraintGradient( self,i ):
        gradient_as_list = self.__GetGradient(self.equalityConstraintIds[i],CONSTRAINT_SENSITIVITY)
        tol = self.__GetEqualityTolerance(i)
        ref = self.Communicator.getReferenceValue(self.equalityConstraintIds[i])
        abs_tolerance = abs(tol*ref)
        return [ u / abs_tolerance for u in gradient_as_list ]

    # --------------------------------------------------------------------------
    def __GetValue( self, functionName ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf( functionName )
        with contextlib.redirect_stdout(None):
            self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )
        value = self.Communicator.getStandardizedValue( functionName )
        return value

    # --------------------------------------------------------------------------
    def __GetGradient( self, functionName, nodal_variable ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestGradientOf( functionName )
        with contextlib.redirect_stdout(None):
            self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )

        gradient_as_array = self.Communicator.getStandardizedGradient( functionName )
        self.__StoreGradientOnNodalVariable( gradient_as_array, nodal_variable )
        # print("Algo.getGradient> gradient {}".format(functionName))
        if functionName == "strain_energy":
            self.GeometryUtilities.ComputeUnitSurfaceNormals()
            self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( nodal_variable )
            # print("Algo.getGradient> project on normals, gradient {}".format(functionName))
        return self.__ReadNodalVariableToList(nodal_variable)

    def __GetEqualityTolerance(self,i):
        tolset = self.OptimizationSettings["constraints"][i]["tolerance"]
        if tolset.IsDouble():
            return tolset.GetDouble()
        tolmat = tolset.GetMatrix()
        for i in range(tolmat.Size1()):
            tol = tolmat[i,0]
            untilIterOuter = tolmat[i,1]
            if self.iterOuter<= untilIterOuter or untilIterOuter==-1:
                return tol
        return tol

    # --------------------------------------------------------------------------
    def __RevertPossibleShapeModificationsDuringAnalysis( self ):
        self.ModelPartController.SetMeshToReferenceMesh()
        self.ModelPartController.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    # get length and directions
    def __ValueFormatToLengthFormat(self,value,gradient):
        if self.dampingIsSpecified:
            gradient = self.__DampVector(gradient)
        gradient = self.__MapVector(gradient)
        if self.dampingIsSpecified:
            gradient = self.__DampVector(gradient)
        ninf = norminf3d(gradient)
        eDir = [-gradient[i]/ninf for i in range(len(gradient))]
        gradL = sum(a*b for a,b in zip(gradient,eDir)) # dot product
        l = -value/ gradL
        return (l,eDir,gradL)

    # --------------------------------------------------------------------------
    # map a vector
    def __MapVector(self, gradient_list):
        # Here we use OBJECTIVE_SENSITIVITY as a temporary variable
        self.__WriteListOfValuesToNodalVariable( gradient_list, OBJECTIVE_SENSITIVITY )
        tmp_timer = timer_factory.CreateTimer()
        # print("Algo.mapVector> map...",end="\r")
        with contextlib.redirect_stdout(None):
            startTime = tmp_timer.StartTimer()
            self.Mapper.MapToDesignSpace( OBJECTIVE_SENSITIVITY, MAPPED_OBJECTIVE_SENSITIVITY )
            self.Mapper.MapToGeometrySpace( MAPPED_OBJECTIVE_SENSITIVITY, OBJECTIVE_SENSITIVITY )
        # print("Algo.mapVector> map in {:.2f} s".format(tmp_timer.GetTotalTime()))
        return self.__ReadNodalVariableToList( OBJECTIVE_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __DampVector(self, gradient_list):
        # Here we use OBJECTIVE_SENSITIVITY as a temporary variable
        self.__WriteListOfValuesToNodalVariable( gradient_list, OBJECTIVE_SENSITIVITY )
        tmp_timer = timer_factory.CreateTimer()
        # print("Algo.mapVector> damp...",end="\r")
        startTime = tmp_timer.StartTimer()
        with contextlib.redirect_stdout(None): # damp twice
            self.DampingUtilities.DampNodalVariable( OBJECTIVE_SENSITIVITY )
        # print("Algo.mapVector> damp in {:.2f} s".format(tmp_timer.GetTotalTime()))
        return self.__ReadNodalVariableToList( OBJECTIVE_SENSITIVITY )

    # --------------------------------------------------------------------------
    # build the call order, 0 is objective, [1...] constraints shifted by 1, [... -1] equality shifted by 1 and take opposite
    def __GetCallOrder(self,timeObjective,timeConstraint,timeEquality):
        timeResp = [(timeObjective,0)] + [(t,i+1) for i,t in enumerate(timeConstraint)] + [(t,-i-1) for i,t in enumerate(timeEquality)]
        sortedTimeResp = sorted(timeResp, key=lambda x: x[0])
        return [r for t,r in sortedTimeResp]

    # --------------------------------------------------------------------------
    def __InInterval( self, value, bound):
        return min(max(value,-bound),bound)


##########################################################################################
# Kratos / Python IO

    # --------------------------------------------------------------------------
    def __StoreGradientOnNodalVariable( self, gradient, variable_name ):
        for nodeId, tmp_gradient in gradient.items():
            self.OptimizationModelPart.Nodes[nodeId].SetSolutionStepValue(variable_name,0,tmp_gradient)

    # --------------------------------------------------------------------------
    def __WriteListOfValuesToNodalVariable(self, list_of_values, nodal_variable):
        k = 0
        for node in self.DesignSurface.Nodes:
            tmp_values = list_of_values[k:k+3]
            k = k+3
            node.SetSolutionStepValue(nodal_variable,tmp_values)

    # --------------------------------------------------------------------------
    def __ReadNodalVariableToList(self, nodal_variable): # e.g.: OBJECTIVE_SENSITIVITY
        list_of_values = []
        for node in self.DesignSurface.Nodes:
            tmp_values = node.GetSolutionStepValue(nodal_variable)
            list_of_values.append(tmp_values[0])
            list_of_values.append(tmp_values[1])
            list_of_values.append(tmp_values[2])
        return list_of_values

    def __ReadDesignSurfaceToList(self):
        list_of_values = []
        for node in self.DesignSurface.Nodes:
            list_of_values.append(node.X)
            list_of_values.append(node.Y)
            list_of_values.append(node.Z)
        return list_of_values