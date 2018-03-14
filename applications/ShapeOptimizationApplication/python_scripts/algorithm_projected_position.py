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

from projected_position_modules.matrix import norminf3d
from projected_position_modules.projection import plusBall

# ==============================================================================
class AlgorithmProjectedPosition( OptimizationAlgorithm ) :

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

        self.projectionOnNormalsIsSpecified = OptimizationSettings["optimization_algorithm"]["project_gradients_on_surface_normals"].GetBool()
        self.dampingIsSpecified = OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()

        self.GeometryUtilities = GeometryUtilities( self.DesignSurface )
        self.OptimizationUtilities = OptimizationUtilities( self.DesignSurface, OptimizationSettings )
        if self.dampingIsSpecified:
            damping_regions = self.ModelPartController.GetDampingRegions()
            self.DampingUtilities = DampingUtilities( self.DesignSurface, damping_regions, self.OptimizationSettings )

        # Parameters for the algorithm
        self.iterOuterMax = OptimizationSettings["optimization_algorithm"]["max_outer_iterations"].GetInt()
        self.iterInnerMax = OptimizationSettings["optimization_algorithm"]["max_inner_iterations"].GetInt()
        self.maxStepLength = OptimizationSettings["optimization_algorithm"]["max_step_length"].GetDouble()

        self.dlInequalityIncrRel = OptimizationSettings["optimization_algorithm"]["delta_l_inequality"].GetDouble()
        self.dlEqualityIncrRel = OptimizationSettings["optimization_algorithm"]["delta_l_equality"].GetDouble()
        self.dlUnfeasibleRel = OptimizationSettings["optimization_algorithm"]["delta_l_infeasible"].GetDouble()

        self.enforceDecreasingObjectiveIfUnfeasible = OptimizationSettings["optimization_algorithm"]["decreasing_objective_if_infeasible"].GetBool()
        self.optimizeInitOfLInner = OptimizationSettings["optimization_algorithm"]["efficient_initialization_of_l"].GetBool()

        self.onlyObjectiveId = OptimizationSettings["objectives"][0]["identifier"].GetString()
        self.inequalityConstraintIds = []
        self.constraintTypes = []
        self.constraintReferences = []
        self.equalityConstraintIds = []
        self.equalityTolerances = []

        # read the constraint and equality response function in the params
        for i in range(self.OptimizationSettings["constraints"].size()):
            cname = OptimizationSettings["constraints"][i]["identifier"].GetString()
            ctype = OptimizationSettings["constraints"][i]["type"].GetString()
            if ctype == "=":
                self.equalityConstraintIds.append(cname)
                self.equalityTolerances.append(OptimizationSettings["constraints"][i]["tolerance"].GetDouble())
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

    # --------------------------------------------------------------------------
    def initializeOptimizationLoop( self ):
        self.Analyzer.initializeBeforeOptimizationLoop()
        self.ModelPartController.InitializeMeshController()
        self.DataLogger.StartTimer()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def runOptimizationLoop( self ):

        print("solve__ppa> begin")
        self.iterCounter = 0 # cumulative number of inner iterations

        xOuter = []
        for node in self.DesignSurface.Nodes:
            xOuter.append(node.X)
            xOuter.append(node.Y)
            xOuter.append(node.Z)
        xInner = xOuter[:]
        xInit = xOuter[:]

        # init variables
        gradInequalityOuter = [0 for _ in range(self.nc)] # item will be replaced by lists
        gradEqualityOuter = [0 for _ in range(self.ne)] # item will be replaced by lists
        eConstraint = [0 for _ in range(self.nc)] # item will be replaced by lists
        eEquality = [0 for _ in range(self.ne)] # item will be replaced by lists
        gradConstraintL = [0 for _ in range(self.nc)]
        gradEqualityL = [0 for _ in range(self.ne)]
        inequalityOuter = [0 for _ in range(self.nc)]
        equalityOuter = [0 for _ in range(self.ne)]
        lConstraint = [0 for _ in range(self.nc)]
        lEquality = [0 for _ in range(self.ne)]
        lInequality0 = [0 for _ in range(self.nc)]
        lEquality0 = [0 for _ in range(self.ne)]

        # get values at xInit
        timer = timer_factory.CreateTimer()
        timeConstraint = [0 for _ in range(self.nc)]
        timeEquality = [0 for _ in range(self.ne)]

        # Compute response for initial design
        timer.StartNewLap()
        objectiveOuter = self.__GetObjective()
        print("Objective_value: ", objectiveOuter)
        timeObjective = timer.GetLapTime()
        for i in range(self.nc):
            timer.StartNewLap()
            inequalityOuter[i] = self.__GetInequalityConstraint()
            print("Inequality_value: ", inequalityOuter[i])
            timeConstraint[i] = timer.GetLapTime()
        for i in range(self.ne):
            timer.StartNewLap()
            equalityOuter[i] = self.__GetEqualityConstraint()
            print("equality_value: ", equalityOuter[i])
            timeEquality[i] = timer.GetLapTime()

        self.callOrder = self.__GetCallOrder(timeObjective,timeConstraint,timeEquality) # build callOrder by increasing computation time
        print("solve_ppa> callOrder = ",self.callOrder)

        # self.DataLogger.ResponseLogger.LogIteration({"iterOuter": 0, "inequality":inequalityOuter,"equality": equalityOuter, "objective": objectiveOuter},p)

        # get the constraints and equality violated by xInit
        iInequalityInfeasible = [i for i,vconstr in enumerate(inequalityOuter) if vconstr>0]
        iEqualityInfeasible = [i for i,veq in enumerate(equalityOuter) if abs(veq)>1]
        print("solve_ppa> iInequalityInfeasible=",iInequalityInfeasible,"   iEqualityInfeasible=",iEqualityInfeasible)

        stepLength = self.maxStepLength


        iterOuter = 0
        while True:

            iterOuter = iterOuter+1

            print("\n>===================================================================")
            print("> ",self.DataLogger.GetTimeStamp(),": Starting outer iteration ", iterOuter)
            print(">===================================================================\n")

            # compute gradients on xOuter
            gradObjectiveOuter = self.__GetObjectiveGradient()
            for i in range(self.nc):
                gradInequalityOuter[i] = self.__GetInequalityConstraintGradient()
            for i in range(self.ne):
                gradEqualityOuter[i] = self.__GetEqualityConstraintGradient()

            self.__RevertPossibleShapeModificationsDuringAnalysis()

            self.GeometryUtilities.ComputeUnitSurfaceNormals()
            # self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( OBJECTIVE_SENSITIVITY )
            # self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( CONSTRAINT_SENSITIVITY )

            # gradObjectiveOuter = []
            # gradInequalityOuter = [[]]
            # self.__ReadNodalVariableToList( OBJECTIVE_SENSITIVITY, gradObjectiveOuter )
            # self.__ReadNodalVariableToList( CONSTRAINT_SENSITIVITY, gradInequalityOuter[0] )







            # convert from value-gradient format to length-direction format
            _,eObjective,_ = self.__valueFormatToLengthFormat(1,gradObjectiveOuter)
            for i in range(self.nc):
                lInequality0[i],eConstraint[i],gradConstraintL[i] = self.__valueFormatToLengthFormat(constraintOuter[i],gradInequalityOuter[i])
            for i in range(self.ne):
                lEquality0[i],eEquality[i],gradEqualityL[i] = self.__valueFormatToLengthFormat(equalityOuter[i],gradEqualityOuter[i])

            ## Initialization for the first inner iteration
            stepLength = min(1.2*stepLength,self.maxStepLength)

            if self.optimizeInitOfLInner: # make a guess on lConstraint and lEquality to reduce the number of inner iterations (might be counterproductive in some case, for instance: strain equality)
                for i in range(self.nc):
                    lConstraint[i] = max( lInequality0[i]+self.dlInequalityIncrRel*stepLength , lConstraint[i]-self.dlInequalityIncrRel*stepLength )
                for i in range(self.ne):
                    lEquality[i] = 0.3*self.__inInterval(lEquality0[i],0.5*stepLength)+0.7*lEquality[i]
            else: # safe assumptions
                lConstraint = lInequality0
                lEquality = lEquality0
            # set l for the unfeasible constraints and equalities
            for i in iInequalityInfeasible:
                lConstraint[i] = min( lInequality0[i] , self.dlUnfeasibleRel*stepLength )
            for i in iEqualityInfeasible:
                lEquality[i] = self.__inInterval(lEquality[i],self.dlUnfeasibleRel*stepLength)

            #######################################################################
            ## Inner iterations
            iterInner = 0
            while True:
                iterInner += 1
                self.iterCounter += 1

                print("\n>===================================================================")
                print("> ",self.DataLogger.GetTimeStamp(),": Starting inner iteration ", iterInner)
                print(">===================================================================\n")

                print("\nsolve_ppa> iterOuter {}: iterInner {}".format(iterOuter,iterInner))

                # convert (stepLength,lConstraint,lEquality) into delta_x
                [dx,lObjective,sInit,sDx,lXConstraint] = plusBall(stepLength,eObjective,lConstraint,eConstraint,lEquality,eEquality)

                # for itr in range(len(dx)):
                #     dx[itr] = 1.0


                dxLastInner = [xOuter[i]+dx[i]-xInner[i] for i in range(self.n)]
                xInner = [xOuter[i]+dx[i] for i in range(self.n)]
                self.__writeListOfValuesToNodalVariable(dxLastInner,SHAPE_UPDATE)

                self.ModelPartController.UpdateMeshAccordingInputVariable( SHAPE_UPDATE )
                self.ModelPartController.SetReferenceMeshToMesh()


                dxAbsolute = [xInner[i] - xInit[i] for i in range(self.n)]
                self.__writeListOfValuesToNodalVariable(dxAbsolute,SHAPE_CHANGE)
                self.DataLogger.DesignLogger.LogCurrentDesign((iterOuter*100+iterInner))


                # # Generate new shape
                # print(self.DesignSurface.Nodes[113].GetSolutionStepValue(SHAPE_UPDATE))
                # self.__writeDeltaValuesToNodalVariable(dx, SHAPE_UPDATE)
                # print(self.DesignSurface.Nodes[113].GetSolutionStepValue(SHAPE_UPDATE))
                # self.ModelPartController.UpdateReferenceMeshAccordingInputVariable( SHAPE_UPDATE )
                # # self.ModelPartController.SetReferenceMeshToMesh()
                # # self.__determineAbsoluteChanges()

                # for itr in range(len(dx)):
                #     dx[itr] = 1.0


                # self.updateReferenceMesh(dxLastInner)
                # self.__writeValuesToNodalVariable(dxLastInner, SHAPE_UPDATE)
                # self.ModelPartController.UpdateReferenceMeshAccordingInputVariable( SHAPE_UPDATE )
                # self.DataLogger.DesignLogger.LogCurrentDesign(self.iterCounter)
                # responseLoggerIterationData = {"iterOuter": iterOuter, "iterInner": iterInner, "stepLength": stepLength, "sDx": sDx, "sInit": sInit, "lConstraint": lConstraint[:], "lEquality": lEquality[:] }

                ## Check response values and update (stepLength,lConstraint,lEquality) if xInner is not feasible
                isFeasible = True
                if sInit>stepLength or sDx>stepLength: # check if is inside ball when projecting, or if projBall could find precisely the right step length
                    isFeasible = False
                    endIter = "sInit" if sInit>stepLength else "sDx"
                    stepLength = stepLength/2
                    for i in range(self.nc):
                        lConstraint[i] = min(lInequality0[i],-self.dlInequalityIncrRel*stepLength)
                    for i in range(self.ne):
                        lEquality[i] = self.__inInterval(lEquality0[i],self.dlEqualityIncrRel*stepLength)

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
                    for i in callOrder:
                        if i>0: # constraint
                            i = i-1 # get the right index to access the lists
                            constraintInner[i] = self.__GetInequalityConstraint()
                            print("inequality_value: ", constraintInner[i])
                            if constraintInner[i]>0:
                                isFeasible = False
                                endIter = "inequality"
                                dlConstraint = -constraintInner[i]/gradConstraintL[i] + self.dlInequalityIncrRel*stepLength
                                lConstraint[i] = lXConstraint[i]+dlConstraint
                                break
                        elif i<0: # equality
                            i = -i-1 # get the right index to access the lists
                            equalityInner[i] = self.__GetEqualityConstraint()
                            print("equality_value: ", equalityInner[i])
                            if abs(equalityInner[i])>1:
                                isFeasible = False
                                endIter = "equality"
                                # dlEquality = self.__inInterval(-equalityInner[i]/gradEqualityL[i],self.dlEqualityIncrRel*stepLength)
                                dlEquality = self.__inInterval(-equalityInner[i]/gradEqualityL[i], 0.2*stepLength)
                                lEquality[i] = lEquality[i]+dlEquality
                                break
                        elif i==0: # objective
                            objectiveInner = self.__GetObjective()
                            print("Objective_value: ", objectiveInner)
                            # may not enforce decreasing objective if there are initially unfeasible constraints, and if said so in the parameters
                            enforceDecreasingObjective = self.enforceDecreasingObjectiveIfUnfeasible or (len(iInequalityInfeasible)+len(iEqualityInfeasible)==0)
                            if objectiveInner>=objectiveOuter and enforceDecreasingObjective:
                                isFeasible = False
                                endIter = "objective"
                                stepLength = stepLength/2
                                for i in range(self.nc):
                                    lConstraint[i] = min(lInequality0[i],-self.dlInequalityIncrRel*stepLength)
                                for i in range(self.ne):
                                    lEquality[i] = self.__inInterval(lEquality0[i],self.dlEqualityIncrRel*stepLength)
                                break
                # End of inner iteration
                if isFeasible or iterInner==self.iterInnerMax or stepLength==1e-5:
                    endIter = "feasible" if isFeasible else ("innerMax" if iterInner==self.iterInnerMax else "s<1e-5" )
                    print("solve_ppa> endIter: {}".format(endIter))
                    break # end of an outer iteration
                else:
                    print("solve_ppa> endIter: {}".format(endIter))
                    # responseLoggerIterationData.update( {"inequality": constraintInner, "equality": equalityInner, "objective": objectiveInner, "endIter": endIter} )
                    # self.DataLogger.ResponseLogger.LogIteration(responseLoggerIterationData,p)

            ## End of outer iteration: prepare next outer iteration
            # check if initially unfeasible constraints and equalities are now satisfied. If so, remove it from i*Unfeasible
            for i in iInequalityInfeasible:
                constraintInner[i] = self.getConstraint(i,False)
                if constraintInner[i]<=0:
                    iInequalityInfeasible.remove(i)
            for i in iEqualityInfeasible:
                equalityInner[i] = self.getEquality(i,False)
                if abs(equalityInner[i])<=1:
                    iEqualityInfeasible.remove(i)

            xOuter = xInner[:]
            # self.ModelPartController.SetReferenceMeshToMesh()
            # VariableUtils().SetToZero_VectorVar(SHAPE_UPDATE,self.DesignSurface.Nodes)

            # for node in self.DesignSurface.Nodes:
            #     zero_update = Vector(3)


            # print(self.DesignSurface.Nodes[1].X0)
            # print(self.DesignSurface.Nodes[1].X)
            # input()

            objectiveOuter = objectiveInner
            constraintOuter = constraintInner
            equalityOuter = equalityInner

            # responseLoggerIterationData.update( {"inequality": constraintInner, "equality": equalityInner, "objective": objectiveInner} )

            if iterOuter==self.iterOuterMax or stepLength<1e-5 or iterInner==self.iterInnerMax:
                # responseLoggerIterationData["endIter"] = "outerMax" if iterOuter==self.iterOuterMax else ("s<1e-5" if stepLength<1e-5 else "innerMax")
                # self.DataLogger.ResponseLogger.LogIteration(responseLoggerIterationData,p)
                break # terminate algorithm
            # else:
                # responseLoggerIterationData["endIter"] = endIter
                # self.DataLogger.ResponseLogger.LogIteration(responseLoggerIterationData,p)






















        # for self.optimizationIteration in range(1,self.iterOuterMax):
        #     print("\n>===================================================================")
        #     print("> ",self.DataLogger.GetTimeStamp(),": Starting optimization iteration ", self.optimizationIteration)
        #     print(">===================================================================\n")

        #     self.__initializeNewShape()

        #     self.__analyzeShape()

        #     if self.projectionOnNormalsIsSpecified:
        #         self.__projectSensitivitiesOnSurfaceNormals()

        #     if self.dampingIsSpecified:
        #         self.__dampSensitivities()

        #     self.__computeShapeUpdate()

        #     if self.dampingIsSpecified:
        #         self.__dampShapeUpdate()

        #     self.__logCurrentOptimizationStep()

        #     self.__timeOptimizationStep()

        #     if self.__isAlgorithmConverged():
        #         break
        #     else:
        #         self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def finalizeOptimizationLoop( self ):
        self.DataLogger.FinalizeDataLogging()
        self.Analyzer.finalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __GetObjective( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf( self.onlyObjectiveId )
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )
        return self.Communicator.getStandardizedValue( self.onlyObjectiveId )

    # --------------------------------------------------------------------------
    def __GetInequalityConstraint( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf( self.inequalityConstraintIds[0] )
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )
        return self.Communicator.getStandardizedValue( self.inequalityConstraintIds[0] )

    # --------------------------------------------------------------------------
    def __GetEqualityConstraint( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf( self.equalityConstraintIds[0] )
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )
        return self.Communicator.getStandardizedValue( self.equalityConstraintIds[0] ) / (self.equalityTolerances[0] * self.Communicator.getReferenceValue(self.equalityConstraintIds[0]))

    # --------------------------------------------------------------------------
    def __GetObjectiveGradient( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestGradientOf( self.onlyObjectiveId )
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )

        gradientOfObjectiveFunction = self.Communicator.getStandardizedGradient( self.onlyObjectiveId )
        self.__storeGradientOnNodalVariable( gradientOfObjectiveFunction, OBJECTIVE_SENSITIVITY )

        gradient_as_list = []
        self.__SortVariableValuesIntoList( OBJECTIVE_SENSITIVITY, gradient_as_list )
        return gradient_as_list

    # --------------------------------------------------------------------------
    def __GetInequalityConstraintGradient( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestGradientOf( self.inequalityConstraintIds[0] )
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )

        gradientOfInequalityConstraint = self.Communicator.getStandardizedGradient( self.inequalityConstraintIds[0] )
        self.__storeGradientOnNodalVariable( gradientOfInequalityConstraint, CONSTRAINT_SENSITIVITY )

        gradient_as_list = []
        self.__SortVariableValuesIntoList( CONSTRAINT_SENSITIVITY, gradient_as_list )
        return gradient_as_list

    # --------------------------------------------------------------------------
    def __GetEqualityConstraintGradient( self ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestGradientOf( self.equalityConstraintIds[0] )
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )

        gradientOfEqualityConstraint = self.Communicator.getStandardizedGradient( self.equalityConstraintIds[0] )
        self.__storeGradientOnNodalVariable( gradientOfEqualityConstraint, CONSTRAINT_SENSITIVITY )

        gradient_as_list = []
        self.__SortVariableValuesIntoList( CONSTRAINT_SENSITIVITY, gradient_as_list )
        return [u/ (self.equalityTolerances[0] * self.Communicator.getReferenceValue(self.equalityConstraintIds[0])) for u in gradient_as_list]

    # --------------------------------------------------------------------------
    def __RevertPossibleShapeModificationsDuringAnalysis( self ):
        self.ModelPartController.SetMeshToReferenceMesh()
        self.ModelPartController.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def __projectSensitivitiesOnSurfaceNormals( self ):
        self.GeometryUtilities.ComputeUnitSurfaceNormals()
        self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( OBJECTIVE_SENSITIVITY )
        self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    # get length and directions
    def __valueFormatToLengthFormat(self,value,gradient):
        gradientMapped = self.__mapVector(gradient)
        ninf = norminf3d(gradientMapped)
        eDir = [-gradientMapped[i]/ninf for i in range(len(gradientMapped))]
        gradL = sum(a*b for a,b in zip(gradientMapped,eDir)) # dot product
        l = -value/ gradL
        return (l,eDir,gradL)

    # --------------------------------------------------------------------------
    # map a vector
    def __mapVector(self, gradient_list):
        self.__writeListOfValuesToNodalVariable( gradient_list, OBJECTIVE_SENSITIVITY )

        tmp_timer = timer_factory.CreateTimer()

        print("Algo.mapVector> map...",end="\r")
        startTime = tmp_timer.StartTimer()
        self.Mapper.MapToDesignSpace( OBJECTIVE_SENSITIVITY, MAPPED_OBJECTIVE_SENSITIVITY )
        self.Mapper.MapToGeometrySpace( MAPPED_OBJECTIVE_SENSITIVITY, OBJECTIVE_SENSITIVITY )
        print("Algo.mapVector> map in {:.2f} s".format(tmp_timer.GetTotalTime()))

        list_of_mapped_gradients = []
        self.__SortVariableValuesIntoList( OBJECTIVE_SENSITIVITY, list_of_mapped_gradients )
        return list_of_mapped_gradients

    # --------------------------------------------------------------------------
    def __storeGradientOnNodalVariable( self, gradient, variable_name ):
        for nodeId, tmp_gradient in gradient.items():
            self.OptimizationModelPart.Nodes[nodeId].SetSolutionStepValue(variable_name,0,tmp_gradient)

    # --------------------------------------------------------------------------
    def __SortVariableValuesIntoList( self, variable_name, list_of_values ):
        for node in self.DesignSurface.Nodes:
            variable_value = node.GetSolutionStepValue(variable_name)
            list_of_values.append(variable_value[0])
            list_of_values.append(variable_value[1])
            list_of_values.append(variable_value[2])

    # --------------------------------------------------------------------------
    # build the call order, 0 is objective, [1...] constraints shifted by 1, [... -1] equality shifted by 1 and take opposite
    def __GetCallOrder(self,timeObjective,timeConstraint,timeEquality):
        timeResp = [(timeObjective,0)] + [(t,i+1) for i,t in enumerate(timeConstraint)] + [(t,-i-1) for i,t in enumerate(timeEquality)]
        sortedTimeResp = sorted(timeResp, key=lambda x: x[0])
        return [r for t,r in sortedTimeResp]

    # --------------------------------------------------------------------------
    def __writeListOfValuesToNodalVariable(self, list_of_values, nodal_variable):
        k = 0
        for node in self.DesignSurface.Nodes:
            tmp_values = list_of_values[k:k+3]
            k = k+3
            node.SetSolutionStepValue(nodal_variable,tmp_values)

    # --------------------------------------------------------------------------
    def __writeValuesToNodalVariable(self, list_of_values, nodal_variable):
        k = 0
        for node in self.DesignSurface.Nodes:
            tmp_values = list_of_values[k:k+3]
            # tmp_values[0] -= node.GetSolutionStepValue(SHAPE_UPDATE_X)
            # tmp_values[1] -= node.GetSolutionStepValue(SHAPE_UPDATE_Y)
            # tmp_values[2] -= node.GetSolutionStepValue(SHAPE_UPDATE_Z)
            k = k+3
            node.SetSolutionStepValue(nodal_variable,tmp_values)

    # --------------------------------------------------------------------------
    def __ReadNodalVariableToList(self, nodal_variable, list_of_values): # e.g.: OBJECTIVE_SENSITIVITY
        for node in self.DesignSurface.Nodes:
            tmp_values = node.GetSolutionStepValue(nodal_variable)
            list_of_values.append(tmp_values[0])
            list_of_values.append(tmp_values[1])
            list_of_values.append(tmp_values[2])

    # --------------------------------------------------------------------------
    def __inInterval( self, value, bound):
        return min(max(value,-bound),bound)

    # # --------------------------------------------------------------------------
    # def __initializeNewShape( self ):
    #     self.ModelPartController.UpdateMeshAccordingInputVariable( SHAPE_UPDATE )
    #     self.ModelPartController.SetReferenceMeshToMesh()

    # # --------------------------------------------------------------------------
    # def __dampSensitivities( self ):
    #     self.DampingUtilities.DampNodalVariable( OBJECTIVE_SENSITIVITY )
    #     self.DampingUtilities.DampNodalVariable( CONSTRAINT_SENSITIVITY )

    # # --------------------------------------------------------------------------
    # def __dampShapeUpdate( self ):
    #     self.DampingUtilities.DampNodalVariable( SHAPE_UPDATE )

    # # --------------------------------------------------------------------------
    # def __logCurrentOptimizationStep( self ):
    #     self.DataLogger.LogCurrentData( self.optimizationIteration )

    # # --------------------------------------------------------------------------
    # def __timeOptimizationStep( self ):
    #     print("\n> Time needed for current optimization step = ", self.DataLogger.GetLapTime(), "s")
    #     print("> Time needed for total optimization so far = ", self.DataLogger.GetTotalTime(), "s")

    # # --------------------------------------------------------------------------
    # def __isAlgorithmConverged( self ):

    #     if self.optimizationIteration > 1 :

    #         # Check if maximum iterations were reached
    #         if self.optimizationIteration == self.maxIterations:
    #             print("\n> Maximal iterations of optimization problem reached!")
    #             return True

    #         relativeChangeOfObjectiveValue = self.DataLogger.GetValue( "RELATIVE_CHANGE_OF_OBJECTIVE_VALUE" )

    #         # Check for relative tolerance
    #         relativeTolerance = self.OptimizationSettings["optimization_algorithm"]["relative_tolerance"].GetDouble()
    #         if abs(relativeChangeOfObjectiveValue) < relativeTolerance:
    #             print("\n> Optimization problem converged within a relative objective tolerance of ",relativeTolerance,"%.")
    #             return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges( self ):
        # self.OptimizationUtilities.AddFirstVariableToSecondVariable( CONTROL_POINT_UPDATE, CONTROL_POINT_CHANGE )
        self.OptimizationUtilities.AddFirstVariableToSecondVariable( SHAPE_UPDATE, SHAPE_CHANGE )

# ==============================================================================
