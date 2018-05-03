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
import data_logger_factory

# Additional imports
import timer_factory as timer_factory
import logging
import contextlib
from copy import deepcopy
import datetime
import os.path

import projected_position_modules.utils as utils
from projected_position_modules.matrix import norminf3d
from projected_position_modules.projection import stepDirectionRule, printListValues

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

        #####################################################################################
        # Utilities

        self.ModelPartController = ModelPartController
        self.Analyzer = Analyzer
        self.Communicator = Communicator
        self.Mapper = Mapper
        self.DataLogger = DataLogger
        self.OptimizationSettings = OptimizationSettings
        self.OptimizationModelPart = ModelPartController.GetOptimizationModelPart()
        self.DesignSurface = ModelPartController.GetDesignSurface()
        self.GeometryUtilities = GeometryUtilities( self.DesignSurface )
        self.OptimizationUtilities = OptimizationUtilities( self.DesignSurface, OptimizationSettings )
        self.dampingIsSpecified = OptimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool()
        if OptimizationSettings["design_variables"]["damping"].Has("damp_only_after_mapping"):
            self.dampOnlyAfterMapping = OptimizationSettings["design_variables"]["damping"]["damp_only_after_mapping"].GetBool()
        else:
            self.dampOnlyAfterMapping = False
        if self.dampingIsSpecified:
            damping_regions = self.ModelPartController.GetDampingRegions()
            self.DampingUtilities = DampingUtilities( self.DesignSurface, damping_regions, self.OptimizationSettings )

        Communicator.p = self # to access p from communicator

        ######################################################################################
        # Algorithm parameters

        algoParamTemplate = [ # attributeName, jsonName, jsontype, is_optional, default_value
            ["maxStepLength","max_step_length","double",False],
            ["sphereNorm","sphere_norm","string",True,"geom_norminf3d"],  # "geom_norminf3d" or "coord_norm1"
            ["sphereSum","sphere_sum","string",True,"projection"], # "addition" or "projection"
            ["stepLengthTermination","step_length_termination","double",True,0],
            ["stepLengthDecreaseFactor","step_length_decrease_factor","double",True,2],
            ["stepLengthIncreaseFactor","step_length_increase_factor","double",True,1.2],
            ["enforceFinalFeasibilityEveryIter","enforce_final_feasibility_every_iterations","int",True,-1], # default deactivate
            ["iterMax","max_iterations","int",True,-1], # deactivate
            ["conjugateGradient","conjugate_gradient","bool",True,False] # for rosenbrock
        ]
        utils.ReadAlgoParameters(self,algoParamTemplate,OptimizationSettings["optimization_algorithm"])

        #############################################################
        # Response functions parameters

        objectiveParamTemplate = [ # propertyName, jsonName, jsontype, optional, default_value
            ["objectiveId","identifier","string",False],
            ["objectiveMinShare","minimum_share","double",True,0], # deactivate
            ["objectiveProjectNormals","project_gradients_on_surface_normals","bool",True,False]
        ]
        inequalityParamTemplate = [
            ["inequalityIds","identifier","string",False],
            ["inequalityTypes","type","string",False],
            ["inequalityReferences","reference",["double","string"],False],
            ["inequalityMaxShare","maximum_share","double",True,9999], # deactivate
            ["inequalityProjectNormals","project_gradients_on_surface_normals","bool",True,False]
        ]
        equalityParamTemplate = [
            ["equalityIds","identifier","string",False],
            ["equalityReferences","reference",["double","string"],False],
            ["equalityMaxShare","maximum_share","double",True,9999], # deactivate
            ["equalityFinalTolerances","final_tolerance","double",True,None],
            ["equalityProjectNormals","project_gradients_on_surface_normals","bool",True,False]
        ]
        utils.ReadResponseParameters(self,objectiveParamTemplate,inequalityParamTemplate,equalityParamTemplate,OptimizationSettings)
        self.projectNormalsIds = utils.BuildProjectNormalsIds(self)

        #####################################################################################
        # Calculated parameters and initialize variables

        self.runId = utils.BuildRunId(self)
        self.resultsFolder = "../results/"+self.runId
        OptimizationSettings.resultsFolder = self.resultsFolder
        OptimizationSettings.runId = self.runId
        self.DataLogger =  data_logger_factory.CreateDataLogger( ModelPartController, Communicator, OptimizationSettings ) # recreate datalogger to have results folder named with runid

        self.n = 3*self.DesignSurface.NumberOfNodes()
        self.ni = len(self.inequalityIds)
        self.ne = len(self.equalityIds)

        self.isRosenbrock = "rosenbrock" in self.objectiveId
        if self.isRosenbrock:
            xi = OptimizationSettings["optimization_algorithm"]["x_init_rosenbrock"].GetVector()
            self.xInit = [xi[0],xi[1]]
            self.xRosenbrock = deepcopy(self.xInit) # to be accessible by the analyser
            self.n = 2
        else:
            self.xInit = self.__ReadDesignSurfaceToList()

        self.isEnforcingFinalFeasibility = False
        self.isTerminationPhase = False
        self.enforceFinalFeasibilityNextIter = self.enforceFinalFeasibilityEveryIter
        self.lastIterUpdateStepLength = 0

        ##################################
        # for logging
        self.__LogParam() # fill self.param
        self.nodeIds, self.idIndex = utils.BuildNodeIdIndex(self.DesignSurface)
        self.r = {} # store results
        if self.isRosenbrock:
            self.param["rosenbrockFactor"] = OptimizationSettings["objectives"][0]["kratos_response_settings"]["factor_b"].GetDouble()
            self.param["rosenbrockInit"] = self.xInit

    # --------------------------------------------------------------------------
    def initializeOptimizationLoop( self ):
        self.Analyzer.initializeBeforeOptimizationLoop()
        self.ModelPartController.InitializeMeshController()
        self.DataLogger.StartTimer()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def runOptimizationLoop( self ):
        # init variables
        stepLength = self.maxStepLength
        self.stepLength = stepLength
        lastIterUpdateStepLength = 0
        iterCounter = 0
        while True:
            iterCounter += 1
            self.iterCounter = iterCounter
            print("\nAlgo> iter {}".format(iterCounter))
            print("Algo> stepLength: {:.3g}".format(stepLength))

            # compute values and gradients and convert to length-directions format
            objective,inequality,equality = self.__GetValues()
            gradObjective,gradInequality,gradEquality = self.__GetGradients()
            eObjective,lInequality,eInequality,lEquality,eEquality = self.__GetLengthsDirections(gradObjective,inequality,gradInequality,equality,gradEquality)

            # update step length
            stepLength = self.__StepLengthRule(stepLength,objective,lInequality,lEquality)

            lBarInequality = [l/stepLength for l in lInequality]
            lBarEquality = [l/stepLength for l in lEquality]

            # manage final feasibility check before calculating next step
            isTerminated = self.__CheckTerminationAndFinalFeasibility(iterCounter,stepLength,inequality,equality)
            if isTerminated: # termination is here, not at the end of the loop (in theory, should be just after GetValues())
                self.__LogVariables(locals())
                break

            # convert (stepLength,lInequality,lEquality) into delta_x
            deltaXBar,lambdaBarObjective,lambdaBarInequality,lambdaBarEquality,sInit,sDx = stepDirectionRule(eObjective,lBarInequality,eInequality,lBarEquality,eEquality,self)
            deltaX = [stepLength*deltaXBar[i] for i in range(self.n)]

            # update mesh
            self.__UpdateMesh(deltaX)

            # log all the variables in the scope except the too long lists (eObjective, dx, etc)
            self.__LogVariables(locals())


    # --------------------------------------------------------------------------
    def finalizeOptimizationLoop( self ):
        self.DataLogger.FinalizeDataLogging()
        self.Analyzer.finalizeAfterOptimizationLoop()


#################################################################################################
# Algorithm steps

    def __GetValues(self):
        objective = self.__GetObjectiveValue()
        inequality = [0 for _ in range(self.ni)]
        equality = [0 for _ in range(self.ne)]
        for i in range(self.ni):
            inequality[i] = self.__GetInequalityValue(i)
        for i in range(self.ne):
            equality[i] = self.__GetEqualityValue(i)
        return objective,inequality,equality

    def __GetGradients(self):
        gradObjective = self.__GetObjectiveGradient()
        gradInequality = [0 for _ in range(self.ni)]
        gradEquality = [0 for _ in range(self.ne)]
        for i in range(self.ni):
            gradInequality[i] = self.__GetInequalityGradient(i)
        for i in range(self.ne):
            gradEquality[i] = self.__GetEqualityGradient(i)
        return gradObjective,gradInequality,gradEquality

    def __GetLengthsDirections(self,gradObjective,inequality,gradInequality,equality,gradEquality):
        self.__RevertPossibleShapeModificationsDuringAnalysis() # must be done before mapping and damping in conversion to length-direction
        lInequality = [0 for _ in inequality]
        eInequality = [0 for _ in inequality]
        gradInequalityL = [0 for _ in inequality]
        lEquality = [0 for _ in equality]
        eEquality = [0 for _ in equality]
        gradEqualityL = [0 for _ in equality]

        _,eObjective,gradObjectiveL = self.__ValueFormatToLengthFormat(1,gradObjective)
        if self.isRosenbrock and self.conjugateGradient and self.iterCounter>1:
            eObjective = [u+v for u,v in zip(eObjective,self.r["eObjective"][-1])]
            eObjective = [u/norminf3d(eObjective) for u in eObjective]
        for i in range(self.ni):
            lInequality[i],eInequality[i],gradInequalityL[i] = self.__ValueFormatToLengthFormat(inequality[i],gradInequality[i])
        for i in range(self.ne):
            lEquality[i],eEquality[i],gradEqualityL[i] = self.__ValueFormatToLengthFormat(equality[i],gradEquality[i])

        # return eObjective,gradObjectiveL,lInequality,eInequality,gradInequalityL,lEquality,eEquality,gradEqualityL
        return eObjective,lInequality,eInequality,lEquality,eEquality

    def __CheckTerminationAndFinalFeasibility(self,iterCounter,stepLength,inequality,equality):
        # condition to enter termination, terminate only after having enforce the final feasibility
        if not self.isTerminationPhase and (iterCounter==self.iterMax or stepLength<self.stepLengthTermination):
            self.isTerminationPhase = True
            self.isEnforcingFinalFeasibility = True
            self.reduceOscillationFactor = 1
            if iterCounter==self.iterMax:
                print("Algo> begin final feasibility for termination (max nb of iter reached)")
            else:
                print("Algo> begin final feasibility for termination (stepLength<epsilon)")

        isFinalFeasibilityEnforced = self.__isFinalFeasibilityEnforced(inequality,equality)

        # if termination phase and final feasible, exit algo
        if self.isTerminationPhase and isFinalFeasibilityEnforced:
            print("Algo> final feasibility enforced. Termination")
            return True

        if self.isEnforcingFinalFeasibility and not isFinalFeasibilityEnforced:
            if sum(self.r["equality"][-1][i]*equality[i]<0 for i in range(len(equality))) >0: # reduce factor when equality overshoot and change sign, avoid oscillation
                self.reduceOscillationFactor = 0.7

        # enter final feasibility if flag file has been created in current folder
        isFlagFeasibility = os.path.isfile("final_feasibility_now.txt")
        if isFlagFeasibility:
            os.rename("final_feasibility_now.txt","final_feasibility_seen.txt")

        # enforce final feasibility before termination
        if isFinalFeasibilityEnforced:
            self.enforceFinalFeasibilityNextIter = iterCounter + self.enforceFinalFeasibilityEveryIter
            if self.isEnforcingFinalFeasibility:
                print("Algo> end enforce final feasibility. Enforce final feasibility next iter:", self.enforceFinalFeasibilityNextIter)
            elif self.enforceFinalFeasibilityEveryIter>0:
                print("Algo> final feasibility spontaneously enforced. Enforce final feasibility next iter:", self.enforceFinalFeasibilityNextIter)
            self.isEnforcingFinalFeasibility = False
        # enter final feasibility every iter
        elif iterCounter==self.enforceFinalFeasibilityNextIter or isFlagFeasibility:
            self.isEnforcingFinalFeasibility = True
            self.reduceOscillationFactor = 1
            if isFlagFeasibility:
                print("Algo> begin enforce final feasibility (flag final_feasibility_now.txt seen)")
            else:
                print("Algo> begin enforce final feasibility (every iter)")

        return False

    def __UpdateMesh(self,dx):
        if self.isRosenbrock:
            self.xRosenbrock = [self.xRosenbrock[i]+dx[i] for i in range(len(dx))]
            return
        # log previous shape
        self.__LogDesign()
        # log shape change of new shape
        x = self.__ReadDesignSurfaceToList()
        dxAbsolute = [x[i]+dx[i] - self.xInit[i] for i in range(self.n)]
        self.__WriteListToNodalVariable(dxAbsolute,SHAPE_CHANGE)
        #update mesh
        self.__WriteListToNodalVariable(dx,SHAPE_UPDATE)
        self.ModelPartController.UpdateMeshAccordingInputVariable( SHAPE_UPDATE )
        self.ModelPartController.SetReferenceMeshToMesh()

    def __StepLengthRule(self,stepLength,objective,lInequality,lEquality):
        # avoid interaction with final feasibility
        if self.isEnforcingFinalFeasibility:
            self.lastIterUpdateStepLength = self.iterCounter
            return stepLength

        # do nothing
        if self.iterCounter<3 or self.iterCounter<self.lastIterUpdateStepLength+2:
            return stepLength

        nSmooth = 5
        histObjSmooth = self.r["objective"][-(nSmooth-1):]+[objective]

        # reduce steplength if objective oscillate and no progress on distance feasible domain
        if self.__ReduceStepLengthCondition(objective,lInequality,lEquality):
            self.lastIterUpdateStepLength = self.iterCounter
            print("Algo> oscillating: decrease step length")
            return stepLength/self.stepLengthDecreaseFactor

        # increase slowly step length if all fine
        if self.iterCounter>=nSmooth and self.__IsSmoothProgress(histObjSmooth):
            self.lastIterUpdateStepLength = self.iterCounter
            print("Algo> smooth progress: increase step length")
            return min(self.maxStepLength, stepLength*self.stepLengthIncreaseFactor)

        #do nothing
        self.stepLength = stepLength
        return stepLength

    # express quantities relatively to step length
    # def __ConvertStepLength(self,stepLength,gradObjectiveL,lInequality0,gradInequalityL,lEquality0,gradEqualityL):
    #     gradObj = gradObjectiveL*stepLength
    #     gradIneq = [g*stepLength for g in gradInequalityL]
    #     gradEq = [g*stepLength for g in gradEqualityL]
    #     lineq = [l/stepLength for l in lInequality0]
    #     leq = [l/stepLength for l in lEquality0]
    #     return gradObj,lineq,gradIneq,leq,gradEq


##########################################################################################################
# Helper functions

    # --------------------------------------------------------------------------
    def __GetObjectiveValue( self ):
        value = self.__GetValue( self.objectiveId )
        if "objective" in self.r and len(self.r["objective"])>0:
            objectiveInit = self.r["objective"][0]
            percentInit = (objectiveInit-value)/objectiveInit*100
        else:
            percentInit = 0
        print("Algo> value objective {} = {:.3g} -> {:.2f}%".format(self.objectiveId,value,percentInit))
        return value

    # --------------------------------------------------------------------------
    def __GetInequalityValue( self,i ):
        value = self.__GetValue( self.inequalityIds[i] )
        self.inequalityReferences[i] = self.Communicator.getReferenceValue(self.inequalityIds[i])
        if self.inequalityReferences[i] != 0:
            percentRef = value/self.inequalityReferences[i]*100
            print("Algo> value inequality {} = {:.3g} -> {:.2f}%".format(self.inequalityIds[i],value,percentRef))
        else:
            print("Algo> value inequality {} = {:.3g}".format(self.inequalityIds[i],value))
        return value

    # --------------------------------------------------------------------------
    def __GetEqualityValue( self,i ):
        value = self.__GetValue( self.equalityIds[i] )
        self.equalityReferences[i] = self.Communicator.getReferenceValue(self.equalityIds[i])
        eqref = self.equalityReferences[i]
        eqref = 1 if eqref==0 else eqref
        percentRef = value/eqref*100
        print("Algo> value equality {} = {:.3g} -> {:.2f}%".format(self.equalityIds[i],value,percentRef))
        return value
    # --------------------------------------------------------------------------
    def __GetObjectiveGradient( self ):
        return self.__GetGradient(self.objectiveId,OBJECTIVE_SENSITIVITY)

    # --------------------------------------------------------------------------
    def __GetInequalityGradient( self,i ):
        return self.__GetGradient(self.inequalityIds[i],CONSTRAINT_SENSITIVITY)

    # --------------------------------------------------------------------------
    def __GetEqualityGradient( self,i ):
        return self.__GetGradient(self.equalityIds[i],CONSTRAINT_SENSITIVITY)

    # --------------------------------------------------------------------------
    def __GetValue( self, responseId ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf( responseId )
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )
        value = self.Communicator.getStandardizedValue( responseId )
        return value

    # --------------------------------------------------------------------------
    def __GetGradient( self, responseId, nodal_variable ):
        self.Communicator.initializeCommunication()
        self.Communicator.requestGradientOf( responseId )
        self.Analyzer.analyzeDesignAndReportToCommunicator( self.DesignSurface, self.iterCounter, self.Communicator )

        gradient = self.Communicator.getStandardizedGradient( responseId ) # kratos array or python list
        if self.isRosenbrock: # rosenbrock case
            return gradient

        self.__WriteGradientToNodalVariable( gradient, nodal_variable )
        # print("Algo.getGradient> gradient {}".format(responseId))
        if responseId in self.projectNormalsIds:
            self.GeometryUtilities.ComputeUnitSurfaceNormals()
            self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals( nodal_variable )
            # print("Algo.getGradient> project on normals, gradient {}".format(responseId))
        return self.__ReadNodalVariableToList(nodal_variable)

    # --------------------------------------------------------------------------
    def __RevertPossibleShapeModificationsDuringAnalysis( self ):
        if not self.isRosenbrock:
            self.ModelPartController.SetMeshToReferenceMesh()
            self.ModelPartController.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    # get length and directions.    project on normals (done in getGradient),damp,map,damp. tolerance
    def __ValueFormatToLengthFormat(self,value,gradient):
        gradientOriginal = deepcopy(gradient)
        if self.dampingIsSpecified:
            if self.dampOnlyAfterMapping:
                print("CAREFUL: DAMP ONLY ONCE")
            else:
                gradient = self.__DampVector(gradient)
                print("damp also before mapping")
        gradient = self.__MapVector(gradient)
        if self.dampingIsSpecified:
            gradient = self.__DampVector(gradient)
        ninf = norminf3d(gradient)
        eDir = [-gradient[i]/ninf for i in range(len(gradient))]
        gradL = sum(a*b for a,b in zip(gradientOriginal,eDir)) # dot product
        l = -value/gradL
        return l,eDir,gradL


    # --------------------------------------------------------------------------
    # map a vector
    def __MapVector(self, gradient_list):
        if self.isRosenbrock:
            return gradient_list
        # Here we use OBJECTIVE_SENSITIVITY as a temporary variable
        self.__WriteListToNodalVariable( gradient_list, CONSTRAINT_SENSITIVITY )
        tmp_timer = timer_factory.CreateTimer()
        # print("Algo.mapVector> map...",end="\r")
        with contextlib.redirect_stdout(None):
            startTime = tmp_timer.StartTimer()
            self.Mapper.MapToDesignSpace( CONSTRAINT_SENSITIVITY, MAPPED_CONSTRAINT_SENSITIVITY )
            self.Mapper.MapToGeometrySpace( MAPPED_CONSTRAINT_SENSITIVITY, CONSTRAINT_SENSITIVITY )
        # print("Algo.mapVector> map in {:.2f} s".format(tmp_timer.GetTotalTime()))
        return self.__ReadNodalVariableToList( CONSTRAINT_SENSITIVITY )

    # --------------------------------------------------------------------------
    def __DampVector(self, gradient_list):
        if self.isRosenbrock:
            return gradient_list
        # Here we use OBJECTIVE_SENSITIVITY as a temporary variable
        self.__WriteListToNodalVariable( gradient_list, CONSTRAINT_SENSITIVITY )
        tmp_timer = timer_factory.CreateTimer()
        # print("Algo.mapVector> damp...",end="\r")
        startTime = tmp_timer.StartTimer()
        with contextlib.redirect_stdout(None):
            self.DampingUtilities.DampNodalVariable( CONSTRAINT_SENSITIVITY )
        # print("Algo.mapVector> damp in {:.2f} s".format(tmp_timer.GetTotalTime()))
        return self.__ReadNodalVariableToList( CONSTRAINT_SENSITIVITY )

    def __isFinalFeasibilityEnforced(self,inequality,equality):
        for ineq in inequality:
            if ineq>0:
                return False
        for i,eq in enumerate(equality):
            if self.equalityFinalTolerances[i] is not None:
                absoluteFinalEqualityTolerance = self.equalityFinalTolerances[i]*self.equalityReferences[i]
                if abs(equality[i])>absoluteFinalEqualityTolerance:
                    return False
        return True

    # objective oscillate and no progress feasible
    def __ReduceStepLengthCondition(self,objective,lInequality,lEquality):
        r = self.r

        # objective oscillate
        v1,v2,v3 = r["objective"][-2:]+[objective]
        deltaObj = v1-v2
        if (v2<v1 and v2<v3) and v3>v1-0.1*deltaObj:
            objectiveOscillate = True
        else:
            objectiveOscillate = False

        # no progress objective
        if v3>v1:
            noProgressObjective = True
        else:
            noProgressObjective = False

        # no progress distance feasible
        distF = utils.simpleDistanceToFeasibleDomain(lInequality,lEquality)
        sOld = r["stepLength"][-2]
        distFOld = utils.simpleDistanceToFeasibleDomain([l*sOld for l in r["lInequality"][-2]] , [l*sOld for l in r["lEquality"][-2]])
        if distF>distFOld - 0.2*sOld:
            noProgressFeasible = True
        else:
            noProgressFeasible = False

        if (noProgressObjective or objectiveOscillate) and noProgressFeasible:
            return True
        else:
            return False

    def __IsSmoothProgress(self,histObj):
        n = len(histObj)
        v = histObj
        oscv = [0 for i in range(n-1)]
        for i in range(n-2):
            oscv[i] = (v[i+1]-v[i])*(v[i+2]-v[i+1])<0
        if sum(oscv)==0: # no oscillation at all
            return True
        else:
            return False

    def __LogVariables(self,locDict):
        if self.isRosenbrock:
            locDict["xRosenbrock"] = self.xRosenbrock
        if self.iterCounter == 2:
            self.__LogParam() # do it once after the beginning to get the reference values of the response functions

        logData = {}
        for name,value in locDict.items():
            if name=="self":
                continue
            if isinstance(value,list) and len(value)>0: # too long
                if len(value)>100:
                    continue
                if isinstance(value[0],list) and len(value[0])>100: # nested list too long
                    continue
            if name not in self.r:
                self.r[name] = []
            self.r[name].append(deepcopy(value))
            logData[name] = deepcopy(value)
        self.DataLogger.ResponseLogger.LogIteration(logData,self)

    def __LogDesign(self):
        # save all nodal variable to the gid file, and create next step
        if not self.isRosenbrock:
            self.DataLogger.DesignLogger.LogCurrentDesign(self.iterCounter)

    def __LogParam(self):
        if not hasattr(self,"param"):
            self.param = {}
        for k,v in vars(self).items():
            if isinstance(v,(str,int,float,bool)):
                self.param[k] = deepcopy(v)
            if isinstance(v,list):
                if len(v)>100: # too long
                    continue
                if len(v)>0 and isinstance(v[0],list) and len(v[0])>100: # nested list too long
                    continue
                self.param[k] = deepcopy(v)

##########################################################################################
# Kratos / Python IO

    # gradient can be a python list or a kratos array
    def __WriteGradientToNodalVariable(self,gradient,nodal_variable):
        if isinstance(gradient,list):
            self.__WriteListToNodalVariable(gradient,nodal_variable)
        else:
            self.__WriteArrayToNodalVariable(gradient,nodal_variable)

    # --------------------------------------------------------------------------
    def __WriteArrayToNodalVariable( self, kratos_array, nodal_variable ):
        for nodeId, tmp_gradient in kratos_array.items():
            self.OptimizationModelPart.Nodes[nodeId].SetSolutionStepValue(nodal_variable,0,tmp_gradient)

    # --------------------------------------------------------------------------
    def __WriteListToNodalVariable(self, pylist, nodal_variable):
        if self.isRosenbrock:
            return
        k = 0
        for node in self.DesignSurface.Nodes:
            tmp_values = pylist[k:k+3]
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

    # non private
    def ReadDesignSurfaceToList(self):
        return self.__ReadDesignSurfaceToList()
