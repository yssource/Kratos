from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from random import uniform

#TODO: test DEM bounding box

import os
import sys
import math
import itertools
from math import fabs
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
################ ADD HERE THE OTHER APPLICATIONS ################
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ErodibleBedApplication import *        #and the application for the soil simulation
from KratosMultiphysics.FSIApplication import *        #and the application for the soil simulation
from KratosMultiphysics.CustomDEMUtilsApplication import *        #and the application for the soil simulation
############### END OF EXTRA APPLICATIONS NEEDED ###############

import CFD_DEM_coupling
import swimming_DEM_procedures as SDP
import swimming_DEM_gid_output
import embedded
import variables_management as vars_man
import time as timer
import os

import bed_convection_solver

import swimming_DEM_algorithm

try:
    import define_output  # MA: some GUI write this file, some others not!
except ImportError:
    pass

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ:
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
    print("Running under MPI...........")
else:
    # DEM Application
    import DEM_procedures
    import DEM_material_test_script
    print("Running under OpenMP........")

sys.path.insert(0,'')

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.console_output_file_name = 'console_output.txt'
        self.path_to_console_out_file = os.getcwd()
        self.path_to_console_out_file += '/' + self.console_output_file_name
        self.log = open(self.path_to_console_out_file, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass


BaseAlgorithm = swimming_DEM_algorithm.Algorithm
class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = dict()):
        super(Algorithm,self).__init__(varying_parameters)

    def Initialize(self):

        print("\nInitializing Problem...")
        sys.stdout.flush()

        self.run_code = self.GetRunCode()

        # Moving to the recently created folder
        os.chdir(self.main_path)
        [self.post_path, data_and_results, self.graphs_path, MPI_results] = self.procedures.CreateDirectories(str(self.main_path), str(self.pp.CFD_DEM["problem_name"].GetString()), self.run_code)
        SDP.CopyInputFilesIntoFolder(self.main_path, self.post_path)

        #self.mixed_model_part = self.all_model_parts.Get('MixedPart')

        self.vars_man = vars_man
        self.vars_man.ConstructListsOfVariables(self.pp)

        self.fluid_algorithm.fluid_model_part.AddNodalSolutionStepVariable(HEIGHT)  #--#
        self.fluid_algorithm.fluid_model_part.AddNodalSolutionStepVariable(SEDIMENT_VELOCITY)  #--#
        self.disperse_phase_algorithm.rigid_face_model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT)  #--#
        self.disperse_phase_algorithm.rigid_face_model_part.AddNodalSolutionStepVariable(THICKNESS)  #--#
        self.CreateMeshSolver()  #--#


        self.FluidInitialize()
        self.DispersePhaseInitialize()

        self.SetAllModelParts()

        self.CreateBedModelPart() #--#
        self.BedSolverInitialize() #--#
        self.InitializeMeshSolver() #--#

        #initializing hydrostatic pressure!
        self.highest_z_coord=-100000.0
        for node in self.fluid_model_part.Nodes:
            if node.Z>self.highest_z_coord:
               self.highest_z_coord=node.Z
        self.SetOutletPressure()      

        self.SetCutsOutput()

        self.swimming_DEM_gid_io = swimming_DEM_gid_output.SwimmingDEMGiDOutput(self.pp.problem_name,
                                                                           self.pp.VolumeOutput,
                                                                           self.pp.GiDPostMode,
                                                                           self.pp.GiDMultiFileFlag,
                                                                           self.pp.GiDWriteMeshFlag,
                                                                           self.pp.GiDWriteConditionsFlag)

        self.swimming_DEM_gid_io.initialize_swimming_DEM_results(self.disperse_phase_algorithm.spheres_model_part,
                                                                 self.disperse_phase_algorithm.cluster_model_part,
                                                                 self.disperse_phase_algorithm.rigid_face_model_part,
                                                                 self.mixed_model_part)

        self.SetDragOutput()

        self.SetPointGraphPrinter()

        self.TransferGravityFromDisperseToFluid()

        # coarse-graining: applying changes to the physical properties of the model to adjust for
        # the similarity transformation if required (fluid effects only).
        SDP.ApplySimilarityTransformations(self.fluid_model_part,
                                           self.pp.CFD_DEM["similarity_transformation_type"].GetInt(),
                                           self.pp.CFD_DEM["model_over_real_diameter_factor"].GetDouble())

        # creating a Post Utils object that executes several post-related tasks
        self.post_utils = SDP.PostUtils(self.swimming_DEM_gid_io,
                                        self.pp,
                                        self.fluid_model_part,
                                        self.disperse_phase_algorithm.spheres_model_part,
                                        self.disperse_phase_algorithm.cluster_model_part,
                                        self.disperse_phase_algorithm.rigid_face_model_part,
                                        self.mixed_model_part)

        #creating an erodiblebed util
        self.erosion_utils = BedErosionUtility(self.bed_model_part) 

        self.deposit_particles_utils = DepositParticlesUtility(self.disperse_phase_algorithm.rigid_face_model_part) 

        # creating an IOTools object to perform other printing tasks
        self.io_tools = SDP.IOTools(self.pp)

        # creating a projection module for the fluid-DEM coupling
        self.h_min = 0.01
        n_balls = 1
        fluid_volume = 10
        self.pp.CFD_DEM.n_particles_in_depth = int(math.sqrt(n_balls / fluid_volume)) # only relevant in 2D problems
        # creating a physical calculations module to analyse the DEM model_part
        dem_physics_calculator = SphericElementGlobalPhysicsCalculator(self.disperse_phase_algorithm.spheres_model_part)

        field_utility = self.GetFieldUtility()

        if self.pp.CFD_DEM["coupling_level_type"].GetInt():

            if self.pp.CFD_DEM["meso_scale_length"].GetDouble() <= 0.0 and self.disperse_phase_algorithm.spheres_model_part.NumberOfElements(0) > 0:
                biggest_size = 2 * dem_physics_calculator.CalculateMaxNodalVariable(self.disperse_phase_algorithm.spheres_model_part, RADIUS)
                self.pp.CFD_DEM.meso_scale_length  = 20 * biggest_size

            elif self.disperse_phase_algorithm.spheres_model_part.NumberOfElements(0) == 0:
                self.pp.CFD_DEM.meso_scale_length  = 1.0

            self.projection_module = CFD_DEM_coupling.ProjectionModule(self.fluid_model_part, self.disperse_phase_algorithm.spheres_model_part, self.disperse_phase_algorithm.rigid_face_model_part, self.pp.domain_size, self.pp, field_utility)
            self.projection_module.UpdateDatabase(self.h_min)

        # creating a custom functions calculator for the implementation of additional custom functions
        self.custom_functions_tool = SDP.FunctionsCalculator(self.pp)

        # creating a derivative recovery tool to calculate the necessary derivatives from the fluid solution (gradient, laplacian, material acceleration...)
        self.derivative_recovery_tool = DerivativeRecoveryTool3D(self.fluid_model_part)

        # creating a stationarity assessment tool
        self.stationarity_tool = SDP.StationarityAssessmentTool(self.pp.CFD_DEM["max_pressure_variation_rate_tol"].GetDouble() , self.custom_functions_tool)

        # creating a debug tool
        self.dem_volume_tool = self.GetVolumeDebugTool()

        #self.SetEmbeddedTools()

        self.KRATOSprint("Initialization Complete" + "\n")
        sys.stdout.flush()

        self.step           = 0
        self.time           = self.pp.Start_time
        self.Dt             = self.pp.Dt
        self.out            = self.Dt
        self.final_time     = self.pp.CFD_DEM["FinalTime"].GetDouble()
        self.output_time    = self.pp.CFD_DEM["OutputTimeStep"].GetDouble()

        self.report.Prepare(self.timer, self.pp.CFD_DEM["ControlTime"].GetDouble())

        #first_print = True; index_5 = 1; index_10 = 1; index_50 = 1; control = 0.0

        if self.pp.CFD_DEM["ModelDataInfo"].GetBool():
            os.chdir(data_and_results)
            if self.pp.CFD_DEM.ContactMeshOption == "ON":
                (coordination_number) = self.procedures.ModelData(self.disperse_phase_algorithm.spheres_model_part, self.solver) # Calculates the mean number of neighbours the mean radius, etc..
                self.KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")
                os.chdir(self.main_path)
            else:
                self.KRATOSprint("Activate Contact Mesh for ModelData information")

        if self.pp.CFD_DEM["flow_in_porous_medium_option"].GetBool():
            fluid_frac_util = SDP.FluidFractionFieldUtility(self.fluid_model_part, self.pp.CFD_DEM.min_fluid_fraction )

            for field in self.pp.fluid_fraction_fields:
                fluid_frac_util.AppendLinearField(field)

            fluid_frac_util.AddFluidFractionField()

        if self.pp.CFD_DEM["flow_in_porous_DEM_medium_option"].GetBool():
            SDP.FixModelPart(self.disperse_phase_algorithm.spheres_model_part)

        # choosing the directory in which we want to work (print to)

        os.chdir(self.post_path)



        ######################################################################################################################################

        #                      I N I T I A L I Z I N G    T I M E    L O O P     ...   ( M I X E D    F L U I D / D E M    B L O C K )

        ######################################################################################################################################

        # setting up loop counters: Counter(steps_per_tick_step, initial_step, active_or_inactive_boolean, dead_or_not)
        self.fluid_solve_counter          = self.GetFluidSolveCounter()
        #self.embedded_counter             = self.GetEmbeddedCounter()
        self.DEM_to_fluid_counter         = self.GetBackwardCouplingCounter()
        self.derivative_recovery_counter  = self.GetRecoveryCounter()
        self.stationarity_counter         = self.GetStationarityCounter()
        self.print_counter                = self.GetPrintCounter()
        self.debug_info_counter           = self.GetDebugInfo()
        self.particles_results_counter    = self.GetParticlesResultsCounter()
        self.quadrature_counter           = self.GetHistoryForceQuadratureCounter()
        self.mat_deriv_averager           = SDP.Averager(1, 3)
        self.laplacian_averager           = SDP.Averager(1, 3)

        ##############################################################################
        #                                                                            #
        #    MAIN LOOP                                                               #
        #                                                                            #
        ##############################################################################

        self.DEM_step     = 0      # necessary to get a good random insertion of particles   # relevant to the stationarity assessment tool
        self.time_dem     = 0.0
        self.Dt_DEM       = self.disperse_phase_algorithm.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.disperse_phase_algorithm.rigid_face_model_part.ProcessInfo[DELTA_TIME] = self.Dt_DEM
        self.disperse_phase_algorithm.cluster_model_part.ProcessInfo[DELTA_TIME] = self.Dt_DEM
        self.stationarity = False

        self.report.total_steps_expected = int(self.pp.CFD_DEM["FinalTime"].GetDouble() / self.Dt_DEM)

        self.KRATOSprint(self.report.BeginReport(self.timer))

        # creating a Post Utils object that executes several post-related tasks
        self.post_utils_DEM = DEM_procedures.PostUtils(self.pp.CFD_DEM, self.disperse_phase_algorithm.spheres_model_part)

        SDP.InitializeVariablesWithNonZeroValues(self.fluid_model_part, self.disperse_phase_algorithm.spheres_model_part, self.pp) # otherwise variables are set to 0 by default

        self.CreateFemFluidCouplingList()

        self.InitializeBedShape()
        self.TransferHeightFromBedModelPartToFluidModelPart()

        (self.erosion_utils).MoveBedMesh(0.5)
        self.DeformFluidAndBedModelPartUsingHeight()
        self.TransferDisplacementsFromFluidPartToWallsFEMPart()
        (self.erosion_utils).MoveMeshUsingUserVariable(self.disperse_phase_algorithm.rigid_face_model_part.Nodes,MESH_DISPLACEMENT)

        self.FindHighestNodeAndElementIds()
        self.creator_destructor=ParticleCreatorDestructor() 
        self.creator_destructor.SetMaxNodeId(self.highest_node_id)

        self.SetUpResultsDatabase()

        # ANALYTICS BEGIN
        self.pp.CFD_DEM.perform_analytics_option = False

        if self.pp.CFD_DEM.perform_analytics_option:
            import analytics
            variables_to_measure = [PRESSURE]
            steps_between_measurements = 100
            gauge = analytics.Gauge(self.fluid_model_part, self.Dt, self.final_time, variables_to_measure, steps_between_measurements)
            point_coors = [0.0, 0.0, 0.01]
            target_node = SDP.FindClosestNode(self.fluid_model_part, point_coors)
            target_id = target_node.Id
            print(target_node.X, target_node.Y, target_node.Z)
            print(target_id)
            def condition(node):
                return node.Id == target_id

            gauge.ConstructArrayOfNodes(condition)
            print(gauge.variables)
            #print_analytics_counter = SDP.Counter( 5 * steps_between_measurements, 1, 1) # MA: not used anywhere?
        # ANALYTICS END

        import derivative_recovery.derivative_recovery_strategy as derivative_recoverer

        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(self.pp, self.fluid_model_part, self.derivative_recovery_tool, self.custom_functions_tool)

        self.FillHistoryForcePrecalculatedVectors()

        self.PerformZeroStepInitializations()

        self.pp.rigid_faces_nodal_results += ["THICKNESS"]
        self.pp.rigid_faces_nodal_results += ["MESH_DISPLACEMENT"]
        self.post_utils.Writeresults(self.time)
 
    def RunMainTemporalLoop(self):

        while self.time <= self.final_time:

            self.time = self.time + self.Dt
            self.step += 1
            self.TransferTimeToFluidSolver()
            self.CloneTimeStep()
            self.TellTime(self.time)


            if self.pp.CFD_DEM["coupling_scheme_type"].GetString() == "UpdatedDEM":
                time_final_DEM_substepping = self.time + self.Dt

            else:
                time_final_DEM_substepping = self.time

            #self.PerformEmbeddedOperations() #TODO: it's crashing

            # solving the fluid part

            if self.step >= 3 and not self.stationarity:
                print("Solving Fluid... (", self.fluid_model_part.NumberOfElements(0), "elements )")
                sys.stdout.flush()

                if self.fluid_solve_counter.Tick():
                    self.SetOutletPressure() ##--##
                    self.FluidSolve(self.time)

            # assessing stationarity

                if self.stationarity_counter.Tick():
                    print("Assessing Stationarity...")
                    sys.stdout.flush()
                    self.stationarity = self.stationarity_tool.Assess(self.fluid_model_part)
                    sys.stdout.flush()

            # printing if required

            if self.particles_results_counter.Tick():
                # eliminating remote balls

                #if self.pp.dem.BoundingBoxOption == "ON":
                #    self.creator_destructor.DestroyParticlesOutsideBoundingBox(self.disperse_phase_algorithm.spheres_model_part)

                self.io_tools.PrintParticlesResults(self.pp.variables_to_print_in_file, self.time, self.disperse_phase_algorithm.spheres_model_part)
                #self.graph_printer.PrintGraphs(self.time) #MA: commented out because the constructor was already commented out
                self.PrintDrag(self.drag_list, self.drag_file_output_list, self.fluid_model_part, self.time)

            if self.output_time <= self.out and self.pp.CFD_DEM["coupling_scheme_type"].GetString() == "UpdatedDEM":

                if self.pp.CFD_DEM["coupling_level_type"].GetInt() > 0:
                    self.projection_module.ComputePostProcessResults(self.disperse_phase_algorithm.spheres_model_part.ProcessInfo)

                self.post_utils.Writeresults(self.time)
                self.out = 0

            # solving the DEM part
            self.derivative_recovery_counter.Switch(self.time > self.pp.CFD_DEM["interaction_start_time"].GetDouble())

            if self.derivative_recovery_counter.Tick():
                self.recovery.Recover()



            self.ComputeNodalFlux()
            self.InjectParticlesIntoStream()
            self.TransferHeightFromBedModelPartToFluidModelPart()
            #self.DeformFluidAndBedModelPartUsingHeight()
            self.TransferDisplacementsFromFluidPartToWallsFEMPart()
            (self.erosion_utils).MoveMeshUsingUserVariable(self.disperse_phase_algorithm.rigid_face_model_part.Nodes,MESH_DISPLACEMENT)

            print("Solving DEM... (", self.disperse_phase_algorithm.spheres_model_part.NumberOfElements(0), "elements )")
            sys.stdout.flush()
            first_dem_iter = True

            for self.time_dem in self.yield_DEM_time(self.time_dem, time_final_DEM_substepping, self.Dt_DEM):
                self.DEM_step += 1   # this variable is necessary to get a good random insertion of particles
                self.disperse_phase_algorithm.spheres_model_part.ProcessInfo[TIME_STEPS]    = self.DEM_step
                self.disperse_phase_algorithm.rigid_face_model_part.ProcessInfo[TIME_STEPS] = self.DEM_step
                self.disperse_phase_algorithm.cluster_model_part.ProcessInfo[TIME_STEPS]    = self.DEM_step

                self.PerformInitialDEMStepOperations(self.time_dem)

                if self.time >= self.pp.CFD_DEM["interaction_start_time"].GetDouble() and self.pp.CFD_DEM["coupling_level_type"].GetInt() and (self.pp.CFD_DEM["project_at_every_substep_option"].GetBool() or first_dem_iter):
                    if self.pp.CFD_DEM["basset_force_type"].GetInt() > 0:
                        print("WARNING. COMPUTING HISTORICAL FORCE. TOO SLOW!")

                    if self.pp.CFD_DEM["coupling_scheme_type"].GetString() == "UpdatedDEM":
                        self.ApplyForwardCoupling()

                    else:
                        self.ApplyForwardCoupling((time_final_DEM_substepping - self.time_dem) / self.Dt)

                        if self.pp.CFD_DEM["IntegrationScheme"].GetString() in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
                            self.DEMSolve(self.time_dem)
                            self.ApplyForwardCouplingOfVelocityOnly(self.time_dem)
                        else:
                            if self.pp.CFD_DEM["basset_force_type"].GetInt() > 0:
                                node.SetSolutionStepValue(SLIP_VELOCITY_X, vx)
                                node.SetSolutionStepValue(SLIP_VELOCITY_Y, vy)

                        if self.quadrature_counter.Tick():
                            self.AppendValuesForTheHistoryForce()

                # performing the time integration of the DEM part

                self.disperse_phase_algorithm.spheres_model_part.ProcessInfo[TIME]    = self.time_dem
                self.disperse_phase_algorithm.rigid_face_model_part.ProcessInfo[TIME] = self.time_dem
                self.disperse_phase_algorithm.cluster_model_part.ProcessInfo[TIME]    = self.time_dem

                if self.pp.do_solve_dem:
                    self.DEMSolve(self.time_dem)

                #self.disperse_phase_algorithm.DEMFEMProcedures.MoveAllMeshes(self.all_model_parts, self.time_dem, self.Dt_DEM)

                #### TIME CONTROL ##################################

                # adding DEM elements by the inlet:
                #--# REPLACING THE ORIGINAL INLET WITH OUR NEW BED INLET - Granada
                #--#if self.pp.CFD_DEM.dem_inlet_option:
                #--#    self.disperse_phase_algorithm.DEM_inlet.CreateElementsFromInletMesh(self.disperse_phase_algorithm.spheres_model_part, self.disperse_phase_algorithm.cluster_model_part, self.disperse_phase_algorithm.creator_destructor)  # After solving, to make sure that neighbours are already set.
                #--# ADD HERE ALL THE THINGS RELATED TO INLET - Granada
                #--# FIRST WE MUST COMPUTE THE NODAL FLUX AT THE BOTTOM.
                #--# REMOVED FROM HERE!!

                if self.output_time <= self.out and self.pp.CFD_DEM["coupling_scheme_type"].GetString() == "UpdatedFluid":

                    if self.pp.CFD_DEM["coupling_level_type"].GetInt():
                        self.projection_module.ComputePostProcessResults(self.disperse_phase_algorithm.spheres_model_part.ProcessInfo)

                    self.post_utils.Writeresults(self.time_dem)
                    self.out = 0

                #updating nodal information on the walls model part using the particles that landed there.
                self.deposit_particles_utils.IncreaseNodalThicknessUsingParticles(first_dem_iter)
		#now we have to transfer this information to the fluid model part, and from there transfer it to the bed model part.

                self.out = self.out + self.Dt_DEM
                first_dem_iter = False

                # applying DEM-to-fluid coupling

                if self.DEM_to_fluid_counter.Tick() and self.time >= self.pp.CFD_DEM["interaction_start_time"].GetDouble():
                    self.projection_module.ProjectFromParticles()



            print("FINISHED DEM")
            #now we can use the new thickness (due to particle deposition) to deform the 3 meshes:
            self.deposit_particles_utils.NormalizeThickness()
            self.TransferThicknessFEMPartToFEMFluidAndBedDisplamecents()
            (self.erosion_utils).MoveMeshUsingUserVariable(self.bed_model_part.Nodes,MESH_DISPLACEMENT)            

            #### PRINTING GRAPHS ####
            os.chdir(self.graphs_path)
            # measuring mean velocities in a certain control volume (the 'velocity trap')
            if self.pp.CFD_DEM["VelocityTrapOption"].GetBool():
                self.post_utils_DEM.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", self.time)

            os.chdir(self.post_path)

            # coupling checks (debugging)
            if self.debug_info_counter.Tick():
                self.dem_volume_tool.UpdateDataAndPrint(self.pp.CFD_DEM["fluid_domain_volume"].GetDouble())

            # printing if required

            if self.particles_results_counter.Tick():
                self.io_tools.PrintParticlesResults(self.pp.variables_to_print_in_file, self.time, self.disperse_phase_algorithm.spheres_model_part)
                #self.graph_printer.PrintGraphs(self.time) #MA: commented out because the constructor was already commented out
                self.PrintDrag(self.drag_list, self.drag_file_output_list, self.fluid_model_part, self.time)

    def ComputeNodalFlux(self): 
        for nodeb in self.bed_model_part.Nodes:
            nodeb.SetSolutionStepValue(SEDIMENT_VELOCITY_X,0.0)
            nodeb.SetSolutionStepValue(SEDIMENT_VELOCITY_Y,0.0)
            nodeb.SetSolutionStepValue(SEDIMENT_VELOCITY_Z,0.0)
            nodeb.SetSolutionStepValue(NODAL_AREA,0.0) #THIS WILL BE THE NUMBER OF ELEMENTS; NOT THE NODAL AREA: ASSUMING SAME-AREA ELEMS!!!!

        i = 0
        for i in range(0,len(self.bed_nodes_in_fluid_model_part)):
            self.bed_nodes_in_fluid_model_part[i].SetSolutionStepValue(SEDIMENT_VELOCITY_X,0.0)
            self.bed_nodes_in_fluid_model_part[i].SetSolutionStepValue(SEDIMENT_VELOCITY_Y,0.0)
            self.bed_nodes_in_fluid_model_part[i].SetSolutionStepValue(SEDIMENT_VELOCITY_Z,0.0)
            self.bed_nodes_in_fluid_model_part[i].SetSolutionStepValue(NODAL_AREA,0.0) #THIS WILL BE THE NUMBER OF ELEMENTS; NOT THE NODAL AREA: ASSUMING SAME-AREA ELEMS!!!!
            i += 1

        ################## ADVECTIVE VELOCITY ####################
        for ie in range(0,len(self.bed_elements_in_fluid_model_part)):
            nnodes = 0
            coorX = [0,0,0,0]
            coorY = [0,0,0,0]
            coorZ = [0,0,0,0]
            id_list = [0,0,0]
            v1    = [0,0,0]
            v2    = [0,0,0]
            v3    = [0,0,0]
            vn    = [0,0,0]
            for nodef in self.bed_elements_in_fluid_model_part[ie].GetNodes():
            #for nodef in elemf.GetNodes():
                if nodef.Z0 < self.lowest_z+1.0e-6:
                   coorX[nnodes] = nodef.X
                   coorY[nnodes] = nodef.Y
                   coorZ[nnodes] = nodef.Z
                   id_list[nnodes] = nodef.Id
                   nnodes = nnodes + 1
                else:
                   coorX[3] = nodef.X
                   coorY[3] = nodef.Y
                   coorZ[3] = nodef.Z
                   ux = nodef.GetSolutionStepValue(VELOCITY_X)
                   uy = nodef.GetSolutionStepValue(VELOCITY_Y)
                   uz = nodef.GetSolutionStepValue(VELOCITY_Z)
            if nnodes == 3:
               v1[0] = coorX[0] - coorX[1]
               v1[1] = coorY[0] - coorY[1]
               v1[2] = coorZ[0] - coorZ[1]
               v2[0] = coorX[2] - coorX[1]
               v2[1] = coorY[2] - coorY[1]
               v2[2] = coorZ[2] - coorZ[1]
               v3[0] = coorX[3] - coorX[1]
               v3[1] = coorY[3] - coorY[1]
               v3[2] = coorZ[3] - coorZ[1]
               vn[0] = v1[1]*v2[2] - v1[2]*v2[1]
               vn[1] = v1[2]*v2[0] - v1[0]*v2[2]
               vn[2] = v1[0]*v2[1] - v1[1]*v2[0]
               modvn = math.sqrt(vn[0]**2 + vn[1]**2 + vn[2]**2)
               vn[0] = vn[0]/modvn
               vn[1] = vn[1]/modvn
               vn[2] = vn[2]/modvn
               dtet = (vn[0]*v3[0] + vn[1]*v3[1] + vn[2]*v3[2])
               if dtet <0:
                  dtet  = -dtet
                  vn[0] = -vn[0]
                  vn[1] = -vn[1]
                  vn[2] = -vn[2]
               #print("vn",elemf.Id,vn[0],vn[1],vn[2])
               auxc = (ux*vn[0] + uy*vn[1] + uz*vn[2])
               utx  = ux - auxc*vn[0]
               uty  = uy - auxc*vn[1]
               utz  = uz - auxc*vn[2]
               #print("ut",utx,uty,utz)
               z0 = 0.0001
               ustx = 0.41*utx/math.log(dtet/z0)
               usty = 0.41*uty/math.log(dtet/z0)
               ustz = 0.41*utz/math.log(dtet/z0)
               rhom = 2650
               rhof = 1000
               lamb = 0.5
               rhos = rhom*(1-lamb)
               grv  = 9.81
               dpart = 0.0004            
               th  = math.acos(vn[2])
               #ut  = math.sqrt(math.sin(th)/0.625 + math.cos(th))*uto
               modust = math.sqrt(ustx**2 + usty**2 + ustz**2)
               modustxy = math.sqrt(ustx**2 + usty**2)+1.0e-9
               #print("qsx",qsx,"ut/mod",ut/modust,"modust",modust,"ut",ut,"ut/modust",(1.0-ut/modust)) 
               #################SEDIMENT MODELS ####################################
               taus   = rhof*modust**2/((rhom-rhof)*grv*dpart)
               tauscr = 0.05
               ################# Engelund and Fredsoe (1976)  ######################
               #phis   = 18.74*(taus-tauscr)*(math.sqrt(taus)-0.7*math.sqrt(tauscr))
               ################# Meyer-Peter and Muler (1948) ######################
               phis   = 8.0*(max((taus-tauscr),0))**1.5
               qs     = 1.0*rhom*phis*math.sqrt((rhom/rhof-1.0)*grv*dpart**3)
               qsx    = qs*ustx/modustxy 
               qsy    = qs*usty/modustxy 
               #print("Ux",ux,"Uy",uy,"Uz",uz,"utx",utx,"uty",uty,"utz",utz)            
               #print("taus",taus,"ustx",ustx,"usty",usty,"ustz",ustz,"mod",modust,"phis",phis,"qs",qs,"qsx",qsx,"qsy",qsy)
               for nodef in self.bed_elements_in_fluid_model_part[ie].GetNodes():
                   if (nodef.Id==id_list[0]) or (nodef.Id==id_list[1]) or (nodef.Id==id_list[2]):
                      height = nodef.GetSolutionStepValue(HEIGHT,0)
                      height = height - 0.48
                      #print("height",height,"modust",modust)       
                      velax = qsx/(rhos*(abs(height)+1.0e-50))
                      velay = qsy/(rhos*(abs(height)+1.0e-50))
                      velxaux = nodef.GetSolutionStepValue(SEDIMENT_VELOCITY_X,0)
                      velyaux = nodef.GetSolutionStepValue(SEDIMENT_VELOCITY_Y,0)
                      nodef.SetSolutionStepValue(SEDIMENT_VELOCITY_X,velax/3.0 + velxaux) 
                      nodef.SetSolutionStepValue(SEDIMENT_VELOCITY_Y,velay/3.0 + velyaux)          
                      nodef.SetSolutionStepValue(SEDIMENT_VELOCITY_Z,0)
                      nodal_areaaux=nodef.GetSolutionStepValue(NODAL_AREA,0)
                      nodef.SetSolutionStepValue(NODAL_AREA,1.0/3.0 + nodal_areaaux)          
        #dividing by nodal area (actually number of contributions, assuming constant area elements.)              
        for nodeb in self.bed_model_part.Nodes:
            nodal_area=nodef.GetSolutionStepValue(NODAL_AREA,0)
            velxaux = nodef.GetSolutionStepValue(SEDIMENT_VELOCITY_X,0)
            velyaux = nodef.GetSolutionStepValue(SEDIMENT_VELOCITY_Y,0)
            nodeb.SetSolutionStepValue(SEDIMENT_VELOCITY_X,velxaux/nodal_area)
            nodeb.SetSolutionStepValue(SEDIMENT_VELOCITY_Y,velyaux/nodal_area)

        i = 0
        for nodeb in self.bed_model_part.Nodes:
           sdx = self.bed_nodes_in_fluid_model_part[i].GetSolutionStepValue(SEDIMENT_VELOCITY_X,0)
           sdy = self.bed_nodes_in_fluid_model_part[i].GetSolutionStepValue(SEDIMENT_VELOCITY_Y,0)      
           sdz = self.bed_nodes_in_fluid_model_part[i].GetSolutionStepValue(SEDIMENT_VELOCITY_Z,0)
           i += 1
           #print("sdx, sdy, sdz", sdx, sdy,sdz)
           nodeb.SetSolutionStepValue(VELOCITY_X,sdx)
           nodeb.SetSolutionStepValue(VELOCITY_Y,sdy)
           nodeb.SetSolutionStepValue(VELOCITY_Z,sdz)

           nodeb.SetSolutionStepValue(SEDIMENT_VELOCITY_X,sdx)
           nodeb.SetSolutionStepValue(SEDIMENT_VELOCITY_Y,sdy)
           nodeb.SetSolutionStepValue(SEDIMENT_VELOCITY_Z,sdz)


    def CreateBedModelPart(self):
       print("Creating Bed Model Part")
       self.bed_model_part = ModelPart("BedModelPart");  #we create a model part  
       bed_convection_solver.AddVariables(self.bed_model_part)

       #finding upper part and lower(soil) of the domain (to fix those node and avoid deforming the mesh there)
       self.highest_z=-1000000.0
       self.lowest_z=10000.0
       for node in self.fluid_model_part.Nodes:
           #right now all the nodes of the bottom are flagged as GetValue(IS_STRUCTURE)=true.
           #but this thing is used to both activate the wall_law drag and the rotation.
           #we DO NOT want these nodes to be rotated since i have fixed Y and Z velocity in these nodes, which is better that the normal procedure.
           #therefore i have modified the condition so that it uses GetSolutionStepValue(IS_STRUCTURE)=true to add the contributions but
           #uses GetValue(IS_STRUCTURE)=true to rotate.
           #right now all the nodes of the bottom have GetValue(IS_STRUCTURE)=true. to begin I copy this to GetSolutionStepValue(IS_STRUCTURE) so that all of them contribute to the drag
           node.SetSolutionStepValue(IS_STRUCTURE,(node.GetValue(IS_STRUCTURE)))
           #having done that, i set GetValue(IS_STRUCTURE)=false in the corner nodes (fixed Y) so that these nodes are not rotated, yet have contributions.
           if node.IsFixed(VELOCITY_Y):
               node.SetValue(IS_STRUCTURE,0.0)

           if node.Z>self.highest_z:
               self.highest_z=node.Z
           if node.Z<self.lowest_z:
                self.lowest_z=node.Z

       #now we can fix the displacement and copy the bottom(soil) nodes into the bed model part.
       print("highest z=",self.highest_z)
       print("lowest z=",self.lowest_z)
       self.bed_nodes_in_fluid_model_part = []
       for node in self.fluid_model_part.Nodes:
           node.Fix(MESH_DISPLACEMENT_X)
           node.Fix(MESH_DISPLACEMENT_Y)
           if node.Z>(self.highest_z-1.0e-6):
               node.Fix(MESH_DISPLACEMENT_X)
               node.Fix(MESH_DISPLACEMENT_Y)
               node.Fix(MESH_DISPLACEMENT_Z)
           if node.Z<(self.lowest_z+1.0e-6):  # a secondary criterion can be added here if there are undeformable bottom zones. Fix displ there!!
               node.Fix(MESH_DISPLACEMENT_X)
               node.Fix(MESH_DISPLACEMENT_Y)
               node.Fix(MESH_DISPLACEMENT_Z)
               self.bed_nodes_in_fluid_model_part.append(node)
               Id = node.Id #the IDs of the fluid and bed model part are the same!
               x = node.X
               y = node.Y
               z = node.Z
               self.bed_model_part.CreateNewNode(Id,x,y,z)
       #we save a list of the fluid elements touching the bed, needed to compute the shear stress (not a model part!!)
       self.bed_elements_in_fluid_model_part = []
       for elemf in self.fluid_model_part.Elements:
           nnodes = 0
           for nodef in elemf.GetNodes():
               if nodef.Z0<self.lowest_z+1e-6:
                  nnodes = nnodes + 1
           if nnodes == 3:
              self.bed_elements_in_fluid_model_part.append(elemf)

       #now the elements of the soil_model_part, which must be a copy of the conditions in the fluid model part.
       elem_Id=1
       for condition in self.fluid_model_part.Conditions:
           is_erosionable=True
           for node in condition.GetNodes():     
              if node.Z>self.lowest_z+1.0e-6:
                   is_erosionable=False
           if is_erosionable: #OK, 3 nodes at the bed (all at Z=lowest_z)
              node_list = [1,1,1]
              position = 0
              for node in condition.GetNodes():     
                  node_Id = node.Id
                  node_list[position] = node_Id
                  position +=1
              #the previous order is only useful if the triangle has positive area.
              area = self.ComputeTriangleAreaInXYPlane(condition)
              #if not. we change the order.
              if area<0.0:
                  temp_int = node_list[2]
                  node_list[2] = node_list[1]
                  node_list[1] = temp_int
              new_elem = self.bed_model_part.CreateNewElement("ErodibleBed2D",elem_Id,node_list,self.fluid_model_part.Properties[0])
              area = self.ComputeTriangleAreaInXYPlane(new_elem)
              if area<0.0:
                  print("ERROR: COULD NOT ORIENT BED ELEMENTS CORRECTLY")
                  quit()
              elem_Id+=1
       #the buffer size should be set up here after the mesh is read for the first time
       self.bed_model_part.SetBufferSize(2)
       #and the ddofs
       bed_convection_solver.AddDofs(self.bed_model_part)


    def InitializeBedShape(self):
       ##### INITIALIZING THE HEIGHT!!!!!!!!! TO BE CHANGED LATER #####
       ########### SINUSOIDAL #########################################
       dsp = 3.000
       lam = 0.740
       hdn = 0.050
       for node in self.bed_model_part.Nodes:
           node.SetSolutionStepValue(HEIGHT,0.5)
           if dsp < node.X < (dsp+3.0*lam):
              node.SetSolutionStepValue(HEIGHT,0.5+hdn*math.sin((node.X-dsp)/lam*math.pi)**2)
           if node.X<0.0:
              node.Fix(HEIGHT)
           #node.SetSolutionStepValue(HEIGHT,0.50)
   
    def TransferHeightFromBedModelPartToFluidModelPart(self):
        i = 0
        for nodeb in self.bed_model_part.Nodes:
            height = nodeb.GetSolutionStepValue(HEIGHT)
            #print("altura",height)
            self.bed_nodes_in_fluid_model_part[i].SetSolutionStepValue(MESH_DISPLACEMENT_Z,height-0.5)
            self.bed_nodes_in_fluid_model_part[i].SetSolutionStepValue(HEIGHT,height)
            i += 1

    def TransferDisplacementsFromFluidPartToWallsFEMPart(self): 
        i = 0
        for nodeb in self.disperse_phase_algorithm.rigid_face_model_part.Nodes:
            displ = self.FEM_nodes_in_fluid_model_part[i].GetSolutionStepValue(MESH_DISPLACEMENT)
            nodeb.SetSolutionStepValue(MESH_DISPLACEMENT, displ)
            i += 1

    def BedSolverInitialize(self):
        self.bed_solver = bed_convection_solver.BedConvectionSolver(self.bed_model_part, self.pp.domain_size-1) #1 dimension smaller. must always be in the xy plane!!!!!
        self.bed_solver.Initialize()
        #to print the bed model part in GID format
        gid_mode = GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need  
        multifile = MultiFileFlag.SingleFile #MultipleFiles
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteElementsOnly
        self.bed_gid_io = GidIO("results_bed",gid_mode,multifile,deformed_mesh_flag,write_conditions)
        bed_mesh_name = 0.0
        self.bed_gid_io.InitializeMesh(0.0);
        self.bed_gid_io.WriteMesh((self.bed_model_part).GetMesh());
        self.bed_gid_io.FinalizeMesh()
        self.bed_gid_io.InitializeResults(bed_mesh_name,(self.bed_model_part).GetMesh())

    def DeformFluidAndBedModelPartUsingHeight(self):
        i = 0
        for nodeb in self.bed_model_part.Nodes:
            #########
            #########
            # I DONT KNOW IF THIS IS ENOUGH TO HAVE A DEFORMED ELEMENT. MAYBE WE SHOULD CHANGE THIS TO : node.Z = node.Z0  + MESH_DISPLACEMENT_Z
            #########
            #########
            #print("altura",height)
            height = nodeb.GetSolutionStepValue(HEIGHT) 
            self.bed_nodes_in_fluid_model_part[i].SetSolutionStepValue(MESH_DISPLACEMENT_Z,height-0.5)
            i += 1
        print("MESH SOLVER")
        self.mesh_sol.Solve();
        print("finished deforming fluid mesh")
        self.erosion_utils.MoveBedMesh(0.5)
        print("finished deforming bed mesh")

        for node in self.fluid_model_part.Nodes:
           node.SetSolutionStepValue(MESH_VELOCITY,[0.0,0.0,0.0])


    def CreateMeshSolver(self):
        os.chdir(self.main_path)
        mesh_solver_parameter_file = open("ProjectParameters_mesh_deformation.json",'r')
        MeshSolverProjectParameters = Parameters( mesh_solver_parameter_file.read())
        mesh_solver_module = __import__(MeshSolverProjectParameters["solver_settings"]["solver_type"].GetString())
        #import mesh_solver_structural_similarity as mesh_solver_module # <------ ALTERNATIVE TO LINE ABOVE
        self.mesh_sol = mesh_solver_module.CreateSolver(self.fluid_algorithm.fluid_model_part, MeshSolverProjectParameters["solver_settings"])
        self.mesh_sol.AddVariables()


    def InitializeMeshSolver(self):
        self.mesh_sol.AddDofs()
        self.mesh_sol.Initialize()
        self.fluid_model_part.CloneTimeStep(0.0000000001)

        #self.InitializeBedShape(self) #just in case it was not done before .WARNING, WE ARE DESTROYING ANY PREVIOUS INFORMATION!!!!!

    def InjectParticlesIntoStream(self):
        #first we reset how much mass we have lost/gained. the mass loss due to particle injection will be substracted from the thickness
        for node in self.bed_model_part.Nodes:
              node.SetSolutionStepValue(THICKNESS,0.0)
              node.SetSolutionStepValue(NODAL_AREA,0.0)

        properties = self.spheres_model_part.GetProperties()[0]
        if (properties.GetValue(PARTICLE_DENSITY)< 0.00001):
             print("error: properties[0] not defined in DEM model part. Cannot inject particles	")
             quit()
        element_name = "SphericSwimmingParticle3D"
        for elem in self.bed_model_part.Elements:
              #first we must get the average flux.
              #FOR THE FOLLOWING LINES USE CalculateVolume2D in C++, but we need the coordinates anyway for the seeding!
              average_flux=0.0;
              x0 = elem.GetNodes()[0].X
              y0 = elem.GetNodes()[0].Y
              z0 = elem.GetNodes()[0].Z
              x1 = elem.GetNodes()[1].X
              y1 = elem.GetNodes()[1].Y
              z1 = elem.GetNodes()[1].Z
              x2 = elem.GetNodes()[2].X
              y2 = elem.GetNodes()[2].Y
              z2 = elem.GetNodes()[2].Z
              
              area = self.ComputeTriangleAreaInXYPlane(elem)
             
              for node in elem.GetNodes():     
                  flux = node.GetSolutionStepValue(SEDIMENT_VELOCITY) * node.GetSolutionStepValue(HEIGHT)
                  average_specific_flux = 0.33333333333333333 * flux
             
              volume_loss_in_this_element = 0.0 
              #finally the flux (per time unit) coming out of this triangle should be:
              total_flux = average_specific_flux * area
              #first we define the number of particles.
              nnew_particles=1
              #now we inject the particles:
              for particle in range(0,nnew_particles):
                  N0= uniform(0.0,1.0)
                  N1 = uniform(0.0,1.0)
                  N2 = uniform(0.0,1.0)
                  Nsum=N0+N1+N2
                  N0= N0/Nsum
                  N1= N1/Nsum
                  N2= N2/Nsum
                  coordinates = Array3()
                  coordinates[0] = x0*N0 + x1*N1 + x2*N2
                  coordinates[1] = y0*N0 + y1*N1 + y2*N2
                  coordinates[2] = z0*N0 + z1*N1 + z2*N2        
                  #coordinates[0]=uniform(0.0,7.6)
                  #coordinates[1]=uniform(0.0,0.6)
                  #coordinates[2] = 0.005  
                  radius = 0.0025
                  coordinates[2] = coordinates[2] +  radius*1.1 #it must not be touching the bed when injected 
                  self.highest_elem_id = self.highest_elem_id + 1
                  created_element = self.creator_destructor.CreateSphericParticle(self.disperse_phase_algorithm.spheres_model_part, self.highest_elem_id , coordinates, properties, radius, element_name)            
                  created_node = created_element.GetNodes()[0]
                  created_node.SetSolutionStepValue(VELOCITY_X,0.0)
                  created_node.SetSolutionStepValue(VELOCITY_Y,0.0)
                  created_node.SetSolutionStepValue(VELOCITY_Z, 0.5)

                  #compute here the thickness loss!
                  #be careful if this is specific (per square meter) or absolute, to decide how to use the area!
                  #volume_loss_in_this_element = volume_loss_in_this_element + 0.00001 # whatever this element represents
                  volume_loss_in_this_element = volume_loss_in_this_element + 0.0

              thickness_loss_in_this_element = volume_loss_in_this_element / area 
              for node in elem.GetNodes():  
                  old_nodal_area = node.GetSolutionStepValue(NODAL_AREA)
                  weight =  area*1.0/3.0
                  node.SetSolutionStepValue(NODAL_AREA,old_nodal_area + weight)
                  old_thickness = node.GetSolutionStepValue(THICKNESS)
                  node.SetSolutionStepValue(THICKNESS,old_thickness - thickness_loss_in_this_element * weight ) #it is a loss, so substracting!

        #dividing by nodal area for weighting
        for node in self.bed_model_part.Nodes:
              old_thx = node.GetSolutionStepValue(THICKNESS)
              area = node.GetSolutionStepValue(NODAL_AREA)
              thickness_change = old_thx/area
              node.SetSolutionStepValue(THICKNESS,thickness_change)
              #we can finally compute the new height
              old_height = node.GetSolutionStepValue(HEIGHT)
              node.SetSolutionStepValue(HEIGHT,old_height + thickness_change )



        
        
    def CreateFemFluidCouplingList(self): #creates a list of pointer nodes of size (wall(aka_fem)_nodes) , in which the ith pointer of this list points to the fluid node that corresponds to the FEM (wall) node in the ith position of bed_model_part.nodes
       #now we can fix the displacement and copy the bottom(soil) nodes into the bed model part.
       print("highest z=",self.highest_z)
       print("lowest z=",self.lowest_z)
       self.FEM_nodes_in_fluid_model_part = []
       for node in self.disperse_phase_algorithm.rigid_face_model_part.Nodes:
           coord_wall_node = [ node.X , node.Y , node.Z ] #the IDs of the fluid and bed model part are the same!
           failed_to_find_couple=True
           #print(node.Id)
           for node_fluid in self.fluid_model_part.Nodes:
                    coord_fluid_node = [ node_fluid.X , node_fluid.Y , node_fluid.Z ]
                    #now we compare:
                    sq_dist = 0.0
                    for iii in range(0,3):
                        sq_dist = sq_dist + (coord_wall_node[iii]-coord_fluid_node[iii])**2
                    if sq_dist<1.0e-9: #we found it!!
                        self.FEM_nodes_in_fluid_model_part.append(node_fluid)
                        failed_to_find_couple=False
                        break
           if (failed_to_find_couple):
               print("I couldnt find a match for node at coordinates ",coord_wall_node)
               print("this function must be called when the coupled nodes are still in the undeformed position!!")
               print(die)
       print("finished creating the coupling list between fluid nodes and wall(fem) nodes")

    def FindHighestNodeAndElementIds(self):
        self.highest_elem_id = 0
        self.highest_node_id = 0
        #we have to check ALL the model parts to find the highest id
        for node in itertools.chain( self.fluid_model_part.Nodes,
                                     self.disperse_phase_algorithm.spheres_model_part.Nodes,
                                     self.disperse_phase_algorithm.cluster_model_part.Nodes,
                                     self.disperse_phase_algorithm.rigid_face_model_part.Nodes,
                                     self.mixed_model_part.Nodes):
            if node.Id>self.highest_node_id:
               self.highest_node_id = node.Id

        for elem in itertools.chain( self.fluid_model_part.Elements,
                                     self.disperse_phase_algorithm.spheres_model_part.Elements,
                                     self.disperse_phase_algorithm.cluster_model_part.Elements,
                                     self.disperse_phase_algorithm.rigid_face_model_part.Elements,
                                     self.mixed_model_part.Elements):
            if elem.Id>self.highest_elem_id:
               self.highest_elem_id = elem.Id

        print("highest_elem_id = ", self.highest_elem_id)
        print("highest_node_id = ", self.highest_node_id)


    def ComputeTriangleAreaInXYPlane(self,elem):
        x0 = elem.GetNodes()[0].X
        y0 = elem.GetNodes()[0].Y
        x1 = elem.GetNodes()[1].X
        y1 = elem.GetNodes()[1].Y
        x2 = elem.GetNodes()[2].X
        y2 = elem.GetNodes()[2].Y
        
        x10 = x1 - x0;
        y10 = y1 - y0;
        x20 = x2 - x0;
        y20 = y2 - y0;
        detJ = x10 * y20-y10 * x20;
        area = (0.5*detJ)
        return area;
        
    def TransferThicknessFEMPartToFEMFluidAndBedDisplamecents(self): 
        #first from the walls to the fluid
        i = 0
        for nodeb in self.disperse_phase_algorithm.rigid_face_model_part.Nodes:
            thickess_in_walls = nodeb.GetSolutionStepValue(THICKNESS)
            mesh_displ_z_old = nodeb.GetSolutionStepValue(MESH_DISPLACEMENT_Z)
            mesh_displ_z = mesh_displ_z_old +  thickess_in_walls
            nodeb.SetSolutionStepValue(MESH_DISPLACEMENT_Z,mesh_displ_z)
            self.FEM_nodes_in_fluid_model_part[i].SetSolutionStepValue(MESH_DISPLACEMENT_Z,mesh_displ_z)
            #FOR POSTPROCESSING ONLY: DISPLACEMENT FROM THE INITIAL CONDITION:
            displ_z_old = nodeb.GetSolutionStepValue(DISPLACEMENT_Z)
            displ_z = displ_z_old +  thickess_in_walls
            nodeb.SetSolutionStepValue(DISPLACEMENT_Z,displ_z)
            i += 1
        #now fro the fluid to the bed
        i = 0
        for nodeb in self.bed_model_part.Nodes:
            mesh_displ_z = self.bed_nodes_in_fluid_model_part[i].GetSolutionStepValue(MESH_DISPLACEMENT_Z)
            nodeb.SetSolutionStepValue(HEIGHT,mesh_displ_z+0.5)
            nodeb.SetSolutionStepValue(MESH_DISPLACEMENT_Z,mesh_displ_z)
            i += 1

    def SetOutletPressure(self):         
        for node in self.fluid_model_part.Nodes:
           if node.IsFixed(PRESSURE):
               rel_z_coord = self.highest_z_coord - node.Z
               node.SetSolutionStepValue(PRESSURE,rel_z_coord * 9.8 * 1000.0 ) #WARNING, ASSUMING Z GRAVITY = 9.8 and water density = 1000 !!!!!!!
