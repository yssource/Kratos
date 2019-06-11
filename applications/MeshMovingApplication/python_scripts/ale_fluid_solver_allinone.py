from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.MeshMovingApplication.python_solvers_wrapper_mesh_motion as mesh_mothion_solvers_wrapper
import KratosMultiphysics.FluidDynamicsApplication.python_solvers_wrapper_fluid as fluid_solvers_wrapper

class AleFluidSolver(PythonSolver):
    def __init__(self, model, solver_settings, parallelism):

        self.settings.ValidateAndAssignDefaults(KratosMultiphysics.Parameters("""{
                                                        "echo_level" : 0,
                                                        "solver_type"                 : "ale_fluid",
                                                        "start_fluid_solution_time"   : 0.0,
                                                        "ale_boundary_parts"          : [ ],
                                                        "mesh_motion_parts"           : [ ],
                                                        "fluid_solver_settings"       : { },
                                                        "mesh_motion_solver_settings" : { },
                                                        "mesh_velocity_calculation"   : { }
                                                        }"""))
        self.start_fluid_solution_time = self.settings["start_fluid_solution_time"].GetDouble()


        fluid_model_part_name = self.settings["fluid_solver_settings"]["model_part_name"].GetString()
        if not self.model.HasModelPart(fluid_model_part_name):
            self.model.CreateModelPart(fluid_model_part_name)

        ## Checking if reactions are being computed in the fluid
        if self.settings["fluid_solver_settings"].Has("compute_reactions"):
            if self.settings["fluid_solver_settings"]["compute_reactions"].GetBool() == False:
                self.settings["fluid_solver_settings"]["compute_reactions"].SetBool(True)
                warn_msg  = '"compute_reactions" is switched off for the fluid-solver, '
                warn_msg += 'switching it on!'
                KratosMultiphysics.Logger.PrintWarning("::[AleFluidSolver]::", warn_msg)
        else:
            self.settings["fluid_solver_settings"].AddEmptyValue("compute_reactions").SetBool(True)
            info_msg = 'Setting "compute_reactions" to true for the fluid-solver'
            KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", info_msg)

        ## Creating the fluid solver
        self.fluid_solver = fluid_solvers_wrapper.CreateSolverByParameters(
            self.model, self.settings["fluid_solver_settings"], self.parallelism)

        # Doing this after the Fluid-solver-settings have been validated to access the settings
        self._SelectMeshVelocityCalculationSettings()

        self.__InitializeMeshVelocityComputation()

        ## Creating the mesh-motion solver
        # Making sure the settings are consistent btw fluid and mesh-motion
        if self.settings["mesh_motion_solver_settings"].Has("model_part_name"):
            if not fluid_model_part_name == self.settings["mesh_motion_solver_settings"]["model_part_name"].GetString():
                err_msg =  'Fluid- and Mesh-Solver have to use the same MainModelPart ("model_part_name")!\n'
                err_msg += 'Use "mesh_motion_parts" for specifying mesh-motion on sub-model-parts'
                raise Exception(err_msg)
        else:
            self.settings["mesh_motion_solver_settings"].AddValue("model_part_name", self.settings["fluid_solver_settings"]["model_part_name"])

        domain_size = self.settings["fluid_solver_settings"]["domain_size"].GetInt()
        if self.settings["mesh_motion_solver_settings"].Has("domain_size"):
            mesh_motion_domain_size = self.settings["mesh_motion_solver_settings"]["domain_size"].GetInt()
            if not domain_size == mesh_motion_domain_size:
                raise Exception('Fluid- and Mesh-Solver have to use the same "domain_size"!')
        else:
            self.settings["mesh_motion_solver_settings"].AddValue("domain_size", self.settings["fluid_solver_settings"]["domain_size"])

        # Constructing the mesh-solver with the entire mesh
        # if no submodelparts are specified then this is used for the computation of the mesh-motion
        # otherwise it only adds the dofs and the variables (to the entire ModelPart!)
        self.mesh_motion_solver_full_mesh = mesh_mothion_solvers_wrapper.CreateSolverByParameters(
            self.model, self.settings["mesh_motion_solver_settings"], self.parallelism)

        # Getting the min_buffer_size from both solvers
        # and assigning it to the fluid_solver, bcs this one handles the model_part
        self.fluid_solver.min_buffer_size = max(
            [ self.fluid_solver.GetMinimumBufferSize(),
              self.mesh_motion_solver_full_mesh.GetMinimumBufferSize(),
              KratosMultiphysics.GetMinimumBufferSize(self.time_int_helper) ] )

        KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "echo_level"                  : 0,
            "solver_type"                 : "ale_fluid",
            "start_fluid_solution_time"   : 0.0,
            "ale_boundary_parts"          : [ ],
            "mesh_motion_parts"           : [ ],
            "fluid_solver_settings"       : { },
            "mesh_motion_solver_settings" : { },
            "mesh_velocity_calculation"   : { }
        }""")
        return this_defaults

    def AddVariables(self):
        self.mesh_motion_solver_full_mesh.AddVariables()
        self.fluid_solver.AddVariables()

        # Adding Variables used for computation of Mesh-Velocity
        time_scheme = self.settings["mesh_velocity_calculation"]["time_scheme"].GetString()
        main_model_part = self.model[self.settings["fluid_solver_settings"]["model_part_name"].GetString()]
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        if not time_scheme.startswith("bdf"): # bdfx does not need MESH_ACCELERATION
            main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_ACCELERATION)

        KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", "Variables Added")

    def AddDofs(self):
        self.mesh_motion_solver_full_mesh.AddDofs()
        self.fluid_solver.AddDofs()
        KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", "DOFs Added")

    def Initialize(self):
        # Saving the ALE-interface-parts for later
        # this can only be done AFTER reading the ModelPart
        main_model_part_name = self.settings["fluid_solver_settings"]["model_part_name"].GetString()

        ale_boundary_parts_params = self.settings["ale_boundary_parts"]
        self.ale_boundary_parts = []
        for i_name in range(ale_boundary_parts_params.size()):
            sub_model_part_name = ale_boundary_parts_params[i_name].GetString()
            full_model_part_name = main_model_part_name + "." + sub_model_part_name
            self.ale_boundary_parts.append(self.model[full_model_part_name])

        mesh_motion_parts_params = self.settings["mesh_motion_parts"]
        self.mesh_motion_solvers = []
        if mesh_motion_parts_params.size() == 0:
            # the entire Fluid-ModelPart is used in the Mesh-Solver
            self.mesh_motion_solvers.append(self.mesh_motion_solver_full_mesh)
        else:
            # SubModelParts of the Fluid-ModelPart are used in the Mesh-Solver
            # each SubModelPart has its own mesh-solver
            # Note that these solvers do NOT need to call AddVariables and AddDofs
            # since this is done already for the MainModelPart
            for i_name in range(mesh_motion_parts_params.size()):
                sub_model_part_name = mesh_motion_parts_params[i_name].GetString()
                if sub_model_part_name == main_model_part_name:
                    err_msg =  'The MainModelPart cannot be used as one of the Sub-Mesh-Solvers!\n'
                    err_msg += 'Remove "mesh_motion_parts" for specifying mesh-motion on the MainModelPart'
                    raise Exception(err_msg)
                full_model_part_name = main_model_part_name + "." + sub_model_part_name
                sub_mesh_solver_settings = self.settings["mesh_motion_solver_settings"].Clone()
                sub_mesh_solver_settings["model_part_name"].SetString(full_model_part_name)

                self.mesh_motion_solvers.append(mesh_mothion_solvers_wrapper.CreateSolverByParameters(
                    self.model, sub_mesh_solver_settings, self.parallelism))

        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Initialize()
        self.fluid_solver.Initialize()

        KratosMultiphysics.Logger.PrintInfo("::[AleFluidSolver]::", "Finished initialization")

    def ImportModelPart(self):
        self.fluid_solver.ImportModelPart() # only ONE mesh_solver imports the ModelPart

    def PrepareModelPart(self):
        # Doing it ONLY for the fluid solver (since this contains filling the buffer)
        self.fluid_solver.PrepareModelPart()

    def AdvanceInTime(self, current_time):
        # Doing it ONLY for the fluid solver
        return self.fluid_solver.AdvanceInTime(current_time)

    def Finalize(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Finalize()
        self.fluid_solver.Finalize()

    def InitializeSolutionStep(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.InitializeSolutionStep()
        self.fluid_solver.InitializeSolutionStep()

    def Predict(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Predict()
        self.fluid_solver.Predict()

    def FinalizeSolutionStep(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.FinalizeSolutionStep()
        self.fluid_solver.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        is_converged = True
        for mesh_solver in self.mesh_motion_solvers:
            is_converged &= mesh_solver.SolveSolutionStep()

        for mesh_solver in self.mesh_motion_solvers:
            KratosMultiphysics.MeshMovingApplication.CalculateMeshVelocities(
                mesh_solver.GetComputingModelPart(),
                self.time_int_helper)

        if self.fluid_solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] >= self.start_fluid_solution_time:
            for mp in self.ale_boundary_parts:
                KratosMultiphysics.VariableUtils().CopyVectorVar(
                    KratosMultiphysics.MESH_VELOCITY,
                    KratosMultiphysics.VELOCITY,
                    mp.GetCommunicator().LocalMesh().Nodes)
                mp.GetCommunicator().SynchronizeVariable(KratosMultiphysics.VELOCITY)
            is_converged &= self.fluid_solver.SolveSolutionStep()

        return is_converged

    def Check(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Check()
        self.fluid_solver.Check()

    def Clear(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.Clear()
        self.fluid_solver.Clear()

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart() # this is the same as the one used in the MeshSolver

    def GetFluidSolver(self):
        return self.fluid_solver

    def GetMeshMotionSolver(self):
        if len(self.mesh_motion_solvers) > 1:
            raise Exception('More than one mesh-motion-solver \
                exists, please use "GetMeshMotionSolvers"')
        return self.mesh_motion_solvers[0]

    def GetMeshMotionSolvers(self):
        return self.mesh_motion_solvers

    def MoveMesh(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.MoveMesh()


    def __InitializeMeshVelocityComputation(self):
        '''Initializing the helper-class for the time-integration
        '''
        time_int_settings = self.settings["mesh_velocity_calculation"]
        time_scheme = time_int_settings["time_scheme"].GetString()

        if time_scheme == "bdf1":
            self.time_int_helper = KratosMultiphysics.BDF1()
        elif time_scheme == "bdf2":
            self.time_int_helper = KratosMultiphysics.BDF2()
        elif time_scheme == "newmark":
            self.time_int_helper = KratosMultiphysics.Newmark()
        elif time_scheme == "bossak":
            if time_int_settings.Has("alpha_m"):
                alpha_m = time_int_settings["alpha_m"].GetDouble()
                self.time_int_helper = KratosMultiphysics.Bossak(alpha_m)
            else:
                self.time_int_helper = KratosMultiphysics.Bossak()
        elif time_scheme == "generalized_alpha":
            alpha_m = time_int_settings["alpha_m"].GetDouble()
            alpha_f = time_int_settings["alpha_f"].GetDouble()
            self.time_int_helper = KratosMultiphysics.GeneralizedAlpha(alpha_m, alpha_f)
        else:
            err_msg =  'The requested time scheme "' + time_scheme + '" is not available!\n'
            err_msg += 'Available options are: "bdf1", "bdf2", '
            err_msg += '"newmark", "bossak", "generalized_alpha"'
            raise Exception(err_msg)

    def _SelectMeshVelocityCalculationSettings(self):
        '''Selecting the time-integration for the MESH_VELOCITY to be consistent
        with the fluid time-integration.
        By now the parameters of the FluidSolver have been validated, which means
        that the time-integration method used by the fluid can be queried
        '''
        mesh_vel_calc_settings = self.settings["mesh_velocity_calculation"]
        fluid_settings = self.settings["fluid_solver_settings"]

        fluid_solver_type = fluid_settings["solver_type"].GetString()
        if fluid_solver_type == "monolithic" or fluid_solver_type == "Monolithic":
            if fluid_settings.Has("time_scheme"):
                time_scheme_fluid = fluid_settings["time_scheme"].GetString()
                alpha_fluid = fluid_settings["alpha"].GetDouble()
                if mesh_vel_calc_settings.Has("time_scheme"):
                    time_scheme_mesh_vel = mesh_vel_calc_settings["time_scheme"].GetString()
                    if time_scheme_fluid != time_scheme_mesh_vel:
                        info_msg  = '"time_scheme" of the fluid (' + time_scheme_fluid
                        info_msg += ') is different from the\n"time_scheme" used for the '
                        info_msg += 'calculation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                        KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                else:
                    mesh_vel_calc_settings.AddValue("time_scheme", fluid_settings["time_scheme"])
                    info_msg  = 'setting "time_scheme" of the mesh-solver for the\ncalculation of the '
                    info_msg += 'mesh-velocity to "' + time_scheme_fluid + '" to be consistent with the\n'
                    info_msg += '"time_scheme" of the fluid'
                    KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                if not mesh_vel_calc_settings["time_scheme"].GetString().startswith("bdf"):
                    if mesh_vel_calc_settings.Has("alpha_m"):
                        alpha_mesh_vel = mesh_vel_calc_settings["alpha_m"].GetDouble()
                        if abs(alpha_fluid-alpha_mesh_vel) > 1e-12:
                            info_msg  = '"alpha" of the fluid (' + str(alpha_fluid)
                            info_msg += ') is different from the\n"alpha_m" used for the '
                            info_msg += 'calculation of the mesh-velocity (' + str(alpha_mesh_vel) + ')'
                            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
                    else:
                        mesh_vel_calc_settings.AddValue("alpha_m", fluid_settings["alpha"])
                        info_msg  = 'setting "alpha_m" of the mesh-solver for the\ncalculation of the '
                        info_msg += 'mesh-velocity to "' + str(alpha_fluid) + '" to be consistent\nwith '
                        info_msg += '"alpha" of the fluid'
                        KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        elif fluid_solver_type == "fractional_step" or fluid_solver_type == "FractionalStep":
            # currently fractional step always uses BDF2
            if mesh_vel_calc_settings.Has("time_scheme"):
                time_scheme_mesh_vel = mesh_vel_calc_settings["time_scheme"].GetString()
                if time_scheme_mesh_vel != "bdf2":
                    info_msg  = '"time_scheme" of the fluid (bdf2) '
                    info_msg += 'is different from the\n"time_scheme" used for the '
                    info_msg += 'calculation of the mesh-velocity (' + time_scheme_mesh_vel + ')'
                    KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)
            else:
                mesh_vel_calc_settings.AddEmptyValue("time_scheme").SetString("bdf2")
                info_msg  = 'setting "time_scheme" of the mesh-solver for the\ncalculation of the '
                info_msg += 'mesh-velocity to "bdf2" to be consistent with the\n'
                info_msg += '"time_scheme" of the fluid'
                KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        if not mesh_vel_calc_settings.Has("time_scheme"):
            mesh_vel_calc_settings.AddEmptyValue("time_scheme").SetString("UNSPECIFIED")
            warn_msg  = 'unknown "solver_type" of the fluid-solver, therefore '
            warn_msg += 'no automatic selection\nof "time_scheme" for the calculation '
            warn_msg += 'of the mesh-velocity performed'
            KratosMultiphysics.Logger.PrintWarning("::[ALEFluidSolver]::", warn_msg)
