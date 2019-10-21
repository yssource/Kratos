from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return CoupledPfemFluidThermalSolver(main_model_part, custom_settings)

class CoupledPfemFluidThermalSolver(PythonSolver):

    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass = True
        
        super(CoupledPfemFluidThermalSolver, self).__init__(model, custom_settings)

        #default_settings = KratosMultiphysics.Parameters()
        #    
        ##self.settings = custom_settings
        ### Overwrite the default settings with user-provided parameters
        #self.settings.ValidateAndAssignDefaults(default_settings)

        ## Get domain size
        self.domain_size = self.settings["fluid_solver_settings"]["domain_size"].GetInt()
        
        from KratosMultiphysics.PfemFluidDynamicsApplication import pfem_fluid_solver
        self.fluid_solver = pfem_fluid_solver.CreateSolver(self.model,self.settings["fluid_solver_settings"]) 

        from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "coupled_pfem_fluid_thermal_solver",
            "model_part_name": "PfemFluidModelPart",
            "echo_level"                         : 1,
            "fluid_solver_settings":{
                "physics_type"   : "fluid",
                "domain_size": 2,
                "time_stepping"               : {
                    "automatic_time_step" : false,
                    "time_step"           : 0.001
                },
                "model_import_settings":{
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                },
                "buffer_size": 3,
                "echo_level": 1,
                "reform_dofs_at_each_step": false,
                "clear_storage": false,
                "compute_reactions": true,
                "move_mesh_flag": true,
                "dofs"                : [],
                "stabilization_factor": 1.0,
                "line_search": false,
                "compute_contact_forces": false,
                "block_builder": false,
                "component_wise": false,
                "predictor_corrector": true,
                "time_order": 2,
                "maximum_velocity_iterations": 1,
                "maximum_pressure_iterations": 7,
                "velocity_tolerance": 1e-5,
                "pressure_tolerance": 1e-5,
                "pressure_linear_solver_settings":  {
                    "solver_type"                    : "amgcl",
                    "max_iteration"                  : 5000,
                    "tolerance"                      : 1e-9,
                    "provide_coordinates"            : false,
                    "scaling"                        : false,
                    "smoother_type"                  : "damped_jacobi",
                    "krylov_type"                    : "cg",
                    "coarsening_type"                : "aggregation",
                    "verbosity"                      : 0
                },
                "velocity_linear_solver_settings": {
                    "solver_type"                    : "bicgstab",
                    "max_iteration"                  : 5000,
                    "tolerance"                      : 1e-9,
                    "preconditioner_type"            : "none",
                    "scaling"                        : false
                },
                "solving_strategy_settings":{
                   "time_step_prediction_level": 0,
                   "max_delta_time": 1.0e-5,
                   "fraction_delta_time": 0.9,
                   "rayleigh_damping": false,
                   "rayleigh_alpha": 0.0,
                   "rayleigh_beta" : 0.0
                },
                "bodies_list": [],
                "problem_domain_sub_model_part_list": [],
                "processes_sub_model_part_list": [],
                "constraints_process_list": [],
                "loads_process_list"       : [],
                "output_process_list"      : [],
                "output_configuration"     : {},
                "problem_process_list"     : [],
                "processes"                : {},
                "output_processes"         : {},
                "check_process_list": []
            },
            "thermal_solver_settings": {
                "solver_type": "Transient",
                "analysis_type": "linear",
                "computing_model_part_name": "fluid_computing_domain",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                        "materials_filename": "ThermalMaterials.json"
                }
            },
            "coupling_settings": {}
        }""")

        this_defaults.AddMissingParameters(super(CoupledPfemFluidThermalSolver, cls).GetDefaultSettings())

        return this_defaults
        
    def AddVariables(self):
        # Import the fluid and thermal solver variables. Then merge them to have them in both fluid and thermal solvers.
        self.fluid_solver.AddVariables()
        self.thermal_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part)
        print("::[Coupled Pfem Fluid Thermal Solver]:: Variables MERGED")

    def ImportModelPart(self):
        # Call the fluid solver to import the model part from the mdpa
        self.fluid_solver.ImportModelPart() # import model fluid model part and call pfem_check_and_prepare_model_process_fluid
        
        self.ImportThermalProperties()
        #self.AddThermalNodes()
        #self.AddThermalElements()
        #self.thermal_solver._assign_nodally_properties()

        if (not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain")): 
            self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")

        #self.thermal_solver._set_and_fill_buffer()
        #self.thermal_solver._execute_after_reading()
        #NOT USED# self.fluid_solver._ImportModelPart(self.fluid_solver.main_model_part,self.settings["fluid_solver_settings"]["model_import_settings"])
        #NOT USED# Save the convection diffusion settings
        #NOT USED#convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        #NOT USED
        #NOT USED# Here the fluid model part is cloned to be thermal model part so that the nodes are shared
        #self.thermal_process_info = self.thermal_solver.main_model_part.ProcessInfo
        #modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        #if self.domain_size == 2:
        #    modeler.GenerateModelPart(self.fluid_solver.main_model_part,
        #    #modeler.GenerateModelPart(self.fluid_solver.GetComputingModelPart(),
        #                              self.thermal_solver.main_model_part,
        #                              "EulerianConvDiff2D")#,
        #                              #"ThermalFace2D2N")
        #else:
        #    modeler.GenerateModelPart(self.fluid_solver.main_model_part,
        #                              self.thermal_solver.main_model_part,
        #                              "EulerianConvDiff3D",
        #                              "ThermalFace3D3N")
        #self.thermal_solver.main_model_part.ProcessInfo = self.thermal_process_info
        #print("::[Coupled Pfem Fluid Thermal Solver]:: Thermal_model_part CLONED")
        #NOT USED
        #NOT USED# self.fluid_solver.ImportModelPart()
        #NOT USED
        #NOT USED# CheckAndPrepareModelProcess creates the thermal_computational model part
        #NOT USED#from KratosMultiphysics.ConvectionDiffusionApplication import check_and_prepare_model_process_convection_diffusion
        #NOT USED#check_and_prepare_model_process_convection_diffusion.CheckAndPrepareModelProcess(self.thermal_solver.main_model_part, params).Execute()
        #NOT USED
        #NOT USED# Set the saved convection diffusion settings to the new thermal model part
        #NOT USED#self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)
        #NOT USED
        #NOT USED# self.settings["thermal_solver_settings"]["computing_model_part_name"].SetString("fluid_computing_domain")
        #NOT USED#print(1)
        #if (not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain")): 
        #    self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")
        #print(1)

    def PrepareModelPart(self):
        self.fluid_solver.PrepareModelPart()
        self.thermal_solver.PrepareModelPart()
        #TODO: it seems to be not used

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.thermal_solver.Initialize()

    def InitializeStrategy(self):
        self.fluid_solver.InitializeStrategy()

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        self.thermal_solver.AdvanceInTime(current_time)
        #NOT USED#new_time_thermal = self.thermal_solver.AdvanceInTime(current_time)
        #NOT USED#print ("Fluid new_time is {} and thermal_new_time is {}".format(new_time, new_time_thermal))
        #NOT USED#TODO: allow different time steppings for the thermal and the fluid part (?)
        return new_time

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        self.thermal_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.thermal_solver.Predict()

    def SolveSolutionStep(self):
        # one_way_coupling --> fluid_to_thermal
        #TODO: density and viscosity temperature dependencies
        fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        
        # one_way_coupling --> fluid_to_thermal
        self.UpdateThermalVelocityField()
        thermal_is_converged = self.thermal_solver.SolveSolutionStep()



        #NOT USED #new_temperature = []
        #NOT USED for node in self.fluid_solver.main_model_part.Nodes:
        #NOT USED #    new_temperature.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        #NOT USED     print("Thermal node Id: {} temp: {}".format(node.Id, node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)))
        #NOT USED #count=0
        #NOT USED ##print(new_temperature)
        #NOT USED #for node in self.fluid_solver.main_model_part.Nodes:
        #NOT USED #    #print(node)
        #NOT USED #    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, new_temperature[count])
        #NOT USED #    #print("Fluid node Id: {} temp: {}".format(node.Id, node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)))
        #NOT USED #    count+=1
        #NOT USED #for elem in self.thermal_solver.main_model_part.Elements:
        #NOT USED #    print(elem.GetProperty)
        return (fluid_is_converged and thermal_is_converged)

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def Solve(self):
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()
    
    #def AddConvDiffElementsAndNodes(self):
        #if (not self.thermal_solver.model.HasModelPart("ThermalModelPart")): 
        #    self.thermal_solver.model.CreateModelPart("ThermalModelPart")
        #    #convection_diffusion_settings = self.fluid_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        #    #self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)
        #    self.thermal_solver.main_model_part.ProcessInfo=self.thermal_process_info
        #
        #if (not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain")): 
        #    convection_diffusion_computational_model_part = self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")
        # copying the nodes
        #thermal_computing_domain = self.thermal_solver.main_model_part.GetSubModelPart("thermal_computing_domain")
        #inlet2_model_part.CreateNewNode(4, 4.00,0.00,0.00)
        #for node in self.fluid_solver.main_model_part.Nodes:
        #    self.thermal_solver.main_model_part.AddNode(node,0)
            #thermal_computing_domain.CreateNewNode(node.Id,node.X, node.Y, node.Z)
            #print(node.Id)
            #print(self.fluid_solver.main_model_part.Nodes[node.Id])
        #for node in self.thermal_solver.main_model_part.Nodes:
        #    print(node.Id)
        #    print(self.thermal_solver.main_model_part.Nodes[node.Id])

        # Import ConvDiff constitutive laws.
        #if not hasattr(self,'ConvDiffProperties'):#(not self.thermal_solver.main_model_part.HasProperties(0)):
        #    materials_imported = self.thermal_solver.import_materials()
        #    if materials_imported:
        #        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "ConvDiff materials were successfully imported.")
        #        self.ConvDiffProperties = self.thermal_solver.main_model_part.Properties[0]
        #    else:
        #        KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "ConvDiff materials were not imported.")
            #TODO: number of properties to be automatic 
            #ConvDiffProperties = self.thermal_solver.main_model_part.Properties[0]

         # creating new thermal elements
        #for elem in self.fluid_solver.main_model_part.Elements:
            #node_ids = []
            #if self.domain_size == 2:
            #    if len(elem.GetNodes()) == 3:
            #        #for i in range(len(elem.GetNodes())):
            #        node_ids = [elem.GetNode(0).Id, elem.GetNode(1).Id, elem.GetNode(2).Id]
            #        #print(elem.Id)
            #        #print(node_ids)
            #        self.thermal_solver.main_model_part.CreateNewElement("EulerianConvDiff2D", elem.Id, node_ids, self.thermal_solver.main_model_part.Properties[0])
            #else:
            #    node_ids = [elem.GetNode(0).Id, elem.GetNode(1).Id, elem.GetNode(2).Id, elem.GetNode(3).Id]
            #    self.thermal_solver.main_model_part.CreateNewElement("EulerianConvDiff3D", elem.Id, node_ids, self.thermal_solver.main_model_part.Properties[0])
        #for elem in self.thermal_solver.main_model_part.Elements:
        #    print(elem)
        #transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(convection_diffusion_computational_model_part, self.thermal_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS)
        #transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.thermal_solver.main_model_part.GetSubModelPart("thermal_computing_domain"), self.thermal_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS)
        #transfer_process.Execute()
        #KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part)
        #print(1)

    #def DeleteConvDiffElementsAndNodes(self):
        #self.DeleteThermalModelPartElements()
        #self.DeleteThermalModelPartNodes()
        #print(1)
    def UpdateThermalVelocityField(self):
        #print("fluid_model_part velocity field")
        #for node in self.fluid_solver.main_model_part.Nodes:
        #    print(node.Id)
        #    print(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))      
        #print("thermal_model_part velocity field")
        #for node in self.thermal_solver.main_model_part.Nodes:
        #    print(node.Id)
        #    print(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
        #
        #print("fluid_model_part temperature field")
        #for node in self.fluid_solver.main_model_part.Nodes:
        #    print(node.Id)
        #    print(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))      
        #print("thermal_model_part temperature field")
        #for node in self.thermal_solver.main_model_part.Nodes:
        #    #print(node.Id)
        #    print("node Id: {} temp: {}".format(node.Id, node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)))
        

        #print("fluid_model_part mesh velocity field")
        #for node in self.fluid_solver.main_model_part.Nodes:
        #    print(node.Id)
        #    print(node.GetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_X))

        #Update mesh velocity
        new_mesh_velocity_x = []
        new_mesh_velocity_y = []
        #count=0
        for node in self.fluid_solver.main_model_part.Nodes:
            new_mesh_velocity_x.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
            #node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_X, new_mesh_velocity_x)
            new_mesh_velocity_y.append(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
            #node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_Y, new_mesh_velocity_y)
            #print("Fluid node Id: {} VEL_X: {}".format(node.Id, node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)))
            #count+=1
        count=0
        for node in self.thermal_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_X, new_mesh_velocity_x[count])
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_Y, new_mesh_velocity_y[count])
            #print("Thermal node Id: {} MESH_VEL_X: {}".format(node.Id, node.GetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_X)))
            count+=1
        #print(new_mesh_velocity_y)
        #print("thermal_model_part UPDATED mesh velocity field")
        #for node in self.thermal_solver.main_model_part.Nodes:
        #    print(node.Id)
        #    print(node.GetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_X))

        print("end")
    #def CloneFluidComputingModelPart(self):

        #if (not self.thermal_solver.model.HasModelPart("ThermalModelPart")): 
        #    self.thermal_solver.model.CreateModelPart("ThermalModelPart")
            #convection_diffusion_settings = self.fluid_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
            #self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)
            #self.fluid_solver.main_model_part.ProcessInfo.SetValue()
            

        #if (not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain")): 
            #convection_diffusion_computational_model_part = self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")
        #self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")
        #
        #modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        #if self.domain_size == 2:
        #    modeler.GenerateModelPart(self.fluid_solver.GetComputingModelPart(),
        #                              self.thermal_solver.main_model_part,
        #                              "EulerianConvDiff2D")
        #else:
        #    modeler.GenerateModelPart(self.fluid_solver.main_model_part,
        #                              self.thermal_solver.main_model_part,
        #                              "EulerianConvDiff3D",
        #                              "ThermalFace3D3N")
        #print("::[Coupled Pfem Fluid Thermal Solver]:: Thermal_model_part CLONED")
        #self.thermal_solver.main_model_part.ProcessInfo=self.thermal_process_info
        #transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.thermal_solver.GetComputingModelPart(), self.thermal_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS)
        #transfer_process.Execute()
        #print(1)
    def DeleteThermalElements(self):
        for elem in self.thermal_solver.main_model_part.Elements:
            elem.Set(KratosMultiphysics.TO_ERASE, True)
        self.thermal_solver.main_model_part.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)
    def DeleteThermalNodes(self):
        for node in self.thermal_solver.main_model_part.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, True)
        self.thermal_solver.main_model_part.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
        #print("end deleting nodes")
        #for node in self.fluid_solver.main_model_part.Nodes:
        #    print(node)
    def ImportThermalProperties(self):
        # Import ConvDiff constitutive laws.
        materials_imported = self.thermal_solver.import_materials()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "ConvDiff materials were successfully imported.")
            #self.ConvDiffProperties = self.thermal_solver.main_model_part.Properties[0]
        else:
            KratosMultiphysics.Logger.PrintInfo("::[ConvectionDiffusionBaseSolver]:: ", "ConvDiff materials were not imported.")
    #def DeleteThermaModelPart(self):    
        #self.thermal_process_info = self.thermal_solver.main_model_part.ProcessInfo
        #self.thermal_solver.model.DeleteModelPart("ThermalModelPart")

    def AddThermalElements(self):
        for elem in self.fluid_solver.main_model_part.Elements:
            node_ids = []
            if self.domain_size == 2:
                if len(elem.GetNodes()) == 3:
                    node_ids = [elem.GetNode(0).Id, elem.GetNode(1).Id, elem.GetNode(2).Id]                    
                    self.thermal_solver.main_model_part.CreateNewElement("EulerianConvDiff2D", elem.Id, node_ids, self.thermal_solver.main_model_part.Properties[0])
                    #self.thermal_solver.main_model_part.AddElement(elem,0)
            else:
                print("3D case to be implemented")
                #TODO: complete 3D case...
        
        # transfer elements to thermal_computing_domain
        #transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.thermal_solver.main_model_part.GetSubModelPart("thermal_computing_domain"), self.thermal_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDELEMENTS)
        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.thermal_solver.main_model_part.GetSubModelPart("thermal_computing_domain"), self.thermal_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.ELEMENTS)
        transfer_process.Execute()

    def AddThermalNodes(self):
        if (not self.thermal_solver.model.HasModelPart("ThermalModelPart")): 
            self.thermal_solver.model.CreateModelPart("ThermalModelPart")

        if (not self.thermal_solver.main_model_part.HasSubModelPart("thermal_computing_domain")): 
            convection_diffusion_computational_model_part = self.thermal_solver.main_model_part.CreateSubModelPart("thermal_computing_domain")
            convection_diffusion_computational_model_part.ProcessInfo = self.thermal_solver.main_model_part.ProcessInfo
            convection_diffusion_computational_model_part.Properties  = self.thermal_solver.main_model_part.Properties
        # copying the nodes
        for node in self.fluid_solver.main_model_part.Nodes:
            self.thermal_solver.main_model_part.AddNode(node,0)
        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.thermal_solver.main_model_part.GetSubModelPart("thermal_computing_domain"), self.thermal_solver.main_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES)
        transfer_process.Execute()
        #print(1)