from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return FluidSolver(model, custom_settings)

class FluidSolver(PythonSolver):

    def __init__(self, model, settings):

        super(FluidSolver,self).__init__(model, settings)

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = None
        self.condition_name = None
        self.min_buffer_size = 3

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    def AddVariables(self):
        raise Exception("Trying to call FluidSolver.AddVariables(). Implement the AddVariables() method in the specific derived solver.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("FluidSolver", "Fluid solver DOFs added correctly.")

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Set fluid properties from materials json file
            materials_imported = self._SetPhysicalProperties()
            if not materials_imported:
                KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "Material properties have not been imported. Check \'material_import_settings\' in your ProjectParameters.json.")
            ## Replace default elements and conditions
            self._ReplaceElementsAndConditions()
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)

        KratosMultiphysics.Logger.PrintInfo("FluidSolver", "Model reading finished.")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("FluidSolver", "Model export finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Initialize(self):
        raise Exception("Calling FluidSolver.Initialize() base method. Please implement a custom Initialize() method for your solver.")

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.InitializeSolutionStep()

    def Predict(self):
        if self._TimeBufferIsInitialized():
            self.solver.Predict()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.solver.SolveSolutionStep()
            if not is_converged:
                msg  = "Fluid solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
                msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
                KratosMultiphysics.Logger.PrintWarning("FluidSolver",msg)
            return is_converged
        else:
            return True

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            (self.solver).FinalizeSolutionStep()

    def Check(self):
        (self.solver).Check()

    def Clear(self):
        (self.solver).Clear()

    def GetComputingModelPart(self):
        if not self.main_model_part.HasSubModelPart("fluid_computational_model_part"):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    ## FluidSolver specific methods.

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.GetMinimumBufferSize()

    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self._GetElementNumNodes()
        cond_num_nodes = self._GetConditionNumNodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## If there are no elements and/or conditions, default to triangles/tetra meshes to avoid breaking the ReplaceElementsAndConditionsProcess
        ## This only affects the input name (if there are no elements or conditions to replace, nothing is replaced).
        if elem_num_nodes == 0:
            elem_num_nodes = domain_size + 1
        if cond_num_nodes == 0:
            cond_num_nodes = domain_size

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D" + str(int(elem_num_nodes)) + "N"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{}""")
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    def _GetElementNumNodes(self):
        if self.main_model_part.NumberOfElements() != 0:
            if sys.version_info[0] >= 3: # python3 syntax
                element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
            else: # python2 syntax
                element_num_nodes = len(self.main_model_part.Elements.__iter__().next().GetNodes())
        else:
            element_num_nodes = 0

        element_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(element_num_nodes)
        return element_num_nodes

    def _GetConditionNumNodes(self):
        if self.main_model_part.NumberOfConditions() != 0:
            if sys.version_info[0] >= 3: # python3 syntax
                condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())
            else: # python2 syntax
                condition_num_nodes = len(self.main_model_part.Conditions.__iter__().next().GetNodes())
        else:
            condition_num_nodes = 0

        condition_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(condition_num_nodes)
        return condition_num_nodes

    def _ExecuteCheckAndPrepare(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

    def _ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

        return delta_time

    def _GetAutomaticTimeSteppingUtility(self):
        if (self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.GetComputingModelPart(),
                                                                     self.settings["time_stepping"])
        else:
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.GetComputingModelPart(),
                                                                     self.settings["time_stepping"])

        return EstimateDeltaTimeUtility

    def _SetPhysicalProperties(self):
        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False

        # If the element uses nodal material properties, transfer them to the nodes
        if self.element_has_nodal_properties:
            self._SetNodalProperties()

        return materials_imported

    def _SetNodalProperties(self):
        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            kin_viscosity = dyn_viscosity / rho
            break
        else:
            raise Exception("No fluid elements found in the main model part.")
        # Transfer the obtained properties to the nodes
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
