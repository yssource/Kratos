from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.CompressiblePotentialFlowApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return AdjointPotentialSolver(main_model_part, custom_settings)

class AdjointPotentialSolver:

    def __init__(self, main_model_part, custom_settings):
        self.MoveMeshFlag = False

        self.main_model_part = main_model_part

        # default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "adjoint_potential_solver",
            "scheme_settings" : {
                "scheme_type" : "steady"
            },
            "response_function_settings" : {
                "response_type" : "lift"
            },
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "linear_solver_settings" : {
                "solver_type" : "AMGCL"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts"  : [""],
            "no_skin_parts":[""],
            "echo_level"  : 0
        }""")

        # overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of AdjointPotentialSolver finished")

    def GetMinimumBufferSize(self):
        return 1

    def AddVariables(self):

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        
        print("variables for the adjoint fluid solver added correctly")

    def AddDofs(self):

        print("ADDING DOFS")
        
        for node in self.main_model_part.Nodes:
            node.AddDof(KratosMultiphysics.ADJOINT_NEGATIVE_FACE_PRESSURE)
            node.AddDof(KratosMultiphysics.ADJOINT_POSITIVE_FACE_PRESSURE)


        print("DOFs for the potential adjoint solver added correctly.")  
    
    
    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # here it would be the place to import restart data if required
            print(self.settings["model_import_settings"]["input_filename"].GetString())
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            # here we shall check that the input read has the shape we like
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            aux_params.AddValue("skin_parts",self.settings["skin_parts"])

            throw_errors = False
            KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part,throw_errors).Execute()
            # here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name": "AdjointCompressiblePotentialFlowElement3D4N",
                    "condition_name": "AdjointPotentialWallCondition3D3N"
                    }
                    """)
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name": "AdjointCompressiblePotentialFlowElement2D3N",
                    "condition_name": "AdjointPotentialWallCondition2D2N"
                    }
                    """)
            else:
                raise Exception("domain size is not 2 or 3")

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

            #import check_and_prepare_model_process_fluid
            #check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

            #here we read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
            #for el in self.main_model_part.Elements:
            #    rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            #    break

            #KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)

        else:
            raise Exception("Other input options are not yet implemented.")

        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("Model reading finished.")


    


    def Initialize(self):

        #self.computing_model_part = self.GetComputingModelPart()

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.settings["response_function_settings"]["response_type"].GetString() == "lift":
            if (domain_size == 2):
                self.response_function = AdjointFluidApplication.PotentialFlowLiftResponseFunction2D(self.main_model_part, self.settings["response_function_settings"])
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        else:
            raise Exception("invalid response_type: " + self.settings["response_function_settings"]["response_type"].GetString())

        if self.settings["scheme_settings"]["scheme_type"].GetString() == "steady":
            self.time_scheme = AdjointFluidApplication.AdjointSteadyPotentialScheme(self.settings["scheme_settings"], self.response_function)
        else:
            raise Exception("invalid scheme_type: " + self.settings["scheme_settings"]["scheme_type"].GetString())

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.main_model_part,
                                                                     self.time_scheme,
                                                                     self.linear_solver,
                                                                     builder_and_solver,
                                                                     False,
                                                                     False,
                                                                     False,
                                                                     False)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Check()       

        print ("Adjoint potential solver initialization finished.")

    def GetComputingModelPart(self):
        # get the submodelpart generated in CheckAndPrepareModelProcess
        #return self.main_model_part.GetSubModelPart("fluid_computational_model_part")
        return self.main_model_part

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def DivergenceClearance(self):
        pass
        
    def SolverInitialize(self):
        self.solver.Initialize()

    def SolverInitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()

    def SolverPredict(self):
        self.solver.Predict()

    def SolverSolveSolutionStep(self):
        self.solver.SolveSolutionStep()

    def SolverFinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

    def Solve(self):
        #for elem in self.model_part.Elements:
        #    if(elem is wake):
        #        elem.Set(ACTIVE,False)
        self.solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()
