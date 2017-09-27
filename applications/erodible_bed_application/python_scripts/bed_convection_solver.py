#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ErodibleBedApplication import *
#from KratosMultiphysics.ExternalSolversApplication import *

def AddVariables(model_part):  #this way er only need one command to add all the variables to our problem 
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(PROJECTED_HEIGHT);
    model_part.AddNodalSolutionStepVariable(DELTA_HEIGHT);
    model_part.AddNodalSolutionStepVariable(YP);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    #model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    #model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(MEAN_SIZE);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(SEDIMENT_VELOCITY);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(THICKNESS);
    model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);

def AddDofs(model_part):
    for node in model_part.Nodes:

        #adding dofs
        node.AddDof(HEIGHT);

    print ("variables for the Poisson solver added correctly")

class BedConvectionSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):  #constructor of the class 

        self.model_part = model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
        self.poisson_linear_solver =  SkylineLUFactorizationSolver()  #we set the type of solver that we want 
        
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-6,1e-9)  #tolerance for the solver 

        self.domain_size = domain_size
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        (self.neighbour_search).Execute()
        self.neighbour_elements_search= FindElementalNeighboursProcess(model_part,domain_size,number_of_avg_elems)
        (self.neighbour_elements_search).Execute()
        ##calculate normals
        self.normal_tools = BodyNormalCalculationUtils()
        self.normal_tools.CalculateBodyNormals(self.model_part,self.domain_size);  

    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = False
        import strategy_python
        self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.poisson_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
      
        self.VariableUtils = VariableUtils()

        maximum_number_of_particles= 8*self.domain_size
        self.moveparticles = MoveParticleUtilityBedTransport2D(self.model_part,maximum_number_of_particles)  
        self.moveparticles.MountBin()      
        self.erosion_utils = BedErosionUtility(self.model_part) 
        self.ExplicitStrategy = Erodible_Bed_Explicit_Strategy(self.model_part,self.domain_size, MoveMeshFlag) 
    #######################################################################   
    def Solve(self):
        (self.erosion_utils).CalculateSedimentVelocity();
        (self.moveparticles).CalculateVelOverElemSize();
        (self.moveparticles).MoveParticles();
        pre_minimum_number_of_particles=self.domain_size;
        (self.moveparticles).PreReseed(pre_minimum_number_of_particles);    
        (self.moveparticles).TransferLagrangianToEulerian();
        (self.VariableUtils).CopyScalarVar(PROJECTED_HEIGHT,HEIGHT,self.model_part.Nodes)   
        (self.moveparticles).ResetBoundaryConditions()
        (self.moveparticles).CopyScalarVarToPreviousTimeStep(HEIGHT,self.model_part.Nodes)  

        self.CalculateVelocityDiveregence()
        delta_time = self.model_part.ProcessInfo.GetValue(DELTA_TIME)
        for node in self.model_part.Nodes:
              node.SetSolutionStepValue(HEIGHT, (node.GetSolutionStepValue(HEIGHT) - delta_time * node.GetSolutionStepValue(HEIGHT) * node.GetSolutionStepValue(DIVPROJ)) ) 
        #(self.solver).Solve()

        (self.moveparticles).CalculateDeltaHeight();        
        (self.moveparticles).CorrectParticlesWithoutMovingUsingDeltaHeight();
        post_minimum_number_of_particles=self.domain_size*2;
        (self.moveparticles).PostReseed(post_minimum_number_of_particles);        

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

    #######################################################################   
    def CalculateVelocityDiveregence(self):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 0)
        (self.ExplicitStrategy).InitializeSolutionStep();
        (self.ExplicitStrategy).AssembleLoop();
        (self.ExplicitStrategy).FinalizeSolutionStep();
