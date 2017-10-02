# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import json as json

# ======================================================================================================================================
# Class definition
# ======================================================================================================================================

class CADReconstrutionUtilities():

    # --------------------------------------------------------------------------
    def __init__( self, fem_filename, cad_geometry_filename, cad_integration_data_filename ):
        self.FEMInputFilename = fem_filename
        self.CADGeometryFilename = cad_geometry_filename
        self.CADIntegrationDataFilename = cad_integration_data_filename

        # Internal parameters to specify reconstruction method
        
        # Gernal strategy parameters
        self.ReconstructionStrategy = "mapping"    # mapping / least_squares
        self.FEMGaussIntegrationDegree = 5

        self.SolutionStrategy = "penalty_approach" # penalty_approach 
        self.SolutionIterations = 1
        self.PenaltyFactorForDisplacementCoupling = 1e3
        self.PenaltyFactorForDirichletConstraints = 1e3
        self.PenaltyIncreaseFactor = 10
       
        # Parameters to edit input data       
        self.FERefinementLevel = 0

        # Projection settings
        self.ParameterResolutionForProjection = [ 100, 100 ]
        self.MaxProjectionIterations = 20
        self.ProjectionTolerance = 1e-5

        # Linear solver
        self.LinearSolver = SuperLUSolver()
        # DiagPrecond = DiagonalPreconditioner()
        # self.LinearSolver =  BICGSTABSolver(1e-9, 5000, DiagPrecond)
        # self.LinearSolver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, 1e-9, 300, 2, 10)        

        # Parameters for output
        self.ParameterResolutionForOutputOfSurfacePoints = [ 50, 50 ]

    # --------------------------------------------------------------------------
    def Initialize( self ):
        self.__ReadFEData()
        self.__RefineFEModel()
        self.__ReadCADData()
        self.__InitializeConstraints()
        self.__CreateReconstructionDataBase()
        self.__CreateReconstructionOutputWriter()

    # --------------------------------------------------------------------------
    def __ReadFEData( self ):
        self.FEModelPart = ModelPart("name_of_empty_mdpa")
        self.FEModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
        model_part_io = ModelPartIO(self.FEMInputFilename)
        model_part_io.ReadModelPart(self.FEModelPart)

    # --------------------------------------------------------------------------
    def __RefineFEModel( self ):       

        # Assign pseudo material to elements (required by refinement)
        prop_id = 1
        prop = self.FEModelPart.Properties[prop_id]
        mat = LinearElasticPlaneStress2DLaw()
        prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

        for refinement_level in range(0,self.FERefinementLevel):

            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            nodal_neighbour_search = FindNodalNeighboursProcess(self.FEModelPart, number_of_avg_elems, number_of_avg_nodes)
            neighbour_calculator = FindElementalNeighboursProcess(self.FEModelPart,2,10)
            nodal_neighbour_search.Execute()
            neighbour_calculator.Execute()

            for elem in self.FEModelPart.Elements:
                elem.SetValue(SPLIT_ELEMENT,True)

            refine_on_reference = False
            interpolate_internal_variables = True
            Refine = LocalRefineTriangleMesh(self.FEModelPart)
            Refine.LocalRefineMesh(refine_on_reference, interpolate_internal_variables)

    # --------------------------------------------------------------------------
    def __ReadCADData( self ):
        self.CADGeometry = {}
        with open(self.CADGeometryFilename) as cad_data1:
            self.CADGeometry = json.load(cad_data1)
        self.CADIntegrationData = {}
        with open(self.CADIntegrationDataFilename) as cad_data2:
            self.CADIntegrationData = json.load(cad_data2)       

    # --------------------------------------------------------------------------
    def __InitializeConstraints( self ):
        self.IsDisplacementCouplingSpecifiedForAllCouplingPoints = False
        self.AreDirichletConstraintsSpecified = False
    
    # --------------------------------------------------------------------------
    def __CreateReconstructionDataBase( self ):
        self.DataBase = ReconstructionDataBase(self.FEModelPart, self.CADGeometry, self.CADIntegrationData)
        self.DataBase.Create()

    # --------------------------------------------------------------------------
    def __CreateReconstructionOutputWriter( self ):
        self.OutputWriter = ReconstructionOutputWriter( self.DataBase )    

    # --------------------------------------------------------------------------
    def SetDisplacementCouplingOnAllCouplingPoints( self ):
        self.IsDisplacementCouplingSpecifiedForAllCouplingPoints = True

    # --------------------------------------------------------------------------
    def SetDirichletConstraints( self, list_of_condition_settings ):
        self.AreDirichletConstraintsSpecified = True
        self.DirichletConditions = list_of_condition_settings

    # --------------------------------------------------------------------------
    def PerformReconstruction( self ):
        self.__CreateReconstructionConditions()
        self.__CreateSolverForReconstruction()
        self.__SolveReconstructionEquation()

    # --------------------------------------------------------------------------
    def __CreateReconstructionConditions( self ):
        self.ConditionsContainer = ReconstructionConditionContainer( self.DataBase )

        if self.ReconstructionStrategy == "mapping":
            self.ConditionsContainer.CreateDisplacementMappingConditions( self.ParameterResolutionForProjection, 
                                                                          self.FEMGaussIntegrationDegree,
                                                                          self.MaxProjectionIterations,
                                                                          self.ProjectionTolerance )
        else:
            raise ValueError( "The following reconstruction strategy does not exist: ", self.ReconstructionStrategy )

        if self.IsDisplacementCouplingSpecifiedForAllCouplingPoints: 
            self.ConditionsContainer.CreateDisplacementCouplingConstraintsOnAllCouplingPoints( self.PenaltyFactorForDisplacementCoupling )

        if self.AreDirichletConstraintsSpecified:
            self.ConditionsContainer.CreateDirichletConstraints( self.DirichletConstraints, self.PenaltyFactorForDirichletConstraints )

    # --------------------------------------------------------------------------
    def __CreateSolverForReconstruction( self ):
        self.ReconstructionSolver = CADReconstructionSolver( self.DataBase, self.ConditionsContainer, self.LinearSolver )

    # --------------------------------------------------------------------------
    def __SolveReconstructionEquation( self ): 
        if(self.SolutionStrategy == "penalty_approach"):
            self.__SolveReconstructionEquationUsingPenaltyApproach()
        else:
            raise ValueError( "The following solution strategy does not exist: ", self.SolutionStrategy )
        
    # --------------------------------------------------------------------------
    def __SolveReconstructionEquationUsingPenaltyApproach( self ):  
        
        self.ReconstructionSolver.InitializeEquationSystem()            
        
        for iteration in range(1,self.SolutionIterations+1): 
            self.ReconstructionSolver.ComputeLHS()
            self.ReconstructionSolver.ComputeRHS()
            self.ReconstructionSolver.RegularizeEquationSystem()
            self.ReconstructionSolver.SolveEquationSystem()
            if self.SolutionIterations>1:  
                self.ReconstructionSolver.IncreaseAllPenaltyFactorsByInputFactor( self.PenaltyIncreaseFactor )

    # --------------------------------------------------------------------------
    def OutputFEData( self ):
        from gid_output import GiDOutput
        fem_output_filename = self.FEMInputFilename+"_as_used_for_reconstruction"
        nodal_results=["SHAPE_CHANGE_ABSOLUTE"]
        gauss_points_results=[]
        VolumeOutput = True
        GiDPostMode = "Binary"
        GiDWriteMeshFlag = False
        GiDWriteConditionsFlag = True
        GiDWriteParticlesFlag = False
        GiDMultiFileFlag = "Single"
        gig_io = GiDOutput(fem_output_filename, VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
        gig_io.initialize_results(self.FEModelPart)
        gig_io.write_results(1, self.FEModelPart, nodal_results, gauss_points_results)
        gig_io.finalize_results()

    # --------------------------------------------------------------------------
    def OutputCADSurfacePoints( self, file_to_write ):
        self.OutputWriter.OutputCADSurfacePoints( file_to_write, self.ParameterResolutionForOutputOfSurfacePoints )

    # --------------------------------------------------------------------------
    def OutputGaussPointsOfFEMesh( self, file_to_write ):
        self.OutputWriter.OutputGaussPointsOfFEMesh( file_to_write, self.FEMGaussIntegrationDegree )

    # --------------------------------------------------------------------------
    def OutputControlPointDisplacementsInRhinoFormat( self, file_to_write ):
        self.OutputWriter.OutputControlPointDisplacementsInRhinoFormat( file_to_write )              

# ======================================================================================================================================
# Reconstruction
# ======================================================================================================================================    

# Input parameters
fem_filename = "tripod"
cad_geometry_filename = "tripod_geometry.json" 
cad_integration_data_filename = "tripod_integration_data.json"

# Initialize Reconstruction
CADReconstructionUtility = CADReconstrutionUtilities( fem_filename, cad_geometry_filename, cad_integration_data_filename )
CADReconstructionUtility.Initialize()

# Set Boundary Conditions
CADReconstructionUtility.SetDisplacementCouplingOnAllCouplingPoints()

# Some output before reconstruction
CADReconstructionUtility.OutputCADSurfacePoints( "surface_points_of_cad_geometry.txt" )

# Perform reconstruction
CADReconstructionUtility.PerformReconstruction()

# Some output
CADReconstructionUtility.OutputFEData()
CADReconstructionUtility.OutputCADSurfacePoints( "surface_points_of_updated_cad_geometry.txt" )
CADReconstructionUtility.OutputGaussPointsOfFEMesh( "gauss_points_of_fe_mesh.txt" )
CADReconstructionUtility.OutputControlPointDisplacementsInRhinoFormat( "tripod.post.res" )




# # ======================================================================================================================================
# # Mapping
# # ======================================================================================================================================    

# # Create CAD-mapper
# linear_solver = SuperLUSolver()
# # DiagPrecond = DiagonalPreconditioner()
# # linear_solver =  BICGSTABSolver(1e-9, 5000, DiagPrecond)
# # linear_solver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, 1e-9, 300, 2, 10)
# mapper = CADMapper(fe_model_part,cad_geometry,cad_integration_data,linear_solver)

# # Compute mapping matrix
# u_resolution = 300
# v_resolution = 300
# mapper.compute_mapping_matrix(u_resolution,v_resolution)

# # Apply boundary conditions
# penalty_factor_displacement_coupling = 1e3
# penalty_factor_rotation_coupling = 1e3
# penalty_factor_dirichlet_condition = 1e3
# edges_with_specific_dirichlet_conditions = [ ]
# edges_with_enforced_tangent_continuity = [ ]
# mapper.apply_boundary_conditions( penalty_factor_displacement_coupling, 
#                                   penalty_factor_rotation_coupling, 
#                                   penalty_factor_dirichlet_condition,
#                                   edges_with_specific_dirichlet_conditions,
#                                   edges_with_enforced_tangent_continuity )

# # Perform mapping
# mapper.map_to_cad_space()