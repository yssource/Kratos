from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMApp

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Other imports
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.KratosUnittest import isclose as t_isclose

import KratosMultiphysics.ExternalSolversApplication as ESA

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

class PreTestMembrane(UnitTest.TestCase):

    def _add_variables(self,mp,explicit_dynamics=False):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(SMApp.POINT_LOAD)

    def _add_dofs(self,mp):
        # Adding the dofs AND their corresponding reaction!
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

    def _create_nodes_3d4n(self,mp):
        mp.CreateNewNode(1,   0.0,  0.0,  0.0)
        mp.CreateNewNode(2,   2.0,  0.0,  0.0)
        mp.CreateNewNode(3,   2.0,  1.0,  0.0)
        mp.CreateNewNode(4,   0.0,  1.0,  0.0)

    def _create_elements_3d4n(self,mp):
        element_name = "PreStressMembraneElement3D4N"
        mp.CreateNewElement(element_name,  1 , [1, 2, 3, 4], mp.GetProperties()[1])

    def _apply_dirichlet_X_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
    def _apply_dirichlet_Y_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
    def _apply_dirichlet_Z_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)

    def _apply_neumann_BCs(self,mp):

        for node in mp.Nodes:
            node.SetSolutionStepValue(SMApp.POINT_LOAD, 0,[500.0, 0.0, 0.0])
            mp.CreateNewCondition("PointLoadCondition3D1N",node.Id,[node.Id],mp.GetProperties()[1])
        print(mp)

    def _apply_material_properties(self,mp):
        #define properties
        mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS,206.9e6)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO,0.0)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS,0.01)
        mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY,7850.0)

        constitutive_law = SMApp.LinearElasticPlaneStress2DLaw()

        local_axis_1 = KratosMultiphysics.Vector(3)
        local_axis_1[0] = 1.0
        local_axis_1[1] = 0.0
        local_axis_1[2] = 0.0

        local_axis_2= KratosMultiphysics.Vector(3)
        local_axis_2[0] = 0.0
        local_axis_2[1] = 1.0
        local_axis_2[2] = 0.0

        mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,constitutive_law)
        mp.GetProperties()[1].SetValue(SMApp.PROJECTION_TYPE_COMBO,"planar")
        mp.GetProperties()[1].SetValue(SMApp.PRESTRESS_AXIS_1_GLOBAL,local_axis_1)
        mp.GetProperties()[1].SetValue(SMApp.PRESTRESS_AXIS_2_GLOBAL,local_axis_2)

        prestress = KratosMultiphysics.Vector(3)
        prestress[0]=1e2        #1e4
        prestress[1]=1e2
        prestress[2]=0.0
        mp.GetProperties()[1].SetValue(SMApp.PRESTRESS_VECTOR,prestress)


    def _solve_static(self,mp):
        linear_solver = ESA.SuperLUSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        # compute_reactions = False
        # reform_step_dofs = False
        # calculate_norm_dx = True
        # move_mesh_flag = True
        # strategy = KratosMultiphysics.ResidualBasedLinearStrategy(mp,scheme,linear_solver,
        #                                                           compute_reactions,reform_step_dofs,calculate_norm_dx,move_mesh_flag)

        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-6,1e-9)
        max_iters = 1000
        compute_reactions = False
        reform_step_dofs = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                        scheme,
                                                                        linear_solver,
                                                                        convergence_criterion,
                                                                        builder_and_solver,
                                                                        max_iters,
                                                                        compute_reactions,
                                                                        reform_step_dofs,
                                                                        move_mesh_flag)
        strategy.SetEchoLevel(3)
        strategy.Check()
        strategy.Solve()

    def _check_static_results(self,node,displacement_results):
        #check that the results are exact on the node
        displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        self.assertAlmostEqual(displacement[0], displacement_results[0], 8)
        self.assertAlmostEqual(displacement[1], displacement_results[1], 8)
        self.assertAlmostEqual(displacement[2], displacement_results[2], 8)

    def _set_up_system_3d4n(self,current_model):
        mp = current_model.CreateModelPart("Structure")
        mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes_3d4n(mp)
        self._add_dofs(mp)
        self._create_elements_3d4n(mp)

        bcs_neumann = mp.CreateSubModelPart("BoundaryCondtionsNeumann")
        bcs_neumann.AddNodes([2,3])

        self._apply_neumann_BCs(bcs_neumann)
        print(bcs_neumann)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet_Z = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_1")
        bcs_dirichlet_Z.AddNodes([1,2,3,4])
        self._apply_dirichlet_Z_BCs(bcs_dirichlet_Z)

        bcs_dirichlet_Y = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_2")
        bcs_dirichlet_Y.AddNodes([1])
        self._apply_dirichlet_Y_BCs(bcs_dirichlet_Y)

        bcs_dirichlet_X = mp.CreateSubModelPart("BoundaryCondtionsDirichlet_3")
        bcs_dirichlet_X.AddNodes([1,4])
        self._apply_dirichlet_X_BCs(bcs_dirichlet_X)
        print(bcs_dirichlet_Z)

        # self._apply_neumann_BCs(mp,3)

        return mp

    def __post_process(self, main_model_part):
        from gid_output_process import GiDOutputProcess
        self.gid_output = GiDOutputProcess(main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["DISPLACEMENT"],
                                                "gauss_point_results" : []
                                            }
                                        }
                                        """)
                                    )

        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()


class StaticTestMembrane(PreTestMembrane):

    def test_membrane_3d4n_static(self):
        displacement_results_1 = [0.0 , 0.0 , 0.0]
        displacement_results_2 = [0.000964981929606997 , 8.392107497161046e-10 , 0.0]
        displacement_results_3 = [0.0009649849296071892 , -4.824857148915281e-07 , 0.0]
        displacement_results_4 = [0.0 , -4.824857148915281e-07 , 0.0]

        current_model = KratosMultiphysics.Model()

        mp = self._set_up_system_3d4n(current_model)

        self._solve_static(mp)

        self._check_static_results(mp.Nodes[1],displacement_results_1)
        self._check_static_results(mp.Nodes[2],displacement_results_2)
        self._check_static_results(mp.Nodes[3],displacement_results_3)
        self._check_static_results(mp.Nodes[4],displacement_results_4)

        #self.__post_process(mp)


if __name__ == '__main__':
    UnitTest.main()
