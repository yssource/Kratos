
from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from new_linear_solver_factory import ConstructSolver

class TestEigenIterativeSolver(KratosUnittest.TestCase):
    def _execute_eigen_iterative_solver_test(self, class_name, solver_type):
        # check if solver is available
        if (not hasattr(EigenSolversApplication, class_name)):
            self.skipTest(class_name + " is not included in the compilation of the EigenSolversApplication")

        space = KratosMultiphysics.UblasSparseSpace()

        settings = KratosMultiphysics.Parameters('{ "solver_type" : "' + solver_type + '''",
                                                    "max_iteration": 1000,
                                                    "tolerance": 1e-8,
                                                    "echo_level": 1 }''')

        solver = ConstructSolver(settings)

        a = KratosMultiphysics.CompressedMatrix()

        KratosMultiphysics.ReadMatrixMarketMatrix('../../../kratos/tests/A.mm', a) # symmetric test matrix

        dimension = a.Size1()

        b_exp = KratosMultiphysics.Vector(dimension) # [1, 2, ..., dimension-1, dimension]

        for i in range(dimension):
           b_exp[i] = i + 1

        x = KratosMultiphysics.Vector(dimension)

        solver.Solve(a, x, b_exp)

        b_act = KratosMultiphysics.Vector(dimension)
        space.Mult(a, x, b_act)

        for i in range(dimension):
           self.assertAlmostEqual(b_act[i], b_exp[i], 4)

    def test_eigen_cg(self):
        self._execute_eigen_iterative_solver_test('ConjugateGradientSolver', 'eigen_cg')

    def test_eigen_bicgstab(self):
        self._execute_eigen_iterative_solver_test('BiCGSTABSolver', 'eigen_bicgstab')

if __name__ == '__main__':
    KratosUnittest.main()
