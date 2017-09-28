from KratosMultiphysics import *
from KratosMultiphysics.CompressiblePotentialFlowApplication import *

import KratosMultiphysics.KratosUnittest as KratosUnittest
import random

class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        self.delta_time = 1.0
        # create test model part
        self.model_part = ModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        prop = self.model_part.GetProperties()[0]
        self.model_part.CreateNewElement("AdjointCompressiblePotentialFlowElement2D3N", 1, [1, 2, 3], prop)
        self.model_part.CreateNewElement("AdjointCompressiblePotentialFlowElement2D3N", 2, [1, 2, 3], prop)
        #self.model_part.SetBufferSize(2)
        self.model_part.ProcessInfo[OSS_SWITCH] = 0

        self.potential_element = self.model_part.GetElement(1)
        self.adjoint_element = self.model_part.GetElement(2)

        self.AssignSolutionStepDataSet1(0)
        self.AssignSolutionStepDataSet2(1)

    def AssignSolutionStepDataSet1(self, step=0):
        # generate nodal solution step test data
        random.seed(1.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,step,random.random())


    def AssignSolutionStepDataSet2(self, step=0):
        # generate nodal solution step test data
        random.seed(2.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,step,random.random())
            

    def zeroVector(self,size):
        v = Vector(size)
        for i in range(size):
            v[i] = 0.0
        return v

    def assertMatrixAlmostEqual(self, matrix1, matrix2, prec=7):
        self.assertEqual(matrix1.Size1(), matrix2.Size1())
        self.assertEqual(matrix1.Size2(), matrix2.Size2())
        for i in range(matrix1.Size1()):
            for j in range(matrix1.Size2()):
                self.assertAlmostEqual(matrix1[i,j], matrix2[i,j], prec)

    def test_CalculateLeftHandSide(self):
        # unperturbed residual
        LHS = Matrix(3,3)
        RHS = self.zeroVector(3)
        Potential = Vector(3)
        self.potential_element.CalculateLocalSystem(LHS,RHS,self.model_part.ProcessInfo)
        self.potential_element.GetFirstDerivativesVector(Potential,0)#TODO: change GetFirstDerivativesVector to GetFirstDerivativesVector
        res0 = LHS * Potential
        # finite difference approximation
        h = 0.0000001
        FDAdjointMatrix = Matrix(3,3)
        row_index = 0
        for node in self.model_part.Nodes:
            # Phi1
            potential = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE,0)
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,potential+h)
            self.potential_element.CalculateLocalSystem(LHS,RHS,self.model_part.ProcessInfo)
            self.potential_element.GetFirstDerivativesVector(Potential,0)
            #print('Potentials =', "{0:.16f}".format(Potential[0]))
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,potential)
            res = LHS * Potential
            for j in range(3):
                FDAdjointMatrix[row_index,j] = (res[j] - res0[j]) / h
            row_index = row_index + 1

        # analytical implementation
        AdjointMatrix = Matrix(3,3)
        self.potential_element.CalculateLeftHandSide(AdjointMatrix,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDAdjointMatrix, AdjointMatrix)

    def test_CalculateSensitivityMatrix(self):
        # unperturbed residual
        LHS = Matrix(3,3)
        RHS = self.zeroVector(3)
        Potential = Vector(3)
        self.potential_element.CalculateLocalSystem(LHS,RHS,self.model_part.ProcessInfo)
        self.potential_element.GetFirstDerivativesVector(Potential,0)#TODO: change GetFirstDerivativesVector to GetFirstDerivativesVector
        res0 = LHS * Potential
        # finite difference approximation
        h = 0.00000001
        FDShapeDerivativeMatrix = Matrix(6,3)
        row_index = 0
        for node in self.model_part.Nodes:
            # X
            x = node.X
            node.X = x+h
            self.potential_element.CalculateLocalSystem(LHS,RHS,self.model_part.ProcessInfo)
            node.X = x
            res = LHS * Potential
            for j in range(3):
                FDShapeDerivativeMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            # Y
            y = node.Y
            node.Y = y+h
            self.potential_element.CalculateLocalSystem(LHS,RHS,self.model_part.ProcessInfo)
            node.Y = y
            res = LHS * Potential
            for j in range(3):
                FDShapeDerivativeMatrix[row_index,j] = -(res[j] - res0[j]) / h
            row_index = row_index + 1
            
        # analytical implementation
        ShapeDerivativeMatrix = Matrix(6,3)
        self.adjoint_element.CalculateSensitivityMatrix(SHAPE_SENSITIVITY,ShapeDerivativeMatrix,self.model_part.ProcessInfo)
        self.assertMatrixAlmostEqual(FDShapeDerivativeMatrix, ShapeDerivativeMatrix)


if __name__ == '__main__':
    KratosUnittest.main()
