from KratosMultiphysics import *
from KratosMultiphysics.CompressiblePotentialFlowApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.AdjointFluidApplication import *


import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication
import random
#from sympy import *


class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        # default settings string in json format
        potential_settings = KratosMultiphysics.Parameters("""
        {
           "response_function_settings" : {
                "response_type"            : "lift",
                "sensitivity_model_part_name"  : "Boundary",
                "nodal_sensitivity_variables" : ["SHAPE_SENSITIVITY"],
                "custom_settings" : {
                    "structure_model_part_name" : "Structure",
                    "lift_direction"            : [0.0, 1.0, 0.0]
                }                
            }
        }""") 
        
        #"structure_model_part_name" : "Structure",
        #        "lift_direction"            : [1.0, 0.0, 0.0],
        
        self.settings = potential_settings
        self.settings.ValidateAndAssignDefaults(potential_settings)
        
        # create test model part
        self.model_part = ModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(NORMAL)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        self.model_part.Nodes[1].Set(KratosMultiphysics.STRUCTURE)#Ask Mike or Riccardo: For some reason I loose this setting in the functions and need to do it again, why? I think before it was working
        self.model_part.Nodes[3].Set(KratosMultiphysics.STRUCTURE)
        prop = self.model_part.GetProperties()[0]
        self.model_part.CreateNewElement("CompressiblePotentialFlowElement2D3N", 1, [1, 2, 3], prop)
        self.model_part.CreateNewCondition("PotentialWallCondition2D2N",1, [1,3], prop)
        
        self.potential_element = self.model_part.GetElement(1)
        self.potential_condition = self.model_part.GetCondition(1)
        
        velocityInf = Vector(3)
        velocityInf = (10.0,0,0)
        self.model_part.ProcessInfo[VELOCITY]=velocityInf
        
        
        self.structure_model_part = self.model_part.CreateSubModelPart("Structure")
        
        for cond in self.model_part.Conditions:
            self.structure_model_part.Conditions.append(cond)
        self.potential_condition.SetValue(PRESSURE, 0.48453305038458894)
        
        self.model_part.ProcessInfo[OSS_SWITCH] = 0

        #self.AssignSolutionStepDataSet1(0)
        self.model_part.Nodes[1].SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,2.0)
        self.model_part.Nodes[2].SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,-18.0)
        self.model_part.Nodes[3].SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,-28.0)
                        
        self.response_function = AdjointFluidApplication.PotentialFlowLiftResponseFunction2D(self.model_part, self.settings["response_function_settings"])
        

    def AssignSolutionStepDataSet1(self, step=0):
        # generate nodal solution step test data
        random.seed(1.0)
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,step,random.random())
            #node.SetSolutionStepValue(VELOCITY,step,random.random())


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

    def assertVectorAlmostEqual(self, vector1, vector2, prec=6):
        self.assertEqual(vector1.Size(), vector2.Size())
        for i in range(vector1.Size()):
            self.assertAlmostEqual(vector1[i], vector2[i], prec)
    
    
    def test_CalculateValue(self):
        print('Testing CalculateValue...')
        for node in self.model_part.Nodes:
            if node.Y < 1 and node.X < 1 or node.Y > 0.5:
                node.Set(KratosMultiphysics.STRUCTURE)
                #print(node)
        liftpython = self.ComputeLift()
        liftc = self.response_function.CalculateValue(self.model_part)
        self.assertAlmostEqual(liftc, liftpython, 16)
        
    def test_CalculateGradient(self):
        print('Testing CalculateGradient...')
        #setting structure nodes
        for node in self.model_part.Nodes:
            if node.Y < 1 and node.X < 1 or node.Y > 0.5:
                node.Set(KratosMultiphysics.STRUCTURE)
        #unperturbed lift
        lift0 = self.ComputeLift() 
        
        # finite difference approximation using python
        h = 0.0000001
        FDResponseGradient = Vector(3)
        row_index = 0
        for node in self.model_part.Nodes:
            # Phi
            potential = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE,0)
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,potential+h)
            lift = self.ComputeLift() 
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,potential)            
            FDResponseGradient[row_index] = (lift - lift0) / h
            row_index = row_index + 1
        
        # finite difference approximation2 using c++
        h = 0.0000001
        FDResponseGradient2 = Vector(3)
        row_index = 0
        for node in self.model_part.Nodes:
            # Phi
            potential = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE,0)
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,potential+h)            
            lift = self.response_function.CalculateValue(self.model_part)
            node.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,potential)
            #lift = cp*An[1]
            FDResponseGradient2[row_index] = (lift - lift0) / h
            row_index = row_index + 1
        
        #print('FDResponseGradient =', FDResponseGradient)
        #print('FDResponseGradient2 =', FDResponseGradient2)
        self.assertVectorAlmostEqual(FDResponseGradient, FDResponseGradient2)
                
        # analytical implementation
        AdjointMatrix = Matrix(3,3)
        self.potential_element.CalculateLeftHandSide(AdjointMatrix,self.model_part.ProcessInfo)
        ResponseGradient = Vector(3)        
        self.response_function.CalculateGradient(self.potential_element,AdjointMatrix,ResponseGradient,self.model_part.ProcessInfo)
        self.assertVectorAlmostEqual(FDResponseGradient, ResponseGradient)
        
    
    def test_CalculateSensitivityGradient(self):
        print('Testing CalculateSensitivityGradient...')
        #setting structure nodes
        for node in self.model_part.Nodes:
            if node.Y < 1 and node.X < 1 or node.Y > 0.5:
                node.Set(KratosMultiphysics.STRUCTURE)
        #unperturbed lift
        lift0 = self.ComputeLift()
        
        # finite difference approximation (wo change of PRESSURE in lift function)
        h = 0.0000001
        FDResponseGradient = Vector(6)
        row_index = 0
        for node in self.model_part.Nodes:
            # X
            x = node.X
            node.X = x+h
            lift = self.response_function.CalculateValue(self.model_part)
            node.X = x
            FDResponseGradient[row_index] = (lift - lift0) / h #minus due to missing minus after
            row_index = row_index + 1
            # Y
            y = node.Y
            node.Y = y+h
            lift = self.response_function.CalculateValue(self.model_part)
            node.Y = y
            FDResponseGradient[row_index] = (lift - lift0) / h #minus due to missing minus after    
            row_index = row_index + 1
        
        # finite difference approximation (w change of PRESSURE)
        h = 0.0000001
        FDResponseGradient2 = Vector(6)
        row_index = 0
        for node in self.model_part.Nodes:
            # X
            x = node.X
            node.X = x+h
            lift = self.ComputeLift()
            node.X = x
            FDResponseGradient2[row_index] = (lift - lift0) / h #minus due to missing minus after
            row_index = row_index + 1
            # Y
            y = node.Y
            node.Y = y+h
            lift = self.ComputeLift()
            node.Y = y
            FDResponseGradient2[row_index] = (lift - lift0) / h #minus due to missing minus after    
            row_index = row_index + 1      
                
        #analytical implementation
        ShapeDerivativeMatrix = Matrix(6,3)
        self.potential_element.CalculateSensitivityMatrix(SHAPE_SENSITIVITY,ShapeDerivativeMatrix,self.model_part.ProcessInfo)
        ResponseGradient = Vector(6)
        self.response_function.CalculateSensitivityGradient(self.potential_element,SHAPE_SENSITIVITY,ShapeDerivativeMatrix,ResponseGradient,self.model_part.ProcessInfo)
        #print('FDResponseGradient =',FDResponseGradient)
        #print('FDResponseGradient2 =',FDResponseGradient2)
        #print('ResponseGradient =',ResponseGradient)
        self.assertVectorAlmostEqual(FDResponseGradient, FDResponseGradient2)
        #self.assertVectorAlmostEqual(FDResponseGradient, ResponseGradient)
       
        
        
    def ComputeLift(self):        
               
        An = self.potential_condition.GetNormal()
        An2 = Vector(2)
        An2[0]= self.model_part.Nodes[3].Y - self.model_part.Nodes[1].Y
        An2[1]= -(self.model_part.Nodes[3].X - self.model_part.Nodes[1].X)
                
        DNDe = Matrix(3,2)
        DNDe[0,0] = -1
        DNDe[1,0] = 1
        DNDe[2,0] = 0        
        DNDe[0,1] = -1
        DNDe[1,1] = 0
        DNDe[2,1] = 1
        
        xt = Matrix(2,3)        
        xt[0,0] = self.model_part.Nodes[1].X
        xt[0,1] = self.model_part.Nodes[2].X
        xt[0,2] = self.model_part.Nodes[3].X        
        xt[1,0] = self.model_part.Nodes[1].Y
        xt[1,1] = self.model_part.Nodes[2].Y
        xt[1,2] = self.model_part.Nodes[3].Y
        
        DxDe = xt * DNDe        
        DxDeDet = DxDe[0,0]*DxDe[1,1] - DxDe[1,0]*DxDe[0,1]
        
        DxDeInverse = Matrix(2,2)        
        DxDeInverse[0,0] =  DxDe[1,1]/DxDeDet
        DxDeInverse[0,1] = -DxDe[0,1]/DxDeDet
        DxDeInverse[1,0] = -DxDe[1,0]/DxDeDet
        DxDeInverse[1,1] =  DxDe[0,0]/DxDeDet
        
        #Checking inverse is computed correctly
        #I = DxDe*DxDeInverse
        #print('I =', I)
        #I2 = DxDeInverse*DxDe
        #print('I2 =', I2)

        DNDx = DNDe*DxDeInverse
                
        Potential = Vector(3)
        self.potential_element.GetFirstDerivativesVector(Potential,0)
        
        LocalVelocity = Potential*DNDx
        
        Uinf = 10.0
        cp = 1 - (LocalVelocity*LocalVelocity)/(Uinf*Uinf)
        lift = cp*An2[1]
        
        return lift 
        

if __name__ == '__main__':
    KratosUnittest.main()
