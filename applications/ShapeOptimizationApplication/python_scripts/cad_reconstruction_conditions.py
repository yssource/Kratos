# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Import ANurbs library
import numpy as np

# For automatic differentiation
from HyperJet import HyperJet, Jet
from numpy.testing import assert_array_almost_equal
import EQlib as eq

# ==============================================================================
class ReconstructionCondition(eq.Element):
    def __init__(self):
        eq.Element.__init__(self)

    # --------------------------------------------------------------------------
    def dofs(self):
        dof_list = []
        for node in self.pole_nodes:
            dof_list += [node.x.dof, node.y.dof, node.z.dof]
        return dof_list

    # --------------------------------------------------------------------------
    def GetPenaltyFactor(self):
        if hasattr(self, 'penalty_factor'):
            return self.penalty_factor
        else:
            raise RuntimeError("Trying to get a penalty factor for a condition that does not have any!!")

    # --------------------------------------------------------------------------
    @staticmethod
    def ComputeActual(pole_nodes, shape):
        value = np.zeros(3)
        for i, pole in enumerate(pole_nodes):
            value += shape[i]*pole.act_location

        return value

# ==============================================================================
class ReconstructionConditionWithAD(ReconstructionCondition):
    # --------------------------------------------------------------------------
    @staticmethod
    def ComputeActualJet(pole_nodes, shape):
        value = np.zeros(3)
        for i, pole in enumerate(pole_nodes):
            value += shape[i]*pole.act_location

        dx = np.zeros(len(shape) * 3)
        dy = np.zeros(len(shape) * 3)
        dz = np.zeros(len(shape) * 3)

        dx[0::3] = shape # [N0 0  0  N1 0  0  N2 0  ...]
        dy[1::3] = shape # [0  N0 0  0  N1 0  0  N2 ...]
        dz[2::3] = shape # [0  0  N0 0  0  N1 0  0  ...]

        return np.array([Jet(value[0], dx),
                         Jet(value[1], dy),
                         Jet(value[2], dz)])

    # --------------------------------------------------------------------------
    @staticmethod
    def ComputeActualHyperJet(pole_nodes, shape):
        value = np.zeros(3)
        for i, pole in enumerate(pole_nodes):
            value += shape[i]*pole.act_location

        dx = np.zeros(len(shape) * 3)
        dy = np.zeros(len(shape) * 3)
        dz = np.zeros(len(shape) * 3)

        dx[0::3] = shape # [N0 0  0  N1 0  0  N2 0  ...]
        dy[1::3] = shape # [0  N0 0  0  N1 0  0  N2 ...]
        dz[2::3] = shape # [0  0  N0 0  0  N1 0  0  ...]

        return np.array([HyperJet(value[0], dx),
                         HyperJet(value[1], dy),
                         HyperJet(value[2], dz)])

# ==============================================================================
class DistanceMinimizationCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, fe_node, pole_nodes, shape_functions, variable_to_map, weight):
        super().__init__()

        self.fe_node = fe_node
        self.pole_nodes = pole_nodes
        self.shape_function_values = np.asarray(shape_functions[0], float)
        self.variable_to_map = variable_to_map
        self.weight = weight

        node_initial_coords = np.array([self.fe_node.X, self.fe_node.Y, self.fe_node.Z])
        nodal_update = np.array(self.fe_node.GetSolutionStepValue(self.variable_to_map))
        self.fe_node_coords = node_initial_coords + nodal_update

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        cad_location = super().ComputeActual(self.pole_nodes, self.shape_function_values)
        distance_vec =  cad_location - self.fe_node_coords
        value = 0.5*self.weight*distance_vec.dot(distance_vec)

        local_rhs = - self.weight * np.outer(self.shape_function_values, distance_vec)

        rhs[0::3] = local_rhs[:,0]
        rhs[1::3] = local_rhs[:,1]
        rhs[2::3] = local_rhs[:,2]

        return value

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        num_dofs = len(self.shape_function_values)*3

        # LHS
        local_lhs_block = self.weight * np.outer(self.shape_function_values, self.shape_function_values)

        lhs[:] = np.zeros([num_dofs,num_dofs])
        lhs[0::3,0::3] = local_lhs_block
        lhs[1::3,1::3] = local_lhs_block
        lhs[2::3,2::3] = local_lhs_block

        # RHS
        value = self.CalculateRHS(rhs)

        return value

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================
class DistanceMinimizationConditionWithAD(ReconstructionConditionWithAD):
    # --------------------------------------------------------------------------
    def __init__(self, fe_node, pole_nodes, shape_functions, variable_to_map, weight):
        super().__init__()

        self.fe_node = fe_node
        self.pole_nodes = pole_nodes
        self.shape_function_values = np.asarray(shape_functions[0], float)
        self.variable_to_map = variable_to_map
        self.weight = weight

        node_initial_coords = np.array([self.fe_node.X, self.fe_node.Y, self.fe_node.Z])
        nodal_update = np.array(self.fe_node.GetSolutionStepValue(self.variable_to_map))
        self.fe_node_coords = node_initial_coords + nodal_update

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        cad_location = self.ComputeActualJet(self.pole_nodes, self.shape_function_values)

        v = cad_location - self.fe_node_coords

        f = 0.5 * self.weight * v.dot(v)

        rhs[:] = -f.g
        return f.f

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        cad_location = self.ComputeActualHyperJet(self.pole_nodes, self.shape_function_values)

        v = cad_location - self.fe_node_coords

        f = 0.5 * self.weight * v.dot(v)

        rhs[:] = -f.g
        lhs[:] = f.h
        return f.f

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================
class PositionEnforcementCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, target_position, pole_nodes, shape_functions, penalty_factor, weight):
        super().__init__()

        self.target_position = target_position
        self.pole_nodes = pole_nodes
        self.shape_function_values = np.asarray(shape_functions[0], float)
        self.penalty_factor = penalty_factor
        self.weight = weight

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        cad_location = super().ComputeActual(self.pole_nodes, self.shape_function_values)
        distance_vec =  cad_location - self.target_position
        value = 0.5 * self.penalty_factor * self.weight * distance_vec.dot(distance_vec)

        local_rhs = - self.penalty_factor * self.weight * np.outer(self.shape_function_values, distance_vec)

        rhs[0::3] = local_rhs[:,0]
        rhs[1::3] = local_rhs[:,1]
        rhs[2::3] = local_rhs[:,2]

        return value

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        num_dofs = len(self.shape_function_values)*3

        # LHS
        local_lhs_block = self.penalty_factor * self.weight * np.outer(self.shape_function_values, self.shape_function_values)

        lhs[:] = np.zeros([num_dofs,num_dofs])
        lhs[0::3,0::3] = local_lhs_block
        lhs[1::3,1::3] = local_lhs_block
        lhs[2::3,2::3] = local_lhs_block

        # RHS
        value = self.CalculateRHS(rhs)

        return value

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================
class TangentEnforcementConditionWithAD( ReconstructionConditionWithAD ):
    # --------------------------------------------------------------------------
    def __init__(self, target_normal, pole_nodes, shape_functions, penalty_factor, weight):
        super().__init__()

        self.target_normal = target_normal
        self.pole_nodes = pole_nodes
        self.shape_function_derivatives_u = np.asarray(shape_functions[1], float)
        self.shape_function_derivatives_v = np.asarray(shape_functions[2], float)
        self.penalty_factor = penalty_factor
        self.weight = weight

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        a1 = super().ComputeActualJet(self.pole_nodes, self.shape_function_derivatives_u)
        a2 = super().ComputeActualJet(self.pole_nodes, self.shape_function_derivatives_v)

        f_1 = self.penalty_factor * self.weight * np.dot(a1,self.target_normal)**2 * 0.5
        f_2 = self.penalty_factor * self.weight * np.dot(a2,self.target_normal)**2 * 0.5
        value = f_1.f + f_2.f

        rhs[:]  = - f_1.g
        rhs[:] -=   f_2.g

        return value

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        a1 = super().ComputeActualHyperJet(self.pole_nodes, self.shape_function_derivatives_u)
        a2 = super().ComputeActualHyperJet(self.pole_nodes, self.shape_function_derivatives_v)

        f_1 = self.penalty_factor * self.weight * np.dot(a1,self.target_normal)**2 * 0.5
        f_2 = self.penalty_factor * self.weight * np.dot(a2,self.target_normal)**2 * 0.5
        value = f_1.f + f_2.f

        # LHS
        lhs[:]  = f_1.h
        lhs[:] += f_2.h

        # RHS
        rhs[:]  = - f_1.g
        rhs[:] -=   f_2.g

        return value

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================
class DisplacementCouplingCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, pole_nodes_a, pole_nodes_b, shape_functions_a, shape_functions_b, geometry_a, geometry_b, parameters_a, parameters_b, penalty_factor, weight):
        super().__init__()

        self.pole_nodes_a = pole_nodes_a
        self.pole_nodes_b = pole_nodes_b
        self.shape_function_values_a = np.asarray(shape_functions_a[0], float)
        self.shape_function_values_b = np.asarray(shape_functions_b[0], float)
        self.geometry_a = geometry_a # only strored for aposteriori refinement
        self.geometry_b = geometry_b # only strored for aposteriori refinement
        self.parameters_a = parameters_a # only strored for aposteriori refinement
        self.parameters_b = parameters_b # only strored for aposteriori refinement
        self.penalty_factor = penalty_factor
        self.weight = weight

        self.pos_0_a = super().ComputeActual(pole_nodes_a, self.shape_function_values_a)
        self.pos_0_b = super().ComputeActual(pole_nodes_b, self.shape_function_values_b)

    # --------------------------------------------------------------------------
    def dofs(self):
        dof_list = []
        for node in self.pole_nodes_a:
            dof_list += [node.x.dof, node.y.dof, node.z.dof]
        for node in self.pole_nodes_b:
            dof_list += [node.x.dof, node.y.dof, node.z.dof]
        return dof_list

    # --------------------------------------------------------------------------
    def CalculateQualityIndicator(self):
        disp_a = super().ComputeActual(self.pole_nodes_a, self.shape_function_values_a) - self.pos_0_a
        disp_b = super().ComputeActual(self.pole_nodes_b, self.shape_function_values_b) - self.pos_0_b
        return disp_a - disp_b

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        num_dofs = len(self.dofs())
        num_pole_nodes_a = len(self.pole_nodes_a)

        N_ab = np.zeros(num_dofs)

        N_ab[0:num_pole_nodes_a*3:3] = self.shape_function_values_a
        N_ab[1:num_pole_nodes_a*3:3] = self.shape_function_values_a
        N_ab[2:num_pole_nodes_a*3:3] = self.shape_function_values_a

        N_ab[3*num_pole_nodes_a+0::3] -= self.shape_function_values_b
        N_ab[3*num_pole_nodes_a+1::3] -= self.shape_function_values_b
        N_ab[3*num_pole_nodes_a+2::3] -= self.shape_function_values_b

        # RHS
        disp_a = super().ComputeActual(self.pole_nodes_a, self.shape_function_values_a) - self.pos_0_a
        disp_b = super().ComputeActual(self.pole_nodes_b, self.shape_function_values_b) - self.pos_0_b
        delta = disp_a - disp_b

        value = 0.5  * self.penalty_factor * self.weight * delta.dot(delta)

        # a
        delta_ab = np.zeros(num_dofs)
        delta_ab[0:num_pole_nodes_a*3:3] = delta[0]
        delta_ab[1:num_pole_nodes_a*3:3] = delta[1]
        delta_ab[2:num_pole_nodes_a*3:3] = delta[2]

        # b
        delta_ab[3*num_pole_nodes_a+0::3] = delta[0]
        delta_ab[3*num_pole_nodes_a+1::3] = delta[1]
        delta_ab[3*num_pole_nodes_a+2::3] = delta[2]

        rhs[:] = - self.penalty_factor * self.weight * np.multiply(N_ab, delta_ab)

        return value

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        num_dofs = len(self.dofs())
        num_pole_nodes_a = len(self.pole_nodes_a)

        # LHS
        lhs[:] = np.zeros([num_dofs,num_dofs])

        # aa
        local_lhs_aa = self.penalty_factor * self.weight * np.outer(self.shape_function_values_a, self.shape_function_values_a)
        lhs[0:num_pole_nodes_a*3:3 , 0:num_pole_nodes_a*3:3] = local_lhs_aa
        lhs[1:num_pole_nodes_a*3:3 , 1:num_pole_nodes_a*3:3] = local_lhs_aa
        lhs[2:num_pole_nodes_a*3:3 , 2:num_pole_nodes_a*3:3] = local_lhs_aa

        # ab
        local_lhs_ab = - self.penalty_factor * self.weight * np.outer(self.shape_function_values_a, self.shape_function_values_b)
        lhs[0:num_pole_nodes_a*3:3 , 3*num_pole_nodes_a+0::3] = local_lhs_ab
        lhs[1:num_pole_nodes_a*3:3 , 3*num_pole_nodes_a+1::3] = local_lhs_ab
        lhs[2:num_pole_nodes_a*3:3 , 3*num_pole_nodes_a+2::3] = local_lhs_ab

        # ba
        local_lsh_ba = - self.penalty_factor * self.weight * np.outer(self.shape_function_values_b, self.shape_function_values_a)
        lhs[3*num_pole_nodes_a+0::3 , 0:num_pole_nodes_a*3:3] = local_lsh_ba
        lhs[3*num_pole_nodes_a+1::3 , 1:num_pole_nodes_a*3:3] = local_lsh_ba
        lhs[3*num_pole_nodes_a+2::3 , 2:num_pole_nodes_a*3:3] = local_lsh_ba

        # bb
        local_lhs_bb = self.penalty_factor * self.weight * np.outer(self.shape_function_values_b, self.shape_function_values_b)
        lhs[3*num_pole_nodes_a+0::3 , 3*num_pole_nodes_a+0::3] = local_lhs_bb
        lhs[3*num_pole_nodes_a+1::3 , 3*num_pole_nodes_a+1::3] = local_lhs_bb
        lhs[3*num_pole_nodes_a+2::3 , 3*num_pole_nodes_a+2::3] = local_lhs_bb

        # RHS
        value = self.CalculateRHS(rhs)

        return value

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================
class DisplacementCouplingConditionWithAD(ReconstructionConditionWithAD):
    # --------------------------------------------------------------------------
    def __init__(self, pole_nodes_a, pole_nodes_b, shape_functions_a, shape_functions_b, geometry_a, geometry_b, parameters_a, parameters_b, penalty_factor, weight):
        super().__init__()

        self.pole_nodes_a = pole_nodes_a
        self.pole_nodes_b = pole_nodes_b
        self.shape_function_values_a = np.asarray(shape_functions_a[0], float)
        self.shape_function_values_b = np.asarray(shape_functions_b[0], float)
        self.geometry_a = geometry_a # only strored for aposteriori refinement
        self.geometry_b = geometry_b # only strored for aposteriori refinement
        self.parameters_a = parameters_a # only strored for aposteriori refinement
        self.parameters_b = parameters_b # only strored for aposteriori refinement
        self.penalty_factor = penalty_factor
        self.weight = weight

        self.pos_0_a = super().ComputeActual(pole_nodes_a, self.shape_function_values_a)
        self.pos_0_b = super().ComputeActual(pole_nodes_b, self.shape_function_values_b)

    # --------------------------------------------------------------------------
    def dofs(self):
        dof_list = []
        for node in self.pole_nodes_a:
            dof_list += [node.x.dof, node.y.dof, node.z.dof]
        for node in self.pole_nodes_b:
            dof_list += [node.x.dof, node.y.dof, node.z.dof]
        return dof_list

    # --------------------------------------------------------------------------
    def CalculateQualityIndicator(self):
        disp_a = super().ComputeActual(self.pole_nodes_a, self.shape_function_values_a) - self.pos_0_a
        disp_b = super().ComputeActual(self.pole_nodes_b, self.shape_function_values_b) - self.pos_0_b
        return disp_a - disp_b

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        pos_a = super().ComputeActualJet(self.pole_nodes_a, self.shape_function_values_a)
        pos_b = super().ComputeActualJet(self.pole_nodes_b, self.shape_function_values_b)

        len_a = len(pos_a[0])
        len_b = len(pos_b[0])

        pos_a = np.array([hj.enlarge(len_b, False) for hj in pos_a])
        pos_b = np.array([hj.enlarge(len_a, True ) for hj in pos_b])

        position_difference = pos_a - pos_b

        f = 0.5  * self.penalty_factor * self.weight * position_difference.dot(position_difference)

        rhs[:] = -f.g
        return f.f

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        pos_a = super().ComputeActualHyperJet(self.pole_nodes_a, self.shape_function_values_a)
        pos_b = super().ComputeActualHyperJet(self.pole_nodes_b, self.shape_function_values_b)

        len_a = len(pos_a[0])
        len_b = len(pos_b[0])

        pos_a = np.array([hj.enlarge(len_b, False) for hj in pos_a])
        pos_b = np.array([hj.enlarge(len_a, True ) for hj in pos_b])

        position_difference = pos_a - pos_b

        f = 0.5 * self.penalty_factor * self.weight * position_difference.dot(position_difference)

        rhs[:] = -f.g
        lhs[:] = f.h
        return f.f

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================
class RotationCouplingConditionWithAD(ReconstructionConditionWithAD):
    # Note that this condition considers the rotation around a common tangent vector from one of the two coupled faces (T2_edge)
    # This is different to the classical IGA coupling where the rotation is calculated on the individual tangents (which requires a direction/sign check)
    # --------------------------------------------------------------------------
    def __init__(self, pole_nodes_a, pole_nodes_b, T2_edge, shape_functions_a, shape_functions_b, geometry_a, geometry_b, parameters_a, parameters_b, penalty_factor, weight):
        super().__init__()

        self.pole_nodes_a = pole_nodes_a
        self.pole_nodes_b = pole_nodes_b
        self.T2_edge = T2_edge
        self.shape_function_derivatives_u_a = np.asarray(shape_functions_a[1], float)
        self.shape_function_derivatives_v_a = np.asarray(shape_functions_a[2], float)
        self.shape_function_derivatives_u_b = np.asarray(shape_functions_b[1], float)
        self.shape_function_derivatives_v_b = np.asarray(shape_functions_b[2], float)
        self.geometry_a = geometry_a # only strored for aposteriori refinement
        self.geometry_b = geometry_b # only strored for aposteriori refinement
        self.parameters_a = parameters_a # only strored for aposteriori refinement
        self.parameters_b = parameters_b # only strored for aposteriori refinement
        self.penalty_factor = penalty_factor
        self.weight = weight

        A1_a = super().ComputeActual(self.pole_nodes_a, self.shape_function_derivatives_u_a)
        A2_a = super().ComputeActual(self.pole_nodes_a, self.shape_function_derivatives_v_a)
        A3_a = np.cross(A1_a, A2_a)
        self.A3_a = A3_a / np.linalg.norm(A3_a)

        A1_b = super().ComputeActual(self.pole_nodes_b, self.shape_function_derivatives_u_b)
        A2_b = super().ComputeActual(self.pole_nodes_b, self.shape_function_derivatives_v_b)
        A3_b = np.cross(A1_b, A2_b)
        self.A3_b = A3_b / np.linalg.norm(A3_b)

    # --------------------------------------------------------------------------
    def dofs(self):
        dof_list = []
        for node in self.pole_nodes_a:
            dof_list += [node.x.dof, node.y.dof, node.z.dof]
        for node in self.pole_nodes_b:
            dof_list += [node.x.dof, node.y.dof, node.z.dof]
        return dof_list

    # --------------------------------------------------------------------------
    def CalculateQualityIndicator(self):
        a1_a = super().ComputeActual(self.pole_nodes_a, self.shape_function_derivatives_u_a)
        a2_a = super().ComputeActual(self.pole_nodes_a, self.shape_function_derivatives_v_a)
        a3_a = np.cross(a1_a, a2_a)
        a3_a = a3_a / np.linalg.norm(a3_a)

        a1_b = super().ComputeActual(self.pole_nodes_b, self.shape_function_derivatives_u_b)
        a2_b = super().ComputeActual(self.pole_nodes_b, self.shape_function_derivatives_v_b)
        a3_b = np.cross(a1_b, a2_b)
        a3_b = a3_b / np.linalg.norm(a3_b)

        w_a = a3_a - self.A3_a
        w_b = a3_b - self.A3_b

        omega_a = np.cross(self.A3_a, w_a)
        omega_b = np.cross(self.A3_b, w_b)

        angle_a = np.arcsin(np.dot(omega_a, self.T2_edge))
        angle_b = np.arcsin(np.dot(omega_b, self.T2_edge))

        return angle_a - angle_b

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        a1_a = super().ComputeActualJet(self.pole_nodes_a, self.shape_function_derivatives_u_a)
        a2_a = super().ComputeActualJet(self.pole_nodes_a, self.shape_function_derivatives_v_a)
        a3_a = np.cross(a1_a, a2_a)
        a3_a = a3_a / np.linalg.norm(a3_a)

        a1_b = super().ComputeActualJet(self.pole_nodes_b, self.shape_function_derivatives_u_b)
        a2_b = super().ComputeActualJet(self.pole_nodes_b, self.shape_function_derivatives_v_b)
        a3_b = np.cross(a1_b, a2_b)
        a3_b = a3_b / np.linalg.norm(a3_b)

        w_a = a3_a - self.A3_a
        w_b = a3_b - self.A3_b

        omega_a = np.cross(self.A3_a, w_a)
        omega_b = np.cross(self.A3_b, w_b)

        angle_a = np.arcsin(np.dot(omega_a, self.T2_edge))
        angle_b = np.arcsin(np.dot(omega_b, self.T2_edge))

        len_a = len(angle_a)
        len_b = len(angle_b)

        angle_a = angle_a.enlarge(len_b, False)
        angle_b = angle_b.enlarge(len_a, True)

        angular_difference = angle_a - angle_b

        f = 0.5 * self.penalty_factor * self.weight * angular_difference**2

        rhs[:] = -f.g
        return f.f

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        a1_a = super().ComputeActualHyperJet(self.pole_nodes_a, self.shape_function_derivatives_u_a)
        a2_a = super().ComputeActualHyperJet(self.pole_nodes_a, self.shape_function_derivatives_v_a)
        a3_a = np.cross(a1_a, a2_a)
        a3_a = a3_a / np.linalg.norm(a3_a)

        a1_b = super().ComputeActualHyperJet(self.pole_nodes_b, self.shape_function_derivatives_u_b)
        a2_b = super().ComputeActualHyperJet(self.pole_nodes_b, self.shape_function_derivatives_v_b)
        a3_b = np.cross(a1_b, a2_b)
        a3_b = a3_b / np.linalg.norm(a3_b)

        w_a = a3_a - self.A3_a
        w_b = a3_b - self.A3_b

        omega_a = np.cross(self.A3_a, w_a)
        omega_b = np.cross(self.A3_b, w_b)

        angle_a = np.arcsin(np.dot(omega_a, self.T2_edge))
        angle_b = np.arcsin(np.dot(omega_b, self.T2_edge))

        len_a = len(angle_a)
        len_b = len(angle_b)

        angle_a = angle_a.enlarge(len_b, False)
        angle_b = angle_b.enlarge(len_a, True)

        angular_difference = angle_a - angle_b

        f = 0.5 * self.penalty_factor * self.weight * angular_difference**2

        rhs[:] = -f.g
        lhs[:] = f.h
        return f.f

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================
class KLShellConditionWithAD( ReconstructionConditionWithAD ):
    # --------------------------------------------------------------------------
    def __init__(self, pole_nodes, shape_function_values_and_dervatives, penalty_factor, weight):
        super().__init__()

        self.pole_nodes = pole_nodes
        self.shape_function_derivatives_u = np.asarray(shape_function_values_and_dervatives[1], float)
        self.shape_function_derivatives_v = np.asarray(shape_function_values_and_dervatives[2], float)
        self.shape_function_derivatives_uu = np.asarray(shape_function_values_and_dervatives[3], float)
        self.shape_function_derivatives_uv = np.asarray(shape_function_values_and_dervatives[4], float)
        self.shape_function_derivatives_vv = np.asarray(shape_function_values_and_dervatives[5], float)
        self.penalty_factor = penalty_factor
        self.weight = weight

        # Reference configuration
        A1 = super().ComputeActual(pole_nodes, self.shape_function_derivatives_u)
        A2 = super().ComputeActual(pole_nodes, self.shape_function_derivatives_v)
        A1_1 = super().ComputeActual(pole_nodes, self.shape_function_derivatives_uu)
        A1_2 = super().ComputeActual(pole_nodes, self.shape_function_derivatives_uv)
        A2_2 = super().ComputeActual(pole_nodes, self.shape_function_derivatives_vv)

        A11, A22, A12 = np.dot(A1, A1), np.dot(A2, A2), np.dot(A1, A2)
        self.A11, self.A22, self.A12 = A11, A22, A12

        A3 = np.cross(A1,A2)
        self.dA = np.linalg.norm(A3)
        A3 = A3/self.dA

        B11, B12, B22 = np.dot([A1_1, A1_2, A2_2], A3)
        self.B11, self.B12, self.B22 = B11, B12, B22

        # Material properties
        thickness = 1.0
        youngsmodulus = 1.0
        poissonsratio = 0.3

        self.Dm = np.array([
            [1.0, poissonsratio, 0],
            [poissonsratio, 1.0, 0],
            [0, 0, (1.0 - poissonsratio) * 0.5],
        ]) * youngsmodulus * thickness / (1.0 - np.power(poissonsratio,2))

        self.Db = np.array([
            [1.0, poissonsratio, 0],
            [poissonsratio, 1.0, 0],
            [0, 0, (1.0 - poissonsratio) * 0.5],
        ]) * youngsmodulus * np.power(thickness, 3) / (12.0 * (1.0 - np.power(poissonsratio,2)))

        # Transformation
        e1 = A1 / np.linalg.norm(A1)
        e2 = A2 - np.dot(A2, e1) * e1
        e2 /= np.linalg.norm(e2)

        det = A11 * A22 - A12 * A12

        g_ab_con = np.array([A22 / det, A11 / det, -A12 / det])

        g_con1 = g_ab_con[0] * A1 + g_ab_con[2] * A2
        g_con2 = g_ab_con[2] * A1 + g_ab_con[1] * A2

        eg11 = np.dot(e1, g_con1)
        eg12 = np.dot(e1, g_con2)
        eg21 = np.dot(e2, g_con1)
        eg22 = np.dot(e2, g_con2)

        self.Tm = np.array([
            [eg11 * eg11, eg12 * eg12, 2 * eg11 * eg12],
            [eg21 * eg21, eg22 * eg22, 2 * eg21 * eg22],
            [2 * eg11 * eg21, 2 * eg12 * eg22, 2 * (eg11 * eg22 + eg12 * eg21)],
        ])

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        B11, B12, B22 = self.B11, self.B12, self.B22
        A11, A22, A12 = self.A11, self.A22, self.A12

        a1 =super().ComputeActualJet(self.pole_nodes, self.shape_function_derivatives_u)
        a2 =super().ComputeActualJet(self.pole_nodes, self.shape_function_derivatives_v)
        a1_1 =super().ComputeActualJet(self.pole_nodes, self.shape_function_derivatives_uu)
        a1_2 =super().ComputeActualJet(self.pole_nodes, self.shape_function_derivatives_uv)
        a2_2 =super().ComputeActualJet(self.pole_nodes, self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        a11, a22, a12 = np.dot(a1, a1), np.dot(a2, a2), np.dot(a1, a2)
        b11, b12, b22 = np.dot([a1_1, a1_2, a2_2], a3)

        eps = np.dot(self.Tm, [a11 - A11, a22 - A22, a12 - A12]) * 0.5
        kap = np.dot(self.Tm, [B11 - b11, B12 - b12, B22 - b22])

        n = np.dot(eps, self.Dm)
        m = np.dot(kap, self.Db)

        f = self.penalty_factor * self.weight * (np.dot(eps, n) + np.dot(kap, m))

        rhs[:] = -f.g
        return f.f

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        B11, B12, B22 = self.B11, self.B12, self.B22
        A11, A22, A12 = self.A11, self.A22, self.A12

        a1 =super().ComputeActualHyperJet(self.pole_nodes, self.shape_function_derivatives_u)
        a2 =super().ComputeActualHyperJet(self.pole_nodes, self.shape_function_derivatives_v)
        a1_1 =super().ComputeActualHyperJet(self.pole_nodes, self.shape_function_derivatives_uu)
        a1_2 =super().ComputeActualHyperJet(self.pole_nodes, self.shape_function_derivatives_uv)
        a2_2 =super().ComputeActualHyperJet(self.pole_nodes, self.shape_function_derivatives_vv)

        a3 = np.cross(a1,a2)
        a3 /= np.linalg.norm(a3)

        a11, a22, a12 = np.dot(a1, a1), np.dot(a2, a2), np.dot(a1, a2)
        b11, b12, b22 = np.dot([a1_1, a1_2, a2_2], a3)

        eps = np.dot(self.Tm, [a11 - A11, a22 - A22, a12 - A12]) * 0.5
        kap = np.dot(self.Tm, [B11 - b11, B12 - b12, B22 - b22])

        n = np.dot(eps, self.Dm)
        m = np.dot(kap, self.Db)

        f = self.penalty_factor * self.weight * (np.dot(eps, n) + np.dot(kap, m))

        rhs[:] = -f.g
        lhs[:] = f.h
        return f.f

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================
class AlphaRegularizationCondition(ReconstructionCondition):
    # --------------------------------------------------------------------------
    def __init__(self, target_pole_rs, greville_parameters, surface_geometry, pole_nodes, pole_ids, shape_functions, alpha):
        super().__init__()

        self.target_pole_rs = target_pole_rs
        self.greville_parameters = greville_parameters
        self.geometry_data = surface_geometry.Data()
        self.pole_nodes = pole_nodes
        self.shape_function_values = np.asarray(shape_functions[0], float)
        self.penalty_factor = alpha

        r = target_pole_rs[0]
        s = target_pole_rs[1]
        id_of_pole_associated_to_grevile_point = pole_ids.index(r*self.geometry_data.NbPolesV() + s)

        I = np.zeros(len(self.shape_function_values))
        I[id_of_pole_associated_to_grevile_point] = 1

        self.IminN =  (I - shape_functions)

    # --------------------------------------------------------------------------
    def CalculateRHS(self, rhs):
        pole_coords = self.geometry_data.Pole(self.target_pole_rs[0],self.target_pole_rs[1])
        greville_point_coords = self.geometry_data.PointAt(self.greville_parameters[0], self.greville_parameters[1])

        distance_vec = pole_coords-greville_point_coords
        distance_vec = np.array(distance_vec)

        value = 0.5 * self.penalty_factor * distance_vec.dot(distance_vec)

        # RHS
        rhs[0::3] = -self.penalty_factor * self.IminN * distance_vec[0]
        rhs[1::3] = -self.penalty_factor * self.IminN * distance_vec[1]
        rhs[2::3] = -self.penalty_factor * self.IminN * distance_vec[2]

        return value

    # --------------------------------------------------------------------------
    def CalculateLocalSystem(self, rhs, lhs):
        num_dofs = len(self.shape_function_values)*3

        # LHS
        local_lhs_block = self.penalty_factor * np.outer(self.IminN, self.IminN)

        lhs[:] = np.zeros([num_dofs,num_dofs])
        lhs[0::3,0::3] = local_lhs_block
        lhs[1::3,1::3] = local_lhs_block
        lhs[2::3,2::3] = local_lhs_block

        # RHS
        value = self.CalculateRHS(rhs)

        return value

    # --------------------------------------------------------------------------
    def compute(self, rhs, lhs):
        if len(lhs) == 0:
            return self.CalculateRHS(rhs)
        return self.CalculateLocalSystem(rhs, lhs)

# ==============================================================================