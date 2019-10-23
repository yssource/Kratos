import os,sys
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.RANSModellingApplication as KratosRANS

from  KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
missing_applications_message = ["Missing required application(s):",]
have_required_applications = CheckIfApplicationsAvailable("HDF5Application")
if have_required_applications:
    import KratosMultiphysics.HDF5Application as kh5
else:
    missing_applications_message.append("HDF5Application")

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_analysis import AdjointFluidAnalysis

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

def str_vector(vector):
    size = vector.Size()
    str_vec = "[" + str(size) + "]("
    for i in range(size-1):
        str_vec += "{0:1.18e},".format(float(vector[i]))
    str_vec += "{0:1.18e})".format(float(vector[size-1]))
    return str_vec

def str_matrix(matrix):
    size1 = matrix.Size1()
    size2 = matrix.Size2()
    str_vec = "[" + str(size1) + "," + str(size2) + "]("
    for i in range(size1):
        str_vec += "("
        for j in range(size2-1):
            str_vec += "{0:1.18e},".format(float(matrix[i,j]))
        str_vec += "{0:1.18e}),".format(float(matrix[i, size2-1]))
    return str_vec[:-1] + ")"

@KratosUnittest.skipUnless(have_required_applications," ".join(missing_applications_message))
class AdjointKEpsilonSensitivity2D(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _removeH5Files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utils.DeleteFileIfExisting(name)

    def _readNodalCoordinates(self,node_id,model_part_file_name):
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        line = lines[node_id] # assumes consecutive node numbering starting with 1
        components = line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        return [float(components[i]) for i in range(1,4)]

    def _writeNodalCoordinates(self,node_id,coords,model_part_file_name):
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        node_lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        old_line = node_lines[node_id] # assumes consecutive node numbering starting with 1
        components = old_line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        new_line = '{:5d}'.format(node_id) + ' ' \
             + '{:19.10f}'.format(coords[0]) + ' ' \
             + '{:19.10f}'.format(coords[1]) + ' ' \
             + '{:19.10f}'.format(coords[2]) + '\n'
        lines[lines.index(old_line)] = new_line
        with open(model_part_file_name + '.mdpa', 'w') as model_part_file:
            model_part_file.writelines(lines)

    def _computeFiniteDifferenceDragSensitivity(self,node_ids,step_size,model_part_file_name,drag_direction,drag_file_name):
        sensitivity = []
        # unperturbed drag
        self.fd_string += "Running unperturbed....\n"
        self.solve(model_part_file_name)
        drag0 = _getTimeAveragedDrag(drag_direction,drag_file_name)
        for node_id in node_ids:
            node_sensitivity = []
            coord = self._readNodalCoordinates(node_id,model_part_file_name)
            # X + h
            perturbed_coord = [coord[0] + step_size, coord[1], coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,model_part_file_name)
            self.fd_string += "Running " + str(node_id) + " x perturbation " + str(step_size) + "\n"
            self.solve(model_part_file_name)
            drag = _getTimeAveragedDrag(drag_direction,drag_file_name)
            node_sensitivity.append((drag - drag0) / step_size)
            # Y + h
            perturbed_coord = [coord[0], coord[1] + step_size, coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,model_part_file_name)
            self.fd_string += "Running " + str(node_id) + " y perturbation " + str(step_size) + "\n"
            self.solve(model_part_file_name)
            drag = _getTimeAveragedDrag(drag_direction,drag_file_name)
            node_sensitivity.append((drag - drag0) / step_size)
            sensitivity.append(node_sensitivity)
            # return mdpa file to unperturbed state
            self._writeNodalCoordinates(node_id,coord,model_part_file_name)
        return sensitivity

    def _readParameters(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        return project_parameters

    def _createFluidTest(self, parameter_file_name):
        test = FluidDynamicsAnalysis(Model(), self._readParameters(parameter_file_name))
        return test

    def _createAdjointTest(self, parameter_file_name):
        test = AdjointFluidAnalysis(Model(), self._readParameters(parameter_file_name))
        return test

    def solve(self, parameter_file_name):
        test = self._createFluidTest(parameter_file_name)
        test.Run()
        model = test._GetSolver().main_model_part.GetModel()
        self.CalculatePerturbedMatrices(model, 1e-8)
        vms_model_part = model["MainModelPart"]
        k_model_part = model["TurbulenceModelPart_RansEvmK"]
        epsilon_model_part = model["TurbulenceModelPart_RansEvmEpsilon"]

        vms_residual = self.CalculateResidual(vms_model_part)
        self.fd_string += "vms_residual : " + str_vector(vms_residual) + "\n"

        k_residual = self.CalculateResidual(k_model_part)
        self.fd_string += "k_residual : " + str_vector(k_residual) + "\n"

        epsilon_residual = self.CalculateResidual(epsilon_model_part)
        self.fd_string += "epsilon_residual : " + str_vector(epsilon_residual) + "\n"

        self.fd_string += "Finished Calculation...\n"

    def RunNutCalculations(self, model):
        nutk_params = Parameters(r"""
                {
                        "model_part_name": "MainModelPart",
                        "echo_level": 0
                }
        """)

        process_nutk = KratosRANS.RansNutKEpsilonHighReCalculationProcess(model, nutk_params)
        process_nutk.Execute()

        nuty_params = Parameters(r"""
                {
                        "model_part_name": "MainModelPart.Boundary",
                        "echo_level": 0
                }
        """)

        process_nuty = KratosRANS.RansNutYPlusWallFunctionProcess(model, nuty_params)
        process_nuty.Execute()

        model_part = model.GetModelPart("MainModelPart")
        for node in model_part.Nodes:
            nu = node.GetSolutionStepValue(KINEMATIC_VISCOSITY)
            nu_t = node.GetSolutionStepValue(TURBULENT_VISCOSITY)
            node.SetSolutionStepValue(VISCOSITY, 0, nu + nu_t)

    def CalculatePerturbedMatrices(self, model, step_size):
        if (hasattr(self, "done")):
            return

        self.done = True

        str_out = ""
        vms_model_part = model["MainModelPart"]
        k_model_part = model["TurbulenceModelPart_RansEvmK"]
        epsilon_model_part = model["TurbulenceModelPart_RansEvmEpsilon"]

        str_out += self.CalculateEquationPerturbation(model, vms_model_part, "vms", -step_size, 3)
        str_out += self.CalculateEquationPerturbation(model, k_model_part, "k", -step_size, 2)
        str_out += self.CalculateEquationPerturbation(model, epsilon_model_part, "epsilon", -step_size, 2)

        with open("/media/suneth/Projects/LRZ Sync+Share/PhD/Workout/k_epsilon_fd_sensitivities/fd_matrices.out", "w") as file_output:
            file_output.write(str_out)

    def CalculateEquationPerturbation(self, model, model_part, cap, step_size, local_size):
        str_out = ""
        # vms_vms
        self.RunNutCalculations(model)
        residual_0 = self.CalculateResidual(model_part)
        vms_vms = Matrix(local_size * 3,residual_0.Size())
        for c in range(3):
            for k in range(2):
                self.PerturbVectorVariable(model_part, c+1, VELOCITY, k, step_size)
                self.RunNutCalculations(model)
                residual = self.CalculateResidual(model_part)
                self.PerturbVectorVariable(model_part, c+1, VELOCITY, k, -step_size)
                self.RunNutCalculations(model)
                for a in range(residual.Size()):
                    vms_vms[c*local_size + k, a] = (residual[a] - residual_0[a]) / step_size

                if local_size==3:
                    self.PerturbScalarVariable(model_part, c+1, PRESSURE, step_size)
                    self.RunNutCalculations(model)
                    residual = self.CalculateResidual(model_part)
                    self.PerturbScalarVariable(model_part, c+1, PRESSURE, -step_size)
                    self.RunNutCalculations(model)
                    for a in range(residual.Size()):
                        vms_vms[c*local_size + 2, a] = (residual[a] - residual_0[a]) / step_size
        str_out += cap + "_vms : " + str_matrix(vms_vms) + "\n"

        # vms_k
        self.RunNutCalculations(model)
        residual_0 = self.CalculateResidual(model_part)
        vms_k = Matrix(3,residual_0.Size())
        for c in range(3):
                self.PerturbScalarVariable(model_part, c+1, KratosRANS.TURBULENT_KINETIC_ENERGY, step_size)
                self.RunNutCalculations(model)
                residual = self.CalculateResidual(model_part)
                self.PerturbScalarVariable(model_part, c+1, KratosRANS.TURBULENT_KINETIC_ENERGY, -step_size)
                self.RunNutCalculations(model)
                for a in range(residual.Size()):
                    vms_k[c, a] = (residual[a] - residual_0[a]) / step_size
        str_out += cap + "_k : " + str_matrix(vms_k) + "\n"

        # vms_epsilon
        self.RunNutCalculations(model)
        residual_0 = self.CalculateResidual(model_part)
        vms_epsilon = Matrix(3,residual_0.Size())
        for c in range(3):
                self.PerturbScalarVariable(model_part, c+1, KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, step_size)
                self.RunNutCalculations(model)
                residual = self.CalculateResidual(model_part)
                self.PerturbScalarVariable(model_part, c+1, KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, -step_size)
                self.RunNutCalculations(model)
                for a in range(residual.Size()):
                    vms_epsilon[c, a] = (residual[a] - residual_0[a]) / step_size
        str_out += cap + "_epsilon : " + str_matrix(vms_epsilon) + "\n"

        return str_out

    def PerturbScalarVariable(self, model_part, node_id, variable, step_size):
        node = model_part.GetNode(node_id)
        scalar_value = node.GetSolutionStepValue(variable)
        node.SetSolutionStepValue(variable, 0, scalar_value + step_size)

    def PerturbVectorVariable(self, model_part, node_id, variable, direction, step_size):
        node = model_part.GetNode(node_id)
        vector_value = node.GetSolutionStepValue(variable)
        vector_value[direction] += step_size
        node.SetSolutionStepValue(variable, 0, vector_value)

    def CalculateResidual(self, model_part):
        element = model_part.GetElement(1)
        process_info = model_part.ProcessInfo
        lhs = Matrix()
        rhs = Vector()

        element.CalculateLocalSystem(lhs, rhs, process_info)
        element.CalculateLocalVelocityContribution(lhs, rhs, process_info)

        return rhs

    def testOneElementSteady(self):
        self.fd_string = ""
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('AdjointKEpsilonSensitivity2DTest/one_element_steady_test')
            # solve adjoint
            test = AdjointFluidAnalysis(Model(), self._readParameters('AdjointKEpsilonSensitivity2DTest/one_element_steady_test_adjoint'))
            test.Run()
            Sensitivity = [[]]
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
            Sensitivity[0].append(test._GetSolver().main_model_part.GetNode(1).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))

            # calculate sensitivity by finite difference
            step_size = 0.00000001
            FDSensitivity = self._computeFiniteDifferenceDragSensitivity([1],step_size,'./AdjointKEpsilonSensitivity2DTest/one_element_steady_test',[1.0,0.0,0.0],'./MainModelPart.Structure_drag.dat')
            sys.stdout.flush()
            with open("fd_out.dat", "w") as file_output:
                file_output.write(self.fd_string)

            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 3)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 3)
            self._removeH5Files("MainModelPart")
            kratos_utils.DeleteFileIfExisting("./AdjointKEpsilonSensitivity2DTest/one_element_test.time")
            kratos_utils.DeleteFileIfExisting("./Structure_drag.dat")
            kratos_utils.DeleteFileIfExisting("./one_element.post.bin")
            kratos_utils.DeleteFileIfExisting("./tests.post.lst")

    def tearDown(self):
        pass

def _getTimeAveragedDrag(direction,drag_file_name):
    time_steps, reactions = _readDrag(drag_file_name)
    total_drag = 0.0
    for reaction in reactions:
        total_drag += reaction[0]*direction[0]+reaction[1]*direction[1]+reaction[2]*direction[2]
    if len(time_steps) > 1:
        delta_time = time_steps[1] - time_steps[0]
        total_drag *= delta_time
    return total_drag

def _readDrag(filename):
        with open(filename, "r") as file_input:
            lines = file_input.readlines()
        time_steps = []
        reaction = []
        for line in lines:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            time_step_data = [float(v) for v in line.split()]
            time, fx, fy, fz = time_step_data
            time_steps.append(time)
            reaction.append([fx, fy, fz])
        return time_steps, reaction

if __name__ == '__main__':
    KratosUnittest.main()
