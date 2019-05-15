import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return ComputeForcesOnNodesProcess(Model, settings["Parameters"])

# all the processes python processes should be derived from "python_process"

class ComputeForcesOnNodesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "create_output_file": false
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.body_model_part = Model[settings["model_part_name"].GetString()]
        self.create_output_file = settings["create_output_file"].GetBool()

    def ExecuteFinalizeSolutionStep(self):
        self.Execute()

    def Execute(self):
        KratosMultiphysics.Logger.PrintInfo('ComputeForcesOnNodesProcess', 'Computing reactions on nodes')

        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.REACTION, self.body_model_part.Nodes)

        free_stream_velocity = self.body_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_VELOCITY)
        free_stream_density = self.body_model_part.ProcessInfo.GetValue(CPFApp.FREE_STREAM_DENSITY)
        free_stream_velocity_norm = free_stream_velocity.norm_2()
        dynamic_pressure = 0.5*free_stream_density*free_stream_velocity_norm**2

        largest_force = 0.0 # comparison value to give warning in the end if no reaction force was computed
        for cond in self.body_model_part.Conditions:
            condition_normal = cond.GetGeometry().Normal()
            pressure_coefficient = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
            for node in cond.GetNodes():
                added_force = condition_normal*(pressure_coefficient/2.0)*dynamic_pressure
                nodal_force = node.GetSolutionStepValue(KratosMultiphysics.REACTION) + added_force
                largest_force = max(largest_force, nodal_force.norm_2())
                node.SetSolutionStepValue(KratosMultiphysics.REACTION, nodal_force)

        if largest_force<1e-300:
            KratosMultiphysics.Logger.PrintWarning('Warning','all reaction forces are zero')
        total_force = KratosMultiphysics.VariableUtils().SumHistoricalNodeVectorVariable(KratosMultiphysics.REACTION, self.body_model_part, 0)

        KratosMultiphysics.Logger.PrintInfo('ComputeForcesOnNodesProcess','Lift Force = ', total_force[1])
        KratosMultiphysics.Logger.PrintInfo('ComputeForcesOnNodesProcess','Drag Force = ', total_force[0])
        KratosMultiphysics.Logger.PrintInfo('ComputeForcesOnNodesProcess','Side Force = ', total_force[2])

        if self.create_output_file:
            with open("cl_points_with_lift.dat", 'w') as cl_file:
                cl_file.write('{0:15.12f}'.format(total_force[1]/dynamic_pressure))
