from KratosMultiphysics import *
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlowApplication
from numpy import *
import itertools
import loads_output
import math

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeLiftProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "far_field_model_part_name"   : "PotentialWallCondition2D_Far_field_Auto1",
                "mesh_id": 0,
                "velocity_infinity": [1.0,0.0,0],
                "angle_of_attack": 0.0,
                "meshsize": 1.0,
                "energy_reference": 1.0,
                "potential_energy_reference": 1.0,
                "reference_area": 1
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.far_field_model_part = Model[settings["far_field_model_part_name"].GetString()]

        self.mesh_size = settings["meshsize"].GetDouble()
        print('mesh size =', self.mesh_size)

        self.energy_reference = 2*1e8/3 #settings["energy_reference"].GetDouble()
        print('self.energy_reference = ', self.energy_reference)

        self.total_potential_energy_reference = -5*1e5 #settings["potential_energy_reference"].GetDouble()
        #print('self.total_potential_energy_reference = ', self.total_potential_energy_reference)

        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part,
                                                                       self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Neigbour search tool instance
        AvgElemNum = 10
        AvgNodeNum = 10
        nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fluid_model_part, AvgElemNum, AvgNodeNum)
        # Find neighbours
        nodal_neighbour_search.Execute()

    def ExecuteFinalizeSolutionStep(self):
        print('COMPUTE LIFT')

        self.work_dir = '/home/inigo/simulations/naca0012/07_salome/06_Rectangle/'
        
    
        #compute the internal energy norm and relative error
        internal_energy_sum = 0.0
        for element in self.fluid_model_part.Elements:
            internal_energy_sum += element.GetValue(KratosMultiphysics.INTERNAL_ENERGY)

        relative_error_energy_norm = math.sqrt(abs(internal_energy_sum - self.energy_reference)/abs(self.energy_reference))

        external_energy_sum = 0.0
        for cond in self.far_field_model_part.Conditions:
            external_energy_sum += cond.GetValue(EXTERNAL_ENERGY)

        total_potential_energy = internal_energy_sum - external_energy_sum

        relative_error_energy_norm_variant = math.sqrt(abs(total_potential_energy - self.total_potential_energy_reference)/abs(self.energy_reference))

        NumberOfNodes = self.fluid_model_part.NumberOfNodes()

        print('\n internal_energy_sum = ', internal_energy_sum)
        print('\n external_energy_sum = ', external_energy_sum)
        print('\n total_potential_energy = ', total_potential_energy)

        print('\n absolute_error = ', abs(internal_energy_sum - self.energy_reference))

        print('\n relative_error_energy_norm = ', relative_error_energy_norm)
        #print('\n relative_error_energy_norm_variant = ', relative_error_energy_norm_variant)
    
        with open (self.work_dir + "mesh_refinement_loads.dat",'a') as loads_file:
            loads_file.write('{0:16.1e} {1:15.1e} \n'.format(NumberOfNodes, relative_error_energy_norm))
            loads_file.flush()

        with open(self.work_dir + "plots/results/all_cases.dat",'a') as aoa_file:
            aoa_file.write('{0:16.1e} {1:15.1e} \n'.format(NumberOfNodes, relative_error_energy_norm))
            aoa_file.flush()

        energy_h_results_file_name = self.work_dir + "plots/relative_error_energy_norm/data/energy/energy_h_results.dat"
        with open(energy_h_results_file_name,'a') as energy_h_file:
            energy_h_file.write('{0:16.6e} {1:16.6e}\n'.format(self.mesh_size, relative_error_energy_norm))
            energy_h_file.flush()

        energy_n_results_file_name = self.work_dir + "plots/relative_error_energy_norm/data/energy/energy_n_results.dat"
        with open(energy_n_results_file_name,'a') as energy_n_file:
            energy_n_file.write('{0:16.6e} {1:16.6e}\n'.format(NumberOfNodes, relative_error_energy_norm))
            energy_n_file.flush()

        energy_variant_h_file_name = self.work_dir + "plots/relative_error_energy_norm/data/energy/energy_variant_h_results.dat"
        with open(energy_variant_h_file_name,'a') as energy_variant_h_file:
            energy_variant_h_file.write('{0:16.6e} {1:16.6e}\n'.format(self.mesh_size, relative_error_energy_norm_variant))
            energy_variant_h_file.flush()

        energy_variant_n_file_name = self.work_dir + "plots/relative_error_energy_norm/data/energy/energy_variant_n_results.dat"
        with open(energy_variant_n_file_name,'a') as energy_variant_n_file:
            energy_variant_n_file.write('{0:16.6e} {1:16.6e}\n'.format(NumberOfNodes, relative_error_energy_norm_variant))
            energy_variant_n_file.flush()