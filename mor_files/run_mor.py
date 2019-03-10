from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

# import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, atan, cos, sin, exp, pi, ceil
from random import randrange
import numpy as np

import sys

def _add_variables(mp):
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)     
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)   
    mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)

def _apply_material_properties(mp,num):

    cl = StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()
    mp.GetProperties()[num].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)
    #aluminum
    # mp.GetProperties()[num].SetValue(KratosMultiphysics.DENSITY,2698.9)
    # mp.GetProperties()[num].SetValue(KratosMultiphysics.YOUNG_MODULUS,0.70e11)
    # mp.GetProperties()[num].SetValue(KratosMultiphysics.THICKNESS,0.002)
    # mp.GetProperties()[num].SetValue(KratosMultiphysics.CROSS_AREA,0.04)
    # mp.GetProperties()[num].SetValue(KratosMultiphysics.POISSON_RATIO,0.34)
    a = 0.1
    mp.GetProperties()[num].SetValue(KratosMultiphysics.DENSITY,7850)
    mp.GetProperties()[num].SetValue(KratosMultiphysics.YOUNG_MODULUS,2.1e11)
    mp.GetProperties()[num].SetValue(KratosMultiphysics.POISSON_RATIO,0.3)
    mp.GetProperties()[num].SetValue(StructuralMechanicsApplication.CROSS_AREA,a**2)
    mp.GetProperties()[num].SetValue(StructuralMechanicsApplication.I22,a**4/12)
    mp.GetProperties()[num].SetValue(StructuralMechanicsApplication.I33,a**4/12)
    mp.GetProperties()[num].SetValue(StructuralMechanicsApplication.TORSIONAL_INERTIA,a**4/12)

def _solve_eigen(model,echo,output=False):
    if output:
        import process_factory
        output_settings = KratosMultiphysics.Parameters("""
            {
                "output_process_list":[{
                "python_module"   : "postprocess_eigenvalues_process",
                "kratos_module"   : "KratosMultiphysics.StructuralMechanicsApplication",
                "help"                  : "This process postprocces the eigen values for GiD",
                "process_name"          : "PostProcessEigenvaluesProcess",
                "Parameters" : {
                    "result_file_name" : "modal",
                    "result_file_format_use_ascii" : false,
                    "computing_model_part_name"   : "mor",
                    "animation_steps"   :  20,
                    "list_of_result_variables" : ["DISPLACEMENT"],
                    "label_type" : "angular_frequency"
                }
                }]
            }
            """)
        list_of_processes = process_factory.KratosProcessFactory(model).ConstructListOfProcesses(output_settings["output_process_list"])
        for process in list_of_processes:
            print(process)
            process.ExecuteInitialize()
            process.ExecuteBeforeSolutionLoop()

    eigensolver_settings = KratosMultiphysics.Parameters("""
        {
            
                "solver_type": "FEAST",
                "print_feast_output": true,
                "perform_stochastic_estimate": false,
                "solve_eigenvalue_problem": true,
                "lambda_min": 1.0e2,
                "lambda_max": 8.0e7,
                "search_dimension": 20,
                "linear_solver_settings": {
                    "solver_type": "pastix"
                }
            
        }
        """)

    feast_system_solver_settings = eigensolver_settings["linear_solver_settings"]
    feast_system_solver = ExternalSolversApplication.PastixComplexSolver(feast_system_solver_settings)
    linear_solver = ExternalSolversApplication.FEASTSolver(eigensolver_settings, feast_system_solver)
    builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)

    eigen_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
    eig_strategy = StructuralMechanicsApplication.EigensolverStrategy(model.GetModelPart('mor'),
                                                                  eigen_scheme,
                                                                  builder_and_solver)

    eig_strategy.SetEchoLevel(echo)
    eig_strategy.Solve()

    if output:
        for process in list_of_processes:
            process.ExecuteFinalize()
                                                

def _set_and_fill_buffer(mp,buffer_size,delta_time):
    # Set buffer size
    mp.SetBufferSize(buffer_size)

    # Fill buffer
    time = mp.ProcessInfo[KratosMultiphysics.TIME]
    time = time - delta_time * (buffer_size)
    mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
    for size in range(0, buffer_size):
        step = size - (buffer_size -1)
        mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
        time = time + delta_time
        #delta_time is computed from previous time in process_info
        mp.CloneTimeStep(time)

    mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

def _add_dofs_mp(mp):
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X,mp)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y,mp)
    KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z,mp)

def _define_gid_io(fname,mesh):
    gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
    #gid_mode = GiDPostMode.GiD_PostAscii
    ##multifile = MultiFileFlag.MultipleFiles
    multifile = KratosMultiphysics.MultiFileFlag.SingleFile
    deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
    #deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
    write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
    #write_conditions = WriteConditionsFlag.WriteConditionsOnly
    #write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

    gid_io = KratosMultiphysics.GidIO(fname,gid_mode,multifile,deformed_mesh_flag, write_conditions)

    gid_io.InitializeMesh(0.0)
    gid_io.WriteMesh(mesh)
    gid_io.FinalizeMesh()

    return gid_io

def run_analysis():
    #####################
    #parameters
    echo_level = 3
    solving_type = 'mor'

    #RAD/s!!!
    exfreq = 1
    max_exfreq = 5000
    df = .5 * 2.0 * pi
    ######################

    Model = KratosMultiphysics.Model()
    mp = Model.CreateModelPart('mor')

    _add_variables(mp)
    _apply_material_properties(mp,0)

    #create geometry
    n01 = mp.CreateNewNode(1,0.0,0.0,0.0)
    n = 100
    for id in range(1, n+1):
        mp.CreateNewNode(id+1, 2.5*id/n, 0.0,0.0)

    end_node = mp.Nodes[len(mp.Nodes)]

    for id in range(1, n+1):
        mp.CreateNewElement("CrBeamElement3D2N", id, [id, id+1], mp.GetProperties()[0])


    _add_dofs_mp(mp)

    n01.Fix(KratosMultiphysics.DISPLACEMENT_X)
    n01.Fix(KratosMultiphysics.DISPLACEMENT_Y)
    n01.Fix(KratosMultiphysics.DISPLACEMENT_Z)
    n01.Fix(KratosMultiphysics.ROTATION_X)
    n01.Fix(KratosMultiphysics.ROTATION_Y)
    n01.Fix(KratosMultiphysics.ROTATION_Z)

    mp.CreateNewCondition("PointLoadCondition3D1N", 1, [end_node.Id], mp.GetProperties()[0])
    end_node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,[0,-10,0])

    ########
    ########

    if solving_type == 'mor':
        linear_solver_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "SuperLUSolver",
            "max_iteration": 500,
            "tolerance": 1e-9,
            "scaling": false,
            "symmetric_scaling": true,
            "verbosity": 1
        }
        """)

        linear_solver = KratosMultiphysics.LinearSolverFactory().Create(linear_solver_settings)
        builder_and_solver = StructuralMechanicsApplication.MassAndStiffnessBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-12,1e-6)
        
        move_mesh_flag = False
        
        off_strategy = StructuralMechanicsApplication.MorOfflineStrategy(mp, 
                                                                        scheme, 
                                                                        linear_solver,
                                                                        [1, 2, 3],
                                                                        move_mesh_flag)
        strategy = StructuralMechanicsApplication.MorOnlineStrategy(mp, 
                                                                    scheme, 
                                                                    linear_solver,
                                                                    builder_and_solver,
                                                                    off_strategy,
                                                                    [80, 500, 1500, 3000],
                                                                    move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Check()

        result = np.zeros((1,2))

        while(exfreq <= max_exfreq):
            exfreq = exfreq + df
            # print(exfreq)
            mp.CloneTimeStep(exfreq)
            strategy.Solve()
            
            disp = mp.Nodes[end_node.Id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0)
            temp_result = np.asarray([exfreq, abs(disp)])
            result = np.vstack([result, temp_result])

        np.savetxt("disp_mor.csv", np.asarray(result), delimiter=',')

    elif solving_type == 'harmonic':
        _solve_eigen(Model,echo_level,False)
        eigenvector = mp.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]
        # print([ev for ev in mp.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]])
        current_vals = [sqrt(ev) for ev in mp.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]]
        eig_result = [n/2/pi for n in current_vals]
        print('ev in rad/s: ' + str(current_vals))
        print('ev in 1/s: ' + str([ev/2/pi for ev in current_vals]) )

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(KratosMultiphysics.LinearSolver())
        eigen_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        harmonic_strategy = StructuralMechanicsApplication.HarmonicAnalysisStrategy(mp, eigen_scheme, builder_and_solver,False)   
        
        # mp.ProcessInfo.SetValue(StructuralMechanicsApplication.RAYLEIGH_ALPHA, 0.01)
        # mp.ProcessInfo.SetValue(StructuralMechanicsApplication.RAYLEIGH_BETA, 0.01)
        # mp.ProcessInfo.SetValue(StructuralMechanicsApplication.SYSTEM_DAMPING_RATIO, 0.05)
        harmonic_strategy.SetEchoLevel(echo_level)
    
        result = np.zeros((1,2))

        while(exfreq <= max_exfreq):
            exfreq = exfreq + df
            # print(exfreq)
            mp.CloneTimeStep(exfreq)
            harmonic_strategy.Solve()
            
            disp = mp.Nodes[end_node.Id].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0)
            temp_result_x = np.asarray([exfreq, abs(disp)])
            result = np.vstack([result, temp_result_x])

        np.savetxt("disp_harm.csv", np.asarray(result), delimiter=',')
            
    else:
        pass

if __name__ == '__main__':
    run_analysis()
