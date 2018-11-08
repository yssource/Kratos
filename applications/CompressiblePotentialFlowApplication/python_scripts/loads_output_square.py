def write_header(work_dir):
    refinement_file = open(work_dir + "mesh_refinement_loads.dat",'w') 
    refinement_file.write("FULL POTENTIAL APPLICATION LOADS FILE\n\n")
    refinement_file.write('%4s %15s %15s %15s %15s \n\n' % 
                          ("Case", "#Segments", "MeshSize", "# Nodes", "Error"))
    refinement_file.flush()

def write_header_all_cases(work_dir):
    aoa_file = open(work_dir + "plots/results/all_cases.dat",'w') 
    aoa_file.write("FULL POTENTIAL APPLICATION ALL CASES LOADS FILE\n\n")
    aoa_file.write('%4s %15s %15s %15s %15s \n\n' % 
                          ("Case", "#Segments", "MeshSize", "# Nodes", "Error"))
    aoa_file.flush()


def write_case(case, NumberOfSegments, MeshSize, work_dir):
    refinement_file = open(work_dir + "mesh_refinement_loads.dat",'a')
    refinement_file.write('{0:4d} {1:15.1e} {2:15.1e}'.format(case, NumberOfSegments, MeshSize))
    refinement_file.flush()

    aoa_file = open(work_dir + "plots/results/all_cases.dat",'a')
    aoa_file.write('{0:4d} {1:15.1e} {2:15.1e}'.format(case, NumberOfSegments, MeshSize))
    aoa_file.flush()
    
def write_figures_energy(energy_data_directory_name, AOA, work_dir):
    with open(work_dir + 'plots/relative_error_energy_norm/figures_energy_h.tex', 'a') as energy_h_figures_file:
        energy_h_figures_file.write('\n\pgfplotsset{table/search path={' + energy_data_directory_name + '},}\n\n' +
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + energy_data_directory_name + '/energy_h.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$}\n' +
                           '\t\label{fig:energy_h_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n'
                           )
        energy_h_figures_file.flush()

    with open(work_dir + 'plots/relative_error_energy_norm/figures_energy_n.tex', 'a') as energy_n_figures_file:
        energy_n_figures_file.write('\n\pgfplotsset{table/search path={' + energy_data_directory_name + '},}\n\n' +
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + energy_data_directory_name + '/energy_n.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$}\n' +
                           '\t\label{fig:energy_n_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n'
                           )
        energy_n_figures_file.flush()

    with open(work_dir + 'plots/relative_error_energy_norm/figures_energy_variant_h.tex', 'a') as energy_variant_h_figures_file:
        energy_variant_h_figures_file.write('\n\pgfplotsset{table/search path={' + energy_data_directory_name + '},}\n\n' +
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + energy_data_directory_name + '/energy_variant_h.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$}\n' +
                           '\t\label{fig:energy_variant_h_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n'
                           )
        energy_variant_h_figures_file.flush()

    with open(work_dir + 'plots/relative_error_energy_norm/figures_energy_variant_n.tex', 'a') as energy_variant_n_figures_file:
        energy_variant_n_figures_file.write('\n\pgfplotsset{table/search path={' + energy_data_directory_name + '},}\n\n' +
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + energy_data_directory_name + '/energy_variant_n.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$}\n' +
                           '\t\label{fig:energy_variant_n_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n'
                           )
        energy_variant_n_figures_file.flush()

def write_figures_condition(condition_data_directory_name, AOA, work_dir):
    with open(work_dir + 'plots/condition_number/figures_condition.tex', 'a') as condition_figures_file:
        condition_figures_file.write('\n\pgfplotsset{table/search path={' + condition_data_directory_name + '},}\n\n' +
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + condition_data_directory_name + '/condition.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$}\n' +
                           '\t\label{fig:condition_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n'
                           )
        condition_figures_file.flush()

def write_figures_far_field(far_field_data_directory_name, AOA, case, Airfoil_MeshSize,  FarField_MeshSize, work_dir):
    with open(work_dir + 'plots/far_field/figures_far_field_x.tex', 'w') as figures_file_x:
        figures_file_x.write('\n\pgfplotsset{table/search path={' + far_field_data_directory_name + '},}\n\n' +
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + far_field_data_directory_name + '/velocity_norm_x.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$, case = ' + str(case) + 
                           ' Far field mesh size = ' + str(FarField_MeshSize) + ' Airfoil mesh size = ' + str(Airfoil_MeshSize)+ '}\n' +
                           '\t\label{fig:velocity_norm_x_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n\n'
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + far_field_data_directory_name + '/velocity_u_x.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$, case = ' + str(case) + 
                           ' Far field mesh size = ' + str(FarField_MeshSize) + ' Airfoil mesh size = ' + str(Airfoil_MeshSize)+ '}\n' +
                           '\t\label{fig:velocity_u_x_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n\n'
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + far_field_data_directory_name + '/velocity_v_x.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$, case = ' + str(case) + 
                           ' Far field mesh size = ' + str(FarField_MeshSize) + ' Airfoil mesh size = ' + str(Airfoil_MeshSize)+ '}\n' +
                           '\t\label{fig:velocity_v_x_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n'
                           )
        figures_file_x.flush()

    with open(work_dir + 'plots/far_field/figures_far_field_y.tex', 'w') as figures_file_y:
        figures_file_y.write('\n\pgfplotsset{table/search path={' + far_field_data_directory_name + '},}\n\n' +
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + far_field_data_directory_name + '/velocity_norm_y.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$, case = ' + str(case) + 
                           ' Far field mesh size = ' + str(FarField_MeshSize) + ' Airfoil mesh size = ' + str(Airfoil_MeshSize)+ '}\n' +
                           '\t\label{fig:velocity_norm_y_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n\n'
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + far_field_data_directory_name + '/velocity_u_y.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$, case = ' + str(case) + 
                           ' Far field mesh size = ' + str(FarField_MeshSize) + ' Airfoil mesh size = ' + str(Airfoil_MeshSize)+ '}\n' +
                           '\t\label{fig:velocity_u_y_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n\n'
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + far_field_data_directory_name + '/velocity_v_y.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$, case = ' + str(case) + 
                           ' Far field mesh size = ' + str(FarField_MeshSize) + ' Airfoil mesh size = ' + str(Airfoil_MeshSize)+ '}\n' +
                           '\t\label{fig:velocity_v_y_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n'
                           )
        figures_file_y.flush()

