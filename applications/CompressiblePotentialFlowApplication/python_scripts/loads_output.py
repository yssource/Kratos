def write_header(work_dir):
    refinement_file = open(work_dir + "mesh_refinement_loads.dat",'w') 
    refinement_file.write("FULL POTENTIAL APPLICATION LOADS FILE\n\n")
    refinement_file.write('%4s %6s %15s %15s %15s %15s %15s %15s %15s %15s\n\n' % 
                          ("Case", "AOA", "FF_MS", "A_MS", "# Nodes", "Cl_u", "Cl_l", "Cd_u", "Cd_l", "Rz"))
    refinement_file.flush()

def write_header_all_cases(work_dir):
    aoa_file = open(work_dir + "plots/results/all_cases.dat",'w') 
    aoa_file.write("FULL POTENTIAL APPLICATION ALL CASES LOADS FILE\n\n")
    aoa_file.write('%4s %6s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n\n' % 
                          ("Case", "AOA", "FF_MS", "A_MS", "Cl_jump", "# Nodes", "Cl_u", "Cl_l", "Cd_u", "Cd_l", "Rz"))
    aoa_file.flush()


def write_case(case, AOA, FarField_MeshSize, Airfoil_MeshSize,work_dir):
    refinement_file = open(work_dir + "mesh_refinement_loads.dat",'a')
    refinement_file.write('{0:4d} {1:6.2f} {2:15.2f} {3:15.2e}'.format(case, AOA, FarField_MeshSize, Airfoil_MeshSize))
    refinement_file.flush()

    aoa_file = open(work_dir + "plots/results/all_cases.dat",'a')
    aoa_file.write('{0:4d} {1:6.2f} {2:15.2f} {3:15.2e}'.format(case, AOA, FarField_MeshSize, Airfoil_MeshSize))
    aoa_file.flush()

def write_figures(cl_data_directory_name, AOA, work_dir):
    figures_file = open(work_dir + 'plots/figures.tex', 'a')
    figures_file.write('\n\pgfplotsset{table/search path={' + cl_data_directory_name + '},}\n\n' +
                       '\\begin{figure}\n' +
                       '\t\centering\n' +
                       '\t\input{' + cl_data_directory_name + '/cl.tikz}\n' +
                       '\t\caption{$\\alpha = ' + str(AOA) + '\degree$}\n' +
                       '\t\label{fig:cl_AOA_' + str(AOA) + '}\n' +
                       '\end{figure}\n'
                       )
    figures_file.flush()

def write_figures_cd(cd_data_directory_name, AOA, work_dir):
    with open(work_dir + 'plots/cd/figures_cd.tex', 'a') as cd_figures_file:
        cd_figures_file.write('\n\pgfplotsset{table/search path={' + cd_data_directory_name + '},}\n\n' +
                           '\\begin{figure}\n' +
                           '\t\centering\n' +
                           '\t\input{' + cd_data_directory_name + '/cd.tikz}\n' +
                           '\t\caption{$\\alpha = ' + str(AOA) + '\degree$}\n' +
                           '\t\label{fig:cd_AOA_' + str(AOA) + '}\n' +
                           '\end{figure}\n'
                           )
        cd_figures_file.flush()

def write_cl(cl,work_dir):
    cl_aoa_file = open(work_dir + 'plots/aoa/cl_aoa.dat','a')
    cl_aoa_file.write('{0:15f}\n'.format(cl))
    cl_aoa_file.flush()

def write_cp_figures(cp_data_directory_name, AOA, case, Airfoil_MeshSize,  FarField_MeshSize, work_dir):
    figures_file = open(work_dir + 'plots/cp/figures.tex', 'w')
    figures_file.write('\n\pgfplotsset{table/search path={' + cp_data_directory_name + '},}\n\n' +
                       '\\begin{figure}\n' +
                       '\t\centering\n' +
                       '\t\input{' + cp_data_directory_name + '/cp.tikz}\n' +
                       '\t\caption{$\\alpha = ' + str(AOA) + '\degree$, case = ' + str(case) + 
                       ' Far field mesh size = ' + str(FarField_MeshSize) + ' Airfoil mesh size = ' + str(Airfoil_MeshSize)+ '}\n' +
                       '\t\label{fig:cp_AOA_' + str(AOA) + '}\n' +
                       '\end{figure}\n'
                       )
    figures_file.flush()

