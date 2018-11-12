#Removing files
#rm $Input_Dir/output_salome/*
#rm $Input_Dir/mdpas/*
rm cp*
rm main*

#Running Latex
cd $Work_Dir/plots/cl
pdflatex -interaction=batchmode main_cl.tex > main_cl_out.txt
cd $Work_Dir/plots/cl_error
pdflatex -interaction=batchmode main_cl_error.tex > main_cl_error_out.txt
cd $Work_Dir/plots/relative_error_energy_norm
pdflatex -interaction=batchmode main_energy_h.tex > main_energy_h_out.txt
pdflatex -interaction=batchmode main_energy_n.tex > main_energy_n_out.txt
pdflatex -interaction=batchmode main_energy_variant_h.tex > main_energy_variant_h_out.txt
pdflatex -interaction=batchmode main_energy_variant_n.tex > main_energy_variannt_n_out.txt
cd $Work_Dir/plots/cd
pdflatex -interaction=batchmode main_cd.tex > main_cd_out.txt
cd $Work_Dir/plots/aoa/
pdflatex -interaction=batchmode cl_aoa.tex > main_aoa_out.txt
cd $Work_Dir/plots/condition_number/
pdflatex -interaction=batchmode main_condition.tex > main_condition_out.txt

#Copying results
DIRECTORY=/media/inigo/10740FB2740F9A1C/Implementations_testing/05_MeshRefinement
mkdir -p ${DIRECTORY}_${DATE}_${GITBRANCH}_${OUTPUTNAME}
mkdir -p ${DIRECTORY}_${DATE}_${GITBRANCH}_${OUTPUTNAME}/output_gid

cp -r $Work_Dir/plots/ ${DIRECTORY}_${DATE}_${GITBRANCH}_${OUTPUTNAME}
cp -r /media/inigo/10740FB2740F9A1C/Outputs/05_MeshRefinement/A* ${DIRECTORY}_${DATE}_${GITBRANCH}_${OUTPUTNAME}/output_gid

