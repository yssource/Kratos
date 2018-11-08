#Removing files
#rm $Input_Dir/output_salome/*
#rm $Input_Dir/mdpas/*
rm main*
rm val*

#Running Latex
cd $Work_Dir/plots/relative_error_energy_norm
pdflatex -interaction=batchmode main_energy.tex > main_energy_out.txt
#pdflatex  main_energy.tex
cd $Work_Dir/plots/condition_number/
pdflatex -interaction=batchmode main_condition.tex > main_condition_out.txt

#Copying results
DIRECTORY=/media/inigo/10740FB2740F9A1C/Old_Outputs/06_Rectangle
mkdir -p ${DIRECTORY}_${DATE}_${GITBRANCH}
mkdir -p ${DIRECTORY}_${DATE}_${GITBRANCH}/output_gid
cd ../..
cp -r $Work_Dir/plots/ ${DIRECTORY}_${DATE}_${GITBRANCH}
cp -r /media/inigo/10740FB2740F9A1C/Outputs/06_Rectangle/A* ${DIRECTORY}_${DATE}_${GITBRANCH}/output_gid