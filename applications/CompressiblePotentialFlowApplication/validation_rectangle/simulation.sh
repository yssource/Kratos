#!/bin/bash
#
#run this file using the command:
# bash simulation.sh

#Setting paths
source scripts/preamble.sh

#Setting the parameters
source scripts/parameters.sh
source generate_mdpas/set_parameters.sh
cd generate_mdpas/

#Run salome: generate geometry and mesh
python3 runSalome.py

#Convert salomes mesh into mdpa
rm $Work_Dir/mdpas/*
python3 use_converter.py

##Removing plots
cd ..
source scripts/removing_before_kratos_run.sh

#Run Kratos
unbuffer python3 MeshRefinement.py 2>&1 | tee $NEWFILE

source generate_mdpas/unset_parameters.sh

##Remove plots after run, run latex and copy results
source scripts/run_latex.sh

