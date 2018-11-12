#!/bin/bash
#
#run this file using the command:
# bash simulation.sh

#Going to the current directory
echo "The previous current working directory: $PWD"

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
cd $SCRIPTPATH

echo "The current working directory: $PWD"

#Setting paths
source scripts/preamble.sh

#Setting the parameters
source scripts/parameters.sh
source generate_mdpas/set_parameters.sh
cd generate_mdpas/

#Run salome: generate geometry and mesh
#python3 runSalome.py

#Convert salomes mesh into mdpa
#rm $Work_Dir/mdpas/*
#python3 use_converter.py

cd ..

#Run Kratos
source scripts/wake_implementation01.sh
source scripts/wake_implementation02.sh

source generate_mdpas/unset_parameters.sh



