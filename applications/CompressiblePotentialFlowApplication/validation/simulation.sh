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
#source scripts/wake_implementation01.sh
#source scripts/wake_implementation02.sh
#source scripts/wake_implementation03.sh
#source scripts/wake_implementation04.sh
#source scripts/wake_implementation05.sh
#source scripts/wake_implementation06.sh
#source scripts/wake_implementation07.sh
#source scripts/wake_implementation08.sh
#source scripts/wake_implementation09.sh
#source scripts/wake_implementation10.sh

#source scripts/wake_implementation15.sh
#source scripts/wake_implementation16.sh
#source scripts/wake_implementation17.sh
#source scripts/wake_implementation18.sh

#TO BE RUN
#source scripts/wake_implementation11.sh
#source scripts/wake_implementation12.sh
#source scripts/wake_implementation13.sh
#source scripts/wake_implementation14.sh

#source scripts/wake_implementation19.sh
#source scripts/wake_implementation23.sh
#source scripts/wake_implementation21.sh

#source scripts/wake_implementation25.sh
source scripts/wake_implementation27.sh

source generate_mdpas/unset_parameters.sh



