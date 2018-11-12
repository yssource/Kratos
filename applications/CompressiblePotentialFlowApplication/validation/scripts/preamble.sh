#Resetting paths
export PYTHONPATH=""
export LD_LIBRARY_PATH=""

#Setting paths
source /home/inigo/Documents/paths/salomeConverter.sh
source /home/inigo/Documents/paths/kratosWake01.sh
source /home/inigo/intel/mkl/bin/mklvars.sh intel64 lp64

GITBRANCH=$(git symbolic-ref HEAD | sed -e 's,.*/\(.*\),\1,')
Input_Dir=/home/inigo/simulations/naca0012/07_salome/07_MeshRefinement

DATE=`date '+%Y%m%d_%H%M%S'`
FILE=${Input_Dir}/plots/output_terminal.txt
NAME=${FILE%.*}
EXT=${FILE#*.}
NEWFILE=${NAME}_${DATE}_${GITBRANCH}.${EXT}