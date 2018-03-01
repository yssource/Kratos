#include "mpi.h"
#include "includes/excpetion.h"

namespace Kratos {

class MpiUtils {

    static void HandleError(int mpi_erno) {
        if(mpi_erno != MPI_SUCCESS) {
            char estring[MPI_MAX_ERROR_STRING];
            int eclass;
            int elen;

            MPI_Error_class(mpi_erno, &eclass);
            MPI_Error_string(mpi_erno, estring, &elen);

            std::string mpiErrorString(estring, elen);

            KRATOS_ERROR << mpiErrorString << std::endl;
        }
    }
};
}