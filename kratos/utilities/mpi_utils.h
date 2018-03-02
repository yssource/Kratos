//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//                    
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes
#include "mpi.h"

// Project includes
#include "includes/exception.h"

namespace Kratos {

class MpiUtils {
public:
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