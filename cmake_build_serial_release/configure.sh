#!/bin/sh
rm CMakeCache.txt
rm *.cmake
rm -rf CMakeFiles

cmake ..                                                                                                \
-DCMAKE_C_COMPILER=gcc                                                                                  \
-DCMAKE_CXX_COMPILER=g++                                                                                \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++11 -march=native -ftree-vectorize -funroll-loops -ffast-math -Wall -Wno-unused-local-typedefs" \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -march=native -funroll-loops -ffast-math -Wall -Wno-unused-local-typedefs" \
-DCMAKE_BUILD_TYPE=Release                                                                       \
-DCMAKE_INSTALL_RPATH="${KRATOS_PATH_ENV}/libs"                                                    \
-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE                                                                \
-DPYTHON_LIBRARY="/usr/lib/python3.6/config-3.6m-x86_64-linux-gnu/libpython3.6m.so"                     \
-DPYTHON_INCLUDE_DIR="/usr/include/python3.6"                                                           \
-DINSTALL_EMBEDDED_PYTHON=ON                                                                            \
-DBOOST_ROOT="${HOME}/software/boost_1_67_0"                                              \
-DBLAS_LIBRARIES="/usr/lib/x86_64-linux-gnu/blas/libblas.so"                                                          \
-DLAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/lapack/liblapack.so"                                                       \
-DMESHING_APPLICATION=ON                                                                                \
-DMESH_MOVING_APPLICATION=ON                                                                                \
-DEIGEN_SOLVERS_APPLICATION=ON                                                                                \
-DEIGEN_ROOT="${HOME}/software/eigen"								\
-DEXTERNAL_SOLVERS_APPLICATION=ON                                                                       \
-DMAPPING_APPLICATION=ON                                                                                \
-DFLUID_DYNAMICS_APPLICATION=ON                                                                        \
-DSTRUCTURAL_MECHANICS_APPLICATION=ON                                                                        \
-DCHIMERA_APPLICATION=ON										\
-DEMPIRE_APPLICATION=ON                                                                                \
-DEXCLUDE_ITSOL=ON                                                                                      \
-DMETIS_APPLICATION=OFF                                                                                  \
-DPARMETIS_ROOT_DIR="${HOME}/software/ParMETIS/ParMetis-3.2.0"                                  \
-DTRILINOS_APPLICATION=OFF                                                                               \
-DTRILINOS_ROOT="${HOME}/software/Trilinos/trilinos-12.10.1-Source_install"                             \
-DINCLUDE_PASTIX=OFF                                                                                     \
-DSCOTCH_INSTALL_DIR="${HOME}/software/Scotch/scotch_6.0.4_install/lib"                                     \
-DPASTIX_INSTALL_DIR="${HOME}/software/PaStiX/pastix_5.2.3/install"                                     \

make -j6 install
