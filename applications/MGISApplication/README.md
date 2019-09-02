 
# MGISApplication
## Dependencies

### TFEL

	git clone https://github.com/thelfer/tfel.git 

with `configure.sh`:

~~~
#!/bin/sh

clear

#you may want to decomment this the first time you compile
rm CMakeCache.txt
rm *.cmake
rm -rf CMakeFiles\

cmake ..                                                                            \
-DCMAKE_C_COMPILER=/usr/bin/gcc                                                     \
-DCMAKE_CXX_COMPILER=/usr/bin/g++                                                   \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 -fPIC "                     \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3  -fPIC "                                   \
-DCMAKE_BUILD_TYPE=Release                                                          \

#decomment this to have it verbose
# make VERBOSE=1 -j4
make -j4
sudo make install

~~~

### MGIS

 	git clone https://github.com/thelfer/MFrontGenericInterfaceSupport.git 

with `configure.sh`:

~~~
#!/bin/sh

clear

#you may want to decomment this the first time you compile
rm CMakeCache.txt
rm *.cmake
rm -rf CMakeFiles\

cmake ..                                                                            \
-DCMAKE_C_COMPILER=/usr/bin/gcc                                                     \
-DCMAKE_CXX_COMPILER=/usr/bin/g++                                                   \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -fPIC "                                \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3  -fPIC "                                   \
-DCMAKE_BUILD_TYPE=Release                                                          \
-Denable-c-bindings=ON                                                              \
-Denable-fast-math=ON                                                               \

#decomment this to have it verbose
# make VERBOSE=1 -j4
make -j4
sudo make install

~~~