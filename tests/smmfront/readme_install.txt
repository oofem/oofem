# INSTALLATION OF ALL PARTS ON UBUNTU 20.04

# needs clang
sudo apt install clang

# install tfel
mkdir TFEL && cd TFEL
git clone https://github.com/thelfer/tfel.git
mkdir build && cd build
cmake ../tfel -DUSE_PYTHON_EXTENSION="OFF" -DCMAKE_CXX_COMPILER=/usr/bin/clang++-10 -DCMAKE_C_COMPILER=/usr/bin/clang-10 -Dlocal-castem-header=ON -DCMAKE_INSTALL_PREFIX=~/.local
make -j 8
sudo make install

# install MGIS
git clone https://github.com/thelfer/MFrontGenericInterfaceSupport.git
cd MFrontGenericInterfaceSupport
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -Denable-c-bindings=ON -Denable-fortran-bindings=ON -Denable-python-bindings=ON -DPython_ADDITIONAL_VERSIONS=3.8
make -j 8
sudo make install

# install oofem
git clone https://github.com/oofem/oofem.git
cd oofem && mkdir build && cd build
cmake -DUSE_MFRONT=ON -DMFrontGenericInterface_DIR=/usr/local/share/mgis/cmake/ ..
make -j 8
# go to the mfront test directory
cd ../tests/smmfront/plasticity
# create the mfront material shared library for generic interface
CC=clang CXX=clang++ mfront --obuild --interface=generic IsotropicLinearHardeningPlasticity.mfront
# go to the build directory
cd ../../../build
# run the tests
ctest