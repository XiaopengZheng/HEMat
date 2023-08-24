# HEMat
## Step 1. Building and installing HElib
### General prerequisites
* Ubuntu 22.04 LTS
* GNU make >= 4.3
* g++ >= 11.3.0
* cmake >= 3.22
* pthreads
* git >=2.36
### Instructions
1. Installing basis tools:
```
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install patchelf
sudo apt-get install m4
sudo apt-get install git
sudo apt-get install cmake
```
2. Clone HElib from gitbub:
```
sudo git clone https://github.com/homenc/HElib.git
```
3. Install HElib
```
cd HElib
sudo mkdir build
cd build
sudo cmake -DPACKAGE_BUILD=ON -DCMAKE_INSTALL_PREFIX=/home/usr/helib_install ..
sudo make -j16 -lpthread
sudo make install
```
Here HElib is installed in `/home/usr/helib_install`. 

## Step 2. Install Open MPI 
### Instructions

1、Download openmpi from `https://www.open-mpi.org/software/ompi/v4.1/`. For example, download `openmpi-4.1.5.tar.gz`.

2. Install Open MPI (Open MPI is installed in `/usr/local/openmpi` )
```
sudo tar zxf openmpi-4.1.5.tar.gz
cd openmpi-4.1.5
./configure --prefix=/usr/local/openmpi
sudo make
sudo make install
```

3. Set the environment variables
```
MPI_HOME=/usr/local/openmpi
export PATH=${MPI_HOME}/bin:$PATH
export LD_LIBRARY_PATH=${MPI_HOME}/lib:$LD_LIBRARY_PATH
export MANPATH=${MPI_HOME}/share/man:$MANPATH
sudo ldconfig
```

## Step 3. Download the code from `HEMat` and build

### Instructions

1. Clone the code
```
git clone http://github.com/XiaopengZheng/HEMat.git
```

2. Build the code and run (HElib has been installed in `/home/usr/helib_install`)
```
cd HEMat
sudo mkdir build
cd build
sudo cmake -DPACKAGE_BUILD=ON -DCMAKE_INSTALL_PREFIX=/home/usr/helib_install ..
sudo make
cd bin
./BGV_matrix_new
```



