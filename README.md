# CT Mathematics/Algebraic Library

Mathematic/Algebraic algorithms for supporting CT data manipulation.

## Cloning repository:

Basic clone can be done via
```
git clone https://gitlab.stimulate.ovgu.de/vojtech.kulvait/CTIOL.git
```
Or over ssh
```
git clone ssh://git@gitlab.stimulate.ovgu.de:2200/vojtech.kulvait/CTIOL.git
```
However to populate submodules directories it is then needed to issue
```
git submodule init
git submodule update
```

## Submodules
The project contains submodules in the subdirectory submodules. 

### [Plog](https://github.com/SergiusTheBest/plog) logger

Logger Plog is used for logging. It is licensed under the Mozilla Public License Version 2.0.

### [Catch2](https://github.com/catchorg/Catch2)

Testing framework. Licensed under Boost Software License 1.0.

### [CTIOL](ssh://git@gitlab.stimulate.ovgu.de:2200/vojtech.kulvait/CTIOL.git)

Input output routines for asynchronous thread safe reading/writing CT data. The DEN format read/write is implemented.

## How to build:
For build process, make and cmake utilities are required. Install them using
```
apt-get install make cmake
```


Please create separate directory and build there
```
mkdir build
cd build
cmake ..
make
```

## Required libraries
Intel MKL

This library is important for LAPACK routines.

## Using IDEs

#Eclipse
For eclipse to work, create out of source project that is not child or ancesor of the parent directory.
Run the following to generate eclipse project
```
cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug path_to_project
```
After each change in CMakeLists.txt remove eclipse project and regenerate it.
See more at https://gitlab.kitware.com/cmake/community/wikis/doc/editors/Eclipse-CDT4-Generator

