# Parallel DTK2 with MOAB for Python

## Introduction

This is a python interface of DataTransferKit-2.0 with MOAB mesh database for multiphysics spacial coupling. Check [here](https://ornl-cees.github.io/DataTransferKit/) for for more information about DataTransferKit-2.0, [here](http://sigma.mcs.anl.gov/moab-library/) for MOAB.

## Installation and Requirements

Currently, this software has the following requirements for installation:

* A complete DataTransferKit-2.0 installation
* A MOAB installation with MPI enabled
* Python3.5 or higher
* MPI and C++1z
* mpi4py
* Cython (optional) if you want to regenerate the C++ source file.

Installation processes: ```[sudo] env CC=mpicxx python3 setup.py install```

## License

MIT License, Copyright (c) 2018 Qiao Chen

## Contact

Qiao Chen, benechiao@gmail.com
