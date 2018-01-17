# Serial DTK2 with MOAB

## Introduction

This is a python interface of DataTransferKit-2.0 with MOAB mesh database for multiphysics spacial couplings. Check [here](https://ornl-cees.github.io/DataTransferKit/) for for more information about DataTransferKit-2.0, [here](http://sigma.mcs.anl.gov/moab-library/) for MOAB.

## Limitations

Only serial solution transfers are supported. Also, we only wrap moving least square method with radial basis functions method.

## Installation and Requirements

This software is developed and tested with (and only with) this [docker image](https://github.com/unifem/coupler-desktop) environment.

Currently, this software has the following Requirements for installation:

* A complete DataTransferKit-2.0 installation

* A MOAB installation with MPI enabled

* PyMOAB 5.0.1 or higher

* Python3.5 or higher

* MPI and C++1z

Here are optional requirements:

* Cython if you want to regenerate the cython-generated C++ source file.

* mpi4py, see known issues (below)

Installation processes (make system):

* ```make```, build in-place

* ```[sudo] make install```, system-wise installation

* ```make clean_cython```, clean the cython-generated C++ source file

* ```make regen_cython```, regenerate C++ source file with Cython

## Jupyter Notebooks

This software comes with several jupyter notebooks as tutorials. Inside ```./notebooks``` directory, copy the directory into your prefered location and enjoy. For instance, ```cp -r notebooks $HOME/```.

## Known Issue(s)

**NOTE, openmpi may have issues with this software, in order to work correctly, you may need to explicit do ```from mpi4py import MPI``` in your python scripts.**

## License

...

## Contact

Qiao Chen, benechiao@gmail.com
