.. include:: links.txt

.. _install:

Installation
============

Installing this package is not a trivial task due to its heavy dependencies.
ParPyDTK2 has the following installation requirements:

.. _dep:

1. C++11 compiler
2. MPI
3. `MOAB`_
4. `DTK2`_ and `Trilinos`_
5. Python >= 3.5
6. `mpi4py`_
7. `NumPy`_
8. setuptools

And both `MOAB`_ and `DTK2`_ have their own dependencies.

In addition, to build the documentation, the following packages are needed:

1. Sphinx
2. Doxygen
3. breathe
4. numpydoc

The good news is you can install these easily through ``pip``.

.. _install_moab:

Install `MOAB`_
---------------

The `MOAB official README <https://bitbucket.org/fathomteam/moab>`_ has a
very clear description of the installation process. Here we take an excerpt
of our `MOAB Docker image building script <https://github.com/unifem/cht-coupler/blob/meshdb-bin/Dockerfile#L13>`_:

.. code-block:: console

    $ git clone --depth=1 https://bitbucket.org/fathomteam/moab.git
    $ cd moab
    $ autoreconf -fi
    $ ./configure \
        --prefix=/usr/local \
        --with-mpi \
        CC=mpicc \
        CXX=mpicxx \
        FC=mpif90 \
        F77=mpif77 \
        --enable-optimize \
        --enable-shared=yes \
        --with-blas=-lopenblas \
        --with-lapack=-lopenblas \
        --with-scotch=/usr/lib \
        --with-metis=/usr/lib/x86_64-linux-gnu \
        --with-eigen3=/usr/include/eigen3 \
        --with-x \
        --with-cgns \
        --with-netcdf \
        --with-hdf5=/usr/lib/hdf5-openmpi \
        --with-hdf5-ldflags="-L/usr/lib/hdf5-openmpi/lib" \
        --enable-ahf=yes \
        --enable-tools=yes
    $ make && sudo make install

Notice that this is for system installation. Install to your preferred
locations if you don't have the root access. Also, turn off those optional
packages if you don't have them, only MPI and HDF5 are necessary.

.. warning:: You must build it into a shared object!

.. note::

    If you use Ubuntu >= 17.10, all those optional packages are likely to be
    available through ``apt``.

.. _install_dtk2:

Install `DTK2`_
---------------

`DTK2`_ is shipped as a sub-module of `Trilinos`_, so building `Trilinos`_ is
needed. For people who are not familiar with `Trilinos`_, this can be tricky.
Therefore, an excerpt of our `DTK2 Docker image building script <https://github.com/unifem/cht-coupler/blob/mapper-bin/Dockerfile#L42>`_
might be helpful:

.. code-block:: console

    $ export TRILINOS_VERSION=12-12-1
    $ git clone --depth 1 --branch trilinos-release-${TRILINOS_VERSION}
    $ cd Trilinos
    $ git clone --depth 1 --branch dtk-2.0 \
        https://github.com/unifem/DataTransferKit.git
    $ mkdir build && cd build
    $ cmake \
        -DCMAKE_INSTALL_PREFIX:PATH=/usr/local \
        -DCMAKE_BUILD_TYPE:STRING=RELEASE \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
        -DCMAKE_SHARED_LIBS:BOOL=ON \
        -DTPL_ENABLE_MPI:BOOL=ON \
        -DTPL_ENABLE_Boost:BOOL=ON \
        -DBoost_INCLUDE_DIRS:PATH=/usr/include/boost \
        -DTPL_ENABLE_Libmesh:BOOL=OFF \
        -DTPL_ENABLE_MOAB:BOOL=ON \
        -DMOAB_INCLUDE_DIRS=$MOAB_ROOT/include \
        -DMOAB_LIBRARY_DIRS=$MOAB_ROOT/lib \
        -DTPL_ENABLE_Netcdf:BOOL=ON \
        -DTPL_ENABLE_BinUtils:BOOL=OFF \
        -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
        -DTrilinos_ENABLE_ALL_PACKAGES=OFF \
        -DTrilinos_EXTRA_REPOSITORIES="DataTransferKit" \
        -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
        -DTrilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
        -DTrilinos_ENABLE_TESTS:BOOL=OFF \
        -DTrilinos_ENABLE_EXAMPLES:BOOL=OFF \
        -DTrilinos_ENABLE_CXX11:BOOL=ON \
        -DTrilinos_ENABLE_Tpetra:BOOL=ON \
        -DTpetra_INST_INT_UNSIGNED_LONG:BOOL=ON \
        -DTPL_ENABLE_BLAS:BOOL=ON \
        -DTPL_BLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/libopenblas.so \
        -DTPL_ENABLE_LAPACK:BOOL=ON \
        -DTPL_LAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/libopenblas.so \
        -DTPL_ENABLE_Eigen:BOOL=ON \
        -DTPL_Eigen_INCLUDE_DIRS=/usr/include/eigen3 \
        -DTrilinos_ENABLE_DataTransferKit=ON \
        -DDataTransferKit_ENABLE_DBC=ON \
        -DDataTransferKit_ENABLE_TESTS=ON \
        -DDataTransferKit_ENABLE_EXAMPLES=OFF \
        -DDataTransferKit_ENABLE_ClangFormat=OFF \
        -DTPL_ENABLE_BoostLib:BOOL=OFF \
        -DBUILD_SHARED_LIBS:BOOL=ON
    $ make && sudo make install

Again, this assumes root access, adjust this based on your situation. `DTK2`_
needs to stay in the root directory of `Trilinos`_ and be turned on through
switches ``DTrilinos_EXTRA_REPOSITORIES`` and
``DTrilinos_ENABLE_DataTransferKit``. The environment var ``MOAB_ROOT`` is
the place where you install `MOAB`_.

.. note::

    We recommend that install `DTK2`_ from `my personal forked repo <https://github.com/chiao45/DataTransferKit>`_
    since we may add/modify the source codes to make `DTK2`_ more advanced.

.. _install_this:

Install ParPyDTK2
-----------------

Once you have the dependencies setup, installing ParPyDTK2 can be very easy.
The easiest way is through PyPI:

.. code-block:: console

    $ sudo pip3 install parpydtk2

However, this assumes that ParPyDTK2 can find `MOAB`_ and `DTK2`_ on the
system. With different specifications of ``install`` command, ParPyDTK2 can
automatically add different paths to search for `MOAB`_ and `DTK2`_.

.. code-block:: console

    $ pip3 install parpydtk2 --user

will assume `MOAB`_ and `DTK2`_ can be found in ``USER_BASE/{include,lib}``.

.. code-block:: console

    $ pip3 install parpydtk2 --prefix=...

will allow ParPyDTK2 to search `MOAB`_ and `DTK2`_ under the ``prefix``
directory.

The preferred way is to define the environment variables
``PARPYDTK2_MOAB_ROOT`` and ``PARPYDTK2_DTK_ROOT`` before you do
:code:`pip install`. For instance,

.. code-block:: console

    $ export PARPYDTK2_MOAB_ROOT=/path/to/moab/root
    $ export PARPYDTK2_DTK_ROOT=/path/to/dtk/root
    $ pip3 install parpydtk2 --user

.. warning::

    We don't mark `mpi4py`_ as installation dependency, so you need to install
    it manually before you install ParPyDTK2. :code:`pip3 install mpi4py` is
    just fine.

Of course, you can install from source, which can be obtained `here`_. Just
make sure you have all Python :ref:`dependencies <dep>` installed.

.. code-block:: console

    $ git clone -b parallel https://github.com/chiao45/parpydtk2.git
    $ cd parpydtk2
    $ python3 setup.py install --user

.. _use_docker:

Using our Docker container
--------------------------

You can try the package through our pre-built `Docker container <https://hub.docker.com/r/chiao/dtk/>`_.
Two driver scripts are provided in order to easily use the container:

.. only:: html

    1. :download:`parpydtk2_desktop.py<../docker/parpydtk2_desktop.py>`
    2. :download:`parpydtk2_jupyter.py<../docker/parpydtk2_jupyter.py>`

.. only:: latex

    1. `parpydtk2_desktop.py <https://github.com/chiao45/parpydtk2/blob/parallel/docker/parpydtk2_desktop.py>`_
    2. `parpydtk2_jupyter.py <https://github.com/chiao45/parpydtk2/blob/parallel/docker/parpydtk2_jupyter.py>`_

The former will launch a desktop environment through VNC, while the latter will
run the container as a Jupyter server.
