.. include:: links.txt

.. _intro:

Introduction
============

Background & Motivation
-----------------------

Multi-physics coupling has become one of the popular topics. People try
to put different models together and simulate the behaviors in a coupled
fashion. With the popularity of `partitioned approach`, i.e. solve different
domains with different solvers and couple the interface conditions as those
solvers boundary conditions, a robust and accurate interface solution
remapping operator is needed. The solution transfer problem on its own is not
an easy task, since it involves the following research aspects:

1.

    **numerical method**, i.e. consistent, conservative, high-order
    convergence.

2.

    **geometry and data structure**, i.e. efficient and robust treatments of
    `mesh association` of two (potentially more) general surfaces that come
    from different discretization methods (`FEM`, `FVM`, `FDM`, etc.) thus
    having different resolutions.

3.

    **parallel rendezvous & HPC**, i.e. handling migrating meshes that have
    different parallel partitions.

`Data Transfer Kit-2.0`_ (DTK2) is a package that is developed at the
`Oak Ridge National Laboratory <https://www.ornl.gov/>`_. `DTK2`_ provides
parallel solution transfer services with `meshless` (a.k.a. `mesh-free`)
methods, which are relatively easy to implement and computational efficient.
Particularly, we are interested in its `modified moving least square` [#]_
method that is an improvement of traditional MLS fitting in terms of robustness
on featured geometries.

.. [#]

    Slattery, S. Hamilton, T. Evans, “A Modified Moving Least Square Algorithm
    for Solution Transfer on a Spacer Grid Surface”, ANS MC2015 - Joint
    International Conference on Mathematics and Computation (M&C),
    Supercomputing in Nuclear Applications (SNA) and the Monte Carlo (MC)
    Method, Nashville, Tennessee · April 19–23, 2015, on CD-ROM, American
    Nuclear Society, LaGrange Park, IL (2015).

`Mesh-Oriented datABase`_ (MOAB) is an array-based general purpose mesh library
with MPI support. Array-based mesh data structure is more efficient in both
computational cost and memory usage compared to traditional pointer-based
data structures. `MOAB`_ has been adapted in `DTK2`_, so we choose to use it
as our mesh database for this work.

In multi-physics coupling, a flexible software framework is must. The fact is:
the physics solvers may be implemented in different programming languages or
shipped as executable binaries (typically commercial codes), this makes using
static languages difficult. Python, on the other hand, can easily glue
different languages together and drive executable binaries smoothly. Its
built-in reference counting, garbage collection, and pass-by-reference make it
as one of the best choices for developing multi-physics coupling frameworks.
In addition, MPI is well supported through the `mpi4py`_ package.
**This motivates us to develop a Python interface for DTK2!**

License
-------

This package is distributed under MIT License. For detailed information, please
take a look at the `LICENSE <https://github.com/chiao45/parpydtk2/blob/parallel/LICENSE>`_ file.

.. _me:

About Me
--------

I am a Ph.D. candidate who work with `Dr. Jim Jiao <http://www.ams.sunysb.edu/~jiao/>`_
on high-order numerical methods. This work is for testing our software
framework of multi-physics coupling in general, `conjugate heat transfer`
(CHT) in particular.

.. note:: Please be aware that I may not have time to maintain this package.
