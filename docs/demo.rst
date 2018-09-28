.. include:: links.txt

.. _demo:

A Demo
======

Here, we show a demo for transferring solutions between two meshes of the
unit square. We let the ``blue`` mesh participant run in parallel with two
cores, while the ``green`` side is treated as a serial mesh.

.. code-block:: python
    :linenos:

    import numpy as np
    from mpi4py import MPI
    from parpydtk2 import *

    comm = MPI.COMM_WORLD

    blue, green = create_imeshdb_pair(comm)

    assert comm.size == 2
    rank = comm.rank

For demo purpose, we construct the meshes globally.

.. code-block:: python
    :lineno-start: 12

    # create blue meshes on all processes
    cob = np.empty(shape=(16, 3), dtype=float, order='C')
    dh = 0.33333333333333333
    index = 0
    x = 0.0
    for i in range(4):
        y = 0.0
        for j in range(4):
            cob[index] = np.asarray([x, y, 0.0])
            index += 1
            y += dh
        x += dh

    # global IDs, one based
    bgids = np.arange(16, dtype='int32')
    bgids += 1

The ``blue`` side has 16 nodes, the following is for the ``green`` side:

.. code-block:: python
    :lineno-start: 29

    # create green meshes on all processes
    cog = np.empty(shape=(36, 3), dtype=float, order='C')
    dh = 0.2
    index = 0
    x = 0.0
    for i in range(6):
        y = 0.0
        for j in range(6):
            cog[index] = np.asarray([x, y, 0.0])
            index += 1
            y += dh
        x += dh

    # global IDs
    ggids = np.arange(36, dtype='int32')
    ggids += 1

The ``green`` participant has 36 nodes. The next step is to put the data in
the two mesh databases.

.. code-block:: python
    :lineno-start: 46

    # creating the mesh database
    blue.begin_create()
    # local vertices and global IDs
    lcob = cob[8 * rank:8 * (rank + 1), :].copy()
    lbgids = bgids[8 * rank:8 * (rank + 1)].copy()
    blue.create_vertices(lcob)
    blue.assign_gids(lbgids)
    # do not use trivial global ID strategy
    blue.finish_create(False)
    blue.create_field('b')

As we can see, we equally distributed the mesh into the two cores as well as
the corresponding global IDs. In line 54, the ``False`` flag indicates that
the mesh database should use the user-provided global IDs.

.. warning::

    Creating vertices and assigning global IDs must be called between
    :py:func:`~parpydtk2.IMeshDB.begin_create` and
    :py:func:`~parpydtk2.IMeshDB.finish_create`! Otherwise, exceptions are
    thrown.

.. warning::

    Creating fields must be done after
    :py:func:`~parpydtk2.IMeshDB.finish_create`!

Here is the treatment for the "serial" participant:

.. code-block:: python
    :lineno-start: 58

    # NOTE that green is assumed to be serial mesh
    green.begin_create()
    # only create on master rank
    if not rank:
        green.create_vertices(cog)
    # since green is serial, we just use the trivial global IDs
    green.finish_create()  # empty partition is resolve here
    green.create_field('g')  # must after finish create

    assert green.has_empty()

As we can see, only the master process has data.

.. note::

    The :ref:`duplicated <empty_part>` node is handled inside
    :py:func:`~parpydtk2.IMeshDB.finish_create`

With the two participants ready, we can now create our
:py:class:`~parpydtk2.Mapper`.

.. code-block:: python
    :lineno-start: 69

    # create our analytical model, i.e. 10+sin(x)*cos(y)
    bf = 10.0 + np.sin(lcob[:, 0]) * np.cos(lcob[:, 1])
    gf = np.sin(cog[:, 0]) * np.cos(cog[:, 1]) + 10.0

    # Construct our mapper
    mapper = Mapper(blue=blue, green=green)

    # using mmls, blue radius 1.0 green radius 0.6
    mapper.method = MMLS
    mapper.radius_b = 1.0
    mapper.radius_g = 0.6
    mapper.dimension = 2

    # Mapper initialization region
    mapper.begin_initialization()
    mapper.register_coupling_fields(bf='b', gf='g', direct=B2G)
    mapper.register_coupling_fields(bf='b', gf='g', direct=G2B)
    mapper.end_initialization()

Line 69-70 just create an analytic model for error analysis. Line 77-80 are
for parameters, for this case, we use MMLS (default) with ``blue`` side
radius 1.0 and ``green`` side radius 0.6 for searching.

The important part is from line 83 to 86. Particularly speaking, the function
:py:func:`~parpydtk2.Mapper.register_coupling_fields`. It takes three
parameters, where the first two are string tokens that represents the data
fields in ``blue`` and ``green``. The ``direct`` is to indicate the transfer
direction, e.g. :py:attr:`~parpydtk2.B2G` stands for blue to green.

.. warning::

    :py:func:`~parpydtk2.Mapper.register_coupling_fields` must be called within
    :py:func:`~parpydtk2.Mapper.begin_initialization` and
    :py:func:`~parpydtk2.Mapper.end_initialization`.

.. code-block:: python
    :lineno-start: 88

    # NOTE that the following only runs on green master mesh
    if not rank:
        green.assign_field('g', gf)
    # Since green is serial and has empty partition, we must call this to
    # resolve asynchronous values
    green.resolve_empty_partitions('g')

The above section is to assign values on the ``green`` participant. Notice that
it is a "serial" mesh, so we only assign values on the master process. But
resolve the :ref:`duplicated <empty_part>` node is needed, this is done in
line 93.

Finally, the solution transfer part is pretty straightforward:

.. code-block:: python
    :lineno-start: 94

    # solution transfer region
    mapper.begin_transfer()
    mapper.transfer_data(bf='b', gf='g', direct=G2B)
    err_b = (bf - blue.extract_field('b'))/10
    mapper.transfer_data(bf='b', gf='g', direct=B2G)
    err_g = (gf - green.extract_field('g'))/10
    mapper.end_transfer()

    comm.barrier()

    print(rank, 'blue L2-error=%.3e' % (np.linalg.norm(err_b)/np.sqrt(err_b.size)))
    if rank == 0:
        print(0, 'green L2-error=%.3e' % (np.linalg.norm(err_g)/np.sqrt(err_g.size)))

.. only:: html

    This code can be obtained :download:`here parallel2serial.py<../examples/parallel2serial.py>`.

.. only:: latex

    This code can be obtained `here parallel2serial.py <https://github.com/chiao45/parpydtk2/blob/parallel/examples/parallel2serial.py>`_.
