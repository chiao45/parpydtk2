.. include:: links.txt

.. _details:

Some Details
============

.. _global_ids:

Global IDs/Handles
------------------

`MOAB`_ uses global IDs in parallel as well as `DTK2`_. Therefore, ParPyDTK2
expects the user to provide this information. For most applications, this can
be obtained offline.

The global IDs are unique handles of the vertices in a point cloud. For
instance, if you want to distribute two triangles, then each of them has
local IDs from 0 to 2, but a unique global ID that ranges from 1 to 4.

.. note:: ParPyDTK2 uses 1-based indexing for global IDs

.. _empty_part:

Treatment of Empty Partitions
-----------------------------

Both `MOAB`_ and `DTK2`_ don't support empty partitions, which occur pretty
frequently in practice. For instance, couple a serial solver with a parallel
solver, or couple a commercial code with an open-sourced one through socket
and the incoming data from the commercial code probably doesn't align with
the open-sourced side.

To support empty partitions, we duplicate the first node in the master process
to the processes that are empty. However, this is not complete yet, because
the user is not aware of this, so that the values between the duplicated nodes
are not synchronized. Of course, we don't want to let the user to explicitly
handle this extra layer of communication. To resolve this, a member function,
called :py:func:`~parpydtk2.IMeshDB.resolve_empty_partitions`, is
added to class :py:class:`~parpydtk2.IMeshDB` and must be called
**collectively** whenever the user updates the field values.
