.. include:: links.txt

.. _details:

Some Details
============

.. _mesh_less:

Meshless/Mesh-free
------------------

Point clouds are easier to use (for the end user) compared to meshes. However,
there are some drawbacks: 1) using meshes can achieve linear time complexity
for mesh associations [#]_; 2) numerical integrations are not easy; 3) handling
cell-averaged data fields, e.g. solutions that come from FVM, potential can
lead to problems.

.. [#]

    Jiao X, Edelsbrunner X, Heath MT. Mesh association: formulation and
    algorithms. In 8th International Meshing Roundtable. Sandia Report, South
    Lake Tahoe, CA, 1999; 75–82.

For 3), typically, people assume cell-averaged data to be cell-centered data,
which limits to 2nd-order of accuracy. (Notice that this is actually not a
big problem in practice.)

.. _global_ids:

Global IDs/Handles
------------------

`MOAB`_ uses global IDs for parallel communications as well as `DTK2`_.
Therefore, ParPyDTK2 expects the user to provide this information. For most
applications, this can be computed offline.

The global IDs are unique handles of the vertices in a point cloud. For
instance, if one wants to distribute two triangles, then each of them has
local IDs from 0 to 2, but unique global IDs that range from 1 to 4.

.. note:: ParPyDTK2 uses 1-based indexing for global IDs

.. _empty_part:

Treatment of Empty Partitions
-----------------------------

Both `MOAB`_ and `DTK2`_ don't natively support empty partitions, which occur
pretty frequently in practice. For instance, couple a serial solver with a
parallel solver, or couple a commercial code with an open-sourced one through
socket and the incoming data partition graph from the commercial code probably
doesn't align with that of the open-sourced side.

To support empty partitions, we duplicate the first node in the master process
to the processes that are empty. However, this is not complete yet, because
the user is not aware of this, so that the values between across duplicated
nodes are not synchronized. Of course, we don't want to let the user to
explicitly handle this extra layer of communication. To resolve this, a member
function, called :py:func:`~parpydtk2.IMeshDB.resolve_empty_partitions`, is
added to class :py:class:`~parpydtk2.IMeshDB` and must be called
**collectively** whenever the user updates the field values through
:py:func:`~parpydtk2.IMeshDB.assign_field`.
