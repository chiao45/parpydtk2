"""
This demo shows how to transfer values between serial and parallel meshes.
The domain is a flat plane in 2D.
"""


import numpy as np
from mpi4py import MPI
from parpydtk2 import *


comm = MPI.COMM_WORLD


# NOTE implementation should be within a try ... except
try:
    # create our blue and green mesh databases
    blue, green = create_imeshdb_pair(comm)
    assert comm.size == 2
    rank = comm.rank

    ################
    # Preprocessing
    ################

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

    ###########################
    # Build blue mesh database
    ###########################

    # begin to create blue
    # NOTE vertices, fields, and global IDs must be construct before finish
    # creating the mesh database
    blue.begin_create()
    # local vertices and global IDs
    lcob = cob[8 * rank:8 * (rank + 1), :].copy()
    lbgids = bgids[8 * rank:8 * (rank + 1)].copy()
    blue.create_vertices(lcob)
    blue.assign_gids(lbgids)
    # do not use trivial global ID strategy
    blue.finish_create(False)
    # NOTE fields must be created after the mesh has been settled
    blue.create_field('b')

    ###########################
    # Build green mesh database
    ###########################

    # begin to create green mesh
    # NOTE that green is assumed to be serial mesh
    green.begin_create()
    # only create on master rank
    if not rank:
        green.create_vertices(cog)
    # since green is serial, we just use the trivial global IDs
    green.finish_create()  # empty partition is resolve here
    # NOTE fields must be created after the mesh has been settled
    green.create_field('g')

    assert green.has_empty()

    ###################
    # setup the problem
    ###################

    # create our analytical model, i.e. 10+sin(x)*cos(y)
    bf = 10.0 + np.sin(lcob[:, 0]) * np.cos(lcob[:, 1])
    gf = np.sin(cog[:, 0]) * np.cos(cog[:, 1]) + 10.0

    ##############
    # Mapper part
    ##############

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

    # NOTE that the following only runs on green master mesh
    if not rank:
        green.assign_field('g', gf)
    # Since green is serial and has empty partition, we must call this to
    # resolve asynchronous values
    green.resolve_empty_partitions('g')

    # solution transfer region
    mapper.begin_transfer()
    mapper.transfer_data(bf='b', gf='g', direct=G2B)
    err_b = (bf - blue.extract_field('b'))/10
    mapper.transfer_data(bf='b', gf='g', direct=B2G)
    err_g = (gf - green.extract_field('g'))/10
    mapper.end_transfer()

    comm.barrier()

    print(rank, 'blue L2-error=%.3e' %
          (np.linalg.norm(err_b)/np.sqrt(err_b.size)))
    if rank == 0:
        print(0, 'green L2-error=%.3e' %
              (np.linalg.norm(err_g)/np.sqrt(err_g.size)))
except Exception:
    # if something goes wrong, set the error code and raise again
    error.ERROR_CODE = 1
    raise
