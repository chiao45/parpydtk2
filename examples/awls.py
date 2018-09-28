from mpi4py import MPI
import numpy as np
from parpydtk2 import create_imeshdb_pair, Mapper, B2G


if MPI.COMM_WORLD.size > 1:
    import sys
    sys.stderr.write('This program only runs in serial\n')
    sys.stderr.flush()
    sys.exit(1)


# define the model
def exp_func(pts):
    return np.exp(pts[:, 0] + pts[:, 1] + pts[:, 2])


def get_parser():
    import argparse
    parser = argparse.ArgumentParser(description=('AWLS parameters'))
    parser.add_argument(
        '-r',
        '--rho',
        type=float,
        help='local column scaling factor, default is 2.3',
        default=2.3
    )
    parser.add_argument(
        '-s',
        '--silent',
        action='store_true',
        help='no information streaming, default is False',
        default=False
    )
    parser.add_argument(
        '-p',
        '--profile',
        type=str,
        help='profiling file name, default is \'\', i.e. no profiling',
        default=''
    )
    return parser


def main():

    args = get_parser().parse_args()
    verbose = not args.silent
    profiling = bool(args.profile)
    stat_file = args.profile + '%i'

    # load the uniformly refined grids of torus
    point_clouds = np.load('torus.npz')

    l2_errors = []

    for level in range(3):
        source_pts = point_clouds['fine{}'.format(level+1)]
        target_pts = point_clouds['coarse{}'.format(level+1)]

        # create two meshdb
        src, tgt = create_imeshdb_pair()

        # initialize and create field
        src.begin_create()
        src.create_vertices(source_pts)
        src.finish_create()

        src.create_field('src')

        # assign values to source
        src.assign_field('src', exp_func(source_pts))

        # initialize and create field for target
        tgt.begin_create()
        tgt.create_vertices(target_pts)
        tgt.finish_create()

        tgt.create_field('tgt')

        # compute the exact function values for target
        tgt_exact = exp_func(target_pts)

        # create the mapper
        mapper = Mapper(blue=src, green=tgt, profiling=profiling,
                        verbose=verbose, stat_file=stat_file % level)

        mapper.dimension = 3  # we know the spatial dimension is 3
        mapper.awls_conf(dim=2, rho=args.rho)

        mapper.begin_initialization()
        mapper.register_coupling_fields('src', 'tgt', direct=B2G)
        mapper.end_initialization()

        mapper.begin_transfer()
        mapper.transfer_data('src', 'tgt', direct=B2G)
        mapper.end_transfer()

        err = tgt.extract_field('tgt')-tgt_exact
        l2_errors.append(np.linalg.norm(err)/np.linalg.norm(tgt_exact))

    print('l2-errors', l2_errors)
    print('convergence rate: %.3e' % abs(np.log2(l2_errors[-1]/l2_errors[0])/2))


if __name__ == '__main__':
    main()
