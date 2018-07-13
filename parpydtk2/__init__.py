import sys
from ._version import __version__


def _excepthook(exctype, value, traceback):
    # override default excepthook
    sys.__excepthook__(exctype, value, traceback)
    import mpi4py
    mpi4py.rc.initialize = False
    mpi4py.rc.finalize = False
    from mpi4py import MPI
    if MPI.Is_initialized() and MPI.COMM_WORLD.size > 1:
        MPI.COMM_WORLD.Abort(1)


sys.excepthook = _excepthook

B2G = True
G2B = False

MMLS = 0
SPLINE = 1
N2N = 2

WENDLAND2 = 0
WENDLAND4 = 1
WENDLAND6 = 2
WU2 = 3
WU4 = 4
WU6 = 5
BUHMANN3 = 6
