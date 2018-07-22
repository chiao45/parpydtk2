import sys
from ._version import __version__
from .mapper import Mapper


def _excepthook(exctype, value, tb):
    # override default excepthook
    sys.__excepthook__(exctype, value, tb)
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
