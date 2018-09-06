"""The error handler module

.. moduleauthor:: Qiao Chen <benechiao@gmail.com>

Attributes
----------
ERROR_CODE : int
    set this to nonzero values one exceptions have been raised

Examples
--------
>>> import parpydtk2 as dtk
>>> try:
...     # your programs here
... except Expection:
...     dtk.error.ERROR_CODE = 1
...     raise
"""

import atexit


ERROR_CODE = 0


def _atexit():
    if ERROR_CODE != 0:
        import traceback
        import sys
        from mpi4py import MPI

        traceback.print_exc()
        sys.stderr.flush()
        if MPI.Is_initialized() and MPI.COMM_WORLD.size > 1:
            MPI.COMM_WORLD.Abort(ERROR_CODE)


atexit.register(_atexit)
