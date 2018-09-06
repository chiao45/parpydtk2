"""Main module interface of ParPyDTK2

.. moduleauthor:: Qiao Chen <benechiao@gmail.com>

Attributes
----------
B2G : bool
    boolean flag of ``True`` denotes transferring direction from ``blue`` to
    ``green``
G2B : bool
    boolean flag of ``False`` denotes transferring direction from ``green`` to
    ``blue``
MMLS : int
    flag (0) represents using `modified moving least square` method
SPLINE : int
    flag (1) represents using `spline interpolation` method
N2N : int
    flag (2) represents using `nearest node projection` method
WENDLAND2 : int
    flag (0) represents using `Wendland 2nd-order` RBF weights
WENDLAND4 : int
    flag (1) represents using `Wendland 4th-order` RBF weights
WENDLAND6 : int
    flag (2) represents using `Wendland 6th-order` RBF weights
WU2 : int
    flag (3) represents using `Wu 2nd-order` RBF weights
WU4 : int
    flag (4) represents using `Wu 4th-order` RBF weights
WU6 : int
    flag (5) represents using `Wu 6th-order` RBF weights
BUHMANN3 : int
    flag (6) represents using `Buhmann 3rd-order` RBF weights
"""

import sys
import parpydtk2.error_handle as error
from ._version import __version__
from .imeshdb import IMeshDB
from .mapper import Mapper

__b2g__ = True
__g2b__ = False
B2G = __b2g__
G2B = __g2b__

__mmls__ = 0
__spline__ = 1
__n2n__ = 2
MMLS = __mmls__
SPLINE = __spline__
N2N = __n2n__

__wendland2__ = 0
__wendland4__ = 1
__wendland6__ = 2
WENDLAND2 = __wendland2__
WENDLAND4 = __wendland4__
WENDLAND6 = __wendland6__

__wu2__ = 3
__wu4__ = 4
__wu6__ = 5
WU2 = __wu2__
WU4 = __wu4__
WU6 = __wu6__

__buhmann3__ = 6
BUHMANN3 = __buhmann3__


def get_include():
    """Get the abs include path"""
    import os
    return os.path.dirname(os.path.abspath(__file__))
