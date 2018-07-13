from distutils.core import setup
from distutils.extension import Extension
import sys
import re
import numpy
import mpi4py

if sys.version_info[:2] < (3, 5):
    print('ParPyDTK2 requires Python 3.5 or higher.')
    sys.exit(-1)

vfile = open('parpydtk2/_version.py', mode='r')
vstr_raw = vfile.read()
vstr_find = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", vstr_raw, re.M)
if vstr_find:
    version = vstr_find.group(1)
else:
    raise RuntimeError('Unable to find __version__ in parpydtk2/_version.py.')
vfile.close()


_libs = [
    'MOAB',
    'dtk_moabadapters',
    'dtk_interface',
    'teuchoscomm',
    'teuchosparameterlist',
    'teuchoscore',
    'dtk_operators',
    'stdc++',
    'mpi',
]

_lib_dirs = [
    '/usr/local/lib',
    '/usr/lib/x86_64-linux-gnu',
]

_inc_dirs = [
    '.',
    '/usr/local/include',
    numpy.get_include(),
    mpi4py.get_include()
]

ext = Extension(
    'parpydtk2.mapper',
    ['parpydtk2/mapper.cpp'],
    include_dirs=_inc_dirs,
    library_dirs=_lib_dirs,
    libraries=_libs,
    runtime_library_dirs=_lib_dirs,
    extra_compile_args=['-w', '-std=c++1z', '-march=native', '-O3'],
    language='c++'
)

classifiers = [
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: C++',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development',
    'Operating System :: Linux',
    'Intended Audience :: Science/Research',
]


def config_mpi():
    import subprocess
    prog = subprocess.Popen(['mpicxx', '-showme'], stdout=subprocess.PIPE)
    while prog.poll() is None:
        continue
    assert prog.poll() == 0
    output = prog.stdout.readline().decode('UTF-8')  # already in py3
    outputs = output.split(' ')
    inc = []
    lib = []
    for arg in outputs:
        if arg.startswith('-I'):
            inc.append(arg[2:])
        elif arg.startswith('-L'):
            lib.append(arg[2:])
    return (inc, lib)


_mpi_conf = config_mpi()
_inc_dirs += _mpi_conf[0]
_lib_dirs += _mpi_conf[1]


setup(
    name='parpydtk2',
    version=version,
    description='Python Wrapper for Parallel DTK2+MOAB',
    author='Qiao Chen',
    author_email='benechiao@gmail.com',
    keywords='Math',
    packages=['pydtk2'],
    ext_modules=[ext],
    classifiers=classifiers
)
