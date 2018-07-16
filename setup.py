from distutils.core import setup
from distutils.extension import Extension
from distutils.command.clean import clean
import sys
import re
import numpy
import mpi4py
import os

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

os.environ['CC'] = 'mpicxx'


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


class MyClean(clean):
    def run(self):
        import subprocess as sp
        super().run()
        sp.Popen(
            'rm -rf parpydtk2/__pycache__ __pycache__',
            executable='/bin/bash',
            shell=True
        )
        if self.all:
            sp.Popen(
                'rm -rf parpydtk2/*.so',
                executable='/bin/bash',
                shell=True
            )


setup(
    name='parpydtk2',
    version=version,
    description='Python Wrapper for Parallel DTK2+MOAB',
    author='Qiao Chen',
    author_email='benechiao@gmail.com',
    keywords='Math',
    packages=['parpydtk2'],
    ext_modules=[ext],
    classifiers=classifiers,
    cmdclass={'clean': MyClean},
)
