from distutils.core import setup
from distutils.extension import Extension
import sys
import re

if sys.version_info[:2] < (3, 5):
    print('PyDTK2 requires Python 3.5 or higher.')
    sys.exit(-1)

vfile = open('pydtk2/_version.py', mode='r')
vstr_raw = vfile.read()
vstr_find = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", vstr_raw, re.M)
if vstr_find:
    version = vstr_find.group(1)
else:
    raise RuntimeError('Unable to find __version__ in pydtk2/_version.py.')
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
]

_lib_dirs = [
    '/usr/local/lib',
    '/usr/lib/x86_64-linux-gnu',
]

_inc_dirs = [
    '.',
]

ext = Extension(
    'pydtk2.dtk2mapper',
    ['pydtk2/dtk2mapper.cpp'],
    include_dirs=_inc_dirs,
    library_dirs=_lib_dirs,
    libraries=_libs,
    runtime_library_dirs=_lib_dirs,
    extra_compile_args=['-w', '-std=c++1z', '-march=native', '-O3'],
    language='c++')

install_req = [
    'pymoab >= 5.0.1',
]

classifiers = [
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: C++',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development',
    'Operating System :: Linux',
    'Intended Audience :: Science/Research',
]

setup(
    name='pydtk2',
    version=version,
    description='Python Wrapper for Serial DTK2+PyMOAB',
    author='Qiao Chen',
    author_email='benechiao@gmail.com',
    keywords='Math',
    packages=['pydtk2'],
    ext_modules=[ext],
    install_requires=install_req,
    classifiers=classifiers
)
