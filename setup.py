from distutils.command.clean import clean
import os
import sys
import re
import codecs
import numpy
from setuptools import setup, Extension, find_packages
try:
    import mpi4py
except ImportError:
    sys.stderr.write('You must install mpi4py first!\n')
    sys.stderr.flush()
    sys.exit(-1)

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

classifiers = [
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: C++',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development',
    'Operating System :: POSIX :: Linux',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
]

install_requires = [
    'numpy',
    'setuptools',
]

package_data = {
    'parpydtk2': [
        '*.pxd',
        '*.pyx',
        '*.cpp',
        '*.hpp',
        os.path.join('src', '*.hpp'),
    ]
}

# set the compiler to mpi
os.environ['CC'] = mpi4py.get_config()['mpicxx']

LEN1, LEN2 = 80, 5


def ensure_cpp11():
    import tempfile
    from distutils.ccompiler import get_default_compiler, new_compiler
    import distutils
    cmp = new_compiler(get_default_compiler())
    tmp_source = os.path.join(tempfile.gettempdir(), 'foo.cpp')
    f = open(tmp_source, 'w')
    f.write('int main(){return 0;}')
    f.close()

    def safe_remove(filename):
        try:
            os.remove(filename)
        except OSError:
            pass
    try:
        cmp.compile([tmp_source], extra_postargs=['-std=c++11'])
        safe_remove(tmp_source)
        safe_remove(tempfile.gettempprefix())
    except distutils.errors.CompileError:
        safe_remove(tmp_source)
        safe_remove(tempfile.gettempprefix())
        sys.stderr.write('You must use a compiler that supports C++11\n')
        sys.stderr.flush()
        sys.exit(-1)


ensure_cpp11()


# configure paths
include_dirs = [
    '.',
    numpy.get_include(),
    mpi4py.get_include()
]

libraries = [
    'MOAB',
    'dtk_moabadapters',
    'dtk_interface',
    'teuchoscomm',
    'teuchosparameterlist',
    'teuchoscore',
    'dtk_operators'
]

library_dirs = None


def configure_lib(env):
    # check env
    root = os.environ.get(env, None)
    if root is not None:
        return ([os.path.join(root, 'include')],
                [os.path.join(root, 'lib')])
    return None


MOAB = configure_lib('PARPYDTK2_MOAB_ROOT')
DTK = configure_lib('PARPYDTK2_DTK_ROOT')


if MOAB is not None:
    if library_dirs is None:
        library_dirs = []
    include_dirs += MOAB[0]
    library_dirs += MOAB[1]
if DTK is not None:
    if library_dirs is None:
        library_dirs = []
    include_dirs += DTK[0]
    library_dirs += DTK[1]


if MOAB is None or DTK is None:
    def configure():
        # check arguments
        prefix = None
        if '--user' in sys.argv:
            from site import USER_BASE
            prefix = USER_BASE
        else:
            for arg in sys.argv[1:]:
                if arg.startswith('--home'):
                    prefix = arg.split('=')[1]
                    if prefix == '~':
                        prefix = os.environ.get('HOME')
                elif arg.startswith('--prefix'):
                    prefix = arg.split('=')[1]
        if prefix is not None:
            return ([os.path.join(prefix, 'include')],
                    [os.path.join(prefix, 'lib')])
        return None
    CONF = configure()
    if CONF is not None:
        if library_dirs is None:
            library_dirs = []
        include_dirs += CONF[0]
        library_dirs += CONF[1]


rpath = library_dirs

exts = [
    Extension(
        'parpydtk2.imeshdb',
        ['parpydtk2/imeshdb.cpp'],
        include_dirs=include_dirs,
        language='c++',
        libraries=libraries,
        library_dirs=library_dirs,
        runtime_library_dirs=rpath,
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        'parpydtk2.mapper',
        ['parpydtk2/mapper.cpp'],
        include_dirs=include_dirs,
        language='c++',
        libraries=libraries,
        library_dirs=library_dirs,
        runtime_library_dirs=rpath,
        define_macros=[('LEN1', LEN1), ('LEN2', LEN2)],
        extra_compile_args=['-std=c++11'],
    )
]


class MyClean(clean):
    def run(self):
        import subprocess as sp
        super().run()
        sp.Popen(
            'rm -rf parpydtk2/__pycache__ __pycache__ dist parpydtk2.egg-info',
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
    long_description=codecs.open('README.rst', encoding='utf-8').read(),
    url='https://github.com/chiao45/parpydtk2',
    packages=find_packages(),
    license='MIT',
    ext_modules=exts,
    classifiers=classifiers,
    cmdclass={'clean': MyClean},
    install_requires=install_requires,
    package_data=package_data,
)
