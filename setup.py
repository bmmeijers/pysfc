import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

setup(
    name = 'pysfc',
    package_dir = {'': 'src'},
    packages = ['pysfc',],
    ext_modules = cythonize(["src/pysfc/speedups/hilbert.pyx", 
                             "src/pysfc/speedups/query_hilbert.pyx", 
                             "src/pysfc/speedups/relate.pyx",]),
)
