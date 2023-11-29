from Cython.Build import cythonize
from setuptools import setup, Extension

setup(  name='pygol',
        version='1.0',
        py_modules=['pygol'],
        ext_modules =
            cythonize(Extension("pygol",
                        # the extension name
                sources=["pygol.c"],
                        # the Cython source and additional C++ source files
                language="c++", # generate and compile C++ code
                                )
                        )
    )