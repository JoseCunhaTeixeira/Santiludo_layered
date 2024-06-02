import subprocess
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize



if not os.path.exists("./lib/"):
   os.makedirs("./lib/")



extensions = [
    Extension("lib.VGfunctions", 
                ["src/VGfunctions.pyx"],
                language="c++",
                extra_compile_args=["-O3"],
                extra_link_args=["-O3"]),

    Extension("lib.RPfunctions", 
                ["src/RPfunctions.pyx", "src/VGfunctions_src.cpp"],
                language="c++",
                extra_compile_args=["-O3"],
                extra_link_args=["-O3"]),

    Extension("lib.TTDSPfunctions", 
                ["src/TTDSPfunctions.pyx"],
                language="c++",
                extra_compile_args=["-O3"],
                extra_link_args=["-O3"]),
]

setup(
    ext_modules = cythonize(extensions)
)



subprocess.run(['mv', './src/RPfunctions.cpp', './lib/'])
subprocess.run(['mv', './src/VGfunctions.cpp', './lib/'])
subprocess.run(['mv', './src/TTDSPfunctions.cpp', './lib/'])
