#!/usr/bin/env python

import os
import sys
import numpy as np

from setuptools import setup, Extension

# Get the version number from ModelInterface.h

__version__ = None

with open("cthreeML/ModelInterface.h") as f:
        
    for line in f:
                
        if line.find("#define INTERFACE_VERSION")==0:
            
            __version__ = "%i.0.3" % int(line.split(" ")[2])
            
            print __version__
            
            break

if __version__ is None:
    
    raise RuntimeError("Could not probe version from ModelInterface.h")

# Now a global __version__ is available

# A list to store final messages to print at the end
final_messages = []

# Probe whether the user has specified its own boost directory through the BOOSTROOT
# environment variable

boost_root = os.environ.get("BOOSTROOT")

if boost_root:

    # Check that the directory provided actually exists

    if not os.path.exists(boost_root):
        print("\nERROR: the directory %s specified in BOOSTROOT does not exist!" % boost_root)
        sys.exit(-1)

    # The user want to override pre-defined location of boost

    print("\n\n **** Using boost.python from the env. variable $BOOSTROOT (%s)" % boost_root)

    include_dirs = [os.path.join(boost_root, 'include')]
    library_dirs = [os.path.join(boost_root, 'lib')]

    # Check that the include and library directories exist

    if not os.path.exists(include_dirs[0]):
        print("\nERROR: the include directory %s for boost.python does not exist!" % include_dirs[0])

        sys.exit(-1)

    if not os.path.exists(library_dirs[0]):
        print("\nERROR: the library directory %s for boost.python does not exist!" % library_dirs[0])

        sys.exit(-1)

    final_messages.append("Used boost.python from the env. variable BOOSTROOT")
    final_messages.append("     Include dir: %s" % include_dirs)
    final_messages.append("     Library dir: %s" % library_dirs)

else:

        include_dirs = []
        library_dirs = []

        final_messages.append("Using boost.python from the system path.")

# Now add the numpy headers
include_dirs.append(np.get_include())

# Configure the variables to build the external module with the C/C++ wrapper

ext_modules_configuration = [

    Extension("cthreeML.pyModelInterfaceCache",

              ["cthreeML/pyToCppModelInterfaceCache.cxx",],

              libraries=["boost_python"],

              include_dirs=include_dirs,

              library_dirs=library_dirs,
              extra_compile_args = [])]

headers_configuration = ["cthreeML/ModelInterface.h",
                         "cthreeML/pyToCppModelInterfaceCache.h"]


setup(

    name="cthreeML",

    packages=["cthreeML"],

    version=__version__,

    description="The C/C++ bridge for the Multi-Mission Maximum Likelihood framework (github.com/giacomov/3ML)",

    long_description="The C/C++ bridge for the Multi-Mission Maximum Likelihood framework (github.com/giacomov/3ML)",
    
    license='BSD-3',

    author='Giacomo Vianello',

    author_email='giacomo.vianello@gmail.com',

    url='https://github.com/giacomov/cthreeML',

    download_url='https://github.com/giacomov/cthreeML/archive/%s' % __version__,

    keywords=['Likelihood', 'Multi-mission', '3ML', 'HAWC', 'Fermi', 'HESS', 'joint', 'fit', 'bayesian',
              'multi-wavelength'],

    classifiers=[],

    ext_modules=ext_modules_configuration,

    headers=headers_configuration,

    install_requires=[])

# Now print the final messages if there are any

if len(final_messages) > 0:
    print("\n#############")
    print("FINAL NOTES:")
    print("#############")

    print("\n".join(final_messages))
