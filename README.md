# dimple
Macromolecular crystallography pipeline for refinement and ligand screening.

#### See docs at http://ccp4.github.io/dimple

## Installation ##

[![Build Status](https://travis-ci.org/ccp4/dimple.svg?branch=master)
](https://travis-ci.org/ccp4/dimple)

DIMPLE is part of the CCP4 suite. All the work is done
by other programs from the suite, which are run underneath.
Most importantly _refmac_ and _phaser_.
So you need to have CCP4 installed.

To test the latest development version, clone this repository
and run `./dimple/dimple`.

Dimple comes with a little `find-blobs` utility.
If you don't compile it, the version from the CCP4 suite will be used.
Compilation requires clipper library which in turn requires mmdb2 and libccp4.
You may either build everything from source or use conda to download binaries
compiled with GCC 4.8. The latter can be done by
[installing miniconda](http://conda.pydata.org/miniconda.html)
and doing:

    conda install -c mw clipper

Than build find-blobs pointing where the required libraries are:

    export CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0 -std=c++98"  # for GCC 5
    cmake -D CMAKE_PREFIX_PATH=$HOME/miniconda2 .
    make

No need for make install - the binary landed in the same directory
as the python scripts and dimple will use programs from here
in preference to the ones from `$CCP4/bin`.

If it doesn't work - see the contact methods below.

## Comments? ##

Any thoughts, comments or feedback is welcome.
If it doesn't work as expected or doesn't work at all - let us know asap.
Use the [issue tracker](https://github.com/ccp4/dimple/issues) or
email CCP4 helpdesk or
email wojdyr@gmail.com or
[![chat at https://gitter.im/ccp4/dimple](https://badges.gitter.im/ccp4/dimple.svg)](https://gitter.im/ccp4/dimple).

---

Made in Diamond Light Source and in CCP4-Harwell.
