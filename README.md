# dimple
Macromolecular crystallography pipeline for refinement and ligand screening.

Requires several programs from the CCP4 suite to run.

## How it works ##

You provide data (merged mtz) and model (pdb) and tell the pipeline
where it should dump the results:

    $ dimple my.mtz apo.pdb output-dir

If you wish you can give more than one model or pdb code
(only one of them will be used - the one with most similar unit cell):

    $ dimple my.mtz apo.pdb my-other.pdb 4uqi output-dir

DIMPLE will do what it can to quickly return refined model,
pointing to unmodelled electron density blobs - potential ligand sites.

Simplifying a bit:
the pipeline runs macromolecular refinement after a few usual
preparatory steps (I to F, choosing Rfree set, reindexing if needed).
Sometimes it needs to run also Molecular Replacement before refinement.
And at the end it checks for unmodelled blobs - suspected ligands.

It's quick. Running time depends of course on data, model and computer,
but about 3 minutes is typical. With MR it is usually 5-10 minutes,
but from time to time much, much longer.

## Options ##

DIMPLE has a lot of options (`dimple -h` lists all of them),
but since the goal of the pipeline is to make things simple,
we'll focus on one:

 `--slow` (or `-s` for short) -- recommended if you are not in hurry.
DIMPLE will take twice more time, spending it mostly on extra cycles
of refinement. If this is still too fast, give this option twice
- to get 100 cycles of jelly-body refinement.

`-M` is also quite popular. It decides when MR should kick in:
`-M0` - always, `-M1` - never, `-M 0.4` - (default) if Rfactor after
rigid-body refinement is above 0.4. Rfactor for this purpose
is calculated only in data up to 3.5A.

## Installation ##

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
(with g++4 ABI). The latter can be done by
[installing miniconda](http://conda.pydata.org/miniconda.html)
and doing:

    conda install -c mw clipper

Than build find-blobs pointing where the required libraries are:

    cmake . -D CMAKE_PREFIX_PATH=$HOME/miniconda2
    make  # no need for make install

If it doesn't work - see the contact methods below.

## Selected Details ##

TODO - a few steps in the pipeline should be explained in details

Generating pictures -

Choosing _free_ set -

Reindexing (Pointless) -

Finding and scoring blobs - it is rather simplistic now, we need to work on it


## Comments? ##

Any thoughts, comments or feedback is welcome.
If it doesn't work as expected or doesn't work at all - let us know asap.
Feel free to use the issue tracker here or email CCP4 helpdesk or
email wojdyr@gmail.com.

---

Made in Diamond Light Source and in CCP4-Harwell.
