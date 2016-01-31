# dimple
MX pipeline for refinement and ligand screening.

License: Apache.

Requires several programs from the CCP4 suite to run.

## What it does ##

You provide data (merged mtz) and model (pdb).
And a directory where the output should go.

    $ dimple my.mtz my.pdb output-dir

If you wish you can give a few potential models or pdb codes:

    $ dimple my.mtz my.pdb my-other.pdb 4uqi output-dir

DIMPLE will do what it can to quickly return refined model,
pointing to unmodelled electron density blobs - potential ligand sites.

## How it works ##

DIMPLE is part of the CCP4 suite. All the work is done
by other programs from the suite, which are run underneath.
Most importantly _refmac_ and _phaser_.

Simplifying a bit:
the pipeline runs macromolecular refinement after a few usual
preparatory steps (I to F, choosing Rfree set, reindexing if needed).
Sometimes it needs to run also Molecular Replacement before refinement.
And at the end it checks for unmodelled blobs - suspected ligands.

It's quick. Running time depends of course on data, model and computer,
but about 3 minutes is typical. With MR it is usually 5-10 minutes,
but from time to time much longer.

## Options ##

It has a lot of options (`dimple -h` lists all of them),
but since the goal of the pipeline is to make things simple,
we'll focus on one:

 `--slow` (or `-s` for short) -- recommended if you are not in hurry.
DIMPLE will take twice more time, spending it mostly on extra cycles
of refinement. If this is still too fast, give this option twice
- to get 100 cycles of jelly-body refinement.

`-M` is also quite popular. It decides when MR should kick in:
`-M0` - always, `-M1` - never, `-M 0.4` - (default) if after rigid-body
refinement Rfactor in data up to 3.5A is above 0.4.

## Selected Details ##

TODO - a few steps in the pipeline should be explained in details

Generating pictures -

Choosing _free_ set -

Reindexing (Pointless) -

Finding and scoring blobs - it is rather simplistic now, we need to work on it


## Comments? ##

Any thoughts, comments or feedback is welcome.
If it doesn't work as expected - let us know asap.
Feel free to use issue tracker or email CCP4 helpdesk or wojdyr@gmail.com.

---

Made in Diamond Light Source and in CCP4-Harwell.
