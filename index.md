---
layout: default
---

# dim<span id="p">p</span>le

#### Macromolecular crystallography pipeline

#### for refinement and ligand screening

#### based on CCP4 programs.


## How it works ##

You provide the reflection data (merged mtz) and Apo model (pdb) and tell the pipeline
where it should dump the results:

    $ dimple my.mtz apo.pdb output-dir

<script type="text/javascript" src="https://asciinema.org/a/awgcb045doxjstods15xfqi1r.js" id="asciicast-awgcb045doxjstods15xfqi1r" async data-size="13"></script>

Actually you can give any number of models or pdb codes, but
only one of them will be used -- the one with most similar unit cell:

    $ dimple my.mtz apo.pdb my-other.pdb 4uqi output-dir

DIMPLE will do what it can to quickly return a refined model and difference density map
pointing to unmodelled electron density blobs -- potential ligand sites.

Simplifying a bit:
the pipeline runs macromolecular refinement after a few usual
preparatory steps (I to F, choosing Rfree set, reindexing if needed).
Sometimes it needs to run Molecular Replacement before refinement.
And at the end it checks for unmodelled blobs -- suspected ligands.

<script type="text/javascript" src="https://asciinema.org/a/awg0n6qr6ez14oe8ugverg4bb.js" id="asciicast-awg0n6qr6ez14oe8ugverg4bb" async data-size="13"></script>

It's quick. Run time depends of course on the data, such as resolution, size of unit cell, etc., the model and computer,
but about 3 minutes is typical. With MR it is usually 3-10 minutes,
but from time to time much, much longer.

In Diamond (UK's national synchrotron) Dimple is automatically run after
auto-processing [if the user provided a PDB file](http://www.diamond.ac.uk/Beamlines/Mx/Common/Common-Manual/Data-Analysis/Automated-Software-Pipeline.html#DIMPLE)
before the experiment.

## Options ##

DIMPLE has a lot of options (`dimple -h` lists all of them),
but since the goal of the pipeline is to make things simple,
we present here only two of them:

 `--slow` (or `-s` for short) -- recommended if you are not in hurry.
DIMPLE will take twice as long, spending more time on extra cycles
of refinement. If this is still too fast, give this option twice --
to get 100 cycles of jelly-body refinement.

`-M` is also quite popular. It decides when MR should kick in.
`-M 0.4` (the default value) runs MR if the R-factor after rigid-body refinement,
in data up to 3.5Å, is above 0.4.
Edge cases: `-M0` -- always run MR, `-M1` -- never.

## Selected Details ##

**MTZ columns** --
by default, Dimple expects either intensities labeled as IMEAN or
amplitudes labeled as F or FP. To use different columns
specify options `--icolumn` / `--fcolumn`.
Intensities are used only when amplitudes are not in the file, or
when the user explicitely sets `--icolumn` or `--ItoF-prog`.

**Reindexing** --
if the data and model are in compatible but different spacegroups,
we use *pointless* to change the spacegroup in the MTZ file.
We also check all possible _settings_ -- pointless calculates structure
factors from the model and compares the CC on E^2 to find the best
matching settings, reindexing data if necessary.
In some cases this steps saves us a couple minutes by avoiding MR.

**Free reflections** --
Rfree statistics depend to some degree on how lucky is the pseudo-random
set of free flags. Therefore one may prefer to use the same free flags
when comparing data collected from different crystals. That is why Dimple
always assigns the same flags when the same pdb file is used.
Unless the data is in different point group than the model
or the resolution is below 1Å (we had to pick an arbitrary limit).
Alternatively, you can pass an external set of reference flags
(option `--free-r-flags`).

**Different asu volume** --
if asu in the data is much larger than in the model,
we search for multiple copies of the model in MR.
In the opposite case, when asu in the data is smaller,
we make a single ensemble from all the chains
(*phaser.ensembler*) before MR.

**NCS** --
The PDB file format stores non-crystallographic symmetry transformations
in the MTRIXn records. The coordinates of NCS copies may or may not be
in the file, as reflected (with some exceptions) by the MTRIXn/iGiven field.
Dimple is generally aware of these things and should work reasonably in all
cases, although as of Q2 2016 it does not, waiting for a bugfix in refmac.

**Rigid-body** refinement --
both *refmac* and *phaser* can do it.
We could simplify our pipeline by always running phaser before
refmac, but it would be slower (on average).
Not by a large margin, though. Phaser checks if the input model is
already placed correctly and skips the search if it is.

Since we only use data up to 3.5Å for rigid-body, in some cases
5% of free reflections was not enough to give reliable statistics.
On the other hand there is no danger of overfitting in rigid-body.
Thus, we stopped using Rfree set in this step at all.

Actual **refinement** --
we run refmac restrained refinement twice.
The first run (labelled as *jelly*) has jelly-body restraints,
no hydrogens and ignores very high resolution reflections.
In the default mode we run only 4 cycles with these settings,
and 8 cycles of the final refinement. Note that the usual
recommendation for jelly-body refinement is 100+ cycles.
We picked the numbers 4 and 8 after testing various combinations
on hundreds of datasets. This split happened to give slightly
better results than other combinations within the same time limit.

**Scoring blobs** --
it is rather simplistic now, we need to work on it

**Generating pictures** --
we have an option (`-f`) to generate static
images (PNG or JPEG) of the blobs. They are used by
[SynchWeb](https://github.com/DiamondLightSource/SynchWeb) in DLS.
Pictures are generated with Coot+Raster3d - this combines
the familiar look and feel of Coot with nicer graphics and headless
rendering.

[<img src="http://ccp4.github.io/img/blob-th.png" width="290px"/>](http://ccp4.github.io/img/blob-th.png)
[<img src="http://ccp4.github.io/img/blob2.png" width="290px"/>](http://ccp4.github.io/img/blob2.png)

Currently we are working on
[interactive web-based viewer](https://github.com/uglymol/uglymol)
that emulates Coot and hopefully works fast enough to replace static images.

## FAQ ##

* _Why is the final R factor in one refmac run different than the initial
  R factor in succeeding run?_

  Because of different refinement options (hydrogens, resolution).

* _Is Dimple based on the EDNA framework?_

  No, but it was at some point, as an experiment. Being a relatively simple
  project, dimple was well-suited for testing new frameworks.
  The current incarnation of Dimple was started in 2013
  preserving the workflow and parameters from the previous version.
  And then gradually evolved.

## Installation ##

DIMPLE is part of the CCP4 suite. All the work is done
by other programs from the suite, which are run underneath.
Most importantly _refmac_ and _phaser_.
So you need to have CCP4 installed, and `dimple` is already there.

But if you'd like to try the very latest version that hasn't filtered
through to the CCP4 suite yet, get it from
[github.com/ccp4/dimple](https://github.com/ccp4/dimple).

## Comments? ##

Any comments and thoughts how to improve this tool are genuinely welcome.
If it doesn't work as expected or doesn't work at all -- let us know asap.
Use the [issue tracker](https://github.com/ccp4/dimple/issues) or
email CCP4 helpdesk or
email wojdyr@gmail.com or
[![chat at https://gitter.im/ccp4/dimple](https://badges.gitter.im/ccp4/dimple.svg)](https://gitter.im/ccp4/dimple).

