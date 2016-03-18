---
layout: default
---

# dimple
Macromolecular crystallography pipeline for refinement and ligand screening.
Based on CCP4 programs.

## How it works ##

You provide the data (merged mtz) and Apo model (pdb) and tell the pipeline
where it should dump the results:

    $ dimple my.mtz apo.pdb output-dir

<script type="text/javascript" src="https://asciinema.org/a/awgcb045doxjstods15xfqi1r.js" id="asciicast-awgcb045doxjstods15xfqi1r" async data-size="13"></script>

Actually you can give any number of models or pdb codes, but
only one of them will be used -- the one with most similar unit cell:

    $ dimple my.mtz apo.pdb my-other.pdb 4uqi output-dir

DIMPLE will do what it can to quickly return a refined model,
pointing to unmodelled electron density blobs -- potential ligand sites.

Simplifying a bit:
the pipeline runs macromolecular refinement after a few usual
preparatory steps (I to F, choosing Rfree set, reindexing if needed).
Sometimes it needs to run Molecular Replacement before refinement.
And at the end it checks for unmodelled blobs -- suspected ligands.

<script type="text/javascript" src="https://asciinema.org/a/awg0n6qr6ez14oe8ugverg4bb.js" id="asciicast-awg0n6qr6ez14oe8ugverg4bb" async data-size="13"></script>

It's quick. Running time depends of course on data, model and computer,
but about 3 minutes is typical. With MR it is usually 3-10 minutes,
but from time to time much, much longer.

### Options ###

DIMPLE has a lot of options (`dimple -h` lists all of them),
but since the goal of the pipeline is to make things simple,
we present here only two of them:

 `--slow` (or `-s` for short) -- recommended if you are not in hurry.
DIMPLE will take twice more time, spending it mostly on extra cycles
of refinement. If this is still too fast, give this option twice --
to get 100 cycles of jelly-body refinement.

`-M` is also quite popular. It decides when MR should kick in.
`-M 0.4` (the default value) runs MR if R-factor after rigid-body refinement,
in data up to 3.5A, is above 0.4.
Edge cases: `-M0` -- always run MR, `-M1` -- never.

### Selected Details ###

**reindexing (Pointless)** - if the data and model are in compatible
but different spacegroups, we change the spacegroup in the MTZ file.
We also check all possible _settings_ - Pointless calculates structure
factors from the model and compares the CC on E^2 to find the best
matching settings, reindexing data if necessary.
In some cases this steps saves as a couple minutes by avoiding MR.

**_free_ reflections** -
Rfree statistic depend to some degree on how lucky is the pseudo-random
set of free flags. To eliminate this luck factor when comparing
different data collection one may want to use the same set of free
flags. It is possible by passing external set of reference flags
(option `--free-r-flags`), but we wanted to do even without the
reference file. This was implemented by generating the same free set
for the same pdb file and should work if the space group is the same
and the resolution is below 1A (we had to pick arbitrary limit).

**different asu volume** -
if asu in the data is much larger than in the model,
we search for multiple copies of the model in MR.
In the opposite case, when asu in the data is smaller,
we make a single ensemble from all the chains
(*phaser.ensembler*) before MR.

**scoring blobs** - it is rather simplistic now, we need to work on it

**generating pictures** - we have an option (`-f`) to generate static
images (PNG or JPEG) of the blobs. They are used by
[SynchWeb](https://github.com/DiamondLightSource/SynchWeb) in DLS.
Pictures are generated with Coot+Raster3d - this combines
the familiar look and feel of Coot with nicer graphics and headless
rendering.

[<img src="http://i.imgur.com/cTDWIzK.png" width="290px"/>](http://i.imgur.com/cTDWIzK.png)
[<img src="http://i.imgur.com/e61cB7G.png" width="290px"/>](http://i.imgur.com/e61cB7G.png)


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
If it doesn't work as expected or doesn't work at all - let us know asap.
Use the [issue tracker](https://github.com/ccp4/dimple/issues) or
email CCP4 helpdesk or
email wojdyr@gmail.com or
[![chat at https://gitter.im/ccp4/dimple](https://badges.gitter.im/ccp4/dimple.svg)](https://gitter.im/ccp4/dimple).

