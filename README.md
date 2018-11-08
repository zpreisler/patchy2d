# Patchy 2d
![license](https://img.shields.io/badge/license-MIT-green.svg)

Simulation software to simulate patchy particles in 2d.

------

This is a Monte Carlo code. Number of ensembles are supported: Canonical (*NVT*), Grand-Canonical (mu VT), isothermal–isobaric ensemble (*NpT*).
The code is also able to simulate more complex particles such as rods composed by number of single discs.

### Prerequisites
opengl, SDL

### Installing

First one needs to compile & install *lib*

```
git clone https://github.com/zpreisler/lib.git
make
make install
```

Get the source code

```
git clone https://github.com/zpreisler/patchy2d.git
make
make install
```

### Examples

From the command line

```
patchy2d  -n three_patch -s 10000000 -N 512 --npatch 3 --sigma_well 1.1 --patch_width 7 -g 0 -b 30 30 -e 7 --pmod 100 -m 1000
```
The above will generate 512 3-patch particles with interaction range 0.1 and the patch width 7

### Patchy particles
![Snapshot](doc/three_patch.gif)

### Rods
![Snapshot](doc/rods.png)

__Usage:__ patchy2d<br>
[-n] [--name] {name:} Name<br>
[-s] [--step] {step:} Number of MC steps<br>
[-N] [--specie] {specie:} Species<br>
[--new_sigma] {(null)} Update sigma<br>
[--new_width] {(null)} Update width<br>
[--new_mu] {(null)} Update mu<br>
[--new_angle] {(null)} Update angles<br>
[-b] [--box] {box:} Simulation box lengths x and y<br>
[-c] [--copy] {copy:} Copy the simulation box in x and y direction<br>
[-e] [--epsilon] {epsilon:} epsilon -- interaction strength<br>
[--end_epsilon] {end_epsilon:} end epsilon -- final interaction strength<br>
[-p] [--pressure] {pressure:} pressure<br>
[-l] [--lambda] {lambda:} lambda<br>
[-L] [--lambda_coupling] {lambda_coupling:} lambda coupling<br>
[--uy] {uy:} uy -- box tilt<br>
[--max_move] {maximum_displacement:} Maximum displacement<br>
[--max_rot] {maximum_rotation:} Maximum rotation<br>
[--max_vol] {maximum_volume:} Maximum volume change<br>
[--max_uy] {maximum_shape:} Maximum shape change<br>
[--max_dsigma] {maximum_dsigma:} Maximum dsigma<br>
[-m] [--mod] {mod:} Write modulus<br>
[--pmod] {pmod:} Print modulus<br>
[-o] [--optimize] {optimize:} Step size optimization on/off<br>
[-v] [--verbose] {verbose:} Verbose mode on/off<br>
[--snapshot] {snapshot:} Snapshot mode on/off<br>
[--seed] {seed:} Seed<br>
