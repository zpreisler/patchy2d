# Patchy 2d

Simulation software to simulate patchy particles in 2d.

### Patchy particles
![Snapshot](doc/b.png)

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
