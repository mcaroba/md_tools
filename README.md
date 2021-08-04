# MD Tools

Copyright (c) 2018-2021 by **Miguel A. Caro**

Mash up of different codes to do analysis of MD data.

See LICENSE.md for license information.

## Usage

Build the codes by typing "make" after editing the Makefile. Add the bin/ directory to your PATH.

### get_steinhardt_ortho

This code only works for orthorhombic unit cells. It requires an XYZ file. You can either pass the
box dimensions (3D array) from the command line:

 get_steinhardt_ortho atoms_filename=traj.xyz box="10. 10. 10." rmax=3.5 atoms_mode=all trajectory_mode=each lmax=8 every=1

or from the XYZ file if it changes dynamically. In that case, omit the Lv keyword in the command line and
provide a box="lx ly lz" keyword in the "comment" line of the XYZ file (the second line for each trajectory frame).

`atoms_mode` and `trajectory_mode` can be `each` or `all` (default is `all` in both cases). `all` averages over
atoms and trajectory frames, respectively, whereas `each` gives individual Q_l per atom and trajectory frame,
respectively. You can combine them at will.
