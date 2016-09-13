This package implements the Approximate Non-Conformal (ANC) potential in LAMMPS.
The ANC potential covers a wide range of interactions, from hard spheres to
a reference potential, with a tunable parameter called softness.

See the doc/USER/anc pages for the anc/cut, anc/cut/coul/long, and anc/cut/dipole/long pair styles.

The ANC potential is described in the following papers:

- _Nonconformal Potentials and Second Virial Coefficients in Molecular Fluids. 1. Theory_,
  F. del Río, J. E. Ramos, and I. A. McLure,
  [J. Phys. Chem. B 102, 10568 (1998)](http://dx.doi.org/10.1021/jp9831684).  
- _Nonconformal Interaction Models and Thermodynamics of Polar Fluids_,
  E. Ávalos, F. del Río, and S. Lago, 
  [J. Phys. Chem. B 109, 508 (2005)](http://pubs.acs.org/doi/abs/10.1021/jp046735y).

It is important to note that, in this package, the reference interaction is the *Lennard-Jones* potential.

**WARNING**: This is not a final version.
There are some issues that need attention:

- How to deal with the hard core
- Optimization regarding potential parameters

The person who created this package is  Tonalli Rodriguez-Lopez <xclassicx at gmail dot com>
while at the University of Waterloo. Contact him directly if you have questions.
