/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(anc/cut/dipole/long,PairANCCutDipoleLong)

#else

#ifndef LMP_PAIR_ANC_CUT_DIPOLE_LONG_H
#define LMP_PAIR_ANC_CUT_DIPOLE_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairANCCutDipoleLong : public Pair {
 public:
  PairANCCutDipoleLong(class LAMMPS *);
  virtual ~PairANCCutDipoleLong();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  void *extract(const char *, int &);

 protected:
  double cut_lj_global;
  double **cut_lj,**cut_ljsq;
  double cut_coul,cut_coulsq;
  double **epsilon,**rmin,**suavi;
  double **anc1,**anc2,**anc3,**anc4,**offset;
  double g_ewald;
  int ewald_order;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args in pair_style command

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair dipole/long requires atom attributes q, mu, torque

The atom style defined does not have these attributes.

E: Cannot (yet) use 'electron' units with dipoles

This feature is not yet supported.

E: Pair style requires a KSpace style

No kspace style is defined.

*/
