/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Tonalli Rodriguez-Lopez (xclassicx@gmail.com)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_anc_cut.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairANCCut::PairANCCut(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairANCCut::~PairANCCut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(rmin);
    memory->destroy(suavi);
    memory->destroy(anc1);
    memory->destroy(anc2);
    memory->destroy(anc3);
    memory->destroy(anc4);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairANCCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rsq,rcu,r2inv,ranc3inv,ranc6inv,forceanc,factor_anc;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_anc = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_anc = special_anc[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        double softness = suavi[itype][jtype];
        double rmn1     = rmin[itype][jtype];
        double rmn3     = rmn1*rmn1*rmn1;
        r     = sqrt(rsq);
        rcu   = r*r*r;
        r2inv = 1.0/rsq;
        ranc3inv = 1.0/((rcu - rmn3) / softness + rmn3);
        ranc6inv = ranc3inv*ranc3inv;

        forceanc = ranc6inv * (anc1[itype][jtype]*ranc6inv - anc2[itype][jtype])*rcu*ranc3inv/softness;
        fpair = factor_anc*forceanc*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = ranc6inv*(anc3[itype][jtype]*ranc6inv-anc4[itype][jtype]) -
            offset[itype][jtype];
          evdwl *= factor_anc;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairANCCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(rmin,n+1,n+1,"pair:rmin");
  memory->create(suavi,n+1,n+1,"pair:suavi");
  memory->create(anc1,n+1,n+1,"pair:anc1");
  memory->create(anc2,n+1,n+1,"pair:anc2");
  memory->create(anc3,n+1,n+1,"pair:anc3");
  memory->create(anc4,n+1,n+1,"pair:anc4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairANCCut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairANCCut::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double rmin_one = force->numeric(FLERR,arg[3]);
  double suavi_one = force->numeric(FLERR,arg[4]);

  double cut_one = cut_global;
  if (narg == 6) cut_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      rmin[i][j] = rmin_one;
      suavi[i][j] = suavi_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairANCCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  anc1[i][j] = 12.0 * epsilon[i][j] * pow(rmin[i][j],12.0);
  anc2[i][j] = 12.0 * epsilon[i][j] * pow(rmin[i][j],6.0);
  anc3[i][j] = epsilon[i][j] * pow(rmin[i][j],12.0);
  anc4[i][j] = 2.0 * epsilon[i][j] * pow(rmin[i][j],6.0);

  if (offset_flag) {
    double softness = suavi[i][j];
    double rmn1     = rmin[i][j];
    double rmn3     = rmn1*rmn1*rmn1;
    double rc1      = cut[i][j];
    double rc3      = rc1*rc1*rc1;
    double ratio = rmn3 / ((rc3 - rmn3) / softness + rmn3);
    offset[i][j] = epsilon[i][j] * (pow(ratio,4.0) - 2.0*pow(ratio,2.0));
  } else offset[i][j] = 0.0;

  anc1[j][i] = anc1[i][j];
  anc2[j][i] = anc2[i][j];
  anc3[j][i] = anc3[i][j];
  anc4[j][i] = anc4[i][j];
  offset[j][i] = offset[i][j];

  epsilon[j][i] = epsilon[i][j];
  rmin[j][i]    = rmin[i][j];
  suavi[j][i]   = suavi[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double softness = suavi[i][j];
    double rmn1     = rmin[i][j];
    double rmn3     = rmn1*rmn1*rmn1;
    double rmn6 = rmn3*rmn3;
    double rmn9 = rmn3*rmn6;
    double rc1      = cut[i][j];
    double rc3      = rc1*rc1*rc1;
    double zct3 = (rc3 - rmn3) / softness + rmn3;
    double zct6 = zct3*zct3;
    double zct9 = zct3*zct6;
    etail_ij = 2.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      rmn6 * softness * (rmn6 - 6.0*zct6) / (9.0*zct9);
    ptail_ij = 2.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      rmn6 * (-3.0*rmn9*(softness - 1.0)/zct3 + 
        4.0*rmn6*softness +
        6.0*rmn3*(softness-1.0)*zct3 -
        12*softness*zct6) / (9.0*zct9);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairANCCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&rmin[i][j],sizeof(double),1,fp);
        fwrite(&suavi[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairANCCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&rmin[i][j],sizeof(double),1,fp);
          fread(&suavi[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&rmin[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&suavi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairANCCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairANCCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairANCCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,epsilon[i][i],rmin[i][i],suavi[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairANCCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,epsilon[i][j],rmin[i][j],suavi[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairANCCut::single(int i, int j, int itype, int jtype,
                      double rsq, double factor_coul, double factor_anc,
                      double &fforce)
{
  double rcu,r2inv,ranc3inv,ranc6inv,forceanc,phianc;
  double softness = suavi[itype][jtype];
  double rmn1     = rmin[itype][jtype];
  double rmn3     = rmn1*rmn1*rmn1;
  double r        = sqrt(rsq);
  rcu   = r*r*r;
  r2inv = 1.0/rsq;
  ranc3inv = 1.0/((rcu - rmn3) / softness + rmn3);
  ranc6inv = ranc3inv*ranc3inv;
  forceanc = ranc6inv * (anc1[itype][jtype]*ranc6inv - anc2[itype][jtype])*rcu*ranc3inv/softness;
  fforce = factor_anc*forceanc*r2inv;

  phianc = ranc6inv*(anc3[itype][jtype]*ranc6inv-anc4[itype][jtype]) -
    offset[itype][jtype];
  return factor_anc*phianc;
}

/* ---------------------------------------------------------------------- */

void *PairANCCut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"rmin") == 0) return (void *) rmin;
  if (strcmp(str,"suavi") == 0) return (void *) suavi;
  return NULL;
}
