/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(twobandl,Pair2BANDL)

#else

#ifndef LMP_PAIR_2BANDL_H
#define LMP_PAIR_2BANDL_H

#include <stdio.h>
#include "pair.h"

namespace LAMMPS_NS {


class Pair2BANDL : public Pair {
 public:
  friend class FixSemiGrandCanonicalMC;   // Alex Stukowski option

  // public variables so USER-ATC package can access them

  double cutmax;

  // potentials as array data

  int nrho,nrhos,nr;
  int nfrho,nfrhos,nrhor,nrhors,nz2r;
  double **frho,**rhor,**frhos,**rhors,**z2r;
  int *type2frho,**type2rhor,*type2frhos,**type2rhors,**type2z2r;

  // potentials in spline form used for force computation

  double dr,rdr,drho,rdrho,rhomax,drhos,rdrhos,rhosmax;   // rdr = 1./dr
  double ***rhor_spline,***frho_spline,***rhors_spline,***frhos_spline,***z2r_spline;

  Pair2BANDL(class LAMMPS *);
  virtual ~Pair2BANDL();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void swap_eam(double *, double **);

 protected:
  int nmax;                   // allocated size of per-atom arrays
  double cutforcesq;
  double **scale;

  // per-atom arrays

  double *rho,*fp,*rhos,*fps;

  // potentials as file data

  int *map;                   // which element each atom type maps to
/* 
  struct Funcfl {
    char *file;
    int nrho,nrhos,nr;
    double drho,drhos,dr,cut,mass;
    double *frho,*frhos,*rhor,*rhors,*zr;
  };
  Funcfl *funcfl;
  int nfuncfl;
 */
  struct Setfl {
    char **elements;
    int nelements,nrho,nrhos,nr;
    double drho,drhos,dr,cut;
    double *mass;
    double **frho,**frhos,**rhor,**rhors,***z2r;
  };
  Setfl *setfl;


  virtual void allocate();
  virtual void array2spline();
  void interpolate(int, double, double *, double **);
  void grab(FILE *, int, double *);

  virtual void read_file(char *);
  virtual void file2array();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot open EAM potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

*/
