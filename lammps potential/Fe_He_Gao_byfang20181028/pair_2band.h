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

PairStyle(twoband,Pair2BAND)

#else

#ifndef LMP_PAIR_2BAND_H
#define LMP_PAIR_2BAND_H

#include <stdio.h>
#include "pair.h"

namespace LAMMPS_NS {


class Pair2BAND : public Pair {
 public:

//fjzzs potential parameter area start
  double mass_w,mass_he;
//fjzzs for tungsten-tungsten
  double w_af1,w_af2,w_rho_rc,w_rho_0rc,w_rpc;
  double *w_ap,*w_rp,*w_ao,*w_ro; 
//fjzzs for helium-helium
  double herm,c6,c8,c10,heaa,hea,heb,hed,hee,he_rpc;
//fjzzs for tungsten-helum
  double whe_rho_rc,whe_rpc,dNs,psi,rb;
  double *whe_ap,*whe_rp,*whe_ao;
//fjzzs for fermi function
  double bf_w,rf_w,bf_he,rf_he,bf_whe,rf_whe;
//fjzzs for ZBL
  double ev,pi,epsil0,bohrads,exn,zeda2,zeda1;
//  double betasa,betasa1,betasa2,rsa,rsa1,rsa2;
//  double *bzaa,*bzba,*bzaa1,*bzba1,*bzaa2,*bzba2;
  double betasa,rsa;
  double *bzaa,*bzba;
//fjzzs potential parameter area end

  // public variables so USER-ATC package can access them

  double cutmax;

  // potentials as array data


  // potentials in spline form used for force computation

//  double dr,rdr,drho,rdrho,rhomax;
//  double ***rhor_spline,***frho_spline,***z2r_spline;

  Pair2BAND(class LAMMPS *);
  virtual ~Pair2BAND();
  virtual void compute(int, int);

  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

  void settings(int, char **);
  
 protected:
  int nmax;                   // allocated size of per-atom arrays
  double cutforcesq;
  double **scale;

  // per-atom arrays

  double *rho,*fp;
  double *rhos,*fps;
  
  // potentials as file data

  int *map;                   // which element each atom type maps to

  struct Setfl {
    char **elements;
    int nelements;   //,nrho,nr;
//    double drho,dr,cut;
    double *mass;
//   double **frho,**rhor,***z2r;
  };
  Setfl *setfl;  

  virtual void allocate();

  // function by fjz

  double xH0(double);
  double vee_pairww(double);
  double vee_pairwhe(double);
  double vee_pairhehe(double);
  double vee_zbl(double,double,double);
  double dvee_pairww(double);
  double dvee_pairwhe(double);
  double dvee_pairhehe(double);
  double dvee_zbl(double,double,double);
  
  double rho_d(int,int,double);
  double rho_s(int,int,double);
  double drho_d(int,int,double);
  double drho_s(int,int,double);
  
  double emb_d(int,double);
  double emb_s(int,double);
  double demb_d(int,double);
  double demb_s(int,double);
  
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
