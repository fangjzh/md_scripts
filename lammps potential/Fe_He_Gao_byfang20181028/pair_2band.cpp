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
   Contributing authors: JZ Fang  (HNU) Sep. 10 2017. 

   Last mod. 2017/09/15 by JZ Fang 
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_2band.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

Pair2BAND::Pair2BAND(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;

  nmax = 0;
  rho = NULL;
  fp = NULL;
  rhos = NULL;
  fps = NULL;
  map = NULL;

//
//
   setfl = NULL ; // 0915 new

   mass_w=0.0;mass_he=0.0;
//fjzzs for tungsten-tungsten
   w_af1=0;w_af2=0;w_rho_rc=0;w_rpc=0;w_rho_0rc=0;
   w_ap=NULL;w_rp=NULL;w_ao=NULL;w_ro=NULL; 
//fjzzs for helium-helium
   herm=0;c6=0;c8=0;c10=0;heaa=0;hea=0;heb=0;hed=0;hee=0;he_rpc=0;
//fjzzs for tungsten-helum
   whe_rho_rc=0;whe_rpc=0;dNs=0;psi=0;rb=0;
   whe_ap=NULL;whe_rp=NULL;whe_ao=NULL;
//fjzzs for fermi function
   bf_w=0;rf_w=0;bf_he=0;rf_he=0;bf_whe=0;rf_whe=0;
//fjzzs for ZBL
   ev=0;pi=0;epsil0=0;bohrads=0;exn=0;zeda2=0;zeda1=0;
   betasa=0;rsa=0;
   bzaa=NULL;bzba=NULL;

  // set comm size needed by this Pair

  comm_forward = 2;
  comm_reverse = 2;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

Pair2BAND::~Pair2BAND()
{
  if (copymode) return;

  memory->destroy(rho);
  memory->destroy(fp);
  memory->destroy(rhos);
  memory->destroy(fps);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
    delete [] map;
    map = NULL;
  }

}

/* ---------------------------------------------------------------------- */

void Pair2BAND::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,rhoip,rhojp,recip,phip,psip,phi,rhoips,rhojps;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {

    memory->destroy(rho);
    memory->destroy(fp);
    memory->destroy(rhos);
    memory->destroy(fps);

    nmax = atom->nmax;

    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
    memory->create(rhos,nmax,"pair:rhos");
    memory->create(fps,nmax,"pair:fps");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density
  
  //temp 
  double vee_pairt,vee_zblt,dvee_pairt,dvee_zblt;
  double qferm,fermi,dfermi;

  if (newton_pair) {
    for (i = 0; i < nall; i++) {rho[i] = 0.0;rhos[i] = 0.0;}
  } else for (i = 0; i < nlocal; i++) {rho[i] = 0.0;rhos[i] = 0.0;}

  // rho = density at each atom
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
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);
	  
      if (rsq < cutforcesq) {
        jtype = type[j];


        rho[i] += rho_d(map[itype],map[jtype],r);
		rhos[i] += rho_s(map[itype],map[jtype],r);
        if (newton_pair || j < nlocal) {
          
          rho[j] += rho_d(map[jtype],map[itype],r);
		  rhos[j] += rho_s(map[jtype],map[itype],r);
        }
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype=type[i];
    fp[i] = demb_d(map[itype],rho[i]);
	fps[i] = demb_s(map[itype],rhos[i]);
    if (eflag) {
      phi = emb_d(map[itype],rho[i])+emb_s(map[itype],rhos[i]);
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // communicate derivative of embedding function

  comm->forward_comm_pair(this);

  // compute forces on each atom
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
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        r = sqrt(rsq);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip


        rhoip = drho_d(map[jtype],map[itype],r);

        rhojp = drho_d(map[itype],map[jtype],r);

        rhoips = drho_s(map[jtype],map[itype],r);

        rhojps = drho_s(map[itype],map[jtype],r);
		
        recip = 1.0/r;
		
		if((map[itype] == 0 )&&(map[jtype] == 0)){
			qferm  = exp( -bf_w*(r - rf_w) );
			fermi= 1.0/(1.0 + qferm);
			dfermi = bf_w*qferm*fermi*fermi;
			
			vee_pairt=vee_pairww(r);
			vee_zblt=vee_zbl(zeda1,zeda1,r);
			phi = fermi*vee_pairt +(1-fermi)*vee_zblt;
			
			dvee_pairt=dvee_pairww(r);
			dvee_zblt=dvee_zbl(zeda1,zeda1,r);
			phip = (vee_pairt-vee_zblt)*dfermi+(dvee_pairt-dvee_zblt)*fermi+dvee_zblt;
			
		}
		else if(((map[itype] == 0 )&&(map[jtype] == 1))||((map[itype] == 1 )&&(map[jtype] == 0))){
			qferm  = exp( -bf_whe*(r - rf_whe) );
			fermi= 1.0/(1.0 + qferm);
			dfermi = bf_whe*qferm*fermi*fermi;
			
			vee_pairt=vee_pairwhe(r);
			vee_zblt=vee_zbl(zeda1,zeda2,r);
			phi = fermi*vee_pairt +(1-fermi)*vee_zblt;
			
			dvee_pairt=dvee_pairwhe(r);
			dvee_zblt=dvee_zbl(zeda1,zeda2,r);
			phip = (vee_pairt-vee_zblt)*dfermi+(dvee_pairt-dvee_zblt)*fermi+dvee_zblt;
		}
		else {
			phi = vee_pairhehe(r);
			phip = dvee_pairhehe(r);
		}
		

		
        psip = fp[i]*rhojp + fp[j]*rhoip + fps[i]*rhojps + fps[j]*rhoips+ phip;
        fpair = -psip*recip;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) evdwl = phi;
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
void Pair2BAND::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  memory->create(scale,n+1,n+1,"pair:scale");

}
/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void Pair2BAND::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void Pair2BAND::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

//  int ilo,ihi,jlo,jhi;
//  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
//  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  // read EAM setfl file

  if (setfl) {
    for (i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
//    memory->destroy(setfl->frho);
//    memory->destroy(setfl->rhor);
//    memory->destroy(setfl->z2r);
    delete setfl;
  }
  setfl = new Setfl();
  
  setfl->nelements = 2;
   setfl->elements = new char *[setfl->nelements];
   setfl->elements[0] = new char [3];
   setfl->elements[1] = new char [3];
  setfl->elements[0]=(char*)"W";
  setfl->elements[1]=(char*)"He";
   setfl->mass = new double[setfl->nelements];
  setfl->mass[0]=183.84;
  setfl->mass[1]=4.0026;
//  read_file(arg[2]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < setfl->nelements; j++)
      if (strcmp(arg[i],setfl->elements[j]) == 0) break;
    if (j < setfl->nelements) map[i-2] = j;
    else error->all(FLERR,"No matching element in 2BAND potential file");
  }

  // clear setflag since coeff() called once with I,J = * *


  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR,i,setfl->mass[map[i]]);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void Pair2BAND::init_style()
{
//  convert read-in files to standard arrays 
     
memory->destroy(w_ap);
memory->destroy(w_rp);
memory->destroy(w_ao);
memory->destroy(w_ro);

memory->destroy(whe_ap);
memory->destroy(whe_rp);
memory->destroy(whe_ao);

memory->destroy(bzaa);
memory->destroy(bzba);


memory->create(w_ap,15,"Pair:wwap");
memory->create(w_rp,15,"Pair:wwrp");
memory->create(w_ao,4,"Pair:wwao");
memory->create(w_ro,4,"Pair:wwro");

memory->create(whe_ap,9,"Pair:wheap");
memory->create(whe_rp,9,"Pair:wherp");
memory->create(whe_ao,3,"Pair:wwheo");

memory->create(bzaa,4,"Pair:zbla");
memory->create(bzba,4,"Pair:zblb");

/////////////
 ev  = 1.602176487e-19;
 pi      = 3.14159265358979324e0;
 epsil0  = 8.854187818e-12;
 bohrads = 0.5291772;
 exn =  0.23;
 zeda2=2;
 zeda1=74;

    //tungsten EAM2 parameters
     w_ap[ 0] =  0.960851701343041e2;
     w_ap[ 1] = -0.184410923895214e3;
     w_ap[ 2] =  0.935784079613550e2;
     w_ap[ 3] = -0.798358265041677e1;
     w_ap[ 4] =  0.747034092936229e1;
     w_ap[ 5] = -0.152756043708453e1;
     w_ap[ 6] =  0.125205932634393e1;
     w_ap[ 7] =  0.163082162159425e1;
     w_ap[ 8] = -0.141854775352260e1;
     w_ap[ 9] = -0.819936046256149e0;
     w_ap[10] =  0.198013514305908e1;
     w_ap[11] = -0.696430179520267e0;
     w_ap[12] =  0.304546909722160e-1;
     w_ap[13] = -0.163131143161660e1;
     w_ap[14] =  0.138409896486177e1;

     w_rp[ 0] =  2.56489750000e0;
     w_rp[ 1] =  2.62979500000e0;
     w_rp[ 2] =  2.69469250000e0;
     w_rp[ 3] =  2.86631750000e0;
     w_rp[ 4] =  2.97304500000e0;
     w_rp[ 5] =  3.07977250000e0;
     w_rp[ 6] =  3.51647250000e0;
     w_rp[ 7] =  3.84644500000e0;
     w_rp[ 8] =  4.17641750000e0;
     w_rp[ 9] =  4.70084500000e0;
     w_rp[10] =  4.89530000000e0;
     w_rp[11] =  5.08975500000e0;
     w_rp[12] =  5.34295250000e0;
     w_rp[13] =  5.40169500000e0;
     w_rp[14] =  5.46043750000e0;

     w_ao[0]  = -0.420429107805055e1;
     w_ao[1]  =  0.518217702261442e0;
     w_ao[2]  =  0.562720834534370e-1;
     w_ao[3]  =  0.344164178842340e-1;

     w_ro[0]  = 2.50000000000e0;
     w_ro[1]  = 3.10000000000e0;
     w_ro[2]  = 3.50000000000e0;
     w_ro[3]  = 4.90000000000e0;

     w_af1 = -5.946454472402710e0;
     w_af2 = -0.049477376935239e0;

     w_rpc = w_rp[14];
     w_rho_rc = w_ro[3];
     w_rho_0rc = 2.002970124727e0;

//tungsten-helium
      whe_ap[0] =  68.16427676 ;
      whe_ap[1] =  14.28131363 ;
      whe_ap[2] =  -2.71270212 ;
      whe_ap[3] = -20.09063955 ;
      whe_ap[4] =   9.89505155 ;
      whe_ap[5] =  -1.39114722 ;
      whe_ap[6] =  -0.02399923 ;
      whe_ap[7] =   0.43256746 ;
      whe_ap[8] =   0.14893734 ;

      whe_rp[0] =  1.36415840;
      whe_rp[1] =  1.47371367;
      whe_rp[2] =  1.70222409;
      whe_rp[3] =  1.84111795;
      whe_rp[4] =  2.32335983;
      whe_rp[5] =  2.81591553;
      whe_rp[6] =  3.01232737;
      whe_rp[7] =  3.17150879;
      whe_rp[8] =  3.63964951;

      whe_ao[0] =     0.02510451;
      whe_ao[1] =     1.18370044;
      whe_ao[2] = 0             ;
	  
      dNs = 20.0;
      rb = 0.529177210818181818;  //! Bohr radius
      psi = 1.712178504/rb;
	  
      whe_rpc  =  whe_rp[8];
      whe_rho_rc = 3.94e0;

//helium
     herm  = 2.9683e0;
     c6    = 1.35186623e0;
     c8    = 0.41495143e0;
     c10   = 0.17151143e0;
     heaa  = 186924.404e0;
     hea   = 10.5717543e0;
     heb   = -2.07758779e0;
     hed   = 1.438e0;
     hee   = 10.956*1.380658e-23/ev;

// ZBL
   bf_w=24.0;
   rf_w=1.11;
   bf_whe=10.0;
   rf_whe=1.11;

//!ZBL End
//////! Initializing parameters end
	
//  int irequest = neighbor->request(this);
  neighbor->request(this,instance_me);   // 0915 new
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
double Pair2BAND::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  scale[j][i] = scale[i][j];  //  0915 new
  
  cutmax = MAX(whe_rho_rc,whe_rpc);
  cutmax = MAX(w_rho_rc,cutmax);
  cutmax = MAX(w_rpc,cutmax);
  cutmax = MAX(he_rpc,cutmax);

  cutforcesq = cutmax*cutmax;

  return cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

int Pair2BAND::pack_forward_comm(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = fp[j];
	buf[m++] = fps[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void Pair2BAND::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {fp[i] = buf[m++];fps[i] = buf[m++];}
}

/* ---------------------------------------------------------------------- */

int Pair2BAND::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {buf[m++] = rho[i];buf[m++] = rhos[i];}
  return m;
}

/* ---------------------------------------------------------------------- */

void Pair2BAND::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
	rhos[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double Pair2BAND::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   swap fp array with one passed in by caller
------------------------------------------------------------------------- */

  double Pair2BAND::xH0(double x){
	double xH0;
	if(x<0.0) xH0=0.0;
    else xH0=1.0;
    return xH0;
  }
  
  double Pair2BAND::vee_pairww(double r)
  {
	double vee_pairww=0.0;
	if(r<w_rpc){
	vee_pairww = w_ap[ 0]*pow((w_rp[ 0]-r),3)*xH0(w_rp[ 0]-r)
				+w_ap[ 1]*pow((w_rp[ 1]-r),3)*xH0(w_rp[ 1]-r)
				+w_ap[ 2]*pow((w_rp[ 2]-r),3)*xH0(w_rp[ 2]-r)
				+w_ap[ 3]*pow((w_rp[ 3]-r),3)*xH0(w_rp[ 3]-r)
				+w_ap[ 4]*pow((w_rp[ 4]-r),3)*xH0(w_rp[ 4]-r)
				+w_ap[ 5]*pow((w_rp[ 5]-r),3)*xH0(w_rp[ 5]-r)
				+w_ap[ 6]*pow((w_rp[ 6]-r),3)*xH0(w_rp[ 6]-r)
				+w_ap[ 7]*pow((w_rp[ 7]-r),3)*xH0(w_rp[ 7]-r)
				+w_ap[ 8]*pow((w_rp[ 8]-r),3)*xH0(w_rp[ 8]-r)
				+w_ap[ 9]*pow((w_rp[ 9]-r),3)*xH0(w_rp[ 9]-r)
				+w_ap[10]*pow((w_rp[10]-r),3)*xH0(w_rp[10]-r)
				+w_ap[11]*pow((w_rp[11]-r),3)*xH0(w_rp[11]-r)
				+w_ap[12]*pow((w_rp[12]-r),3)*xH0(w_rp[12]-r)
				+w_ap[13]*pow((w_rp[13]-r),3)*xH0(w_rp[13]-r)
				+w_ap[14]*pow((w_rp[14]-r),3)*xH0(w_rp[14]-r) ;
	}
	return vee_pairww;
  }
  
  
  double Pair2BAND::vee_pairwhe(double r)
  {
	  double vee_pairwhe=0.0;
	  if(r<whe_rpc){
	vee_pairwhe =whe_ap[1]*pow((whe_rp[1]-r),3)*xH0(whe_rp[1]-r)
				+whe_ap[2]*pow((whe_rp[2]-r),3)*xH0(whe_rp[2]-r)
				+whe_ap[3]*pow((whe_rp[3]-r),3)*xH0(whe_rp[3]-r)
				+whe_ap[4]*pow((whe_rp[4]-r),3)*xH0(whe_rp[4]-r)
				+whe_ap[5]*pow((whe_rp[5]-r),3)*xH0(whe_rp[5]-r)
				+whe_ap[6]*pow((whe_rp[6]-r),3)*xH0(whe_rp[6]-r)
				+whe_ap[7]*pow((whe_rp[7]-r),3)*xH0(whe_rp[7]-r)
				+whe_ap[8]*pow((whe_rp[8]-r),3)*xH0(whe_rp[8]-r)
				+whe_ap[9]*pow((whe_rp[9]-r),3)*xH0(whe_rp[9]-r);
	  }
	  return vee_pairwhe;
  }
  
  
  double Pair2BAND::vee_pairhehe(double r)
  {
	double vee_pairhehe=0.0;
	double xrm,sumc,term,FFx;
	
	if(r<he_rpc)
	{
	xrm=r/herm;
	sumc=c6/pow(xrm,6)+c8/pow(xrm,8)+c10/pow(xrm,10);
	term=heaa*exp(-hea*xrm+heb*xrm*xrm);
	
	if(xrm<hed)
	{
		FFx=exp(-pow((hed/xrm-1.),2));
		vee_pairhehe = hee*(term-FFx*sumc);
		}
	else
	{
		FFx=1.0;
		vee_pairhehe = hee*(term-FFx*sumc);
    }
	}
	return vee_pairhehe;
  }
  
  double Pair2BAND::vee_zbl(double zeda1,double zeda2,double r)
  {
	  double vee_zbl=0.0;
	  if(r<cutmax){
    betasa = (zeda1*zeda2*ev*ev)/(4.0e0*pi*epsil0);
    rsa = 0.4685017e0/(pow(zeda1,exn) +pow(zeda2,exn));
    bzaa[0] = 0.1818*betasa*1.0e10/ev;
    bzaa[1] = 0.5099*betasa*1.0e10/ev;
    bzaa[2] = 0.2802*betasa*1.0e10/ev;
    bzaa[3] = 0.02817*betasa*1.0e10/ev;
    bzba[0] = -3.2/rsa;
    bzba[1] = -0.9423/rsa;
    bzba[2] = -0.4029/rsa;
    bzba[3] = -0.2016/rsa;
	
	vee_zbl = 1./r*( bzaa[0]*exp(bzba[0]*r)+ 
                     bzaa[1]*exp(bzba[1]*r)+
                     bzaa[2]*exp(bzba[2]*r)+
                     bzaa[3]*exp(bzba[3]*r));
	  }
	  
	  return vee_zbl;
  }
  
  
  double Pair2BAND::dvee_pairww(double r)
    {
	double dvee_pairww=0.0;
	if(r<w_rpc){
      dvee_pairww =  -3.*w_ap[ 0]*pow((w_rp[ 0]-r),2)*xH0(w_rp[ 0]-r)-    
                      3.*w_ap[ 1]*pow((w_rp[ 1]-r),2)*xH0(w_rp[ 1]-r)-    
                      3.*w_ap[ 2]*pow((w_rp[ 2]-r),2)*xH0(w_rp[ 2]-r)-    
                      3.*w_ap[ 3]*pow((w_rp[ 3]-r),2)*xH0(w_rp[ 3]-r)-    
                      3.*w_ap[ 4]*pow((w_rp[ 4]-r),2)*xH0(w_rp[ 4]-r)-    
                      3.*w_ap[ 5]*pow((w_rp[ 5]-r),2)*xH0(w_rp[ 5]-r)-    
                      3.*w_ap[ 6]*pow((w_rp[ 6]-r),2)*xH0(w_rp[ 6]-r)-    
                      3.*w_ap[ 7]*pow((w_rp[ 7]-r),2)*xH0(w_rp[ 7]-r)-    
                      3.*w_ap[ 8]*pow((w_rp[ 8]-r),2)*xH0(w_rp[ 8]-r)-    
                      3.*w_ap[ 9]*pow((w_rp[ 9]-r),2)*xH0(w_rp[ 9]-r)-  
                      3.*w_ap[10]*pow((w_rp[10]-r),2)*xH0(w_rp[10]-r)-  
                      3.*w_ap[11]*pow((w_rp[11]-r),2)*xH0(w_rp[11]-r)-  
                      3.*w_ap[12]*pow((w_rp[12]-r),2)*xH0(w_rp[12]-r)-  
                      3.*w_ap[13]*pow((w_rp[13]-r),2)*xH0(w_rp[13]-r)-  
                      3.*w_ap[14]*pow((w_rp[14]-r),2)*xH0(w_rp[14]-r);
	}
	return dvee_pairww;
  }
  
  double Pair2BAND::dvee_pairwhe(double r)
  {
	  double dvee_pairwhe=0.0;
	  if(r<whe_rpc){
     dvee_pairwhe =  -3.*whe_ap[0]*pow((whe_rp[0]-r),2)*xH0(whe_rp[0]-r)-
                      3.*whe_ap[1]*pow((whe_rp[1]-r),2)*xH0(whe_rp[1]-r)-
                      3.*whe_ap[2]*pow((whe_rp[2]-r),2)*xH0(whe_rp[2]-r)-
                      3.*whe_ap[3]*pow((whe_rp[3]-r),2)*xH0(whe_rp[3]-r)-
                      3.*whe_ap[4]*pow((whe_rp[4]-r),2)*xH0(whe_rp[4]-r)-
                      3.*whe_ap[5]*pow((whe_rp[5]-r),2)*xH0(whe_rp[5]-r)-
                      3.*whe_ap[6]*pow((whe_rp[6]-r),2)*xH0(whe_rp[6]-r)-
                      3.*whe_ap[7]*pow((whe_rp[7]-r),2)*xH0(whe_rp[7]-r)-
                      3.*whe_ap[8]*pow((whe_rp[8]-r),2)*xH0(whe_rp[8]-r);
	  }
	  return dvee_pairwhe;  
  }
  
  
  double Pair2BAND::dvee_pairhehe(double r)
  {
	double dvee_pairhehe=0.0;
	double xrm,sumc,term,FFx;
	
	if(r<he_rpc)
	{
	xrm=r/herm;
	sumc=c6/pow(xrm,6)+c8/pow(xrm,8)+c10/pow(xrm,10);
	term=heaa*exp(-hea*xrm+heb*xrm*xrm);
	
	if(xrm<hed)
	{
		FFx=exp(-pow((hed/xrm-1.),2));
		dvee_pairhehe =  hee*(term*(-hea+2.*heb*xrm)-FFx*sumc*
                        (2.*hed*(hed/xrm-1.)/(xrm*xrm))+FFx*(6.*c6/pow(xrm,7)
                         +8.*c8/pow(xrm,9)+10.*c10/pow(xrm,11)))/herm;
		}
	else
	{
		FFx=1.0;
		dvee_pairhehe = hee*(term*(-hea+2.*heb*xrm)+FFx*
                       (6.*c6/pow(xrm,7)+8.*c8/pow(xrm,9)+10.*c10/pow(xrm,11)))/herm;
    }
	}
	return dvee_pairhehe;  
  }
  
  double Pair2BAND::dvee_zbl(double zeda1,double zeda2,double r)
  {
	  double vee_zbl=0.0,dvee_zbl=0.0;
	  if(r<cutmax){
    betasa = (zeda1*zeda2*ev*ev)/(4.0e0*pi*epsil0);
    rsa = 0.4685017e0/(pow(zeda1,exn) +pow(zeda2,exn));
    bzaa[0] = 0.1818*betasa*1.0e10/ev;
    bzaa[1] = 0.5099*betasa*1.0e10/ev;
    bzaa[2] = 0.2802*betasa*1.0e10/ev;
    bzaa[3] = 0.02817*betasa*1.0e10/ev;
    bzba[0] = -3.2/rsa;
    bzba[1] = -0.9423/rsa;
    bzba[2] = -0.4029/rsa;
    bzba[3] = -0.2016/rsa;
	
	vee_zbl = 1./r*( bzaa[0]*exp(bzba[0]*r)+ 
                     bzaa[1]*exp(bzba[1]*r)+
                     bzaa[2]*exp(bzba[2]*r)+
                     bzaa[3]*exp(bzba[3]*r));
    dvee_zbl = 1./r*(  bzaa[0]*bzba[0]*exp(bzba[0]*r)+   
                       bzaa[1]*bzba[1]*exp(bzba[1]*r)+   
                       bzaa[2]*bzba[2]*exp(bzba[2]*r)+   
                       bzaa[3]*bzba[3]*exp(bzba[3]*r))  
                     -1./r*vee_zbl;
	  }
	  
	  return dvee_zbl; 
  }

  double Pair2BAND::rho_d(int i,int j,double r)
  {
	  double rho_d=0.0;
	  if((i==0)&&(j==0)){
		if(r<w_rho_rc){
			if(r<w_rho_0rc){
				rho_d =  w_ao[0]*pow((w_ro[0]-w_rho_rc),3)*xH0(w_ro[0]-w_rho_rc)+
                         w_ao[1]*pow((w_ro[1]-w_rho_rc),3)*xH0(w_ro[1]-w_rho_rc)+
                         w_ao[2]*pow((w_ro[2]-w_rho_rc),3)*xH0(w_ro[2]-w_rho_rc)+
                         w_ao[3]*pow((w_ro[3]-w_rho_rc),3)*xH0(w_ro[3]-w_rho_rc);
			}
			else
			{
				rho_d = w_ao[0]*pow((w_ro[0]-r),3)*xH0(w_ro[0]-r)+
                        w_ao[1]*pow((w_ro[1]-r),3)*xH0(w_ro[1]-r)+
                        w_ao[2]*pow((w_ro[2]-r),3)*xH0(w_ro[2]-r)+
                        w_ao[3]*pow((w_ro[3]-r),3)*xH0(w_ro[3]-r);
			}
		}
	  }
	  return rho_d;
  }
  
  double Pair2BAND::rho_s(int i,int j,double r)
  {
	double rho_s=0.0;
	if (r<whe_rho_rc){
		if(((i==1)&&(j==0))||((i==0)&&(j==1))){
			rho_s = dNs*pow(r,5)*exp(-2.0*psi*r);
		}
	}
	return rho_s;
  }
  
  
  double Pair2BAND::drho_d(int i,int j,double r)
  {
	  double drho_d=0.0;
	  if((i==0)&&(j==0)){
		if(r<w_rho_rc){
			if(r<w_rho_0rc){
				drho_d =  -3.*w_ao[0]*pow((w_ro[0]-w_rho_rc),2)*xH0(w_ro[0]-w_rho_rc)-
                           3.*w_ao[1]*pow((w_ro[1]-w_rho_rc),2)*xH0(w_ro[1]-w_rho_rc)-
                           3.*w_ao[2]*pow((w_ro[2]-w_rho_rc),2)*xH0(w_ro[2]-w_rho_rc)-
                           3.*w_ao[3]*pow((w_ro[3]-w_rho_rc),2)*xH0(w_ro[3]-w_rho_rc);
			}
			else
			{
				drho_d = -3.*w_ao[0]*pow((w_ro[0]-r),2)*xH0(w_ro[0]-r)-
                          3.*w_ao[1]*pow((w_ro[1]-r),2)*xH0(w_ro[1]-r)-
                          3.*w_ao[2]*pow((w_ro[2]-r),2)*xH0(w_ro[2]-r)-
                          3.*w_ao[3]*pow((w_ro[3]-r),2)*xH0(w_ro[3]-r);
			}
		}
	  }
	  return drho_d;  
  }
  
  double Pair2BAND::drho_s(int i,int j,double r)
  {
	double drho_s=0.0;
	if (r<whe_rho_rc){
		if(((i==1)&&(j==0))||((i==0)&&(j==1))){
			drho_s = dNs*pow(r,4)*exp(-2.0*psi*r)*(5.-2.0*r*psi);
		}
	}
	return drho_s;
  }

  double Pair2BAND::emb_d(int i,double rho)
  {
	  double emb_d=0.0;
	  if(i==0){
		emb_d = w_af1*sqrt(rho) + w_af2*rho*rho;
	  }
	  return emb_d;
  }
  
  double Pair2BAND::emb_s(int i,double rhos)
  {
	  double emb_s=0.0;
	  emb_s = whe_ao[0]*sqrt(rhos) +whe_ao[1]*(rhos*rhos) +whe_ao[2]*pow(rhos,4);
	  return emb_s;
  }
  
   double Pair2BAND::demb_d(int i,double rho)
  {
	  double demb_d=0.0;
	  if(i==0){
		  if(rho>10.0e-15){
		demb_d = 0.5*w_af1/sqrt(rho) + 2.*w_af2*rho;
		  }
	  }
	  return demb_d;
  }
  double Pair2BAND::demb_s(int i,double rhos)
   {
	  double demb_s=0.0;
	  if(rhos>10.e-15){
	  demb_s = (whe_ao[0]*0.5)/sqrt(rhos) +
	            whe_ao[1]*(rhos*2) +
	            whe_ao[2]*pow(rhos,3)*4.0;
	  }
	  return demb_s;
  } 
