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

   Last mod. 2017/09/17 by JZ Fang 
------------------------------------------------------------------------- */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_2bandl.h"
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

Pair2BANDL::Pair2BANDL(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;

  nmax = 0;
  rho = NULL;
  fp = NULL;
  rhos = NULL;
  fps = NULL;
  map = NULL;
  type2frho = NULL;
  type2frhos = NULL;
  
//  nfuncfl = 0;
//  funcfl = NULL;

  setfl = NULL;
 // fs = NULL;

  frho = NULL;
  rhor = NULL;
  frhos = NULL;
  rhors = NULL;
  
  z2r = NULL;
  scale = NULL;

  frho_spline = NULL;
  rhor_spline = NULL;
  
  frhos_spline = NULL;
  rhors_spline = NULL;
  
  z2r_spline = NULL;

  // set comm size needed by this Pair

  comm_forward = 2;
  comm_reverse = 2;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

Pair2BANDL::~Pair2BANDL()
{
  if (copymode) return;

  memory->destroy(rho);
  memory->destroy(fp);
  
  memory->destroy(rhos);
  memory->destroy(fps);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] type2frho;
    delete [] type2frhos;
    map = NULL;
    type2frho = NULL;
    type2frhos = NULL;
    memory->destroy(type2rhor);
    memory->destroy(type2rhors);
    memory->destroy(type2z2r);
    memory->destroy(scale);
  }


  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->frhos);
    memory->destroy(setfl->rhors);
    memory->destroy(setfl->z2r);
    delete setfl;
    setfl = NULL;
  }


  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(frhos);
  memory->destroy(rhors);
  memory->destroy(z2r);

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(frhos_spline);
  memory->destroy(rhors_spline);
  memory->destroy(z2r_spline);
}

/* ---------------------------------------------------------------------- */

void Pair2BANDL::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,rhoips,rhojps,z2,z2p,recip,phip,psip,phi;
  double *coeff,*coeffs;
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

  if (newton_pair) {
    for (i = 0; i < nall; i++) {rho[i] = 0.0;rhos[i] = 0.0;}
  } else for (i = 0; i < nlocal; i++){ rho[i] = 0.0; rhos[i] = 0.0;}

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

      if (rsq < cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
		
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		
        coeffs = rhors_spline[type2rhors[jtype][itype]][m];
        rhos[i] += ((coeffs[3]*p + coeffs[4])*p + coeffs[5])*p + coeffs[6];
		
        if (newton_pair || j < nlocal) {
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
		  
          coeffs = rhors_spline[type2rhors[itype][jtype]][m];
          rhos[j] += ((coeffs[3]*p + coeffs[4])*p + coeffs[5])*p + coeffs[6];
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
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (eflag) {
      phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      if (rho[i] > rhomax) phi += fp[i] * (rho[i]-rhomax);
      phi *= scale[type[i]][type[i]];   //??????
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = rhos[i]*rdrhos + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrhos-1));
    p -= m;
    p = MIN(p,1.0);
    coeffs = frhos_spline[type2frhos[type[i]]][m];
    fps[i] = (coeffs[0]*p + coeffs[1])*p + coeffs[2];
    if (eflag) {
	  phi = ((coeffs[3]*p + coeffs[4])*p + coeffs[5])*p + coeffs[6];
	  if (rhos[i] > rhosmax) phi += fps[i] * (rhos[i]-rhosmax);
      phi *= scale[type[i]][type[i]];   //??????
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
        p = r*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
        // scale factor can be applied by thermodynamic integration

        coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
		
        coeffs = rhors_spline[type2rhors[itype][jtype]][m];
        rhoips = (coeffs[0]*p + coeffs[1])*p + coeffs[2];
        coeffs = rhors_spline[type2rhors[jtype][itype]][m];
        rhojps = (coeffs[0]*p + coeffs[1])*p + coeffs[2];
		
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r;
        phi = z2*recip;
        phip = z2p*recip - phi*recip;
        psip = fp[i]*rhojp + fp[j]*rhoip + fps[i]*rhojps + fps[j]*rhoips + phip;
        fpair = -scale[itype][jtype]*psip*recip;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) evdwl = scale[itype][jtype]*phi;
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

void Pair2BANDL::allocate()
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

  type2frho = new int[n+1];
  type2frhos = new int[n+1];
  
  memory->create(type2rhor,n+1,n+1,"pair:type2rhor");
  memory->create(type2rhors,n+1,n+1,"pair:type2rhors");
  
  memory->create(type2z2r,n+1,n+1,"pair:type2z2r");
  memory->create(scale,n+1,n+1,"pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void Pair2BANDL::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void Pair2BANDL::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read EAM setfl file

  if (setfl) {
    for (i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->frhos);
    memory->destroy(setfl->rhors);
    memory->destroy(setfl->z2r);
    delete setfl;
  }
  setfl = new Setfl();
  read_file(arg[2]);

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
    else error->all(FLERR,"No matching element in EAM potential file");
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

void Pair2BANDL::init_style()
{
  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double Pair2BANDL::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  scale[j][i] = scale[i][j];

   
  if (setfl) cutmax = setfl->cut;
  

  cutforcesq = cutmax*cutmax;

  return cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */

void Pair2BANDL::read_file(char *filename)
{
  Setfl *file = setfl;

  // open potential file

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = force->open_potential(filename);
    if (fptr == NULL) {
      char str[128];
      sprintf(str,"Cannot open EAM potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read and broadcast header
  // extract element names from nelements line

  int n;
  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  sscanf(line,"%d",&file->nelements);
  int nwords = atom->count_words(line);
  if (nwords != file->nelements + 1)
    error->all(FLERR,"Incorrect element names in EAM potential file");

  char **words = new char*[file->nelements+1];
  nwords = 0;
  strtok(line," \t\n\r\f");
  while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

  file->elements = new char*[file->nelements];
  for (int i = 0; i < file->nelements; i++) {
    n = strlen(words[i]) + 1;
    file->elements[i] = new char[n];
    strcpy(file->elements[i],words[i]);
  }
  delete [] words;

  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg %d %lg %d %lg %lg",
           &file->nrho,&file->drho,&file->nrhos,&file->drhos,&file->nr,&file->dr,&file->cut);
  }

  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&file->nrhos,1,MPI_INT,0,world);
  MPI_Bcast(&file->drhos,1,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  file->mass = new double[file->nelements];
  memory->create(file->frho,file->nelements,file->nrho+1,"pair:frho");
  memory->create(file->rhor,file->nelements,file->nr+1,"pair:rhor");
  
  memory->create(file->frhos,file->nelements,file->nrhos+1,"pair:frhos");
  memory->create(file->rhors,file->nelements,file->nr+1,"pair:rhors");
  
  memory->create(file->z2r,file->nelements,file->nelements,file->nr+1,
                 "pair:z2r");

  int i,j,tmp;
  for (i = 0; i < file->nelements; i++) {
    if (me == 0) {
      fgets(line,MAXLINE,fptr);
      sscanf(line,"%d %lg",&tmp,&file->mass[i]);
    }
    MPI_Bcast(&file->mass[i],1,MPI_DOUBLE,0,world);

    if (me == 0) grab(fptr,file->nrho,&file->frho[i][1]);
    MPI_Bcast(&file->frho[i][1],file->nrho,MPI_DOUBLE,0,world);
    if (me == 0) grab(fptr,file->nr,&file->rhor[i][1]);
    MPI_Bcast(&file->rhor[i][1],file->nr,MPI_DOUBLE,0,world);
	
    if (me == 0) grab(fptr,file->nrhos,&file->frhos[i][1]);
    MPI_Bcast(&file->frhos[i][1],file->nrhos,MPI_DOUBLE,0,world);
    if (me == 0) grab(fptr,file->nr,&file->rhors[i][1]);
    MPI_Bcast(&file->rhors[i][1],file->nr,MPI_DOUBLE,0,world);
	
  }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fptr,file->nr,&file->z2r[i][j][1]);
      MPI_Bcast(&file->z2r[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }

  // close the potential file

  if (me == 0) fclose(fptr);
}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void Pair2BANDL::file2array()
{
  int i,j,m,n;
  int ntypes = atom->ntypes;

  // set function params directly from setfl file

  nrho = setfl->nrho;
  nrhos = setfl->nrhos;
  nr = setfl->nr;
  drho = setfl->drho;
  drhos = setfl->drhos;
  dr = setfl->dr;
  rhomax = (nrho-1) * drho;
  rhosmax = (nrhos-1) * drhos;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of setfl elements + 1 for zero array

  nfrho = setfl->nelements + 1;
  memory->destroy(frho);
  memory->create(frho,nfrho,nrho+1,"pair:frho");
  
  nfrhos = setfl->nelements + 1;
  memory->destroy(frhos);
  memory->create(frhos,nfrhos,nrhos+1,"pair:frhos");

  // copy each element's frho to global frho

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nrho; m++) frho[i][m] = setfl->frho[i][m];

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nrhos; m++) frhos[i][m] = setfl->frhos[i][m];

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  for (m = 1; m <= nrhos; m++) frhos[nfrhos-1][m] = 0.0;
  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;
	
  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frhos[i] = map[i];
    else type2frhos[i] = nfrhos-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of setfl elements

  nrhor = setfl->nelements;
  memory->destroy(rhor);
  memory->create(rhor,nrhor,nr+1,"pair:rhor");
  
  nrhors = setfl->nelements;
  memory->destroy(rhors);
  memory->create(rhors,nrhors,nr+1,"pair:rhors");

  // copy each element's rhor to global rhor

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nr; m++) rhor[i][m] = setfl->rhor[i][m];

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nr; m++) rhors[i][m] = setfl->rhors[i][m];

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for setfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
	{
	if((map[i]==1)||(map[j]==1)) 
	 {type2rhor[i][j] = 1; }
	else
	 {type2rhor[i][j] = 0 ;}
	}
  
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
	{
		if( ((map[i]==1)&&(map[j]==0)) || ((map[i]==0)&&(map[j]==1)) )
		{type2rhors[i][j] = 1;}
		else
		{type2rhors[i][j] = 0;}
	}

      

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = setfl->nelements * (setfl->nelements+1) / 2;
  memory->destroy(z2r);
  memory->create(z2r,nz2r,nr+1,"pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = setfl->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Pair2BANDL::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;
  rdrhos = 1.0/drhos;
  
  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(frhos_spline);
  memory->destroy(rhors_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline,nfrho,nrho+1,7,"pair:frho");
  memory->create(rhor_spline,nrhor,nr+1,7,"pair:rhor");
  memory->create(frhos_spline,nfrhos,nrhos+1,7,"pair:frhos");
  memory->create(rhors_spline,nrhors,nr+1,7,"pair:rhors");
  memory->create(z2r_spline,nz2r,nr+1,7,"pair:z2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho,drho,frho[i],frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr,dr,rhor[i],rhor_spline[i]);

  for (int i = 0; i < nfrhos; i++)
    interpolate(nrhos,drhos,frhos[i],frhos_spline[i]);

  for (int i = 0; i < nrhors; i++)
    interpolate(nr,dr,rhors[i],rhors_spline[i]);

  for (int i = 0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],z2r_spline[i]);
}

/* ---------------------------------------------------------------------- */

void Pair2BANDL::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                    8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0*(spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void Pair2BANDL::grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fptr);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL," \t\n\r\f"))) list[i++] = atof(ptr);
  }
}

/* ---------------------------------------------------------------------- */

double Pair2BANDL::single(int i, int j, int itype, int jtype,
                       double rsq, double factor_coul, double factor_lj,
                       double &fforce)
{
  int m;
  double r,p,rhoip,rhojp,rhoips,rhojps,z2,z2p,recip,phi,phip,psip;
  double *coeff,*coeffs;

  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);

  coeff = rhor_spline[type2rhor[itype][jtype]][m];
  rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = rhor_spline[type2rhor[jtype][itype]][m];
  rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];

  coeffs = rhors_spline[type2rhors[itype][jtype]][m];
  rhoips = (coeffs[0]*p + coeffs[1])*p + coeffs[2];
  coeffs = rhors_spline[type2rhors[jtype][itype]][m];
  rhojps = (coeffs[0]*p + coeffs[1])*p + coeffs[2];

  coeff = z2r_spline[type2z2r[itype][jtype]][m];
  z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
  z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = fp[i]*rhojp + fp[j]*rhoip + fps[i]*rhojps + fps[j]*rhoips + phip;
  fforce = -psip*recip;

  return phi;
}

/* ---------------------------------------------------------------------- */

int Pair2BANDL::pack_forward_comm(int n, int *list, double *buf,
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

void Pair2BANDL::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {fp[i] = buf[m++];fps[i] = buf[m++];}
}

/* ---------------------------------------------------------------------- */

int Pair2BANDL::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {buf[m++] = rho[i];buf[m++] = rhos[i];}
  return m;
}

/* ---------------------------------------------------------------------- */

void Pair2BANDL::unpack_reverse_comm(int n, int *list, double *buf)
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

double Pair2BANDL::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   swap fp array with one passed in by caller
------------------------------------------------------------------------- */

void Pair2BANDL::swap_eam(double *fp_caller, double **fp_caller_hold)
{
  double *tmp = fp;
  fp = fp_caller;
  *fp_caller_hold = tmp;
}

/* ---------------------------------------------------------------------- */

void *Pair2BANDL::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return NULL;
}
