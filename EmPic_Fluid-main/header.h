#include <math.h>
#include <iostream>
#include <fstream>
#include "threevector.h"
#include "threematrix.h"
#include "stdio.h"
#include "stdlib.h"
#include <vector>
#include "mpi.h"

#ifdef M_PI
  #define PI M_PI
#else
  #define PI 3.14159265358979
#endif

#define PI0 3.14159265359
#define m_e 1.0// 9.10938356e-31
#define mass_p 1.0 //1.6726219e-27
#define eps  1.0   //8.85418782e-12
#define ech  1.0   // 1.60217662e-19
#define cph  1.0 // 1.0
#define kb   1.0   // 1.38064852e-23
#define B_max 1.0  //151.0e-4
#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))
#define  nbuffer 2 

//#include "mpi.h"
using namespace std;

double dt,tmax;
double dt_c,factor_c,factor_FDTD,Delta_r;
double dx,xmin,xmax;
double partweight; // density weight of macroparticles  (originally init_n*Nx*Ny*Nz/Np)

double dy,ymin,ymax;
double dz,zmin,zmax;
int Nx, Ny, Nz; // global (entire domain)
int Nx_local, Ny_local, Nz_local; // local (only proc)
int Nx_nodes_local, Ny_nodes_local, Nz_nodes_local;
int Npart_cell; // initial targed number of particles in cell
int nrel,nmax,nout,nout_print,nemax;
int Npi, Npe;
int Npi_local, Npe_local;
int solver,solver_Maxwell,solver_Boris,bc_x,bc_y,SIR,limiter,variabletype,limitertype;
double gam,epsil,mi_me,vth_i,vth_e;
double init_n,del_n,del_v,vp;
double mass_i,mass_e;
// Added for restart function
int restart;
int RestartFreq;
int initC; // new nt
int naveglobal; // now global headerfile
// subtimestepping
int ncycle;
//int np_emiss = 0;
FILE *fb,*fb1,*fb2,*fb3,*fb4;
FILE *rhof,*phif, *Exf, *Eyf;
char filenam[50], filenam1[50],filenam2[50], filenam3[50],filenam4[50];
char filenamew[50], filenamew1[50], filenamew2[50],filenamew3[50];

double Te_init, Ti_init;
int nrec, number_record;

// For processors
// defined in header_mpi.h 
// =======================
int nproc_x, nproc_y; // how many number of proces in x-y directions
int pid_x, pid_y; // index of local proc in global partitioning
int *ndist_x,  *ndist_y;   // distance (local)
int *ndispl_x, *ndispl_y;  // displacement from (0,0)
int *N_per_proc_x;   // ndlocal's in each proc
int *N_per_proc_y;   // ndlocal's in each proc

double lambda_d, w_pe,w_ce,w_pi,w_ci;
double u_Bohm,cs_ion,cs_electron;
double n0;
double B_0;
double flux_cathode_tot;
vector<double> ener_k,ener_p;


// CELL CENTERS
struct GridPoints{
    double x;
    double y;
    double z;
    double phi;
    double den_i;
    double ux_i;
    double uy_i;
    double uz_i;
    double pxx_i;
    double pxy_i;
    double pyy_i;
    double pxz_i;
    double pyz_i;
    double pzz_i;
    double ex;
    double ey;
    double ez;
    double bx;
    double by;
    double bz;
    double jx;
    double jy;
    double jz;
    double fvpape;
    double fxvx;
};

// NODES
struct FieldPoints{
    double x;
    double y;
    double z;
    double ex;
    double ey;
    double ez;
    double bx; 
    double by;
    double bz;
    double jx;
    double jy;
    double jz;
};


//FLUID VARIABLES
struct GridFluid{
    double n_i;
    double u_ix;
    double u_iy;
    double u_iz; //new
    double m_x_i;
    double m_y_i;
    double m_z_i; //new
    double P_i;
    double e_i;
    double E_i;
    double c_i;
};

struct RKsteps{
    double U_rk_1_0_i;
    double U_rk_2_0_i;
    double U_rk_3_0_i;
    double U_rk_4_0_i;
    double U_rk_5_0_i;

    double U_rk_1_1_i;
    double U_rk_2_1_i;
    double U_rk_3_1_i;
    double U_rk_4_1_i;
    double U_rk_5_1_i;

    double U_rk_1_2_i;
    double U_rk_2_2_i;
    double U_rk_3_2_i;
    double U_rk_4_2_i;
    double U_rk_5_2_i;

    double U_rk_1_3_i;
    double U_rk_2_3_i;
    double U_rk_3_3_i;
    double U_rk_4_3_i;
    double U_rk_5_3_i;
};

struct Coordinate{
    int y;
    int vxe;
    int x;
};

struct Values{
    double x;
    double y;
    double vxe;
};

struct Proc_id{
    int x;
    int vx;
};

struct NODE{
    double x;
    double v;
};

struct CONNECT{
    int sw;
    int nw;
    int ne;
    int se;
};

Proc_id pid;
Coordinate nlocal,nstart,ntot, nlocal_nodes, ntot_nodes;
Values Min, Max, dd;

int p_real,px,pvx;
int nx_global,ny_global,nvxe_global, nx_nodes_global, ny_nodes_global;
double *x;                /* local index*/
double *x_global;         /* global index*/
double *y, *vxe;        /* local index */
int *ind_x, *ind_y, *ind_vxe;  /* convert local to global index */
double *vdf_in_i, *vdf_in_e, *vdf_in_see;  /* global index (VDF for BC) */
double *flux_cathode;
double Ub,Nb,Pt_inflow;
double  neMax;
double CFL=0.8;
double qii = ech;
double qee = -1.0*ech;
int nsteps_output;
int ic_left, ic_right;
double atmp_cathode;

int chooseParticle
(int iflag)
{
    int nv;
    if(iflag==0)
        nv= nlocal.y;
    else
        nv= nlocal.vxe;

    return nv;
}

bool checkEW
(int my_id, int ptgt)
{
   bool tmp;

   if(my_id%px!=ptgt)
     tmp = true;
   else
     tmp = false;

   return tmp;
}

bool checkN
(int my_id, int pn)
{
   bool tmp;

   if(my_id<pn)
     tmp = true;
   else
     tmp = false;

   return tmp;
}

bool checkS
(int my_id, int ps)
{
   bool tmp;

   if(my_id>=ps)
     tmp = true;
   else
     tmp = false;

   return tmp;
}

void checkPID
(int my_id, int p)
{
   if(p!=p_real)
   {
      if(my_id==0) cout << "STOP SIMULATION WRONG P!" <<endl;
      if(my_id==0) cout << "=========== " <<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
   }
}
//==========================================
// Initialize 
//==========================================
void initialize
(int my_id, int p)
{
    
  nmax = (int)(tmax/dt+0.01);
  // ion sub time stepping
  ncycle = 1;
  
  dx = (xmax-xmin)/((double)Nx);
  dy = (ymax-ymin)/((double)Ny);
  dz = (zmax-zmin)/((double)Nz);

  B_0=1.0;
  init_n = 1.0;
  n0 = init_n;

  //Should remove variable (flux_cathode) and also on the HLLC where it is called
  int nlocal_mod, nlocal; 
  flux_cathode  = (double *) malloc(Ny+4*sizeof(int));

  Npi = 10; //Npart_cell*Nx*Ny*Nz;
  Npe = 0.0;//Npart_cell*Nx*Ny*Nz;//Npart_cell * Nx * Ny*Nz;

  partweight =  init_n*Nx*Ny*Nz/Npi;
  double partweight_init = init_n / ((double)Npart_cell);

  
  nout = neMax/dt;
  nrec = (int)(nout);
  naveglobal = 100;  // frequency of sampling global quantities

  mass_e = m_e;
  mass_i = mass_e*mi_me;
  Ti_init = (vth_i*vth_i)*mass_i;
  Te_init = (vth_e*vth_e)*mass_e;
  lambda_d = sqrt(eps*Te_init/(n0*ech));
  w_pe = sqrt((n0*ech*ech)/(eps*mass_e));//PI0;
  w_ce = ((ech*B_0)/mass_e);///PI0;
  w_pi = sqrt((n0*ech*ech)/(eps*mass_i));//PI0;
  w_ci = ((ech*B_0)/mass_i);///PI0;
  u_Bohm = sqrt(ech*gam*(Ti_init+Te_init)/mass_i);
  cs_ion = sqrt(gam*ech*Ti_init/mass_i);
  cs_electron = sqrt(gam*ech*Te_init/mass_e);

   nmax = (int)(tmax/dt+0.01);
   nsteps_output = int(neMax/dt);

    
//    Delta_r = min(dx,dy);
//    Delta_r =min(Delta_r ,dz);
//    dt_c=factor_c*Delta_r ;                    // Courant Limit (normalized)
//    dt_c = 1./sqrt(pow(dx,2)+pow(dx,2)+pow(dx,2));
//    dt=0.5*dt_c;
    factor_c = 1./sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
    factor_FDTD = 0.5*factor_c;

  // 0: relativistic
  // 1: non-relativistic
  nrel = 1;

  // save the input parameters in a file
  //if((my_id==0) && (restart == 1))
  if(my_id==0)
  {
	cout << "... Saving Input Parameters in a file ... " << endl;
	ofstream inputSave;
	char fileS[150];
	sprintf(fileS,"./SimulationsSetuo.dat");
	inputSave.open(fileS);
	inputSave << "=== Summary of Input Parameters ===" << endl;
 	inputSave << " dt = " << dt << endl;
  	inputSave << " tmax = " << tmax << endl;
  	inputSave << " nmax " << nmax << endl;
	inputSave << " Restart = " << restart << endl;
	inputSave << " Nx = " << Nx << endl;
    inputSave << " Ny = " << Ny << endl;
    inputSave << " Nz = " << Nz << endl;
    inputSave << " xmin = " << xmin << endl;
    inputSave << " xmax = " << xmax << endl;
    inputSave << " ymin = " << ymin << endl;
    inputSave << " ymax = " << ymax << endl;
    inputSave << " zmin = " << zmin << endl;
    inputSave << " zmax = " << zmax << endl;
    inputSave << " dx = " << dx << endl;
    inputSave << " dy = " << dy << endl;
    inputSave << " dz = " << dz << endl;
	inputSave << " init_n = " << init_n << endl;
	inputSave << " # part per cell = " << Npart_cell << endl;
	inputSave << " Total # ions = " << Npi <<endl;
    inputSave << " Total # electrons = " << Npe <<endl;
    inputSave << " # of particles per cell (ions) = " << Npi/Nx/Ny/Nz <<endl;
    inputSave << " # of particles per cell (electrons) = " << Npe/Nx/Ny/Nz <<endl;
    inputSave << " particle weight = " << partweight <<endl;
    inputSave << " modified init_n = " << init_n <<endl;
	inputSave << " me = " << mass_e << endl;
    inputSave << " mi = " << mass_i << endl;
    inputSave << " Ti = " << Ti_init << endl;
	inputSave << " Te = " << Te_init << endl;
	inputSave << " Vthe = " << vth_e << endl;
	inputSave << " Vthi = " << vth_i << endl;
    inputSave << " lambda_d = " << lambda_d << endl;
    inputSave << " w_pe = " << w_pe << endl;
    inputSave << " w_ce = " << w_ce << endl;
    inputSave << " w_pi = " << w_pi << endl;
    inputSave << " w_ci = " << w_ci << endl;
	inputSave << " Relativistic = " << nrel << endl;
	inputSave << " nout = " << nout << endl;
	inputSave << " nout print = " << nout_print << endl;
	inputSave << " # records = " << number_record << endl;
	inputSave << " nrec = " << nrec << endl;
	inputSave.close();
  }
}

//==========================================
// Gaussian distribution
//    Maxwellian speed
//==========================================
double gaussian_speed(double vth)
{
    double xx, yy;
    double ranf;

    do
    {
       ranf = ((double) rand()) / ((double)RAND_MAX);
       xx = sqrt(-2.0*log(ranf))*vth; 
    }
    while (ranf<=0.0 || ranf>1.0) ;
    // while (1.0*rand()/RAND_MAX > yy);
    
    return xx;
}

//==========================================
// Cosine distribution
//==========================================
double cosine()
{
    double xx, yy;
    double ranf;

    do
    {
       ranf = ((double) rand() / (RAND_MAX));
       xx = ranf; 
       yy = sin(PI*xx); 
    }
    while ( ((double)rand()/(RAND_MAX)) > yy);
    
    return xx;
}

//==========================================
// Cosine distribution k modes
//==========================================
double cosine_k(int kwave)
{
    double xx, yy;
    double ranf;

    do
    {
       ranf = ((double) rand() / (RAND_MAX));
       xx = ranf; 
       yy = sin(2.0*PI*((double)kwave)*xx); 
    }
    while ( ((double)rand()/(RAND_MAX)) > yy);
    
    return xx;
}


//==========================================
// Gaussian distribution
//    Maxwellian 
//==========================================
double gaussian(double vth)
{
    double gaussian_width=5.0; // +/-2 sigma range
    double xx, yy;
    double ranf;

    do
    {
       ranf = ((double) rand() / (RAND_MAX));
       xx = (2.0*ranf-1.0)*gaussian_width*vth;
       yy = exp(-0.5*xx*xx/vth/vth);
    }
    while ( ((double)rand()/(RAND_MAX)) > yy);
    // while (1.0*rand()/RAND_MAX > yy);
    
    return xx;
}

//==========================================
// Initialize Arrays 
//==========================================
void initializeArray
(GridPoints *** &gc, FieldPoints *** &field)
{
  gc  = (GridPoints***) malloc(Nx*sizeof(GridPoints));
  for(int i=0; i<Nx; i++)
  {
     gc[i]  = (GridPoints**) malloc(Ny*sizeof(GridPoints));
     for(int j=0; j<Ny; j++)
     {
          gc[i][j]  = (GridPoints*) malloc(Nz*sizeof(GridPoints));

      }
  }

  field  = (FieldPoints***) malloc((Nx+1)*sizeof(FieldPoints));
  for(int i=0; i<(Nx+1); i++)
  {
     field[i]  = (FieldPoints**) malloc((Ny+1)*sizeof(FieldPoints));
     for(int j=0; j<(Ny+1); j++)
     {
          field[i][j]  = (FieldPoints*) malloc((Nz+1)*sizeof(FieldPoints));

      }
  }
}
//==========================================
// Initialize fluid  Arrays
//==========================================
void initializeArrayFluid
(GridFluid *** &gf, RKsteps *** &rk, int nxt,int nyt)
{
  gf  = (GridFluid***) malloc(nxt*sizeof(GridFluid));
  for(int i=0; i<nxt; i++)
  {
     gf[i]  = (GridFluid**) malloc(nyt*sizeof(GridFluid));
     for(int j=0; j<nyt; j++)
     {
          gf[i][j]  = (GridFluid*) malloc(1*sizeof(GridFluid));

      }
  }
  rk  = (RKsteps***) malloc(nxt*sizeof(RKsteps));
  for(int i=0; i<nxt; i++)
  {
     rk[i]  = (RKsteps**) malloc(nyt*sizeof(RKsteps));
     for(int j=0; j<nyt; j++)
     {
          rk[i][j]  = (RKsteps*) malloc(1*sizeof(RKsteps));

      }
  }
    
}
