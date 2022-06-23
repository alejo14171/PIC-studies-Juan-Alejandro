#include "header.h"
#include "header_mpi.h"
#include "createFolders.h"
#include "setup.h"
#include "definition.h"
#include "setupGrid.h"
#include "functions_PIC.h"
#include "functions_FDTD.h"
#include "functions_HLLC.h"
#include "output.h"
#include <iostream>
#include <fstream>

int main
(int argc, char *argv[])
{
  int my_id, p;
  int ierr=0;
  int istop_final; // error message from hypre solver
  int ilower[2],iupper[2]; // required for hypre domain

  int nt=0;
  double t;
    

  GridPoints ***gc = NULL;     // Eulerian
  FieldPoints ***field = NULL; // Eulerian
  FieldPoints ***field_copy = NULL; // Eulerian
  vector<ThreeVec> xi, pi;     // particles (Ions)
  vector<double>   qi;         // particles (Ions)
  vector<double>   Vpar;         // particles (Ions)
  vector<double>   Vperp;         // particles (Ions)
  vector<ThreeVec> x_old;
  vector<ThreeVec> x_new;
  vector<ThreeVec> E_p;        // Efield on particles
  vector<ThreeVec> B_p;        // Bfield on particles
        
  int nparts_i_local;
  int nparts_i_total;
   // averaging global parameters (Boeuf)
   // ====================================
   // (xisize,xesize,anode_i,anode_e,cathode_i,cathode_e)
   int numlocal[6];   // local
   int numglobal[6];  // global
    
    //domain decomposition information
    bool notEast, notWest;
    int peast,pwest;
    bool notNorth, notSouth;
    int  psouth,pnorth;

  // ofstream outtime,outspace,outrho,outne,outni,outey;
  ofstream outt; // (for 1D azimuthal)
  ofstream outte, outfe, outfi;
  ofstream outtec; // for 2D defined in output.h (each output = nt)
  ofstream outtecfluid;
  ofstream outtecfield; // for field 
  ofstream outglobal; // for boeuf
 

  // start MPI
  // =========
  ierr = MPI_Init(&argc,&argv);
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &p);
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &my_id);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("%d my_id, %d procs, %d error message\n",my_id,p,ierr);
  // Root should have MPI_COMM_WORLD
  MPI_Bcast(&p,1,MPI_INT,0,MPI_COMM_WORLD);
  if(my_id == 0 ) printf("\n\n main # procs = %d \n\n",p);
  setProperties();
  // create folders
  // ==============
  if(my_id==0) createFolders();
  MPI_Barrier(MPI_COMM_WORLD);
  if(my_id == 0 ) printf("create folder complete \n");

  srand(time(NULL)+my_id);

  // ============
  // Initialize
  // ============
  initialize(my_id, p);

  // check processor information
  // ===========================
  checkProcs(my_id, p);
  assignProcs(my_id, p);
  distributeProcs(MPI_COMM_WORLD,my_id, p);
  getCellSize(my_id,p);
  getCoordinate(my_id);

  if(my_id == 0 ) printf("proc info complete \n");
  MPI_Barrier(MPI_COMM_WORLD);
    
  // Initialize (continue)
  // =====================
  initializeArray(gc,field);
  initializeArray(gc,field_copy);
  MPI_Barrier(MPI_COMM_WORLD);

    
  if(my_id==0) cout << " === Initialize Done (header.h, header_mpi.h) " <<endl;
    
  //FIELDS & PARTICLE SETUP
  // Setup arrays (restart included here)
  // ====================================
    setup(my_id,p,gc,field,field_copy,xi,pi,x_old,x_new,E_p,B_p,qi,Vpar,Vperp);
    //backward
    leapfrog(xi,pi,x_old,x_new,qi,0);
    //forward
    leapfrog(xi,pi,x_old,x_new,qi,1);
    //moments_VDF(gc,xi,pi,qi,0,1,1);
    //current_density_zigzag(gc,x_old,x_new);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_id==0) cout << " === Setup Done (setup.h) " <<endl;

    if(my_id==0)
    {
     cout << endl << " Before Iteration " <<endl;
     cout << " nmax: " << nmax << endl;
     cout << " nrec: " << nrec << endl << endl; 
    }

    MPI_Barrier(MPI_COMM_WORLD);
  // =====================
  // Global quantities
  // =====================
    if(my_id==0)
    {
        outt.open("./global.dat"); //
    }
 
    int nttmp0 = initC;

    int idt_final = 10;
    double _timer_start, _timer_check;
    double _timer[idt_final];
    for(int idt=0;idt<idt_final;++idt)
       _timer[idt] = 0.0;

   _timer_start = MPI_Wtime();
   if(restart == 0) 
   {
	t = 0;
   }
   else
   { 
	t = initC*dt; 
   }

  if(my_id==0)
  {
     cout << endl << "  Iteration START " <<endl << endl;
     cout << " initC is " << initC << " (if restart not zero); t is " << t << endl <<endl;
  } 

  // ==========
  // Iteration
  // ==========
  for(nt=initC;nt<=initC+nmax;++nt) // added nmax for output
  {
      //Start the loop checking the number of particles
      nparts_i_local = (int)xi.size();
      MPI_Reduce(&nparts_i_local,&nparts_i_total,1,MPI_INT,MPI_SUM,0, MPI_COMM_WORLD);
      if(my_id==0) cout << "Time passed:" << dt*nt << " Timestep:"<< nt << " ion count: "<< nparts_i_total << endl ;
      //output ot global quantities
      Temporal_diagnostics(my_id,p,nt,dt*nt,outt,xi,pi,qi,gc);
      //output time
      if(nt%nout==0)
      {
         //initial conditions
         if(nt==0)
         {
             if(nt%nrec==0)
             {
                if(my_id==0) _timer_check = MPI_Wtime();
                if(my_id==0) cout << " entering particle writing ... " <<endl;
                write_particles(my_id,p,nt,xi,pi,qi,E_p,B_p,Vpar,Vperp,0);
                if(my_id==0) cout << " wrote particle  ... [Time: " << MPI_Wtime() - _timer_check << " s]" << endl;
             }
             // velocity move (half timestep): v^{n}-> v^{n+1/2}
             boris(field,xi,pi,E_p,B_p,0,0,0.5);
         }
         else
         {
             if(nt%nrec==0)
             {
                 //velocity move (half timestep): v^{n+1/2}-> v^{n}
                 boris(field,xi,pi,E_p,B_p,0,0,-0.5);
                 if(my_id==0) _timer_check = MPI_Wtime();
                 if(my_id==0) cout << " entering particle writing ... " <<endl;
                 write_particles(my_id,p,nt,xi,pi,qi,E_p,B_p,Vpar,Vperp,0);
                 if(my_id==0) cout << " wrote particle  ... [Time: " << MPI_Wtime() - _timer_check << " s]" << endl;
                 //velocity move (half timestep): v^{n}-> v^{n+1/2}
                 boris(field,xi,pi,E_p,B_p,0,0,0.5);
             }
         }
      }
     else
     {
             // velocity update (one timestep): v^{n+1/2}->v^{n+3/2}
             // ===================================================
             boris(field,xi,pi,E_p,B_p,0,0,1.0);
     }
    // Advance position: r^{n} --> r^{n+1}
     // ================
      leapfrog(xi,pi,x_old,x_new,qi,1);
      //moments_VDF(gc,xi,pi,qi,0,1,1);
      //current_density_zigzag(gc,x_old,x_new);
      //check conservation of number of particles
      if((int)xi.size() != (int)pi.size() || (int)xi.size() != (int)qi.size())
         cout << "ion size error" <<endl;
     t+=dt;
  }
  // ===============
  // After iteration
  // ===============
    outfi.close();
    
  if(my_id==0) 
  {
      cout << nmax << " " <<nt << endl;
  }
   MPI_Finalize();
   return 0;
}

