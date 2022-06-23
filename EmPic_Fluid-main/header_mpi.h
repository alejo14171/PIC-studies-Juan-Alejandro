
// ====================================================================
//
//    check the number of processors as designed
//    
//     RULE:
//     p = 1,2,3 OKAY
//               = 4,6,8,10,12,14 OKAY
//               = 16,20,24,28,32,36,40,44,48,52,56,60 OKAY
//               = 64,72,80,88,96,104,112,120 OKAY
//               = otherwise NG
//
// ====================================================================

void checkProcs
(int my_id, int p)
{

  int check = 0;

  // Currently up to 64 procs (designed)
  if(p > 4 && p%2!=0)  check = 1;  
  if(p > 16 && p%4!=0) check = 1;  
  if(p > 64 && p%8!=0) check = 1;  

  if(check==1)
  {
     if(my_id==0) printf(" STOP simulation wrong p (must be even number \n");
     MPI_Abort(MPI_COMM_WORLD,1);
  } 
}

// ====================================================================
//
//   ASSIGN processors 
//      - # of procs in x & y (nproc_x & nproc_y)
//      - Global location of proc: pid
//
// ====================================================================
void assignProcs
(int my_id, int p)
{
  // Supported up to 256 procs
  //if(p<4)
  //{
    nproc_y = 1;
    nproc_x = p;
  //}
  /*else if(p<16)
  {
    nproc_y = 2;
    nproc_x = p/2;
  }
  else if(p<64)
  {
    nproc_y = 4;
    nproc_x = p/4;
  }
  else if(p<256)
  {
    nproc_y = 8;
    nproc_x = p/8;
  }

  // set at least 2 points
  if(Ny / nproc_y < 3)
  {
     --nproc_y;
     nproc_x= p/nproc_y;
  }

  if(Nx / nproc_x < 3)
  {
     if(my_id==0) printf(" STOP: make sure there are at least 3 cells in each procs for x and y \n");  
     if(my_id==0) printf(" change (1) p, (2) ndglobal.x, or (3) ndglobal.y \n"); 
     MPI_Abort(MPI_COMM_WORLD,1);
  }
 */
}

// ====================================================================
//
//   DISTRIBUTE  processors 
//    - based on hypre assigned hypre_p
//
// ====================================================================
void distributeProcs
(MPI_Comm hypre_comm, 
 int hypre_id, int hypre_p)
{
   int i,j;
   int k=0;

   if(MPI_COMM_NULL!=hypre_comm)
   {
      for(j=0; j<nproc_y; j++)
      {
         for(i=0; i<nproc_x; i++)
         {
             if(hypre_id==k)
             {
                pid_x  = i; 
                pid_y  = j; 
             }
             ++k;
         }
      }
   }
   else
   {
       pid_x  = -1;
       pid_y  = -1;
   }
}
