#include "Alloc.h"
void EMfield_halftime
(int my_id, int p,GridPoints ***xcg,FieldPoints ***field, FieldPoints ***field_tmp,int step, int solver_Maxwell)
{
    double sig_res;
     double dt_tmp;
     double **gz_ly  = newArr(double,Nx+1,Ny+1);
     double **gz_ry  = newArr(double,Nx+1,Ny+1);
     double **gz_lx  = newArr(double,Nx+1,Ny+1);
     double **gz_rx  = newArr(double,Nx+1,Ny+1);
     double **fx_ly  = newArr(double,Nx+1,Ny+1);
     double **fx_ry  = newArr(double,Nx+1,Ny+1);
     double **fy_lx  = newArr(double,Nx+1,Ny+1);
     double **fy_rx  = newArr(double,Nx+1,Ny+1);

     double **gz_t  = newArr(double,Nx+1,Ny+1);
     double **gz_r  = newArr(double,Nx+1,Ny+1);
     double **fy_r  = newArr(double,Nx+1,Ny+1);
     double **fx_t  = newArr(double,Nx+1,Ny+1);
    
    for(int iy=0; iy<=Ny; iy++)
     {
         for(int ix=0; ix<=Nx; ix++)
         {
             gz_ly[ix][iy]  =  0.0;
             gz_ry[ix][iy]  =  0.0;
             gz_lx[ix][iy]  =  0.0;
             gz_rx[ix][iy]  =  0.0;
             fx_ly[ix][iy]  =  0.0;
             fx_ry[ix][iy]  =  0.0;
             fy_lx[ix][iy]  =  0.0;
             fy_rx[ix][iy]  =  0.0;
             gz_t[ix][iy]  =  0.0;
             gz_r[ix][iy]  =  0.0;
             fy_r[ix][iy]  =  0.0;
             fx_t[ix][iy]  =  0.0;
         }
     }
    
    
     //step = 0 backward half time step for initial condition
    if(step==0)
    {
        dt_tmp = -0.5*dt;
    }
    else
    {
        dt_tmp = 0.5*dt;
    }
    
    if(solver_Maxwell == 0)
    {

        //Evolve Bx,By & Bz----> dt
           for(int iy=1; iy<=Ny; iy++)
            {
                for(int ix=0; ix<=Nx; ix++)
                {
                        field[ix][iy][0].bx  = field[ix][iy][0].bx - (((dt_tmp)/dy)*(field[ix][iy][0].ez - field[ix][iy-1][0].ez));
                }
            }

            //Bondary condition Bx: Periodic BC
            for(int ix=0;ix<=Nx;ix++)
            {
                field[ix][0][0].bx = field[ix][Ny][0].bx;
            }

            for(int iy=0; iy<=Ny; iy++)
            {
                for(int ix=1; ix<=Nx; ix++)
                {
                        field[ix][iy][0].by  = field[ix][iy][0].by + (((dt_tmp)/dx)*(field[ix][iy][0].ez - field[ix-1][iy][0].ez));
                }
            }
            
            //Bondary condition By: Periodic BC
            for(int iy=0;iy<=Ny;iy++)
            {
                field[0][iy][0].by = field[Nx][iy][0].by;
            }
            
            for(int iy=1; iy<=Ny; iy++)
            {
                for(int ix=1; ix<=Nx; ix++)
                {
                    field[ix][iy][0].bz  = field[ix][iy][0].bz - ( (((dt_tmp)/dx) *(field[ix][iy][0].ey - field[ix-1][iy][0].ey)) -  (((dt_tmp)/dy) *(field[ix][iy][0].ex - field[ix][iy-1][0].ex)) );
                }
            }
          
             //Bondary condition Bz: Periodic BC
               for(int ix=1;ix<=Nx;ix++)
                {
                  field[ix][0][0].bz = field[ix][Ny][0].bz ;
                }
                for(int iy=1;iy<=Ny;iy++)
                {
                  field[0][iy][0].bz = field[Nx][iy][0].bz ;
                }
                //corner
                field[0][0][0].bz = field[Nx][Ny][0].bz ;
    }
    else if(solver_Maxwell == 1)
    {
        
        for(int iy=1; iy<Ny; iy++)
         {
             for(int ix=1; ix<Nx; ix++)
             {
                 gz_lx[ix][iy]  =  xcg[ix-1][iy][0].ez;
                 gz_rx[ix][iy]  =  xcg[ix][iy][0].ez;
                 gz_ly[ix][iy]  =  xcg[ix][iy-1][0].ez;
                 gz_ry[ix][iy]  =  xcg[ix][iy][0].ez;
                 //Bx on x cell-interface
                 fx_ly[ix][iy]  =  xcg[ix][iy-1][0].bx;
                 fx_ry[ix][iy]  =  xcg[ix][iy][0].bx;
                 //By on y cell-interface
                 fy_lx[ix][iy]  =  xcg[ix-1][iy][0].by;
                 fy_rx[ix][iy]  =  xcg[ix][iy][0].by;
             }
         }
        //first column and first row
        for(int iy=1; iy<Ny; iy++)
        {
            gz_lx[0][iy]  =  xcg[Nx][iy][0].ez;
            gz_rx[0][iy]  =  xcg[0][iy][0].ez;
            gz_ly[0][iy]  =  xcg[0][iy-1][0].ez;
            gz_ry[0][iy]  =  xcg[0][iy][0].ez;
            //Bx on x cell-interface
            fx_ly[0][iy]  =  xcg[0][iy-1][0].bx;
            fx_ry[0][iy]  =  xcg[0][iy][0].bx;
            //By on y cell-interface
            fy_lx[0][iy]  =  xcg[0][iy][0].by;
            fy_rx[0][iy]  =  xcg[0][iy][0].by;
        }
        for(int ix=1; ix<Nx; ix++)
        {
            gz_lx[ix][0]  =  xcg[ix-1][0][0].ez;
            gz_rx[ix][0]  =  xcg[ix][0][0].ez;
            gz_ly[ix][0]  =  xcg[ix][Ny][0].ez;
            gz_ry[ix][0]  =  xcg[ix][0][0].ez;
            //Bx on x cell-interface
            fx_ly[ix][0]  =  xcg[ix][Ny][0].bx;
            fx_ry[ix][0]  =  xcg[ix][0][0].bx;
            //By on y cell-interface
            fy_lx[ix][0]  =  xcg[ix-1][0][0].by;
            fy_rx[ix][0]  =  xcg[ix][0][0].by;
        }
        //nodes
        gz_lx[0][0]  =  xcg[Nx][0][0].ez;
        gz_rx[0][0]  =  xcg[0][0][0].ez;
        gz_ly[0][0]  =  xcg[0][Ny][0].ez;
        gz_ry[0][0]  =  xcg[0][0][0].ez;
        //Bx on x cell-interface
        fx_ly[0][0]  =  xcg[0][Ny][0].bx;
        fx_ry[0][0]  =  xcg[0][0][0].bx;
        //By on y cell-interface
        fy_lx[0][0]  =  xcg[Nx][0][0].by;
        fy_rx[0][0]  =  xcg[0][0][0].by;
        
        //( periodic boundary conditions)
        for(int ix=0; ix<=Nx; ix++)
        {
            gz_lx[ix][Ny]  =  gz_lx[ix][0];
            gz_rx[ix][Ny]  =  gz_rx[ix][0];
            gz_ly[ix][Ny]  =  gz_ly[ix][0];
            gz_ry[ix][Ny]  =  gz_ry[ix][0];
            //Bx on x cell-interface
            fx_ly[ix][Ny]  =  fx_ly[ix][0] ;
            fx_ry[ix][Ny]  =  fx_ry[ix][0];
            //By on y cell-interface
            fy_lx[ix][Ny]  =  fy_lx[ix][0];
            fy_rx[ix][Ny]  =  fy_rx[ix][0];
        }
        for(int iy=0; iy<=Ny; iy++)
        {
            gz_lx[Nx][iy]  =  gz_lx[0][iy];
            gz_rx[Nx][iy]  =  gz_rx[0][iy];
            gz_ly[Nx][iy]  =  gz_ly[0][iy];
            gz_ry[Nx][iy]  =  gz_ry[0][iy];
            //Bx on x cell-interface
            fx_ly[Nx][iy]  =  fx_ly[0][iy] ;
            fx_ry[Nx][iy]  =  fx_ry[0][iy];
            //By on y cell-interface
            fy_lx[Nx][iy]  =  fy_lx[0][iy];
            fy_rx[Nx][iy]  =  fy_rx[0][iy];
        }
        
        for(int iy=0; iy<Ny+1; iy++)
         {
             for(int ix=0; ix<Nx+1; ix++)
             {
                 sig_res = copysign(1.0, 1.0); //last 1.0 is c
                 gz_t[ix][iy]  = ((0.5*(gz_ly[ix][iy] +  gz_ry[ix][iy])) - (sig_res*(fx_ly[ix][iy] -  fx_ry[ix][iy])));
                 gz_r[ix][iy]  =  ((0.5*(gz_lx[ix][iy] +  gz_rx[ix][iy])) + (sig_res*(fy_lx[ix][iy] -  fy_rx[ix][iy])));
                 fy_r[ix][iy]  =  ((0.5*(fy_lx[ix][iy] +  fy_rx[ix][iy])) + (sig_res*(gz_lx[ix][iy] -  gz_rx[ix][iy])));
                 fx_t[ix][iy]  =  ((0.5*(fx_ly[ix][iy] +  fx_ry[ix][iy])) - (sig_res*(gz_ly[ix][iy] -  gz_ry[ix][iy])));
             }
         }
        //gz on nodes
        for(int iy=1; iy<=Ny; iy++)
         {
             for(int ix=1; ix<=Nx; ix++)
             {
                 sig_res = copysign(1.0, 1.0); //last 1.0 is c
                 field[ix][iy][0].ez =  0.25*(gz_r[ix][iy-1] + gz_r[ix][iy] + gz_t[ix-1][iy] + gz_t[ix][iy] )  ;
             }
         }
        //node
        field[Nx][0][0].ez = 0.25*(gz_r[0][Ny] + gz_r[0][0] + gz_t[Nx][0] + gz_t[0][0] );
        
        //periodic boundary for bz on nodes
        for(int iy=0; iy<=Ny; iy++)
        {
            field[0][iy][0].ez = field[Nx][iy][0].ez;
        }
        for(int ix=0; ix<=Nx; ix++)
        {
            field[ix][0][0].ez = field[ix][Ny][0].ez;
        }
        
        //Evolve Bx,By & Bz----> dt
           for(int iy=1; iy<=Ny; iy++)
            {
                for(int ix=0; ix<=Nx; ix++)
                {
                        field[ix][iy][0].bx  = field[ix][iy][0].bx - (((dt_tmp)/dy)*(field[ix][iy][0].ez - field[ix][iy-1][0].ez));
                }
            }
            //Bondary condition Bx: Periodic BC
            for(int ix=0;ix<=Nx;ix++)
            {
                field[ix][0][0].bx = field[ix][Ny][0].bx;
            }
            for(int iy=0; iy<=Ny; iy++)
            {
                for(int ix=1; ix<=Nx; ix++)
                {
                        field[ix][iy][0].by  = field[ix][iy][0].by + (((dt_tmp)/dx)*(field[ix][iy][0].ez - field[ix-1][iy][0].ez));
                }
            }
            //Bondary condition By: Periodic BC
            for(int iy=0;iy<=Ny;iy++)
            {
                field[0][iy][0].by = field[Nx][iy][0].by;
            }
            
            for(int iy=1; iy<=Ny; iy++)
            {
                for(int ix=1; ix<=Nx; ix++)
                {
                    field[ix][iy][0].bz  = field[ix][iy][0].bz - ( (((dt_tmp)/dx) *(field[ix][iy][0].ey - field[ix-1][iy][0].ey)) -  (((dt_tmp)/dy) *(field[ix][iy][0].ex - field[ix][iy-1][0].ex)) );
                }
            }
          
             //Bondary condition Bz: Periodic BC
               for(int ix=1;ix<=Nx;ix++)
                {
                  field[ix][0][0].bz = field[ix][Ny][0].bz ;
                }
                for(int iy=1;iy<=Ny;iy++)
                {
                  field[0][iy][0].bz = field[Nx][iy][0].bz ;
                }
                //corner
                field[0][0][0].bz = field[Nx][Ny][0].bz ;
        
    }
    else
    {
         cout << "wrong Maxwell solver" <<endl;
    }
    
    delArr(gz_ly,Nx+1);
    delArr(gz_ry,Nx+1);
    delArr(gz_lx,Nx+1);
    delArr(gz_rx,Nx+1);
    delArr(fx_ly,Nx+1);
    delArr(fx_ry,Nx+1);
    delArr(fy_lx,Nx+1);
    delArr(fy_rx,Nx+1);
    delArr(gz_t,Nx+1);
    delArr(gz_r,Nx+1);
    delArr(fy_r,Nx+1);
    delArr(fx_t,Nx+1);
    
}

//magentic field
void EMfield
(int my_id, int p,GridPoints ***xcg,
 FieldPoints ***field,FieldPoints ***field_tmp,int solver_Maxwell)
{
    if(solver_Maxwell == 0)
    {
    
    //Evolve Bx,By & Bz----> dt
        for(int iy=1; iy<=Ny; iy++)
         {
             for(int ix=0; ix<=Nx; ix++)
             {
                     field[ix][iy][0].bx  = field[ix][iy][0].bx - (((dt)/dy)*(field[ix][iy][0].ez - field[ix][iy-1][0].ez));
             }
         }

         //Bondary condition Bx: Periodic BC
         for(int ix=0;ix<=Nx;ix++)
         {
             field[ix][0][0].bx = field[ix][Ny][0].bx;
         }

         for(int iy=0; iy<=Ny; iy++)
         {
             for(int ix=1; ix<=Nx; ix++)
             {
                     field[ix][iy][0].by  = field[ix][iy][0].by + (((dt)/dx)*(field[ix][iy][0].ez - field[ix-1][iy][0].ez));
             }
         }
         
         //Bondary condition By: Periodic BC
         for(int iy=0;iy<=Ny;iy++)
         {
             field[0][iy][0].by = field[Nx][iy][0].by;
         }
         
         for(int iy=1; iy<=Ny; iy++)
         {
             for(int ix=1; ix<=Nx; ix++)
             {
                 field[ix][iy][0].bz  = field[ix][iy][0].bz - ( (((dt)/dx) *(field[ix][iy][0].ey - field[ix-1][iy][0].ey)) -  (((dt)/dy) *(field[ix][iy][0].ex - field[ix][iy-1][0].ex)) );
             }
         }
       
          //Bondary condition Bz: Periodic BC
            for(int ix=1;ix<=Nx;ix++)
             {
               field[ix][0][0].bz = field[ix][Ny][0].bz ;
             }
             for(int iy=1;iy<=Ny;iy++)
             {
               field[0][iy][0].bz = field[Nx][iy][0].bz ;
             }
             //corner
             field[0][0][0].bz = field[Nx][Ny][0].bz ;
    
   
        //Evolve Ex,Ey & Ez----> dt
        for(int iy=0; iy<Ny; iy++)
         {
             for(int ix=1; ix<Nx+1; ix++)
             {
                     field[ix][iy][0].ex  = field[ix][iy][0].ex + (((dt)/dy)*(field[ix][iy+1][0].bz - field[ix][iy][0].bz));
             }
         }

         //Bondary condition Ex: Periodic BC
         for(int ix=0;ix<Nx+1;ix++)
         {
             field[ix][Ny][0].ex = field[ix][0][0].ex;
         }
          for(int iy=0;iy<Ny+1;iy++)
          {
              field[0][iy][0].ex = field[Nx][iy][0].ex;
          }
    

         for(int iy=1; iy<Ny+1; iy++)
         {
             for(int ix=0; ix<Nx; ix++)
             {
                     field[ix][iy][0].ey  = field[ix][iy][0].ey - (((dt)/dx)*(field[ix+1][iy][0].bz - field[ix][iy][0].bz));
             }
         }
         
         //Bondary condition Ey: Periodic BC
         for(int iy=0;iy<Ny+1;iy++)
         {
             field[Nx][iy][0].ey = field[0][iy][0].ey;
         }
        for(int ix=0;ix<Nx+1;ix++)
        {
             field[ix][0][0].ey = field[ix][Ny][0].ey;
        }
         
         
         for(int iy=0; iy<Ny; iy++)
         {
             for(int ix=0; ix<Nx; ix++)
             {
                     field[ix][iy][0].ez  = field[ix][iy][0].ez + ( (((dt)/dx)*( field[ix+1][iy][0].by - field[ix][iy][0].by )) - (((dt)/dy)*( field[ix][iy+1][0].bx - field[ix][iy][0].bx )) );
             }
         }
       
          //Bondary condition Ez: Periodic BC
            for(int ix=0;ix<=Nx;ix++)
             {
               field[ix][Ny][0].ez = field[ix][0][0].ez ;
             }
             for(int iy=0;iy<=Ny;iy++)
             {
               field[Nx][iy][0].ez = field[0][iy][0].bz ;
             }
    }
    else if(solver_Maxwell == 1)
    {
        
    }
    else
    {
         cout << "wrong Maxwell solver" <<endl;
    }
    
}
void nodes_to_cc
(int my_id, int p,
 FieldPoints ***field,GridPoints ***xcg)

{
    //Bx,Ex
    for(int iy=0; iy<Ny; iy++)
    {
        for(int ix=0; ix<Nx; ix++)
        {
            xcg[ix][iy][0].bx = 0.25*(field[ix][iy][0].bx + field[ix+1][iy][0].bx+ field[ix][iy+1][0].bx + field[ix+1][iy+1][0].bx);
            xcg[ix][iy][0].ex = 0.25*(field[ix][iy][0].ex + field[ix+1][iy][0].ex + field[ix][iy+1][0].ex + field[ix+1][iy+1][0].ex);
            
            xcg[ix][iy][0].by = 0.25*(field[ix][iy][0].by + field[ix+1][iy][0].by + field[ix][iy+1][0].by + field[ix+1][iy+1][0].by);
            xcg[ix][iy][0].ey = 0.25*(field[ix][iy][0].ey + field[ix+1][iy][0].ey + field[ix][iy+1][0].ey + field[ix+1][iy+1][0].ey);
            
            xcg[ix][iy][0].bz = field[ix][iy][0].bz;
            xcg[ix][iy][0].ez = 0.25*(field[ix][iy][0].ez + field[ix+1][iy][0].ez + field[ix][iy+1][0].ez + field[ix+1][iy+1][0].ez);
	}
    }
    
}

void cc_to_nodes
(int my_id, int p,
 GridPoints ***xgc,FieldPoints ***field,FieldPoints ***field_copy)
{
   //From cc to nodes: B & E
    for(int iy=1;iy<Ny;iy++)
    {
     for(int ix=1;ix<Nx;ix++)
      {
          field[ix][iy][0].ex = 0.25*(xgc[ix-1][iy-1][0].ex + xgc[ix][iy-1][0].ex + xgc[ix-1][iy][0].ex + xgc[ix][iy][0].ex);
          
          field[ix][iy][0].bx = 0.25*(xgc[ix-1][iy-1][0].bx + xgc[ix][iy-1][0].bx + xgc[ix-1][iy][0].bx + xgc[ix][iy][0].bx);
          
          field[ix][iy][0].ey = 0.25*(xgc[ix-1][iy-1][0].ey + xgc[ix][iy-1][0].ey + xgc[ix-1][iy][0].ey + xgc[ix][iy][0].ey);
          
          field[ix][iy][0].by = 0.25*(xgc[ix-1][iy-1][0].by + xgc[ix][iy-1][0].by + xgc[ix-1][iy][0].by + xgc[ix][iy][0].by);
          
          field[ix][iy][0].ez = 0.25*(xgc[ix-1][iy-1][0].ez + xgc[ix][iy-1][0].ez + xgc[ix-1][iy][0].ez + xgc[ix][iy][0].ez);
      }
    }
    
    //for Ez 2D setup
    for(int iy=1;iy<Ny;iy++)
    {
     for(int ix=1;ix<Nx;ix++)
      {
            field[ix][iy][0].bz = xgc[ix-1][iy-1][0].bz;
      }
    }
    

    //boundary condition: periodic
    for(int iy=1;iy<Ny;iy++)
    {
            field[0][iy][0].ex = 0.25*(xgc[0][iy-1][0].ex + xgc[0][iy][0].ex + xgc[Nx-1][iy][0].ex + xgc[Nx-1][iy-1][0].ex);
            field[0][iy][0].bx = 0.25*(xgc[0][iy-1][0].bx + xgc[0][iy][0].bx + xgc[Nx-1][iy][0].bx + xgc[Nx-1][iy-1][0].bx);
        
            field[0][iy][0].ey = 0.25*(xgc[0][iy-1][0].ey + xgc[0][iy][0].ey + xgc[Nx-1][iy][0].ey + xgc[Nx-1][iy-1][0].ey);
            field[0][iy][0].by = 0.25*(xgc[0][iy-1][0].by + xgc[0][iy][0].by + xgc[Nx-1][iy][0].by + xgc[Nx-1][iy-1][0].by);
        
            field[0][iy][0].ez = 0.25*(xgc[0][iy-1][0].ez + xgc[0][iy][0].ez + xgc[Nx-1][iy][0].ez + xgc[Nx-1][iy-1][0].ez);
    
    }
    
    for(int ix=1;ix<Nx;ix++)
    {
        field[ix][0][0].ex = 0.25*(xgc[ix][0][0].ex + xgc[ix-1][0][0].ex + xgc[ix][Ny-1][0].ex + xgc[ix-1][Ny-1][0].ex);
        
        field[ix][0][0].bx = 0.25*(xgc[ix][0][0].bx + xgc[ix-1][0][0].bx + xgc[ix][Ny-1][0].bx + xgc[ix-1][Ny-1][0].bx);
        
        field[ix][0][0].ey = 0.25*(xgc[ix][0][0].ey + xgc[ix-1][0][0].ey + xgc[ix][Ny-1][0].ey + xgc[ix-1][Ny-1][0].ey);
        
        field[ix][0][0].by = 0.25*(xgc[ix][0][0].by + xgc[ix-1][0][0].by + xgc[ix][Ny-1][0].by + xgc[ix-1][Ny-1][0].by);
        
        field[ix][0][0].ez = 0.25*(xgc[ix][0][0].ez + xgc[ix-1][0][0].ez + xgc[ix][Ny-1][0].ez + xgc[ix-1][Ny-1][0].ez);
    }
    
    //corners
    field[0][0][0].ex = 0.25*(xgc[0][0][0].ex + xgc[0][Ny-1][0].ex + xgc[Nx-1][0][0].ex + xgc[Nx-1][Ny-1][0].ex);
    
    field[0][0][0].bx = 0.25*(xgc[0][0][0].bx + xgc[0][Ny-1][0].bx + xgc[Nx-1][0][0].bx + xgc[Nx-1][Ny-1][0].bx);
    
    field[0][0][0].ey = 0.25*(xgc[0][0][0].ey + xgc[0][Ny-1][0].ey + xgc[Nx-1][0][0].ey + xgc[Nx-1][Ny-1][0].ey);
    
    field[0][0][0].by = 0.25*(xgc[0][0][0].by + xgc[0][Ny-1][0].by + xgc[Nx-1][0][0].by + xgc[Nx-1][Ny-1][0].by);
    
    field[0][0][0].ez = 0.25*(xgc[0][0][0].ez + xgc[0][Ny-1][0].ez + xgc[Nx-1][0][0].ez + xgc[Nx-1][Ny-1][0].ez);

    
    for(int ix=0;ix<Nx;ix++)
    {
        field[ix][Ny][0].ex = field[ix][0][0].ex;
        field[ix][Ny][0].bx = field[ix][0][0].bx;
        field[ix][Ny][0].ey = field[ix][0][0].ey;
        field[ix][Ny][0].by = field[ix][0][0].by;
        field[ix][Ny][0].ez = field[ix][0][0].ez;
    }
    
    for(int iy=0;iy<Ny+1;iy++)
    {
        field[Nx][iy][0].ex =  field[0][iy][0].ex;
        field[Nx][iy][0].bx =  field[0][iy][0].bx;
        field[Nx][iy][0].ey =  field[0][iy][0].ey;
        field[Nx][iy][0].by =  field[0][iy][0].by;
        field[Nx][iy][0].ez =  field[0][iy][0].ez;
    }
    
    
       for(int ix=1;ix<=Nx;ix++)
       {
         field[ix][0][0].bz = field[ix][Ny][0].bz ;
       }
       for(int iy=1;iy<=Ny;iy++)
       {
         field[0][iy][0].bz = field[Nx][iy][0].bz ;
       }
       //corner
       field[0][0][0].bz = field[Nx][Ny][0].bz ;
    
    
       for(int iy=0;iy<=Ny;iy++)
       {
        for(int ix=0;ix<=Nx;ix++)
         {
             field_copy[ix][iy][0].bx = field[ix][iy][0].bx;
             field_copy[ix][iy][0].by = field[ix][iy][0].by;
             field_copy[ix][iy][0].bz = field[ix][iy][0].bz;
             field_copy[ix][iy][0].ex = field[ix][iy][0].ex;
             field_copy[ix][iy][0].ey = field[ix][iy][0].ey;
             field_copy[ix][iy][0].ez = field[ix][iy][0].ez;
         }
       }
}
     

void nodes_to_interc
(int my_id, int p,FieldPoints ***field, FieldPoints ***field_copy, int solver_Maxwell)
{
    //FDTF method (different grid for B and E)
    if(solver_Maxwell == 0)
    {
    for(int iy=0;iy<=Ny;iy++)
    {
     for(int ix=0;ix<=Nx;ix++)
      {
          field_copy[ix][iy][0].bx = field[ix][iy][0].bx;
          field_copy[ix][iy][0].by = field[ix][iy][0].by;
          field_copy[ix][iy][0].bz = field[ix][iy][0].bz;
          field_copy[ix][iy][0].ex = field[ix][iy][0].ex;
          field_copy[ix][iy][0].ey = field[ix][iy][0].ey;
          field_copy[ix][iy][0].ez = field[ix][iy][0].ez;
      }
    }
    
    //yellow points
    for(int iy=0;iy<Ny+1;iy++)
       {
        for(int ix=1;ix<Nx+1;ix++)
         {
             field[ix][iy][0].ex = 0.5*(field_copy[ix][iy][0].ex + field_copy[ix-1][iy][0].ex);
             field[ix][iy][0].by = 0.5*(field_copy[ix][iy][0].by + field_copy[ix-1][iy][0].by);
             //field[ix][iy][0].jx = 0.5*(field[ix][iy][0].jx + field[ix][iy-1][0].jx);
         }
       }
        //Boundary condition: Periodic
    for(int iy=0;iy<Ny+1;iy++)
       {
               field[0][iy][0].ex = field[Nx][0][0].ex;
               field[0][iy][0].by = field[Nx][0][0].by;
               //field[0][iy][0].jy = field[Nx][0][0].jy;
       }
    //blue points
    for(int iy=1;iy<Ny+1;iy++)
       {
        for(int ix=0;ix<Nx+1;ix++)
         {
             field[ix][iy][0].ey = 0.5*(field_copy[ix][iy][0].ey + field_copy[ix][iy-1][0].ey);
             field[ix][iy][0].bx = 0.5*(field_copy[ix][iy][0].bx + field_copy[ix][iy-1][0].bx);
             //field[ix][iy][0].jy = 0.5*(field[ix][iy][0].jy + field[ix-1][iy][0].jy);
         }
       }
    //Boundary condition: Periodic
    for(int ix=0;ix<Nx+1;ix++)
    {
            field[ix][0][0].ey = field[ix][Ny][0].ey;
            field[ix][0][0].bx = field[ix][Ny][0].bx;
            //field[0][iy][0].jy = field[Nx][0][0].jy;
    }
  }
    //CT method (same grid for B and E)
    else if(solver_Maxwell == 1)
    {
        
        for(int iy=0;iy<=Ny;iy++)
        {
         for(int ix=0;ix<=Nx;ix++)
          {
              field_copy[ix][iy][0].bx = field[ix][iy][0].bx;
              field_copy[ix][iy][0].by = field[ix][iy][0].by;
              field_copy[ix][iy][0].bz = field[ix][iy][0].bz;
              field_copy[ix][iy][0].ex = field[ix][iy][0].ex;
              field_copy[ix][iy][0].ey = field[ix][iy][0].ey;
              field_copy[ix][iy][0].ez = field[ix][iy][0].ez;
          }
        }
        
        //yellow points
        for(int iy=0;iy<Ny+1;iy++)
           {
            for(int ix=1;ix<Nx+1;ix++)
             {
                 field[ix][iy][0].ex = 0.5*(field_copy[ix][iy][0].ex + field_copy[ix-1][iy][0].ex);
                 field[ix][iy][0].bx = 0.5*(field_copy[ix][iy][0].bx + field_copy[ix-1][iy][0].bx);
                 //field[ix][iy][0].jx = 0.5*(field[ix][iy][0].jx + field[ix][iy-1][0].jx);
             }
           }
            //Boundary condition: Periodic
        for(int iy=0;iy<Ny+1;iy++)
           {
                   field[0][iy][0].ex = field[Nx][0][0].ex;
                   field[0][iy][0].bx = field[Nx][0][0].bx;
                   //field[0][iy][0].jy = field[Nx][0][0].jy;
           }
        //blue points
        for(int iy=1;iy<Ny+1;iy++)
           {
            for(int ix=0;ix<Nx+1;ix++)
             {
                 field[ix][iy][0].ey = 0.5*(field_copy[ix][iy][0].ey + field_copy[ix][iy-1][0].ey);
                 field[ix][iy][0].by = 0.5*(field_copy[ix][iy][0].by + field_copy[ix][iy-1][0].by);
                 //field[ix][iy][0].jy = 0.5*(field[ix][iy][0].jy + field[ix-1][iy][0].jy);
             }
           }
        //Boundary condition: Periodic
        for(int ix=0;ix<Nx+1;ix++)
        {
                field[ix][0][0].ey = field[ix][Ny][0].ey;
                field[ix][0][0].by = field[ix][Ny][0].by;
                //field[0][iy][0].jy = field[Nx][0][0].jy;
        }
        
    }
    else
    {
         cout << "wrong Maxwell solver" <<endl;
    }
}

void interc_to_nodes
(int my_id, int p,FieldPoints ***field, FieldPoints ***field_copy,int solver_Maxwell)
{
    if(solver_Maxwell == 0)
    {

    for(int iy=0;iy<=Ny;iy++)
    {
     for(int ix=0;ix<=Nx;ix++)
      {
          field_copy[ix][iy][0].bx = field[ix][iy][0].bx;
          field_copy[ix][iy][0].by = field[ix][iy][0].by;
          field_copy[ix][iy][0].bz = field[ix][iy][0].bz;
          field_copy[ix][iy][0].ex = field[ix][iy][0].ex;
          field_copy[ix][iy][0].ey = field[ix][iy][0].ey;
          field_copy[ix][iy][0].ez = field[ix][iy][0].ez;
          field_copy[ix][iy][0].jx = field[ix][iy][0].ex;
          field_copy[ix][iy][0].jy = field[ix][iy][0].ey;
          field_copy[ix][iy][0].jz = field[ix][iy][0].ez;
      }
    }
    
    for(int iy=0;iy<Ny;iy++)
       {
        for(int ix=0;ix<Nx;ix++)
         {
             field[ix][iy][0].ex = 0.5*(field_copy[ix+1][iy][0].ex + field_copy[ix][iy][0].ex);
             field[ix][iy][0].by = 0.5*(field_copy[ix+1][iy][0].by + field_copy[ix][iy][0].by);
            // field[ix][iy][0].jx = 0.5*(field[ix][iy][0].jx + field[ix][iy+1][0].jx);
             
         }
       }
    
    for(int iy=0;iy<Ny+1;iy++)
    {
            field[Nx][iy][0].ex = field[0][0][0].ex;
            field[Nx][iy][0].by = field[0][0][0].by;
            //field[Nx][iy][0].jy = field[0][0][0].jy;
    }
    
    for(int iy=0;iy<Ny;iy++)
       {
        for(int ix=0;ix<Nx;ix++)
         {
             field[ix][iy][0].ey = 0.5*(field_copy[ix][iy+1][0].ey + field_copy[ix][iy][0].ey);
             field[ix][iy][0].bx = 0.5*(field_copy[ix][iy+1][0].bx + field_copy[ix][iy][0].bx);
             //field[ix][iy][0].jy = 0.5*(field[ix][iy][0].jy + field[ix+1][iy][0].jy);
         }
       }
    for(int ix=0;ix<Nx+1;ix++)
       {
           field[ix][Ny][0].ey = field[ix][0][0].ey;
           field[ix][Ny][0].bx = field[ix][0][0].bx;
           //field[ix][Ny][0].jx = field[ix][0][0].jx;
       }
   
    //Bz re location
    for(int iy=0;iy<Ny;iy++)
       {
        for(int ix=0;ix<Nx;ix++)
         {
             field[ix][iy][0].bz = 0.25*(field_copy[ix][iy][0].bz + field_copy[ix+1][iy][0].bz + field_copy[ix][iy+1][0].bz + field_copy[ix+1][iy+1][0].bz);
         }
       }
        for(int ix=0;ix<Nx+1;ix++)
       {
            field[ix][Ny][0].bz = field[ix][0][0].bz ;
       }
        for(int iy=0;iy<Ny+1;iy++)
       {
            field[Nx][iy][0].bz =  field[0][iy][0].bz ;
       }
    
    //save fields to fields copy on the node
       for(int iy=0;iy<=Ny;iy++)
             {
              for(int ix=0;ix<=Nx;ix++)
               {
                   field_copy[ix][iy][0].bx = field[ix][iy][0].bx;
                   field_copy[ix][iy][0].by = field[ix][iy][0].by;
                   field_copy[ix][iy][0].bz = field[ix][iy][0].bz;
                   field_copy[ix][iy][0].ex = field[ix][iy][0].ex;
                   field_copy[ix][iy][0].ey = field[ix][iy][0].ey;
                   field_copy[ix][iy][0].ez = field[ix][iy][0].ez;
               }
             }
    }
    else if(solver_Maxwell == 1)
    {
        
         for(int iy=0;iy<=Ny;iy++)
         {
          for(int ix=0;ix<=Nx;ix++)
           {
               field_copy[ix][iy][0].bx = field[ix][iy][0].bx;
               field_copy[ix][iy][0].by = field[ix][iy][0].by;
               field_copy[ix][iy][0].bz = field[ix][iy][0].bz;
               field_copy[ix][iy][0].ex = field[ix][iy][0].ex;
               field_copy[ix][iy][0].ey = field[ix][iy][0].ey;
               field_copy[ix][iy][0].ez = field[ix][iy][0].ez;
               //field_copy[ix][iy][0].jx = field[ix][iy][0].jx;
               //field_copy[ix][iy][0].jy = field[ix][iy][0].jy;
               //field_copy[ix][iy][0].jz = field[ix][iy][0].jz;
           }
         }
         
         for(int iy=0;iy<Ny;iy++)
            {
             for(int ix=0;ix<Nx;ix++)
              {
                  field[ix][iy][0].ex = 0.5*(field_copy[ix+1][iy][0].ex + field_copy[ix][iy][0].ex);
                  field[ix][iy][0].bx = 0.5*(field_copy[ix+1][iy][0].bx + field_copy[ix][iy][0].bx);
                 // field[ix][iy][0].jx = 0.5*(field[ix][iy][0].jx + field[ix][iy+1][0].jx);
                  
              }
            }
         
         for(int iy=0;iy<Ny+1;iy++)
         {
                 field[Nx][iy][0].ex = field[0][0][0].ex;
                 field[Nx][iy][0].bx = field[0][0][0].bx;
                 //field[Nx][iy][0].jy = field[0][0][0].jy;
         }
         
         for(int iy=0;iy<Ny;iy++)
            {
             for(int ix=0;ix<Nx;ix++)
              {
                  field[ix][iy][0].ey = 0.5*(field_copy[ix][iy+1][0].ey + field_copy[ix][iy][0].ey);
                  field[ix][iy][0].by = 0.5*(field_copy[ix][iy+1][0].by + field_copy[ix][iy][0].by);
                  //field[ix][iy][0].jy = 0.5*(field[ix][iy][0].jy + field[ix+1][iy][0].jy);
              }
            }
         for(int ix=0;ix<Nx+1;ix++)
            {
                field[ix][Ny][0].ey = field[ix][0][0].ey;
                field[ix][Ny][0].by = field[ix][0][0].by;
                //field[ix][Ny][0].jx = field[ix][0][0].jx;
            }
        
         //Bz re location
         for(int iy=0;iy<Ny;iy++)
            {
             for(int ix=0;ix<Nx;ix++)
              {
                  field[ix][iy][0].bz = 0.25*(field_copy[ix][iy][0].bz + field_copy[ix+1][iy][0].bz + field_copy[ix][iy+1][0].bz + field_copy[ix+1][iy+1][0].bz);
              }
            }
             for(int ix=0;ix<Nx+1;ix++)
            {
                 field[ix][Ny][0].bz = field[ix][0][0].bz ;
            }
             for(int iy=0;iy<Ny+1;iy++)
            {
                 field[Nx][iy][0].bz =  field[0][iy][0].bz ;
            }
         
         //save fields to fields copy on the node
            for(int iy=0;iy<=Ny;iy++)
                  {
                   for(int ix=0;ix<=Nx;ix++)
                    {
                        field_copy[ix][iy][0].bx = field[ix][iy][0].bx;
                        field_copy[ix][iy][0].by = field[ix][iy][0].by;
                        field_copy[ix][iy][0].bz = field[ix][iy][0].bz;
                        field_copy[ix][iy][0].ex = field[ix][iy][0].ex;
                        field_copy[ix][iy][0].ey = field[ix][iy][0].ey;
                        field_copy[ix][iy][0].ez = field[ix][iy][0].ez;
                    }
                  }
        
    }
    else
    {
         cout << "wrong Maxwell solver" <<endl;
    }
    
 
}
