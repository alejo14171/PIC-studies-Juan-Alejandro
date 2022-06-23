#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//==========================================
//  Setup all initial conditions 
//==========================================
void setup
(int my_id ,int p, 
 GridPoints ***xgc,
 FieldPoints ***field,
 FieldPoints ***field_copy,
 vector<ThreeVec> &xi, vector<ThreeVec> &pi,vector<ThreeVec> &x_old,vector<ThreeVec> &x_new,vector<ThreeVec> &E_p, vector<ThreeVec> &B_p, vector<double> &qi,vector<double> &Vpar,vector<double> &Vperp)
{
    double dxp;
    double dyp;
    double del_x,kp,xa;
    double ranf;
    double xcosine;
    double tmp_far;
    ThreeVec tmp;

    for(int ix=0;ix<Nx;ix++)
        {
            for(int iy=0;iy<Ny;iy++)
            {
                for(int iz=0;iz<Nz;iz++)
                {
                    xgc[ix][iy][iz].x = xmin + dx*0.5 + double(ix)*dx;
                    xgc[ix][iy][iz].y = ymin + dy*0.5 + double(iy)*dy;
                    xgc[ix][iy][iz].z = zmin + dz*0.5 + double(iz)*dz;
                }
            }
        }

    // Initialize cells
    for(int ix=0;ix<Nx+1;ix++)
    {
        for(int iy=0;iy<Ny+1;iy++)
        {
            for(int iz=0;iz<Nz+1;iz++)
            {
                    field[ix][iy][iz].x = xmin + double(ix)*dx;
                    field[ix][iy][iz].y = ymin + double(iy)*dy;
                    field[ix][iy][iz].z = zmin + double(iz)*dz;
                    field[ix][iy][iz].ex = 0.0;
                    field[ix][iy][iz].ey = 0.0;
                    field[ix][iy][iz].ez = 0.0;
                    field[ix][iy][iz].bx = 0.0;
                    field[ix][iy][iz].by = 0.0;
                    field[ix][iy][iz].bz = 0.0;
                    field[ix][iy][iz].jx = 0.0;
                    field[ix][iy][iz].jy = 0.0;
                    field[ix][iy][iz].jz = 0.0;
                    field_copy[ix][iy][iz].ex = 0.0;
                    field_copy[ix][iy][iz].ey = 0.0;
                    field_copy[ix][iy][iz].ez = 0.0;
                    field_copy[ix][iy][iz].bx = 0.0;
                    field_copy[ix][iy][iz].by = 0.0;
                    field_copy[ix][iy][iz].bz = 0.0;
                    field_copy[ix][iy][iz].jx = 0.0;
                    field_copy[ix][iy][iz].jy = 0.0;
                    field_copy[ix][iy][iz].jz = 0.0;
            }
        }
    }

/*
//CC initial condition
    for(int iy=0;iy<Ny;iy++)
    {
     for(int ix=0;ix<Nx;ix++)
      {
          xgc[ix][iy][0].ex = 0.0;
          xgc[ix][iy][0].ez = 0.0;
          xgc[ix][iy][0].bx = 1.0;
          xgc[ix][iy][0].bz = 0.0;
          
          if(xgc[ix][iy][0].x <= 0.5)
          {
              xgc[ix][iy][0].ey = 1.0;
              xgc[ix][iy][0].by = -0.75;
          }
          else
          {
              xgc[ix][iy][0].ey = -1.0;
              xgc[ix][iy][0].by = 0.75;
          }
        }
    }
 */
  
    //nodes initial condition
    for(int iy=0;iy<Ny+1;iy++)
    {
     for(int ix=0;ix<Nx+1;ix++)
       {
           tmp_far = pow((field[ix][iy][0].x),2) + pow((field[ix][iy][0].y),2);
           if(tmp_far==0)
           {
               field[ix][iy][0].bz=0.0;
               field[ix][iy][0].ex=0.0;
               field[ix][iy][0].ey=0.0;
           }
           else
           {
               field[ix][iy][0].bz = sqrt(tmp_far);
               field[ix][iy][0].ex = (-0.01*(field[ix][iy][0].x))/pow(tmp_far,1.5);
               field[ix][iy][0].ey = (-0.01*(field[ix][iy][0].y))/pow(tmp_far,1.5);
           }
       }
     }

    double ranf1,ranf2;
    double vel_0, del_v_x,del_v_y,del_v_z; 


    // No-restart
    if (restart == 0)
    {
	initC = 0;
 
  	for(int ip=0;ip<Npi;ip++)
    	{
      	  if(ip%p==my_id)
      	  {
            qi.push_back(partweight); 
            // qi.push_back(init_n*(double)Nx*(double)Ny*(double)Nz/(double)Npi);

              //ranf = ((double) rand()) / ((double) RAND_MAX);
              tmp.setX(0.5);
              //tmp.setX(xmin + (xmax-xmin-1.0e-15)*ranf);
              //ranf1 = ((double) rand()) / ((double) RAND_MAX);
              tmp.setY(0.0);
              //tmp.setY(ymin + (ymax-ymin-1.0e-15)*ranf1);
              //ranf2 = ((double) rand()) / ((double) RAND_MAX);
              tmp.setZ(0.0);
              //tmp.setZ(zmin+(zmax-zmin)*ranf2);
              xi.push_back(tmp);

              //vel_0 = 0.0;
              vel_0 = gaussian(vth_i);
              del_v_x = 0.1;//vel_0;
              //vel_0 = 0.0;
              vel_0 = gaussian_speed(vth_i);
              ranf = ((double) rand() / (RAND_MAX));
              del_v_y = 0.0;//vel_0 * sin(2.0*PI0*ranf);
              del_v_z = 0.0;//vel_0 * cos(2.0*PI0*ranf);
              tmp.setX(del_v_x);
              tmp.setY(del_v_y);
              tmp.setZ(del_v_z);
              pi.push_back(tmp);
              
//For output porpuses
              tmp.setX(0.0);
              tmp.setY(0.0);
              tmp.setZ(0.0);
              x_old.push_back(tmp);
              tmp.setX(0.0);
              tmp.setY(0.0);
              tmp.setZ(0.0);
              x_new.push_back(tmp);
	          tmp.setX(0.0);
              tmp.setY(0.0);
              tmp.setZ(0.0);
              E_p.push_back(tmp);
              tmp.setX(0.0);
              tmp.setY(0.0);
              tmp.setZ(0.0);
              B_p.push_back(tmp);
              Vpar.push_back(0.0);
              Vperp.push_back(0.0);
              
          }
        }

        cout << " -- my_id " << my_id << "parts total (ppc) :" << (int)pi.size() << " " << (double)((int)pi.size())/((double)Nx*Ny*Nz) <<endl;
        cout << " particle weight " << partweight <<endl;

    }
    // Restart function
    else
    { 
	if (my_id == 0)	cout << " -- loading from files -- " << endl;

	double xr, yr, ur, vr, qr;
	int rnt; // time of last restart cycle

	DIR *dir;
	struct dirent *ent;
	int scounter = -2; // for some reason it counts two more
	if ((dir = opendir ("./output/electron/")) != NULL) 
	{
		// counts the number of files
		while ((ent = readdir(dir)) != NULL) 
		{
		        ++scounter;
		}

		// last saved time cycle will be the number of files divided by number of procs divided by saving frequency (i.e. nrec) -- since it saves at any nrec-1 cycle, we still need to take off -1 
		rnt = (scounter/p-1)*nrec;
		initC = rnt; // because is the next cycle we want to restart the simulation from 
  		closedir (dir);

		if(my_id==0) 
			cout << endl << endl << " Old nt to retrieve the particle files is: " << rnt << " new nt is "  << initC <<endl<<endl;  
	}
       	else 
	{
  		cout <<  " ERROR -  Could not open the directory" << endl;
	}

	
        char filename[80];
	fstream  afile;
        sprintf(filename,"./output/ion/parts_%08d_%02d.dat",rnt,my_id);
        cout << " proc " << my_id << " is " << "reading filename : \t" << filename  <<  endl;
  	afile.open(filename, ios::in);

	//  file format
        //   outf << xx[np].getX() << " \t"<< xx[np].getY() << " \t" <<  px[np].getX() << " \t" << px[np].getY() << " \t" << qx[np] <<  endl;
   	while(afile >> xr >> yr >> ur >> vr >> qr)
	{
		tmp.setX(xr); 
		tmp.setY(yr); 
		ranf2 = ((double) rand() / (RAND_MAX));
        	tmp.setZ(zmin+(zmax-zmin)*ranf2);
       		xi.push_back(tmp);

   		tmp.setX(ur);
                tmp.setY(vr);
        	vel_0 = gaussian_speed(vth_i);
        	ranf = ((double) rand() / (RAND_MAX));
        	del_v_z = vel_0 * sin(2.0*PI0*ranf);
        	tmp.setZ(del_v_z);
        	pi.push_back(tmp);
         
		qi.push_back(qr); 
	}
	afile.close();
    }
}
void read_initial_fields
(int my_id ,int p,
 GridPoints ***xgc,int nt)
{
    double *array_tmp;
    int ntt;
    ntt = nt/nout;
    
    ifstream fin;
    char filename[80];
    array_tmp =  (double*)malloc(((Nx*Ny)+1)*sizeof(double));
    double ddd;
    int counter;
    
 //   if(my_id==0)
    {
        sprintf(filename,"./output/Ex_%d.out",ntt);
        fin.open(filename);
        counter = 0;
        while (counter < Nx*Ny)
        {
            fin >> ddd;
    
            array_tmp[counter] = ddd;
            counter++;
        }
        fin.close();
           //share the array_tmp to all the processors
        //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
   //     MPI_Barrier(MPI_COMM_WORLD);
        counter=0;
        for(int ix=0;ix<Nx;ix++)
        {
            for(int iy=0;iy<Ny;iy++)
            {
                xgc[ix][iy][0].ex = array_tmp[counter];
                counter++;
            }
        }
    
//      if(my_id==0)
      {
          sprintf(filename,"./output/Ey_%d.out",ntt);
          fin.open(filename);
          counter = 0;
          while (counter < Nx*Ny)
          {
              fin >> ddd;
      
              array_tmp[counter] = ddd;
              counter++;
          }
          fin.close();
             //share the array_tmp to all the processors
          //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
     
     //share the array_tmp to all the processors
 //   MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             xgc[ix][iy][0].ey = array_tmp[counter];
             counter++;
         }
    }
    
 //     if(my_id==0)
       {
           sprintf(filename,"./output/Ez_%d.out",ntt);
           fin.open(filename);
           counter = 0;
           while (counter < Nx*Ny)
           {
               fin >> ddd;
       
               array_tmp[counter] = ddd;
               counter++;
           }
           fin.close();
              //share the array_tmp to all the processors
           //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       }
      
      //share the array_tmp to all the processors
 //   MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             xgc[ix][iy][0].ez = array_tmp[counter];
             counter++;
         }
    }
    
 //   if(my_id==0)
       {
           sprintf(filename,"./output/Bx_%d.out",ntt);
           fin.open(filename);
           counter = 0;
           while (counter < Nx*Ny)
           {
               fin >> ddd;
       
               array_tmp[counter] = ddd;
               counter++;
           }
           fin.close();
              //share the array_tmp to all the processors
           //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       }
      
      //share the array_tmp to all the processors
  //  MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             xgc[ix][iy][0].bx = array_tmp[counter];
             counter++;
         }
    }
    
//    if(my_id==0)
    {
        sprintf(filename,"./output/By_%d.out",ntt);
        fin.open(filename);
        counter = 0;
          while (counter < Nx*Ny)
          {
              fin >> ddd;
      
              array_tmp[counter] = ddd;
              counter++;
          }
          fin.close();
           //share the array_tmp to all the processors
        //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
     
     //share the array_tmp to all the processors
  //  MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             xgc[ix][iy][0].by = array_tmp[counter];
             counter++;
         }
    }

 //      if(my_id==0)
       {
           
           sprintf(filename,"./output/Bz_%d.out",ntt);
           fin.open(filename);
           counter = 0;
           while (counter < Nx*Ny)
           {
               fin >> ddd;
       
               array_tmp[counter] = ddd;
               counter++;
           }
           fin.close();
              //share the array_tmp to all the processors
           //MPI_Bcast(&array_tmp,(Nx*Ny)+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       }
      
      //share the array_tmp to all the processors
  // MPI_Barrier(MPI_COMM_WORLD);
    counter=0;
    for(int ix=0;ix<Nx;ix++)
     {
         for(int iy=0;iy<Ny;iy++)
         {
             //cout << array_tmp[counter];
             //cout << endl;
             xgc[ix][iy][0].bz = array_tmp[counter];
             counter++;
         }
    }
free(array_tmp);
}
