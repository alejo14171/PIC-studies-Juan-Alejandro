void getCellSize
(int my_id, int p) 
{
   int ny_mod, nvxe_mod, nx_mod;
   int ny_div, nvxe_div, nx_div;
   int ny_tmp, nvxe_tmp, nx_tmp;

   pvx = nproc_y;
   px = nproc_x; ;

   p_real = px*pvx;
   
   // pvx is the number of procs in y direction
   // px  is the number of procs in x direction
   // mod is the number of procs that have extra cell
   // div is the number of cells (dividable)
   // tmp is the number of cells in each proc
  
   ny_mod = ny_nodes_global%pvx;
   nvxe_mod = nvxe_global%pvx;
   nx_mod = nx_nodes_global%px;

   ny_div = (ny_nodes_global-ny_mod);
   nvxe_div = (nvxe_global-nvxe_mod);
   nx_div = (nx_nodes_global-nx_mod);

   ny_tmp = ny_div/pvx;
   nvxe_tmp = nvxe_div/pvx;
   nx_tmp = nx_div/px;

   if(my_id<ny_mod*px)
         nlocal_nodes.y = ny_tmp+1;
   else
         nlocal_nodes.y = ny_tmp;

   if(my_id<nvxe_mod*px)
         nlocal.vxe = nvxe_tmp+1;
   else
         nlocal.vxe = nvxe_tmp;

   if(my_id%px<nx_mod)
         nlocal_nodes.x = nx_tmp+1;
   else
         nlocal_nodes.x = nx_tmp;

   ny_mod = ny_global%pvx;  
   nvxe_mod = nvxe_global%pvx;  
   nx_mod = nx_global%px;
  
   ny_div = (ny_global-ny_mod);
   nvxe_div = (nvxe_global-nvxe_mod);
   nx_div = (nx_global-nx_mod);

   ny_tmp = ny_div/pvx;
   nvxe_tmp = nvxe_div/pvx;
   nx_tmp = nx_div/px;
   
   if(my_id<ny_mod*px)
         nlocal.y = ny_tmp+1; 
   else
         nlocal.y = ny_tmp;

   if(my_id<nvxe_mod*px)
         nlocal.vxe = nvxe_tmp+1; 
   else
         nlocal.vxe = nvxe_tmp;

   if(my_id%px<nx_mod)
         nlocal.x = nx_tmp+1;
   else
         nlocal.x = nx_tmp;


   // pid is the ID number of proc
   pid.x  = my_id%px;
   pid.vx = (my_id - pid.x)/pvx;

  
   // Start cell in global 
   if(pid.x<nx_mod) 
     nstart.x = pid.x*(nx_tmp+1);
   else
     nstart.x = pid.x*nx_tmp + nx_mod;

   if(pid.vx<ny_mod) 
     nstart.y = pid.vx*(ny_tmp+1);
   else
     nstart.y = pid.vx*ny_tmp + ny_mod;

   if(pid.vx<nvxe_mod) 
     nstart.vxe = pid.vx*(nvxe_tmp+1);
   else
     nstart.vxe = pid.vx*nvxe_tmp + nvxe_mod;

   //nlocal_nodes.x = nlocal.x + 1;
   //nlocal_nodes.y = nlocal.y + 1;

   // Total size  
   ntot.x   = nlocal.x   +2*nbuffer; 
   ntot.y = nlocal.y +2*nbuffer; 
   ntot.vxe = nlocal.vxe +2*nbuffer; 
  
   ntot_nodes.x = nlocal_nodes.x +2*nbuffer;
   ntot_nodes.y = nlocal_nodes.y +2*nbuffer;

printf("---------------------------------------\n");
printf("my_id=%d,p=%d\n",my_id,p);
printf("nlocal.x=%d,nlocal.y=%d\n",nlocal.x,nlocal.y);
printf("nlocal_nodes.x=%d,nlocal_nodes.y=%d\n",nlocal_nodes.x,nlocal_nodes.y);
printf("ntot.x=%d, ntot.y=%d\n",ntot.x,ntot.y);
printf("ntot_nodes.x=%d,ntot_nodes.y=%d\n",ntot_nodes.x,ntot_nodes.y);
printf("---------------------------------------\n");

}


void getCoordinate
/*(Coordinate ntot, Coordinate nstart, 
 Values Min, Values dd,
 double *x, double *y, double *vxe,
 int *ind_x, int *ind_y, int *ind_vxe)
*/
(int my_id)
{
   ind_x = (int*)malloc(ntot.x*sizeof(int));
   ind_y = (int*)malloc(ntot.y*sizeof(int));
   ind_vxe = (int*)malloc(ntot.vxe*sizeof(int));

   x = (double*)malloc(ntot.x*sizeof(double));
   y = (double*)malloc(ntot.y*sizeof(double));
   vxe = (double*)malloc(ntot.vxe*sizeof(double));

   // set index notation the same as vdf_i and vdf_e
   // takes the buffer cells into account 
 for(int i=0; i<ntot.x; i++)
 {
   ind_x[i] = nstart.x + (i-nbuffer);
   x[i] = Min.x + dd.x*ind_x[i] + dd.x/2.0;
 }

 for(int i=0; i<ntot.y; i++)
 {
   ind_y[i] = nstart.y + (i-nbuffer);
   y[i] = Min.y + dd.y*ind_y[i] + dd.y/2.0;
 }
   for(int i=0; i<ntot.vxe; i++)
   {
     ind_vxe[i] = nstart.vxe + (i-nbuffer);
     vxe[i] = Min.vxe + dd.vxe*ind_vxe[i] + dd.vxe/2.0;
   }

   if(my_id==0) 
   {
      cout << "x,y,vxe: " << ntot.x << " "<< ntot.y << " "<< ntot.vxe <<endl;
   }

   // global index
   x_global = (double*)malloc(nx_global*sizeof(double));
   for(int i=0; i<nx_global;i++)
   x_global[i] = Min.x +dd.x*(double)(i);
}
