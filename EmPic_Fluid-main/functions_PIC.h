#include "extractPositions.h"
//==========================================
// Weighting density on cells 
//   irec = 0 (only density) 1 (save mode)
//   bcflag = 0 (wall) , 1 (periodic)
//==========================================
void reduced_VDF
(GridPoints ***gc,
 vector<ThreeVec> &xx, vector<ThreeVec> &px, vector<double> &qx, vector<ThreeVec> &B_p,vector<double> &Vpar, vector<double> &Vperp,
 int iflag,  // ion - electron
 int irec,   // 0- only density, 1-everything
 int bcflag) // bc
{
    double xa,hxleft;
    double ya,hyleft;
    double za,hzleft;
    double px_tmp, py_tmp,pz_tmp;
    double bx_tmp, by_tmp,bz_tmp;
    double bi;
    double d2,dp,dd,v;
    int ix,iy;
    
    int nvpa=511;
    int nvpe=511;
    double xamin=-2.0;
    double xamax=2.0;
    double yamin=0.0;
    double yamax=2.0;
    double da=(xamax-xamin)/(nvpa-1);
    double de=(yamax-yamin)/(nvpe-1);
    double ax,ax1,ay,ay1;
    
    for(int ix=0;ix<Nx;ix++)
    {
        for(int iy=0;iy<Ny;iy++)
        {
            gc[ix][iy][0].fvpape = 0.0;
            gc[ix][iy][0].fxvx = 0.0;
        }
    }

    for(int ip=0;ip<(int)xx.size();ip++)
    {

       px_tmp = px[ip].getX();
       py_tmp = px[ip].getY();
       pz_tmp = px[ip].getZ();
       
       bx_tmp = B_p[ip].getX();
       by_tmp = B_p[ip].getY();
       bz_tmp = B_p[ip].getZ();
       bi = pow((pow(bx_tmp,2) + pow(by_tmp,2) + pow(bz_tmp,2)),0.5);

        if( bi > 1.e-3)
        {
            bx_tmp=bx_tmp/bi;
            by_tmp=by_tmp/bi;
            bz_tmp=bz_tmp/bi;
        }
        d2= pow(px_tmp,2) + pow(py_tmp,2) + pow(pz_tmp,2);
        v = pow(d2,0.5);
        dp=bx_tmp*px_tmp+by_tmp*py_tmp+bz_tmp*pz_tmp;
        dd=d2-(dp*dp);
        d2=0.0;
        if (dd>0) d2=pow(dd,0.5);
        
        ix=int(((dp-xamin)/da)+1);
        ax1=((dp-xamin)/da-double(ix)+1.0);
        ax=1.-ax1;
        iy=int(((d2-yamin)/de)+1);
        ay1=(d2-yamin)/de-double(iy)+1.0;
        ay=1.-ay1;
        if(ix > 0 && iy > 0 && ix < nvpa && iy < nvpe)
        {
            gc[ix][iy][0].fvpape= gc[ix][iy][0].fvpape +(ax*ay);
            gc[ix+1][iy][0].fvpape = gc[ix+1][iy][0].fvpape + (ax1*ay);
            gc[ix][iy+1][0].fvpape =  gc[ix][iy+1][0].fvpape + (ax*ay1);
            gc[ix+1][iy+1][0].fvpape =  gc[ix+1][iy+1][0].fvpape + (ax1*ay1);
        }
        Vpar[ip] = dp;
        Vperp[ip] = d2;
    }
    
    xamin=xmin;
    xamax=xmax;
    yamin=-2.0;
    yamax=2.0;
    da=(xamax-xamin)/(nvpa-1);
    de=(yamax-yamin)/(nvpe-1);
    for(int ip=0;ip<(int)xx.size();ip++)
    {
       px_tmp = xx[ip].getX();
       py_tmp = px[ip].getX();

       ix=int(((px_tmp -xamin)/da)+1);
       ax1=((px_tmp-xamin)/da-double(ix)+1.0);
       ax=1.-ax1;
       iy=int(((py_tmp-yamin)/de)+1);
       ay1=(py_tmp-yamin)/de-double(iy)+1.0;
       ay=1.-ay1;
       if(ix > 0 && iy > 0 && ix < nvpa && iy < nvpe)
       {
            gc[ix][iy][0].fxvx = gc[ix][iy][0].fxvx +(ax*ay);
            gc[ix+1][iy][0].fxvx = gc[ix+1][iy][0].fxvx + (ax1*ay);
            gc[ix][iy+1][0].fxvx =  gc[ix][iy+1][0].fxvx + (ax*ay1);
            gc[ix+1][iy+1][0].fxvx =  gc[ix+1][iy+1][0].fxvx + (ax1*ay1);
        }
    }
    
}

void moments_VDF
(GridPoints ***gc,
 vector<ThreeVec> &xx, vector<ThreeVec> &px, vector<double> &qx,
 int iflag,  // ion - electron 
 int irec,   // 0- only density, 1-everything
 int bcflag) // bc
{
    double xa,hxleft;
    double ya,hyleft;
    double za,hzleft;
    double px_tmp, py_tmp,pz_tmp;
    int ixleft,ixright;
    int iyleft,iyright;
    int izleft,izright;

    for(int ix=0;ix<Nx;ix++)
    {
        for(int iy=0;iy<Ny;iy++)
        {
            for(int iz=0;iz<Nz;iz++)
            {
                    gc[ix][iy][iz].den_i = 0.0;
                    gc[ix][iy][iz].ux_i = 0.0;
                    gc[ix][iy][iz].uy_i = 0.0;
                    gc[ix][iy][iz].uz_i = 0.0;
                    gc[ix][iy][iz].pxx_i = 0.0;
                    gc[ix][iy][iz].pxy_i = 0.0;
                    gc[ix][iy][iz].pyy_i = 0.0;
            }
        }
    }

   // double count= 0.0;
    for(int ip=0;ip<(int)xx.size();ip++)
    {
       xa = xx[ip].getX();
       ya = xx[ip].getY();
//       za = xx[ip].getZ();
       px_tmp = px[ip].getX();
       py_tmp = px[ip].getY();
       pz_tmp = px[ip].getZ();

       extractPositions_x_cell(gc,xa,hxleft,ixright,ixleft,bcflag); // 0: non-periodic, 1: periodic
       extractPositions_y_cell(gc,ya,hyleft,iyright,iyleft,bcflag); // 0: non-periodic, 1: periodic

       // count += (hxright+hxleft)*(hyleft+hyright);
       if(iflag==0)
       {
           gc[ixright][iyright][0].den_i += hxleft*hyleft*qx[ip];
           gc[ixright][iyleft][0].den_i  += hxleft*(1.0-hyleft)*qx[ip];
           gc[ixleft][iyright][0].den_i  += (1.0-hxleft)*hyleft*qx[ip];
           gc[ixleft][iyleft][0].den_i   += (1.0-hxleft)*(1.0-hyleft)*qx[ip];
           
           gc[ixright][iyright][0].ux_i += hxleft*hyleft*px_tmp*qx[ip];
           gc[ixright][iyleft][0].ux_i  += hxleft*(1.0-hyleft)*px_tmp*qx[ip];
           gc[ixleft][iyright][0].ux_i  += (1.0-hxleft)*hyleft*px_tmp*qx[ip];
           gc[ixleft][iyleft][0].ux_i   += (1.0-hxleft)*(1.0-hyleft)*px_tmp*qx[ip];
            
           gc[ixright][iyright][0].uy_i += hxleft*hyleft*py_tmp*qx[ip];
           gc[ixright][iyleft][0].uy_i  += hxleft*(1.0-hyleft)*py_tmp*qx[ip];
           gc[ixleft][iyright][0].uy_i  += (1.0-hxleft)*hyleft*py_tmp*qx[ip];
           gc[ixleft][iyleft][0].uy_i   += (1.0-hxleft)*(1.0-hyleft)*py_tmp*qx[ip];
            
           gc[ixright][iyright][0].uz_i += hxleft*hyleft*pz_tmp*qx[ip];
           gc[ixright][iyleft][0].uz_i  += hxleft*(1.0-hyleft)*pz_tmp*qx[ip];
           gc[ixleft][iyright][0].uz_i  += (1.0-hxleft)*hyleft*pz_tmp*qx[ip];
           gc[ixleft][iyleft][0].uz_i   += (1.0-hxleft)*(1.0-hyleft)*pz_tmp*qx[ip];
            
           gc[ixright][iyright][0].pxx_i += hxleft*hyleft*px_tmp*px_tmp*qx[ip];
           gc[ixright][iyleft][0].pxx_i  += hxleft*(1.0-hyleft)*px_tmp*px_tmp*qx[ip];
           gc[ixleft][iyright][0].pxx_i  += (1.0-hxleft)*hyleft*px_tmp*px_tmp*qx[ip];
           gc[ixleft][iyleft][0].pxx_i   += (1.0-hxleft)*(1.0-hyleft)*px_tmp*px_tmp*qx[ip];
            
           gc[ixright][iyright][0].pxy_i += hxleft*hyleft*px_tmp*py_tmp*qx[ip];
           gc[ixright][iyleft][0].pxy_i  += hxleft*(1.0-hyleft)*px_tmp*py_tmp*qx[ip];
           gc[ixleft][iyright][0].pxy_i  += (1.0-hxleft)*hyleft*px_tmp*py_tmp*qx[ip];
           gc[ixleft][iyleft][0].pxy_i   += (1.0-hxleft)*(1.0-hyleft)*px_tmp*py_tmp*qx[ip];
            
           gc[ixright][iyright][0].pxz_i += hxleft*hyleft*px_tmp*pz_tmp*qx[ip];
           gc[ixright][iyleft][0].pxz_i  += hxleft*(1.0-hyleft)*px_tmp*pz_tmp*qx[ip];
           gc[ixleft][iyright][0].pxz_i  += (1.0-hxleft)*hyleft*px_tmp*pz_tmp*qx[ip];
           gc[ixleft][iyleft][0].pxz_i   += (1.0-hxleft)*(1.0-hyleft)*px_tmp*pz_tmp*qx[ip];
           
           gc[ixright][iyright][0].pyz_i += hxleft*hyleft*py_tmp*pz_tmp*qx[ip];
           gc[ixright][iyleft][0].pyz_i  += hxleft*(1.0-hyleft)*py_tmp*pz_tmp*qx[ip];
           gc[ixleft][iyright][0].pyz_i  += (1.0-hxleft)*hyleft*py_tmp*pz_tmp*qx[ip];
           gc[ixleft][iyleft][0].pyz_i   += (1.0-hxleft)*(1.0-hyleft)*py_tmp*pz_tmp*qx[ip];
           
           gc[ixright][iyright][0].pzz_i += hxleft*hyleft*pz_tmp*pz_tmp*qx[ip];
           gc[ixright][iyleft][0].pzz_i  += hxleft*(1.0-hyleft)*pz_tmp*pz_tmp*qx[ip];
           gc[ixleft][iyright][0].pzz_i  += (1.0-hxleft)*hyleft*pz_tmp*pz_tmp*qx[ip];
           gc[ixleft][iyleft][0].pzz_i   += (1.0-hxleft)*(1.0-hyleft)*pz_tmp*pz_tmp*qx[ip];
           
       }
    }
    
    for(int ix=0;ix<Nx;ix++)
    {
        for(int iy=0;iy<Ny;iy++)
        {
                if(gc[ix][iy][0].den_i<1e-3)
                {
                    gc[ix][iy][0].den_i = 1e-3;
                }
                else
                {
                    gc[ix][iy][0].ux_i = gc[ix][iy][0].ux_i/gc[ix][iy][0].den_i;
                    gc[ix][iy][0].uy_i = gc[ix][iy][0].uy_i/gc[ix][iy][0].den_i;
                    gc[ix][iy][0].uz_i = gc[ix][iy][0].uz_i/gc[ix][iy][0].den_i;
                }
        }
    }
}
//==========================================
// Weighting for electric field 
//  nflag: 0 = pre, 1 = main
//  iflag: 0 = ion, 1 = electron
//==========================================
void boris
(FieldPoints ***field,
 vector<ThreeVec> &xx, vector<ThreeVec> &px,
 vector<ThreeVec> &E_p, vector<ThreeVec> &B_p,
 int iflag,int solver_Boris,
 double dtweight)
{
    double xa,hxleft;
    double ya,hyleft;
    double za,hzleft;
    int ixleft,ixright;
    int iyleft,iyright;
    int izleft,izright;
    double ex_tmp, ey_tmp, ez_tmp;
    double bx_tmp, by_tmp, bz_tmp;

    double px_pre, px_post,energy_ke;

    ThreeVec accel;
    ThreeVec B0;
    ThreeVec tmp;

    double gamma;
    double B0_z;
    double pix_s, piy_s, piz_s;
    double ttx,tty,ttz;
    double ssx,ssy,ssz;
    double vf_x,vf_y,vf_z;
    double vprime_x,vprime_y,vprime_z;
    double vminus_x,vminus_y,vminus_z;
    double vpar,theta_n;
    double bxn,byn,bzn,Bn;
    double vplus_x,vplus_y,vplus_z;
    double tz,ss;
    double q_m;
    double dteff; 

    if(iflag==0)      {q_m = ech/mass_i;  dteff= dt;}
    else if(iflag==1) {q_m = -ech/mass_e; dteff= dt;}
    else   cout << " Not designed iflag =  " << iflag <<endl;

    if(dtweight!=1.0) dteff *= dtweight; // 0.5 or 1 

    double atmp,atmp2;
    atmp  = q_m/2.0*dteff;
    atmp2 = q_m*dteff;
    for(int ip=0;ip<(int)xx.size();ip++)
    {
         vminus_x = 0.0;
         vminus_y = 0.0;
         vminus_z = 0.0;
         vprime_x = 0.0;
         vprime_y = 0.0;
         vprime_z = 0.0;
         vplus_x = 0.0;
         vplus_y = 0.0;
         vplus_z = 0.0;
         vf_x = 0.0;
         vf_y = 0.0;
         vf_z = 0.0;
       
        xa = xx[ip].getX();
        ya = xx[ip].getY();
//       za = xx[ip].getZ();

       extractPositions_x_node(field,xa,hxleft,ixright,ixleft); 
       extractPositions_y_node(field,ya,hyleft,iyright,iyleft);

       ex_tmp =  hxleft*hyleft  *field[ixright][iyright][0].ex
               + hxleft*(1.0-hyleft) *field[ixright][iyleft][0].ex
               + (1.0-hxleft)*hyleft *field[ixleft][iyright][0].ex
               + (1.0-hxleft)*(1.0-hyleft)*field[ixleft][iyleft][0].ex;
       ey_tmp =  hxleft*hyleft  *field[ixright][iyright][0].ey
               + hxleft*(1.0-hyleft) *field[ixright][iyleft][0].ey
               + (1.0-hxleft)*hyleft *field[ixleft][iyright][0].ey
               + (1.0-hxleft)*(1.0-hyleft)*field[ixleft][iyleft][0].ey;
       ez_tmp =  hxleft*hyleft  *field[ixright][iyright][0].ez
               + hxleft*(1.0-hyleft) *field[ixright][iyleft][0].ez
               + (1.0-hxleft)*hyleft *field[ixleft][iyright][0].ez
               + (1.0-hxleft)*(1.0-hyleft)*field[ixleft][iyleft][0].ez;
       bx_tmp =  hxleft*hyleft  *field[ixright][iyright][0].bx
               + hxleft*(1.0-hyleft) *field[ixright][iyleft][0].bx
               + (1.0-hxleft)*hyleft *field[ixleft][iyright][0].bx
               + (1.0-hxleft)*(1.0-hyleft)*field[ixleft][iyleft][0].bx;
       by_tmp =  hxleft*hyleft  *field[ixright][iyright][0].by
               + hxleft*(1.0-hyleft) *field[ixright][iyleft][0].by
               + (1.0-hxleft)*hyleft *field[ixleft][iyright][0].by
               + (1.0-hxleft)*(1.0-hyleft)*field[ixleft][iyleft][0].by;
       bz_tmp =  hxleft*hyleft  *field[ixright][iyright][0].bz
               + hxleft*(1.0-hyleft) *field[ixright][iyleft][0].bz
               + (1.0-hxleft)*hyleft *field[ixleft][iyright][0].bz
               + (1.0-hxleft)*(1.0-hyleft)*field[ixleft][iyleft][0].bz;

       accel.setX(ex_tmp);
       accel.setY(ey_tmp);
       accel.setZ(ez_tmp);
       B0.setX(bx_tmp);
       B0.setY(by_tmp);
       B0.setZ(bz_tmp);

       // Boris push
       gamma = 1.0;
       //gamma = sqrt(1.0+px[ip].square()/cph/cph);
       if(solver_Boris == 0)
       {
            ttx = q_m*B0.getX()*dteff/2.0/gamma;
            tty = q_m*B0.getY()*dteff/2.0/gamma;
            ttz = q_m*B0.getZ()*dteff/2.0/gamma;

            ssx =  2.0*ttx/(1.0+(ttx*ttx));
            ssy =  2.0*tty/(1.0+(tty*tty));
            ssz =  2.0*ttz/(1.0+(ttz*ttz));
                  
            vminus_x = px[ip].getX() + accel.getX()*atmp;
            vminus_y = px[ip].getY() + accel.getY()*atmp;
            vminus_z = px[ip].getZ() + accel.getZ()*atmp;

            //step2a
            vprime_x = vminus_x + ((vminus_y*ttz)-(vminus_z*tty));
            vprime_y = vminus_y + ((vminus_z*ttx)-(vminus_x*ttz));
            vprime_z = vminus_z + ((vminus_x*tty)-(vminus_y*ttx));

            //step2b
            vplus_x = vminus_x + ((vprime_y*ssz)-(vprime_z*ssy));
            vplus_y = vminus_y + ((vprime_z*ssx)-(vprime_x*ssz));
            vplus_z = vminus_z + ((vprime_x*ssy)-(vprime_y*ssx));
             
            vf_x = vplus_x + accel.getX()*atmp;
            vf_y = vplus_y + accel.getY()*atmp;
            vf_z = vplus_z + accel.getZ()*atmp;
            px[ip].setX(vf_x);
            px[ip].setY(vf_y);
            px[ip].setZ(vf_z);
        }
        else if(solver_Boris == 1)
        {
            vminus_x = px[ip].getX() + accel.getX()*atmp;
            vminus_y = px[ip].getY() + accel.getY()*atmp;
            vminus_z = px[ip].getZ() + accel.getZ()*atmp;
                  
            Bn = sqrt(pow(B0.getX(),2) + pow(B0.getY(),2) + pow(B0.getZ(),2));
            if( Bn > 1.e-3)
            {
                bxn = B0.getX()/Bn;
                byn = B0.getY()/Bn;
                bzn = B0.getZ()/Bn;
            }
            theta_n = (q_m*Bn*dteff/gamma);
            vpar =  (bxn*vminus_x)+(byn*vminus_y)+(bzn*vminus_z);
            ssx = vpar*bxn;
            ssy = vpar*byn;
            ssz = vpar*bzn;
                  
            vprime_x = ((vminus_y*bzn)-(vminus_z*byn));
            vprime_y = ((vminus_z*bxn)-(vminus_x*bzn));
            vprime_z = ((vminus_x*byn)-(vminus_y*bxn));
                  
            vplus_x = ssx + ((vminus_x - ssx)*cos(theta_n)) + (vprime_x*sin(theta_n));
            vplus_y = ssy + ((vminus_y - ssy)*cos(theta_n)) + (vprime_y*sin(theta_n));
            vplus_z = ssz + ((vminus_z - ssz)*cos(theta_n)) + (vprime_z*sin(theta_n));
                
            vf_x = vplus_x + accel.getX()*atmp;
            vf_y = vplus_y + accel.getY()*atmp;
            vf_z = vplus_z + accel.getZ()*atmp;
                  
            px[ip].setX(vf_x);
            px[ip].setY(vf_y);
            px[ip].setZ(vf_z);
            
        }
        else
            cout << "error wrong solver_Boris" <<endl;
        
              
        E_p[ip].setX(ex_tmp);
        E_p[ip].setY(ey_tmp);
        E_p[ip].setZ(ez_tmp);
        B_p[ip].setX(bx_tmp);
        B_p[ip].setY(by_tmp);
        B_p[ip].setZ(bz_tmp);
 

    }
}


//==========================================
// Leap Frog 
//
// Periodic boundary condition
// 
//==========================================
void leapfrog
(vector<ThreeVec> &xx, vector<ThreeVec> &px, vector<ThreeVec> &x_old,vector<ThreeVec> &x_new ,vector<double> &qx, int iflag)
{
   double dteff;
   int iy;
//   double xtmp,ytmp,ztmp; 
   double ranf,ranf1;
   double vel_0,del_v_x,del_v_y,del_v_z;
   int count;  // for erasing (electrons)

    if(iflag==0)      { dteff= -1.0*dt;}
    else if(iflag==1) { dteff= dt;}
    else   cout << " Not designed iflag =  " << iflag <<endl;

   for(int ip=0;ip<(int)xx.size();ip++)
   {
      count =0;
        
      //save the old particle position
      x_old[ip].setX(xx[ip].getX());
      x_old[ip].setY(xx[ip].getY());
      x_old[ip].setZ(xx[ip].getZ());

      if(nrel==0) xx[ip] += px[ip]/sqrt(1.0+px[ip].square()/cph/cph)*dteff;
      else        xx[ip] += px[ip]*dteff;
       
       //save the new particle position before BC
       x_new[ip].setX(xx[ip].getX());
       x_new[ip].setY(xx[ip].getY());
       x_new[ip].setZ(xx[ip].getZ());
       

      // Periodic boundary (x)
      if(xx[ip].getX()>=xmax)
      {
           xx[ip].inc(0,-(xmax-xmin));
      }
      else if(xx[ip].getX()<xmin)
      {
           xx[ip].inc(0,+(xmax-xmin));
      }
      // Periodic boundary (y)
      if(xx[ip].getY()>=ymax)
      {
           xx[ip].inc(1,-(ymax-ymin));
      }
      else if(xx[ip].getY()<ymin)
      {
           xx[ip].inc(1,+(ymax-ymin));
      }

      // Periodic boundary (z)
      if(xx[ip].getZ()>=zmax)
      {
           xx[ip].inc(2,-(zmax-zmin));
      }
      else if(xx[ip].getZ()<zmin)
      {
           xx[ip].inc(2,+(zmax-zmin));
      }
       //both x_old,x_new are inside the domain       
      }
}
void current_density_zigzag
(GridPoints ***xgc, vector<ThreeVec> &x_old,vector<ThreeVec> &x_new)
{
   int num_save = 3;
   double *Atemp_local;
   double *Atemp_total;
   int ntotal = Nx*Ny*Nz*num_save;
   Atemp_local = (double*)malloc(Nx*Ny*num_save*sizeof(double));
   Atemp_total = (double*)malloc(Nx*Ny*num_save*sizeof(double));
   
    double xi,xf;
    double yi,yf;
    double zi,zf;
    int i1,j1,k1,i2,j2,k2,il,jl,kl,ill,jll,kll,contsu,ilm,jlm;//,ic=(npx-1)/2,jc=(npy-1)/2;
    double val1x,val2x,val1y,val2y,val1z,val2z,xr,yr,zr;
    double Fx[2],Fy[2],Fz[2],Wx[2],Wy[2],Wz[2];
    int ii[2],jj[2],kk[2];
    double dV=dx*dy;

// Jx(i,j,k), Jy(i,j,k)

    for(int i=0;i<ntotal;i++)
    {
        Atemp_local[i] = 0.0;
        Atemp_total[i] = 0.0;
    }
   for(int i=0;i<Nx;i++)
   {
     for(int j=0;j<Ny;j++)
       {
         for(int k=0;k<Nz;k++)
          {
            xgc[i][j][k].jx=0.0;
            xgc[i][j][k].jy=0.0;
            xgc[i][j][k].jz=0.0;
          }
       }
   }
    
   for(int ip=0;ip<(int)x_old.size();ip++)
    {

        xi = x_old[ip].getX();
        yi = x_old[ip].getY();
        zi = x_old[ip].getZ();
        
        xf = x_new[ip].getX();
        yf = x_new[ip].getY();
        zf = x_new[ip].getZ();
        
        ii[0]=floor(xi/dx); // i1
        jj[0]=floor(yi/dy); // j1
        kk[0]=floor(zi/dz); // k1

        ii[1]=floor(xf/dx);  // i2
        jj[1]=floor(yf/dy);  // j2
        kk[1]=floor(zf/dz);  // k2
        
        val1x=min(ii[0]*dx,ii[1]*dx)+dx;
        val2x=max(max(ii[0]*dx,ii[1]*dx),(xi+xf)/2.0);
        xr=min(val1x,val2x);

        val1y=min(jj[0]*dy,jj[1]*dy)+dy;
        val2y=max(max(jj[0]*dy,jj[1]*dy),(yi+yf)/2.0);
        yr=min(val1y,val2y);

        val1z=min(kk[0]*dz,kk[1]*dz)+dz;
        val2z=max(max(kk[0]*dz,kk[1]*dz),(zi+zf)/2.0);
        zr=min(val1z,val2z);
                 
        Fx[0]=(xr-xi)/dt;   // Fx1
        Fy[0]=(yr-yi)/dt;   // Fy1
        Fz[0]=(zr-zi)/dt;   // Fz1

        Fx[1]=(xf-xr)/dt;   // Fx2
        Fy[1]=(yf-yr)/dt;   // Fy2
        Fz[1]=(zf-zr)/dt;   // Fz2
        

        Wx[0]=((xi+xr)/(2.*dx))-ii[0]; // Wx1
        Wy[0]=((yi+yr)/(2.*dy))-jj[0]; // Wy1
        Wz[0]=((zi+zr)/(2.*dz))-kk[0]; // Wz1

        Wx[1]=((xr+xf)/(2.*dx))-ii[1]; // Wx2
        Wy[1]=((yr+yf)/(2.*dy))-jj[1]; // Wy2
        Wz[1]=((zr+zf)/(2.*dz))-kk[1]; // Wz2
        

        //Periodic Boundary Condition
           il=ii[0];
           jl=jj[0];
           kl=kk[0];
           ill=ii[1];
           jll=jj[1];
           kll=kk[1];
      

            if(ill<0){ill=Nx-1;}
            else if(ill>Nx-1){ill=0;}
            //else{printf("il=%d,jl=%d,ill=%d,jll=%d\n",il,jl,ill,jll);}
        
            if(jll<0){jll=Ny-1;}
            else if(jll>Ny-1){jll=0;}
            //else{printf("il=%d,jl=%d,ill=%d,jll=%d\n",il,jl,ill,jll);}
       
        //Jx
        jlm=jl+1;
        if(jl==Ny-1){jlm = 0;}
           xgc[il][jl][0].jx = xgc[il][jl][0].jx + (Fx[0]*(1.-Wy[0])/dV);
           xgc[il][jlm][0].jx = xgc[il][jlm][0].jx + (Fx[0]*Wy[0]/dV);
        jlm=jll+1;
        if(jll==Ny-1){jlm = 0;}
           xgc[ill][jll][0].jx = xgc[ill][jll][0].jx + (Fx[1]*(1.-Wy[1])/dV);
           xgc[ill][jlm][0].jx = xgc[ill][jlm][0].jx + (Fx[1]*Wy[1]/dV);
        
        //Jy
        ilm=il+1;
        if(il==Nx-1){ilm = 0;}
           xgc[il][jl][0].jy = xgc[il][jl][0].jy + (Fy[0]*(1.-Wx[0])/dV);
           xgc[ilm][jl][0].jy = xgc[ilm][jl][0].jy + (Fy[0]*Wx[0]/dV);

        ilm=ill+1;
        if(ill==Nx-1){ilm = 0;}
           xgc[ill][jll][0].jy = xgc[ill][jll][0].jy + (Fy[1]*(1.-Wx[1])/dV);
           xgc[ilm][jll][0].jy = xgc[ilm][jll][0].jy + (Fy[1]*Wx[1]/dV);
    }

    int in =0;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             Atemp_local[in*Nx*Ny+iy*Nx+ix] = xgc[ix][iy][0].jx;
         }
     }
     ++in;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             Atemp_local[in*Nx*Ny+iy*Nx+ix] = xgc[ix][iy][0].jy;
             
         }
     }
    ++in;
     for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             Atemp_local[in*Nx*Ny+iy*Nx+ix] = xgc[ix][iy][0].jz;
         }
     }
    ++in;
    
    // check consistency
    if(in!=num_save) cout << " wrong num_save " << num_save << " != " << in <<endl;

    // Reduce to root
    MPI_Allreduce(Atemp_local,Atemp_total,ntotal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    in = 0;
    for(int iy=0;iy<Ny;iy++)
    {
        for(int ix=0;ix<Nx;ix++)
        {
             xgc[ix][iy][0].jx =  Atemp_total[in*Nx*Ny+iy*Nx+ix];
        }
    }
    ++in;
    for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
              xgc[ix][iy][0].jy =  Atemp_total[in*Nx*Ny+iy*Nx+ix];
         }
     }
    ++in;
    for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
              xgc[ix][iy][0].jz =  Atemp_total[in*Nx*Ny+iy*Nx+ix];
         }
     }
    ++in;
    
    for(int iy=0;iy<Ny;iy++)
     {
         for(int ix=0;ix<Nx;ix++)
         {
             xgc[ix][iy][0].jz = xgc[ix][iy][0].uz_i*xgc[ix][iy][0].den_i;
         }
     }
        

    free(Atemp_local);
    free(Atemp_total);

 }
