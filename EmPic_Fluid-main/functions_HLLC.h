#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Alloc.h"
//******************************************************************************************************************************************************
//
//                                                                   HLLC FUNCTIONS
//
//******************************************************************************************************************************************************
void ascii_outputfluid(GridFluid ***gf,int tmp, int my_id, int specie)
{
        FILE *fb,*fb1,*fb2,*fb3,*fb4;
        FILE *rhof,*phif, *Exf, *Eyf, *Ezf, *Bxf, *Byf, *Bzf ;
        char filenam[50], filenam1[50],filenam2[50], filenam3[50],filenam4[50], filenam5[50];
        char filenamew[50], filenamew1[50], filenamew2[50],filenamew3[50],filenamew4[50];
/*
sprintf(filenamew, "./results/fields/phi_%4.4d_%d.txt",tmp,my_id);
                       phif = fopen(filenamew, "wb");
                       if (phif == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
sprintf(filenamew1, "./results/fields/Ex_%4.4d_%d.txt",tmp,my_id);
                       Exf = fopen(filenamew1, "wb");
                       if (Exf == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
sprintf(filenamew2, "./results/fields/Ey_%4.4d_%d.txt",tmp,my_id);
                       Eyf = fopen(filenamew2, "wb");
                       if (Eyf == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
sprintf(filenamew3, "./results/fields/rho_%4.4d_%d.txt",tmp,my_id);
                       rhof = fopen(filenamew3, "wb");
                       if (rhof == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
sprintf(filenamew4, "./results/fields/Br_%4.4d_%d.txt",tmp,my_id);
                       Bf = fopen(filenamew4, "wb");
                       if (Bf == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
sprintf(filenam5, "./results/fields/S_%4.4d_%d.txt",tmp,my_id);
                       Sf = fopen(filenam5, "wb");
                       if (Sf == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
*/
    if(specie==0)
    {
                       sprintf(filenam, "./output/results/ion/density_%4.4d.txt",tmp);
                       fb = fopen(filenam, "wb");
                       if (fb == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
                        sprintf(filenam1, "./output/results/ion/velocity_x_%4.4d.txt",tmp);
                       fb1 = fopen(filenam1, "wb");
                       if (fb1 == NULL)
                        {
                          printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
                        sprintf(filenam4, "./output/results/ion/velocity_y_%4.4d.txt",tmp);
                       fb4 = fopen(filenam4, "wb");
                       if (fb4 == NULL)
                        {
                          printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }

                         sprintf(filenam2, "./output/results/ion/pressure_%4.4d.txt",tmp);
                       fb2 = fopen(filenam2, "wb");
                       if (fb2 == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
    }
    else if(specie==1)
    {
                       sprintf(filenam, "./output/results/electron/density_%4.4d.txt",tmp);
                       fb = fopen(filenam, "wb");
                       if (fb == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
                        sprintf(filenam1, "./output/results/electron/velocity_x_%4.4d.txt",tmp);
                       fb1 = fopen(filenam1, "wb");
                       if (fb1 == NULL)
                        {
                          printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }
                        sprintf(filenam4, "./output/results/electron/velocity_y_%4.4d.txt",tmp);
                       fb4 = fopen(filenam4, "wb");
                       if (fb4 == NULL)
                        {
                          printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }

                         sprintf(filenam2, "./output/results/electron/pressure_%4.4d.txt",tmp);
                       fb2 = fopen(filenam2, "wb");
                       if (fb2 == NULL)
                        {
                           printf("I couldn't open fluid_variables_...._.txt for writing.\n");
                           exit(0);
                        }

    }
             for(int i=nbuffer;i<nlocal.x+(1*nbuffer);i++)
                        {
                        for(int j=nbuffer;j<nlocal.y+(1*nbuffer);j++)
                        {
                               fprintf(fb, "%e\t",gf[i][j][0].n_i);
                        }//
                        }
                              fclose(fb);

                        for(int i=nbuffer;i<nlocal.x+(1*nbuffer);i++)
                        {
                        for(int j=nbuffer;j<nlocal.y+(1*nbuffer);j++)
                               {
                                        fprintf(fb1, "%e\t",gf[i][j][0].u_ix);
                                }//
                                }
                              fclose(fb1);

                        for(int i=nbuffer;i<nlocal.x+(1*nbuffer);i++)
                        {
                        for(int j=nbuffer;j<nlocal.y+(1*nbuffer);j++)
                               {

                                fprintf(fb4, "%e\t",gf[i][j][0].u_iy);
                                }//
                               }
                              fclose(fb4);
                        for(int i=nbuffer;i<nlocal.x+(1*nbuffer);i++)
                        {
                        for(int j=nbuffer;j<nlocal.y+(1*nbuffer);j++)
                               {
                                    fprintf(fb2, "%e\t",gf[i][j][0].P_i);
                                }//
                                }
                              fclose(fb2);
/*
                      for(int i=0;i<nlocal_nodes.x+(2*nbuffer);i++)
                        {
                        for(int j=0;j<nlocal_nodes.y+(2*nbuffer);j++)

                            {
                                fprintf(phif, "%e\t",phi[i][j]);
                            }//
                         }
                         fclose(phif);


                 for(int i=nbuffer;i<nlocal.x+(1*nbuffer);i++)
                 {
                   for(int j=nbuffer;j<nlocal.y+(1*nbuffer);j++)
                             {
                                fprintf(Exf, "%e\t",Ex[i][j]);
                                }//
                               }
                              fclose(Exf);


                          for(int i=nbuffer;i<nlocal.x+(1*nbuffer);i++)
                          {
                            for(int j=nbuffer;j<nlocal.y+(1*nbuffer);j++)
                              {
                                fprintf(Eyf, "%e\t",Ey[i][j]);
                               }//
                          }
                           fclose(Eyf);


                         for(int i=nbuffer;i<nlocal.x+(1*nbuffer);i++)
                          {
                            for(int j=nbuffer;j<nlocal.y+(1*nbuffer);j++)
                               {
                                fprintf(Bf, "%e\t",Br[i][j]);
                                }//
                                }
                              fclose(Bf);

                           for(int i=nbuffer;i<nlocal.x+(1*nbuffer);i++)
                          {
                            for(int j=nbuffer;j<nlocal.y+(1*nbuffer);j++)
                            {
                                fprintf(Sf, "%e\t",S[i][j]);
                                }//
                                }
                              fclose(Sf);


                        for(int i=0;i<nlocal_nodes.x+(2*nbuffer);i++)
                          {
                        for(int j=0;j<nlocal_nodes.y+(2*nbuffer);j++)
                              {
                                fprintf(rhof, "%e\t",rho[i][j]);
                                }//
                                }
                              fclose(rhof);
*/
}
double power(double x,double a)
{
    return(exp(a*log(x)));
}
double minn(double m, double n)
{
    double k=0;
    if (m>n)
    {
        k=n;
    }
    else {k=m;}
    return k;
}
double maxn(double m, double n)
{
    double k=0;
    if (m>n)
    {
        k=m;
    }
    else {k=n;}
    return k;
    
}
double Wave_speed_isentropic_estimates(double U_l, double U_r, double A_l, double A_r, int lor)
{
    double a_star, u_star;
    double S_l, S_r;
    
    a_star = 0.5*(A_r+A_l) + ((1.0/4.0)*(gam-1.0)*(U_l-U_r));
    u_star = 0.5*(U_l-U_r) + ((A_l-A_r)/(gam-1.0));
    S_l = min(U_l - A_l, (u_star - a_star));
    S_r = max(U_l + A_l, (u_star + a_star));
    
    if(lor==0){return(S_l);}
    else if(lor==1){return(S_r);}
    else {return -1;}
}
double pressure_estimate_acustic(double n_l,double n_r,double U_l, double U_r,double Pr_l,double Pr_r, double A_l, double A_r, double mtmp)
{
    double p_star, u_star, n_bar, a_bar;
    double S_l, S_r,Ppvrs;
    
    n_bar = 0.5*(n_r+n_l)*mtmp;
    a_bar = 0.5*(A_r+A_l);
    Ppvrs = 0.5*(Pr_l+Pr_r) - (0.5*(U_r-U_l)*a_bar*n_bar);
    p_star = max(0.0,Ppvrs);
    return(p_star);
}
double Average_pressure(double n_l,double n_r,double U_l, double U_r,double Pr_l,double Pr_r, double S_l, double S_r,double S_m,double mtmp)
{
    double p_lr;
    p_lr = 0.5*(Pr_l + Pr_r + mtmp*(n_l*(S_l - U_l)*(S_m - U_l)) + mtmp*(n_r*(S_r-U_r)*(S_m-U_r)));
    return(p_lr);
}

double Wave_speed_linear_estimates(double Rho_l,double Rho_r,double U_l, double U_r,double Pr_l,double Pr_r, double A_l, double A_r, int lor)
{
    double p_star=0, u_star=0, rho_l_star=0, rho_r_star=0, a_bar=0, rho_bar=0, a_l_star=0, a_r_star=0;
    double S_l=0, S_r=0;
    
    rho_bar = 0.5*(Rho_r+Rho_l);
    a_bar = 0.5*(A_r+A_l);
    p_star = 0.5*(Pr_l+Pr_r) - (0.5*(U_r-U_l)*rho_bar*a_bar);
    if(a_bar==0){p_star = 0.5*(Pr_l+Pr_r);}
    u_star = 0.5*(U_r+U_l) - (0.5*(Pr_r-Pr_l)/(rho_bar*a_bar));
    if(a_bar==0.0){ u_star = 0.5*(U_r+U_l);}
    rho_l_star = Rho_l + (((U_l - u_star)*rho_bar)/a_bar);
    if(a_bar==0.0){rho_l_star = Rho_l;}
    rho_r_star = Rho_r + (((u_star - U_r)*rho_bar)/a_bar);
    if(a_bar==0.0){rho_r_star = Rho_r;}
    a_l_star = pow(gam*p_star/rho_l_star,0.5);
    if(p_star==0.0){a_l_star =0.0;}
    a_r_star = pow(gam*p_star/rho_r_star,0.5);
    if(p_star==0.0){a_r_star =0.0;}
    S_l = minn((U_l - A_l), (u_star - a_l_star));
    S_r = max((U_r + A_r), (u_star + a_r_star));
    
    if(lor==0.0){return(S_l);}
    if(lor==1.0){return(S_r);}
}
double Wave_speed_Pressure(double n_l, double n_r,double U_l, double U_r,double Pr_l,double Pr_r, double A_l,double A_r, double Prl,int lor, double mtmp)
{
    double S_l, S_r,S_m,q_l,q_r, H_l,H_r;
    double p_star, u_star, num,den;;
    double z = 0.5*((gam+1.0)/gam);
    
    H_l = Prl/Pr_l;
    H_r = Prl/Pr_r;
    if(Prl <=Pr_l)
    {
        q_l = 1.0;
    }
    else
        q_l = pow((1.0+ (z*(H_l-1.0))),0.5);
    if(Prl <= Pr_r)
    {
        q_r = 1.0;
    }
    else
    {
        q_r = pow((1.0+ (z*(H_r-1.0))),0.5);
    }
    S_l = U_l - (A_l*q_l);
    S_r = U_r + (A_r*q_r);
    num = (Pr_r - Pr_l + mtmp*(n_l*U_l*(S_l - U_l)) - mtmp*(n_r*U_r*(S_r - U_r)));
    den = (mtmp*n_l*(S_l - U_l)) - (mtmp*n_r*(S_r - U_r));
    u_star = num/den;
    S_m = u_star;
    if(lor==0){return(S_l);}
    else if(lor==1){return(S_r);}
    else if(lor==2){return(S_m);}
    else {return -999999;}
}
double U_two_fluid_2D(double n,double Ux,double Uy, double Pr, double m, double q,int nmoment)
{
    double Fm[5];
    double U1, U2, U3,U4, F1, F2, F3,F4;

    U1 = n; //n   
    U2 = n*Ux; //nux 
    U3 = n*Uy; //nuy
    U4 = (Pr/(gam-1.0) + 0.5*m*n*((Ux*Ux)+(Uy*Uy)));  //e
    Fm[0] = U1; //n
    Fm[1] = U2; //nux
    Fm[2] = U3; //nuy
    Fm[3] = U4; //e
    return (Fm[nmoment]);
}
double F_two_fluid_2D(double q, double n,double Ux,double Uy, double Pr, double m ,int nmoment)
{
    double Fm[5];
    double  U2, U3, U4, F1, F2, F3, F4;

    U2 = (n*Ux);  //nux
    U3 = (Uy);    //nuy
    U4 = (Pr/(gam-1.0) + 0.5*m*n*((Ux*Ux)+(Uy*Uy))); //e 

    F1 = U2;    //nux
    F2 = (((U2*U2)/n) + (Pr/m)); //nux^2 + P/m  
    F3= (U2*U3); // n ux uy
    F4 = (((U4 + Pr)*U2)/(n));// (e + P) vx        
    Fm[0] = F1;
    Fm[1] = F2;
    Fm[2] = F3;
    Fm[3] = F4;
     return (Fm[nmoment]);
}
double G_two_fluid_2D(double q, double n,double Ux,double Uy, double Pr, double m ,int nmoment)
{
    double Fm[5];
    double  U2, U3, U4, F1, F2, F3, F4;

    U2 = (n*Uy); 
    U3 = (Ux); 
    U4 = (Pr/(gam-1.0) + 0.5*m*n*((Ux*Ux)+(Uy*Uy)));  

    F1 = U2;  //n uy   
    F2= (U2*U3);  //n ux uy
    F3 = (((U2*U2)/n) + ((Pr/m)));//  nuyuy + P/m 
    F4 = (((U4 + Pr)*U2)/(n));     // (e + P)uy   
    Fm[0] = F1;
    Fm[1] = F2;
    Fm[2] = F3;
    Fm[3] = F4;
    return (Fm[nmoment]);
}
double U_inter(double U_l, double F_l, double S_l,double S_m, double Pr_inter,int  nmoment,double mtmp, double dir)
{
    double U,D;
    
    if(dir==0)
    {
        
        if(nmoment==0){
            D = 0.0;
        }
        if(nmoment==1){
            D = 1.0/mtmp;
        }
        if(nmoment==2){
            D = 0.0;
        }
        if(nmoment==3){
            D = S_m;
        }
    }
    if(dir==1)
    {
        if(nmoment==0){
            D = 0.0;
        }
        if(nmoment==1){
            D = 0.0;
        }
        if(nmoment==2){
            D = 1.0/mtmp;
        }
        if(nmoment==3){
            D = S_m;
        }
    }
    
    U = (((S_l*U_l) - F_l + (Pr_inter*D))/(S_l - S_m));
    return (U);
}
double Flux_inter(double U_l, double F_l, double S_l,double S_m, double Pr_inter,int  nmoment, double mtmp,int dir)
{
    double F,D;
    
    if(dir==0)
    {
        if(nmoment==0){
            D = 0.0;
        }
        if(nmoment==1){
            D = 1.0/mtmp;
        }
        if(nmoment==2){
            D = 0.0;
        }
        if(nmoment==3){
            D = S_m;
        }
    }
    if(dir==1)
    {
        if(nmoment==0){
            D = 0.0;
        }
        if(nmoment==1){
            D = 0.0;
        }
        if(nmoment==2){
            D = 1.0/mtmp;
        }
        if(nmoment==3){
            D = S_m;
        }
    }
    F = ( ((S_m*(S_l*U_l - F_l)) + (S_l*Pr_inter*D))/(S_l - S_m));
    return (F);

}
double Flux_HLLC(double F_r, double F_l,double F_star_r, double F_star_l, double S_l, double S_r, double S_m)
{
    double F_HLLC=0.0;
    
    if(S_l>=0.0)
    {
        F_HLLC = F_l;
    }
    if(S_r<=0.0)
    {
        F_HLLC = F_r;
    }
    if(S_l<=0.0 && S_m>=0.0)
    {
        F_HLLC = F_star_l;
    }
    if(S_m<=0.0 && S_r>=0.0)
    {
        F_HLLC = F_star_r;
    }
    return(F_HLLC);
    
    
}
double Flux_HLL(double U_r, double U_l, double F_r, double F_l, double S_r, double S_l)
{
    double F_HLL=0.0, F_star=0.0;
    
    F_star = (S_r*F_l - S_l*F_r + (S_l*S_r*(U_r - U_l)))/(S_r - S_l);
    
    if(S_l < 0.0 && S_r>0.0)
    {
        F_HLL = F_star;
    }
    if(S_l>=0.0)
    {
        F_HLL = F_l;
    }
    if(S_r<=0.0)
    {
        F_HLL = F_r;
    }
    return(F_HLL);
}
double U_HLL(double ** U_r, double ** U_l, double ** F_r, double ** F_l, double ** S_r, double ** S_l, double **& U_HLL)
{

double ** U_tmp = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_tmp[i]= new double[ntot.y];

double ** U_star = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star[i]= new double[ntot.y];

 
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            U_tmp[i][j] = 0.0;
            U_star[i][j] = 0.0;
        }
    }
    
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            U_star[i][j] = (S_r[i][j]*U_r[i][j] - S_l[i][j]*U_l[i][j] + ((F_l[i][j] - F_r[i][j])))/(S_r[i][j] - S_l[i][j]);
            
            if(S_l[i][j] < 0.0 && S_r[i][j]>0.0)
            {
                U_tmp[i][j] = U_star[i][j];
            }
            if(S_l[i][j]>=0.0)
            {
                U_tmp[i][j] = U_l[i][j];
            }
            if(S_r[i][j]<=0.0)
            {
                U_tmp[i][j] = U_r[i][j];
            }
        }
    }
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            U_HLL[i][j] = U_tmp[i][j];
        }
    }

for (int i = 0; i < ntot.x; i++)
delete[] U_tmp[i] ;
delete[] U_tmp;

for (int i = 0; i < ntot.x; i++)
delete[] U_star[i] ;
delete[] U_star;
}
double F_trans_HLL(double ** U_1, double ** U_2, double ** U_3, double ** U_4, double ** G_1, double ** G_2, double ** G_3,double ** G_4, double **& F_trans,int mtmp)
{

double ** F_tmp = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_tmp[i]= new double[ntot.y];

    int i,j;
    
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            F_tmp[i][j] = 0.0;
        }
    }
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            if(mtmp == 0)
            {
                F_tmp[i][j] = U_2[i][j];
            }
            if(mtmp == 1)
            {
                F_tmp[i][j] = G_3[i][j] + (((U_2[i][j]*U_2[i][j]) - (U_3[i][j]*U_3[i][j]))/U_1[i][j]);
                if(U_1[i][j] == 0)
                {
                    F_tmp[i][j] = G_3[i][j];
                }
            }
            if(mtmp == 2)
            {
                F_tmp[i][j] = (U_3[i][j]*U_2[i][j])/U_1[i][j];
                if(U_1[i][j] == 0)
                {
                    F_tmp[i][j] = 0.0;
                }
            }
            if(mtmp == 3)
            {
                F_tmp[i][j] = (U_2[i][j]*G_4[i][j])/U_3[i][j];
                if(U_3[i][j] == 0)
                {
                    F_tmp[i][j] = 0.0;
                }
            }
            F_trans[i][j] = F_tmp[i][j];
        }
    }
for (int i = 0; i < ntot.x; i++)
delete[] F_tmp[i] ;
delete[] F_tmp;

}
double G_trans_HLL(double ** U_1, double ** U_2, double ** U_3, double ** U_4, double ** F_1, double ** F_2, double ** F_3,double ** F_4, double **& G_trans,int mtmp)
{
    double ** G_tmp = new double *[ntot.x];
    for (int i = 0; i < ntot.x; i++)
    G_tmp[i]= new double[ntot.y];

    int i,j;
    
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            G_tmp[i][j] = 0.0;
        }
    }
    
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            if(mtmp == 0)
            {
                G_tmp[i][j] = U_3[i][j];
            }
            if(mtmp == 1)
            {
                G_tmp[i][j] = U_2[i][j]*U_3[i][j]/U_1[i][j];
                if(U_1[i][j] == 0)
                {
                    G_tmp[i][j] = 0.0;
                }
            }
            if(mtmp == 2)
            {
                G_tmp[i][j] = F_2[i][j] + (((U_3[i][j]*U_3[i][j]) - (U_2[i][j]*U_2[i][j]))/U_1[i][j]);
                if(U_1[i][j] == 0)
                {
                    G_tmp[i][j] = F_2[i][j];
                }
            }
            if(mtmp == 3)
            {
                G_tmp[i][j] = U_3[i][j]*F_4[i][j]/U_2[i][j];
                if(U_2[i][j] == 0)
                {
                    G_tmp[i][j] = 0.0;
                }
            }
            G_trans[i][j] = G_tmp[i][j];
        }
    }

	for (int i = 0; i < ntot.x; i++)
	delete[] G_tmp[i] ;
	delete[] G_tmp;
}
double H_fun(double delta1, double delta2, double & h)
{
    h = (1.0/3.0)*(2.0*delta2 + delta1);
}

double HL_fun(double delta1, double delta2, double & hl)
{
    double signo, h3,m1,m2,m3,aux1,aux2,aux3,aux4;
    signo = copysign(1.0,delta2);
    H_fun(delta1,delta2,h3);
    m1 = signo*h3;
    m2 = signo*delta1;
    m3 = 1.5*fabs(delta2);
    aux1 = minn(2.0*m2,m1);
    aux2 = max(aux1,m3);
    aux3 = minn(m1,aux2);
    aux4 = max(0.0,aux3);
    hl = signo*aux4;
}
double limiter_2D(double ** n0,double ** ux0,double ** uy0, double ** P0, double **& np,double **& nm, double **& uxp,double **& uxm,double **& uyp,double **& uym,double **& Pp,double **& Pm, int variable_type, int limiter_type,int dir, int bc_x, int bc_y, int my_id, bool notWest, bool notEast,bool notSouth,bool notNorth,int specie,double *flux_cath)
{
   int i,j; 
    double ***q_p = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        q_p[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            q_p[i][j] = new double  [4]();
    }
    double ***q_m = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        q_m[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            q_m[i][j] = new double  [4]();
    }
    double ***auxpxp = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        auxpxp[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            auxpxp[i][j] = new double  [4]();
    }
    double ***auxpxm = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        auxpxm[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            auxpxm[i][j] = new double  [4]();
    }
    double ***auxpyp = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        auxpyp[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            auxpyp[i][j] = new double  [4]();
    }
    double ***auxpym = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        auxpym[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            auxpym[i][j] = new double  [4]();
    }
    
     double ***q = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        q[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            q[i][j] = new double  [4]();
    }
    double ***sx = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        sx[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            sx[i][j] = new double  [4]();
    }
    double ***sy = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        sy[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            sy[i][j] = new double  [4]();
    }
    double ***sx1 = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        sx1[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            sx1[i][j] = new double  [4]();
    }
    double ***sy1 = new double **[ntot.x]();
    for (i = 0; i < ntot.x; i++)
    {
        sy1[i] = new double  *[ntot.y]();
        for (j = 0; j < ntot.y; j++)
            sy1[i][j] = new double  [4]();
    }






    double rx=0.0, ry=0.0,rx1=0.0, ry1=0.0, s_a=0.0,s_b=0.0;
    double delta_m_x=0,delta_p_x=0.0, delta_m_y=0.0,delta_p_y=0.0;
    double da=0.0,db=0.0,dc=0.0,dd=0.0;
    double p,px1=0.0,py1=0.0,theta;
    double ql =0.0;
 
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            for(int k=0;k<3;k++)
            {
                q_p[i][j][k] = 0.0;
                q_m[i][j][k] = 0.0;
                auxpxp[i][j][k] = 0.0;
                auxpxm[i][j][k] = 0.0;
                auxpyp[i][j][k] = 0.0;
                auxpym[i][j][k] = 0.0;
                q[i][j][k] = 0.0;
                sx[i][j][k] = 0.0;
                sy[i][j][k] = 0.0;
                sx1[i][j][k] = 0.0;
                sy1[i][j][k] = 0.0;
            }
        }
    }
    
    // define conserved variables in each slice of 3D matrix q
    
    for(int j=0;j<ntot.y;j++)
    {
        for(int i=0;i<ntot.x;i++)
        {
            q[i][j][0] = n0[i][j];
            q[i][j][1] = ux0[i][j];
            q[i][j][2] = uy0[i][j];
            q[i][j][3] = P0[i][j];
        }
    }
    //.....limiters for TVD scheme in both directions.......//
         for (int j=nbuffer-1;j<=ntot.y-nbuffer; j++)
            {
              for (int i=nbuffer-1;i<=ntot.x-nbuffer; i++)
	       {
		for(int k=0;k<4;k++)
 	 	          {
                delta_m_x = q[i][j][k] - q[i-1][j][k];
                delta_p_x = q[i+1][j][k] - q[i][j][k];
                delta_m_y = q[i][j][k] - q[i][j-1][k];
                delta_p_y = q[i][j+1][k] - q[i][j][k];
                
                if (fabs(q[i][j][k]-q[i][j+1][k])<epsil) ry=0.0;
                else
                {
                    ry=delta_m_y/delta_p_y;
                }
                if (fabs(q[i][j][k]-q[i+1][j][k])<epsil) rx=0.0;
                else
                {
                    rx=delta_m_x/delta_p_x;
                }
                if (fabs(q[i][j-1][k]-q[i][j][k])<epsil) ry=0.0;
                else
                {
                    ry1=delta_p_y/delta_m_y;
                }
                if (fabs(q[i-1][j][k]-q[i][j][k])<epsil) rx=0.0;
                else
                {
                    rx1=delta_p_x/delta_m_x;
                }
                
                if (limiter_type==0) //superbee
                {
                    s_a = min(2.0*rx,1.0);
                    s_b = minn(2.0,rx);
                    sx[i][j][k] = max(s_a,s_b);
                    sx[i][j][k] = max(0,sx[i][j][k]);
                    s_a = min(2.0*ry,1.0);
                    s_b = min(2.0,ry);
                    sy[i][j][k] = max(s_a,s_b);
                    sy[i][j][k] = max(0,sy[i][j][k]);
                }
                else if (limiter_type == 1) //van-leer
                {
                    sx[i][j][k] = (rx+fabs(rx))/(1.0+fabs(rx));
                    sy[i][j][k] = (ry+fabs(ry))/(1.0+fabs(ry));
                }
                else if (limiter_type == 2)  // symmetric Koren
                {
                    if (ux0[i][j] < 0 && rx!=0) rx=1/rx;
                    s_a = min((2.0*rx+1.0)/3.0,2.0);
                    s_a = min(s_a,2.0*rx);
                    sx[i][j][k] = max(0.0,s_a);
                    if (uy0[i][j] < 0 && ry!=0) ry=1/ry;
                    s_a = min((2.0*ry+1.0)/3.0,2.0);
                    s_a = min(s_a,2.0*ry);
                    sy[i][j][k] = max(0.0,s_a);
                }
                else if (limiter_type == 3)                    // asymmetric Koren
                {
                    s_a = min((2.0*rx+1.0)/3.0,2.0);
                    s_a = min(s_a,2.0*rx);
                    sx[i][j][k] = max(0.0,s_a);
                    s_a = min((2.0*ry+1.0)/3.0,2.0);
                    s_a = min(s_a,2.0*ry);
                    sy[i][j][k] = max(0.0,s_a);
                }
                else if (limiter_type == 5)                                          //Logaritic
                {
                    px1 = 2.0*((pow(fabs(rx),ql))/(1.0+(pow(fabs(rx),2.0*ql))));
                    s_a = (((pow(px1,2))-1.0)*(pow((px1-1.0),2)));
                    sx[i][j][k] = ((2.0*px1)*( (((pow(px1,2) - (2.0*px1*rx) + 1.0))*log(px1)) - ((1-rx)*((pow(px1,2))-1.0))))/s_a;
                    px1 = 2.0*((pow(fabs(rx1),ql))/(1.0+(pow(fabs(rx1),2.0*ql))));
                    s_a = (((pow(px1,2))-1.0)*(pow((px1-1.0),2)));
                    sx1[i][j][k] = ((2.0*px1)*( (((pow(px1,2) - (2.0*px1*rx1) + 1.0))*log(px1)) - ((1-rx1)*((pow(px1,2))-1.0))))/s_a;
                    py1 = 2.0*((pow(fabs(ry),ql))/(1.0+(pow(fabs(ry),2.0*ql))));
                    s_a = (((pow(py1,2))-1.0)*(pow((py1-1.0),2)));
                    sy[i][j][k] = ((2.0*py1)*( (((pow(py1,2) - (2.0*py1*ry) + 1.0))*log(py1)) - ((1-ry)*((pow(py1,2))-1.0))))/s_a;
                    py1 = 2.0*((pow(fabs(ry1),ql))/(1.0+(pow(fabs(ry1),2.0*ql))));
                    s_a = (((pow(py1,2))-1.0)*(pow((py1-1.0),2)));
                    sy1[i][j][k] = ((2.0*py1)*( (((pow(py1,2) - (2.0*py1*ry1) + 1.0))*log(py1)) - ((1-ry1)*((pow(py1,2))-1.0))))/s_a;
                }
                
            }
        }
    }
    
    
    
    //.....Reconstruct variables at cell edges......//
    //   q_p[0] is not defined yet
    //   q_m[Ntot-nbuffer-1] is not defined yet
    if(limiter_type != 4 || limiter_type != 5)
    {
	 for (int j=nbuffer-1;j<=ntot.y-nbuffer; j++)
        {
            for (int i=nbuffer-1;i<=ntot.x-nbuffer; i++)
            {
                for (int k=0;k<4;k++)
                {
                    if(dir == 0)
                    {
                        if ((ux0[i][j]<0) && (limiter_type == 2))
                        {
                            q_p[i][j][k]  =   q[i][j][k] + 0.5*sx[i][j][k]*(q[i][j][k]-q[i-1][j][k]);
                            q_m[i-1][j][k]  =   q[i][j][k] - 0.5*sx[i][j][k]*(q[i][j][k]-q[i-1][j][k]);
                        }
                        else
                        {
                            q_p[i][j][k]  =   q[i][j][k] + 0.5*sx[i][j][k]*(q[i+1][j][k]-q[i][j][k]);
                            q_m[i-1][j][k]  =   q[i][j][k] - 0.5*sx[i][j][k]*(q[i+1][j][k]-q[i][j][k]);
                        }
                    }
                    if(dir==1)
                    {
                        if ((uy0[i][j]<0) && (limiter_type == 2))
                        {
                            q_p[i][j][k]  =   q[i][j][k] + 0.5*sy[i][j][k]*(q[i][j][k]-q[i][j-1][k]);
                            q_m[i][j-1][k]  =   q[i][j][k] - 0.5*sy[i][j][k]*(q[i][j][k]-q[i][j-1][k]);
                        }
                        else
                        {
                            q_p[i][j][k]  =   q[i][j][k] + 0.5*sy[i][j][k]*(q[i][j+1][k]-q[i][j][k]);
                            q_m[i][j-1][k]  =   q[i][j][k] - 0.5*sy[i][j][k]*(q[i][j+1][k]-q[i][j][k]);
                        }
                    }
                }
            }
        }
    }
    if(limiter_type == 4 || limiter_type == 5)
    {
	 for (int j=nbuffer-1;j<=ntot.y-nbuffer; j++)
        {
            for (int i=nbuffer-1;i<=ntot.x-nbuffer; i++)
            {
                for (int k=0;k<4;k++)
                {
                    if(limiter_type == 4)
                    {
                        if(dir == 0)
                        {
                            q_p[i][j][k]  =   (5.0/6.0)*q[i][j][k] -  (1.0/6.0)*q[i-1][j][k] + (1.0/3.0)*q[i+1][j][k];
                            q_m[i-1][j][k]  =   (5.0/6.0)*q[i][j][k] -  (1.0/6.0)*q[i+1][j][k] + (1.0/3.0)*q[i-1][j][k];
                        }
                        if(dir == 1)
                        {
                            q_p[i][j][k]  =   (5.0/6.0)*q[i][j][k] -  (1.0/6.0)*q[i][j-1][k] + (1.0/3.0)*q[i][j+1][k];
                            q_m[i][j-1][k]  =   (5.0/6.0)*q[i][j][k] -  (1.0/6.0)*q[i][j+1][k] + (1.0/3.0)*q[i][j-1][k];
                        }
                    }
                    if(limiter_type == 5)
                    {
                        if(dir == 0)
                        {
                            if ((ux0[i][j]<0) && (limiter_type == 2))
                            {
                                q_p[i][j][k]  =   q[i][j][k] + 0.5*sx[i][j][k]*(q[i][j][k]-q[i-1][j][k]);
                                q_m[i-1][j][k]  =   q[i][j][k] - 0.5*sx1[i][j][k]*(q[i+1][j][k]-q[i][j][k]);
                            }
                            else
                            {
                                q_p[i][j][k]  =   q[i][j][k] + 0.5*sx[i][j][k]*(q[i+1][j][k]-q[i][j][k]);
                                q_m[i-1][j][k]  =   q[i][j][k] - 0.5*sx1[i][j][k]*(q[i][j][k]-q[i-1][j][k]);
                            }
                        }
                        if(dir==1)
                        {
                            if ((uy0[i][j]<0) && (limiter_type == 2))
                            {
                                q_p[i][j][k]  =   q[i][j][k] + 0.5*sy[i][j][k]*(q[i][j][k]-q[i][j-1][k]);
                                q_m[i][j-1][k]  =   q[i][j][k] - 0.5*sy1[i][j][k]*(q[i][j+1][k]-q[i][j][k]);
                            }
                            else
                            {
                                q_p[i][j][k]  =   q[i][j][k] + 0.5*sy[i][j][k]*(q[i][j+1][k]-q[i][j][k]);
                                q_m[i][j-1][k]  =   q[i][j][k] - 0.5*sy1[i][j][k]*(q[i][j][k]-q[i][j-1][k]);
                            }
                        }
                    }
                }
            }
        }
    }
    //.......Calculating slope limited primitive variables at cell edges.....//
    for(int j=0;j<ntot.y;j++)
    {
        for(int i=0;i<ntot.x;i++)
        {
            np[i][j]=max(epsil,q_p[i][j][0]);
            nm[i][j]=max(epsil,q_m[i][j][0]);
        }
    }
    
    for(int j=0;j<ntot.y;j++)
    {
        for(int i=0;i<ntot.x;i++)
        {
            uxp[i][j] = q_p[i][j][1];
            uxm[i][j] = q_m[i][j][1];
            uyp[i][j] = q_p[i][j][2];
            uym[i][j] = q_m[i][j][2];
            Pp[i][j]=max(0,q_p[i][j][3]);
            Pm[i][j]=max(0,q_m[i][j][3]);
        }
    }


for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] q_p[i][j];
    delete[] q_p[i];
}
delete[] q_p;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] q_m[i][j];
    delete[] q_m[i];
}
delete[] q_m;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] auxpxp[i][j];
    delete[] auxpxp[i];
}
delete[] auxpxp;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] auxpxm[i][j];
    delete[] auxpxm[i];
}
delete[] auxpxm;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] auxpyp[i][j];
    delete[] auxpyp[i];
}
delete[] auxpyp;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] auxpym[i][j];
    delete[] auxpym[i];
}
delete[] auxpym;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] q[i][j];
    delete[] q[i];
}
delete[] q;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] sx[i][j];
    delete[] sx[i];
}
delete[] sx;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] sy[i][j];
    delete[] sy[i];
}
delete[] sy;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] sx1[i][j];
    delete[] sx1[i];
}
delete[] sx1;

for (i = 0; i < ntot.x; i++)
{
    for (j = 0; j < ntot.y; j++)
    delete[] sy1[i][j];
    delete[] sy1[i];
}
delete[] sy1;
}
double Wave_speed_linea_stimates_2D(double n_m, double u_m, double p_m,double n_p, double u_p, double p_p, double & Sl, double & Sr, double & Sm)
{
    double c_p=0, c_m=0;
    c_m = pow(gam*p_m/n_m,0.5);
    c_p = pow(gam*p_p/n_p,0.5);
    Sl = Wave_speed_linear_estimates(n_m,n_p,u_m,u_p,p_m,p_p,c_m,c_p,0);
    Sr = Wave_speed_linear_estimates(n_m,n_p,u_m,u_p,p_m,p_p,c_m,c_p,1);
    Sm = 0.0;
}
double Wave_speed_pressure_2D(double n_m, double u_m, double p_m,double n_p, double u_p, double p_p, double & Sl, double & Sr, double &Sm, double mtmp)
{
    double c_p=0, c_m=0, P_rl;
    if(p_m<=epsil || n_m<=epsil)
    {
     c_m = 0.0;
    }
    else
    {
    c_m = pow(gam*p_m/n_m,0.5);
    }
    if(p_p<=epsil || n_p<=epsil)
    {
    c_p = 0.0;
    }
    else
    {
    c_p = pow(gam*p_p/n_p,0.5);
    }
    P_rl =  pressure_estimate_acustic(n_m,n_p,u_m,u_p,p_m,p_p,c_m,c_p,mtmp);
    Sl = Wave_speed_Pressure(n_m,n_p,u_m,u_p,p_m,p_p,c_m,c_p,P_rl,0,mtmp);
    Sr = Wave_speed_Pressure(n_m,n_p,u_m,u_p,p_m,p_p,c_m,c_p,P_rl,1,mtmp);
    Sm = Wave_speed_Pressure(n_m,n_p,u_m,u_p,p_m,p_p,c_m,c_p,P_rl,2,mtmp);
}
double Wave_speed_maximal_2D(double ** s_l_d,double ** s_r_d,double ** s_l_u, double ** s_r_u, double ** s_d_l,double ** s_u_l,double ** s_d_r, double ** s_u_r, double **& Sl, double **& Sr,double **& Sd,double **& Su)
{
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            Sr[i][j] = max(s_r_u[i][j],s_r_d[i][j]);
            Sl[i][j] = min(s_l_u[i][j],s_l_d[i][j]);
            Su[i][j] = max(s_u_r[i][j],s_u_l[i][j]);
            Sd[i][j] = min(s_d_r[i][j],s_d_l[i][j]);
        }
        
    }
}
double U_F_G_2D(double n_m, double u_px_m, double u_py_m, double p_p_m, double n_p, double u_px_p, double u_py_p, double p_p_p, double & Ua1_m, double & Ua2_m, double & Ua3_m, double & Ua4_m, double & Ua1_p, double & Ua2_p, double & Ua3_p, double & Ua4_p, double & Fa1_m, double & Fa2_m, double & Fa3_m, double & Fa4_m, double & Fa1_p, double & Fa2_p, double & Fa3_p, double & Fa4_p, double & Ga1_m, double & Ga2_m, double & Ga3_m, double & Ga4_m, double & Ga1_p, double & Ga2_p, double & Ga3_p, double & Ga4_p, double qtmp, double mtmp)
{
    Ua1_m = U_two_fluid_2D(n_m,u_px_m,u_py_m,p_p_m,mtmp,qtmp,0);         // n
    Ua1_p = U_two_fluid_2D(n_p,u_px_p,u_py_p,p_p_p,mtmp,qtmp,0);
    Ua2_m = U_two_fluid_2D(n_m,u_px_m,u_py_m,p_p_m,mtmp,qtmp,1);        // nu
    Ua2_p = U_two_fluid_2D(n_p,u_px_p,u_py_p,p_p_p,mtmp,qtmp,1);
    Ua3_m = U_two_fluid_2D(n_m,u_px_m,u_py_m,p_p_m,mtmp,qtmp,2);        // nu
    Ua3_p = U_two_fluid_2D(n_p,u_px_p,u_py_p,p_p_p,mtmp,qtmp,2);
    Ua4_m = U_two_fluid_2D(n_m,u_px_m,u_py_m,p_p_m,mtmp,qtmp,3);
    Ua4_p = U_two_fluid_2D(n_p,u_px_p,u_py_p,p_p_p,mtmp,qtmp,3);  // e
    Fa1_m = F_two_fluid_2D(qtmp,n_m,u_px_m,u_py_m,p_p_m,mtmp,0);
    Fa1_p = F_two_fluid_2D(qtmp,n_p,u_px_p,u_py_p,p_p_p,mtmp,0);
    Fa2_m = F_two_fluid_2D(qtmp,n_m,u_px_m,u_py_m,p_p_m,mtmp,1);
    Fa2_p = F_two_fluid_2D(qtmp,n_p,u_px_p,u_py_p,p_p_p,mtmp,1);
    Fa3_m = F_two_fluid_2D(qtmp,n_m,u_px_m,u_py_m,p_p_m,mtmp,2);
    Fa3_p = F_two_fluid_2D(qtmp,n_p,u_px_p,u_py_p,p_p_p,mtmp,2);
    Fa4_m = F_two_fluid_2D(qtmp,n_m,u_px_m,u_py_m,p_p_m,mtmp,3);
    Fa4_p = F_two_fluid_2D(qtmp,n_p,u_px_p,u_py_p,p_p_p,mtmp,3);
    Ga1_m = G_two_fluid_2D(qtmp,n_m,u_px_m,u_py_m,p_p_m,mtmp,0);
    Ga1_p = G_two_fluid_2D(qtmp,n_p,u_px_p,u_py_p,p_p_p,mtmp,0);
    Ga2_m = G_two_fluid_2D(qtmp,n_m,u_px_m,u_py_m,p_p_m,mtmp,1);
    Ga2_p = G_two_fluid_2D(qtmp,n_p,u_px_p,u_py_p,p_p_p,mtmp,1);
    Ga3_m = G_two_fluid_2D(qtmp,n_m,u_px_m,u_py_m,p_p_m,mtmp,2);
    Ga3_p = G_two_fluid_2D(qtmp,n_p,u_px_p,u_py_p,p_p_p,mtmp,2);
    Ga4_m = G_two_fluid_2D(qtmp,n_m,u_px_m,u_py_m,p_p_m,mtmp,3);
    Ga4_p = G_two_fluid_2D(qtmp,n_p,u_px_p,u_py_p,p_p_p,mtmp,3);
}
double Flux_2D_HLL(double ** Ua1_m, double ** Ua2_m,double ** Ua3_m,double ** Ua4_m,double ** Ua1_p, double ** Ua2_p,double ** Ua3_p,double ** Ua4_p,double ** Fa1_m, double ** Fa2_m,double ** Fa3_m,double ** Fa4_m,double ** Fa1_p, double ** Fa2_p,double ** Fa3_p,double ** Fa4_p,double ** Ga1_m, double ** Ga2_m,double ** Ga3_m,double ** Ga4_m,double ** Ga1_p, double ** Ga2_p,double ** Ga3_p,double ** Ga4_p,double ** S_r_x,double ** S_l_x,double ** S_r_y,double ** S_l_y,double **& F_rho_hll,double **& G_rho_hll,double **& F_mom_x_hll,double **& G_mom_x_hll,double **& F_mom_y_hll,double **& G_mom_y_hll,double **& F_e_hll,double **& G_e_hll)
{
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            F_rho_hll[i][j]=Flux_HLL(Ua1_p[i][j],Ua1_m[i][j],Fa1_p[i][j],Fa1_m[i][j],S_r_x[i][j],S_l_x[i][j]);
            G_rho_hll[i][j]=Flux_HLL(Ua1_p[i][j],Ua1_m[i][j],Ga1_p[i][j],Ga1_m[i][j],S_r_y[i][j],S_l_y[i][j]);
            F_mom_x_hll[i][j]=Flux_HLL(Ua2_p[i][j],Ua2_m[i][j],Fa2_p[i][j],Fa2_m[i][j],S_r_x[i][j],S_l_x[i][j]);
            G_mom_x_hll[i][j]=Flux_HLL(Ua2_p[i][j],Ua2_m[i][j],Ga2_p[i][j],Ga2_m[i][j],S_r_y[i][j],S_l_y[i][j]);
            F_mom_y_hll[i][j]=Flux_HLL(Ua3_p[i][j],Ua3_m[i][j],Fa3_p[i][j],Fa3_m[i][j],S_r_x[i][j],S_l_x[i][j]);
            G_mom_y_hll[i][j]=Flux_HLL(Ua3_p[i][j],Ua3_m[i][j],Ga3_p[i][j],Ga3_m[i][j],S_r_y[i][j],S_l_y[i][j]);
            F_e_hll[i][j]=Flux_HLL(Ua4_p[i][j],Ua4_m[i][j],Fa4_p[i][j],Fa4_m[i][j],S_r_x[i][j],S_l_x[i][j]);
            G_e_hll[i][j]=Flux_HLL(Ua4_p[i][j],Ua4_m[i][j],Ga4_p[i][j],Ga4_m[i][j],S_r_y[i][j],S_l_y[i][j]);
        }
    }
}
double flux_2D_HLLC(double ** Ua1_m, double ** Ua2_m,double ** Ua3_m,double ** Ua4_m,double ** Ua1_p, double ** Ua2_p,double ** Ua3_p,double ** Ua4_p,double ** Fa1_m, double ** Fa2_m,double ** Fa3_m,double ** Fa4_m,double ** Fa1_p, double ** Fa2_p,double ** Fa3_p,double ** Fa4_p,double ** Ga1_m, double ** Ga2_m,double ** Ga3_m,double ** Ga4_m,double ** Ga1_p, double ** Ga2_p,double ** Ga3_p,double ** Ga4_p,double ** S_r,double ** S_l,double ** S_m,double **& F_rho_hll,double **& G_rho_hll,double **& F_mom_x_hll,double **& G_mom_x_hll,double **& F_mom_y_hll,double **& G_mom_y_hll,double **& F_e_hll,double **& G_e_hll,double ** prla,int dir, double mtmp)
{
    
    double s_l_x, s_r_x,s_m_x,s_l_y,s_r_y,s_m_y;
    double U1_star_m, U2_star_m, U3_star_m, U4_star_m, U1_star_p, U2_star_p, U3_star_p,U4_star_p;
    double F1_star_m, F2_star_m, F3_star_m, F4_star_m, F1_star_p, F2_star_p, F3_star_p,F4_star_p;
    double G1_star_m, G2_star_m, G3_star_m, G4_star_m, G1_star_p, G2_star_p, G3_star_p,G4_star_p;
    
    
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            if(dir==0)
            {
                U1_star_m = U_inter(Ua1_m[i][j],Fa1_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],0,mtmp,dir);
                U2_star_m = U_inter(Ua2_m[i][j],Fa2_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],1,mtmp,dir);
                U3_star_m = U_inter(Ua3_m[i][j],Fa3_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],2,mtmp,dir);
                U4_star_m = U_inter(Ua4_m[i][j],Fa4_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],3,mtmp,dir);
                U1_star_p = U_inter(Ua1_p[i][j],Fa1_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],0,mtmp,dir);
                U2_star_p = U_inter(Ua2_p[i][j],Fa2_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],1,mtmp,dir);
                U3_star_p = U_inter(Ua3_p[i][j],Fa3_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],2,mtmp,dir);
                U4_star_p = U_inter(Ua4_p[i][j],Fa4_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],3,mtmp,dir);
            }
            if(dir==1)
            {
                U1_star_m = U_inter(Ua1_m[i][j],Ga1_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],0,mtmp,dir);
                U2_star_m = U_inter(Ua2_m[i][j],Ga2_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],1,mtmp,dir);
                U3_star_m = U_inter(Ua3_m[i][j],Ga3_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],2,mtmp,dir);
                U4_star_m = U_inter(Ua4_m[i][j],Ga4_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],3,mtmp,dir);
                U1_star_p = U_inter(Ua1_p[i][j],Ga1_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],0,mtmp,dir);
                U2_star_p = U_inter(Ua2_p[i][j],Ga2_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],1,mtmp,dir);
                U3_star_p = U_inter(Ua3_p[i][j],Ga3_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],2,mtmp,dir);
                U4_star_p = U_inter(Ua4_p[i][j],Ga4_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],3,mtmp,dir);
            }
            
            F1_star_m = Flux_inter(Ua1_m[i][j],Fa1_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],0,mtmp,dir);
            F2_star_m = Flux_inter(Ua2_m[i][j],Fa2_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],1,mtmp,dir);
            F3_star_m = Flux_inter(Ua3_m[i][j],Fa3_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],2,mtmp,dir);
            F4_star_m = Flux_inter(Ua4_m[i][j],Fa4_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],3,mtmp,dir);
            F1_star_p = Flux_inter(Ua1_p[i][j],Fa1_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],0,mtmp,dir);
            F2_star_p = Flux_inter(Ua2_p[i][j],Fa2_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],1,mtmp,dir);
            F3_star_p = Flux_inter(Ua3_p[i][j],Fa3_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],2,mtmp,dir);
            F4_star_p = Flux_inter(Ua4_p[i][j],Fa4_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],3,mtmp,dir);
            G1_star_m = Flux_inter(Ua1_m[i][j],Ga1_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],0,mtmp,dir);
            G2_star_m = Flux_inter(Ua2_m[i][j],Ga2_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],1,mtmp,dir);
            G3_star_m = Flux_inter(Ua3_m[i][j],Ga3_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],2,mtmp,dir);
            G4_star_m = Flux_inter(Ua4_m[i][j],Ga4_m[i][j],S_l[i][j],S_m[i][j],prla[i][j],3,mtmp,dir);
            G1_star_p = Flux_inter(Ua1_p[i][j],Ga1_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],0,mtmp,dir);
            G2_star_p = Flux_inter(Ua2_p[i][j],Ga2_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],1,mtmp,dir);
            G3_star_p = Flux_inter(Ua3_p[i][j],Ga3_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],2,mtmp,dir);
            G4_star_p = Flux_inter(Ua4_p[i][j],Ga4_p[i][j],S_r[i][j],S_m[i][j],prla[i][j],3,mtmp,dir);
            
            F_rho_hll[i][j]=Flux_HLLC(Fa1_p[i][j],Fa1_m[i][j],F1_star_p,F1_star_m,S_l[i][j],S_r[i][j],S_m[i][j]);
            G_rho_hll[i][j]=Flux_HLLC(Ga1_p[i][j],Ga1_m[i][j],G1_star_p,G1_star_m,S_l[i][j],S_r[i][j],S_m[i][j]);
            F_mom_x_hll[i][j]=Flux_HLLC(Fa2_p[i][j],Fa2_m[i][j],F2_star_p,F2_star_m,S_l[i][j],S_r[i][j],S_m[i][j]);
            G_mom_x_hll[i][j]=Flux_HLLC(Ga2_p[i][j],Ga2_m[i][j],G2_star_p,G2_star_m,S_l[i][j],S_r[i][j],S_m[i][j]);
            F_mom_y_hll[i][j]=Flux_HLLC(Fa3_p[i][j],Fa3_m[i][j],F3_star_p,F3_star_m,S_l[i][j],S_r[i][j],S_m[i][j]);
            G_mom_y_hll[i][j]=Flux_HLLC(Ga3_p[i][j],Ga3_m[i][j],G3_star_p,G3_star_m,S_l[i][j],S_r[i][j],S_m[i][j]);
            F_e_hll[i][j]=Flux_HLLC(Fa4_p[i][j],Fa4_m[i][j],F4_star_p,F4_star_m,S_l[i][j],S_r[i][j],S_m[i][j]);
            G_e_hll[i][j]=Flux_HLLC(Ga4_p[i][j],Ga4_m[i][j],G4_star_p,G4_star_m,S_l[i][j],S_r[i][j],S_m[i][j]);
        }
    }
}
double  U_star_2D(double ** Uru,double ** Uld,double ** Urd,double ** Ulu,double ** Fru,double ** Flu,double ** Frd,double ** Fld,double ** Gru,double ** Grd,double ** Glu,double ** Gld,double ** F_star_R,double ** F_star_L,double ** G_star_U,double ** G_star_D,double ** s_l_d,double ** s_r_d,double ** s_l_u,double ** s_r_u,double ** s_d_l,double **  s_u_l,double ** s_d_r,double ** s_u_r,double ** S_l,double ** S_r,double ** S_d,double ** S_u,double **& U_star)
{
   double ** U_star_tmp = newArr(double, ntot.x,ntot.y);
    double num1= 0.0,num2= 0.0,num3= 0.0,num4= 0.0, den= 0.0;
    
    for(int j=0;j<ntot.y;j++)
    {
        for(int i=0;i<ntot.x;i++)
        {
            U_star_tmp[i][j] = 0.0;
        }
    }
    
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            den = ((S_r[i][j] - S_l[i][j])*(S_u[i][j] - S_d[i][j]));
            num1 = ((S_u[i][j]*S_r[i][j]*Uru[i][j]) +  (S_d[i][j]*S_l[i][j]*Uld[i][j]) - (S_d[i][j]*S_r[i][j]*Urd[i][j]) - (S_u[i][j]*S_l[i][j]*Ulu[i][j]));
            num2 = -0.5*((S_u[i][j]*(Fru[i][j]-Flu[i][j])) + (S_d[i][j]*(Frd[i][j]-Fld[i][j])) + (S_r[i][j]*(Gru[i][j]-Grd[i][j])) + (S_l[i][j]*(Glu[i][j]-Gld[i][j])));
            num3 = -0.5*(F_star_R[i][j] - F_star_L[i][j])*(S_u[i][j] - S_d[i][j]);
            num4 = -0.5*(G_star_U[i][j] - G_star_D[i][j])* (S_r[i][j] - S_l[i][j]);
            U_star_tmp[i][j] = (num1/den) + (num2/den) + (num3/(den)) + (num4/(den));
        }
    }
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            U_star[i][j] = U_star_tmp[i][j];
        }
    } 
    delArr(U_star_tmp,ntot.x);
}
double  F_star_2D(double Uru,double Urd,double Fru,double Flu ,double Frd,double Fld,double Gru,double Grd,double Glu,double Gld,double F_star_L,double F_star_R,double G_star_U,double G_star_D,double F_hll_u,double F_hll_d,double U_star,double S_l_d,double S_r_d,double S_l_u,double S_r_u,double S_d_l,double S_u_l,double S_d_r,double S_u_r,double S_l,double S_r,double S_d,double S_u,double & F_star)
{
    double num1= 0.0,num2= 0.0,num3= 0.0,num4= 0.0, den= 0.0;
    den = ((S_u -S_d)*(S_r-S_l));
    num1 = ((S_u*F_hll_u)/(S_u-S_d)) -( (S_d*F_hll_d)/(S_u-S_d));
    num2 = (S_l*S_r*(Gru-Glu+Gld-Grd))/den;
    num3 = ((S_l*S_u*(Fru-F_star_R)) + (S_r*S_u*(Flu-F_star_L)))/den;
    num4 = ((S_r*S_d*(Fld-F_star_L)) + (S_l*S_d*(Frd-F_star_R)))/den;
    F_star = (num1) - (num2)  - (num3) - (num4);
}
double  G_star_2D(double Uru,double Ulu,double Fru,double Flu,double Frd,double Fld,double Gld,double Gru,double Glu,double Grd,double F_star_R,double F_star_L,double G_star_D,double G_star_U,double G_hll_r,double G_hll_l,double U_star,double S_l_d,double S_r_d,double S_l_u,double S_r_u,double S_d_l,double  S_u_l,double S_d_r,double S_u_r,double S_l,double S_r,double S_d,double S_u,double & G_star)
{
    double num1= 0.0,num2= 0.0,num3= 0.0,num4= 0.0, den= 0.0;
    
    den = ((S_u -S_d)*(S_r-S_l));
    num1 = ((S_r*G_hll_r)/(S_r-S_l)) -( (S_l*G_hll_l)/(S_r-S_l));
    num2 = (S_d*S_u*(Fru-Flu+Fld-Frd))/den;
    num3 = ((S_r*S_d*(Gru-G_star_U)) + (S_l*S_d*(Glu-G_star_U)))/den;
    num4 = ((S_l*S_u*(Gld-G_star_D)) + (S_r*S_u*(Grd-G_star_D)))/den;
    G_star = (num1) - (num2)  - (num3) - (num4);
}
double Flux_HLL_2D(double ** Uru,double ** Uld,double ** Urd,double ** Ulu,double ** Fru,double ** Flu,double ** Frd,double ** Fld,double ** Gru,double ** Grd,double ** Glu,double ** Gld,double ** F_hll_d,double ** F_hll_u,double ** G_hll_l,double ** G_hll_r,double ** F_star_L,double ** F_star_R,double ** G_star_U,double ** G_star_D,double ** S_l_d,double ** S_r_d,double ** S_l_u,double ** S_r_u,double ** S_d_l,double **  S_u_l,double ** S_d_r,double ** S_u_r,double ** S_l,double ** S_r,double ** S_d,double ** S_u,double ** U_star, double **& F_hll_star, double **& G_hll_star)
{
    double **F_star_tmp  = newArr(double, ntot.x,ntot.y);
    double  **G_star_tmp = newArr(double, ntot.x,ntot.y); 
    
    for(int j=0;j<ntot.y;j++)
    {
        for(int i=0;i<ntot.x;i++)
        {
            F_star_tmp[i][j] = 0.0;
            G_star_tmp[i][j] = 0.0;
        }
    }
 
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            //both sub-sonic
            if(S_l[i][j] <= 0.0 && S_r[i][j]>=0.0)
            {
                if( S_d[i][j] <= 0.0 && S_u[i][j]>=0.0 )
                {
                    F_star_2D(Uru[i][j],Urd[i][j],Fru[i][j],Flu[i][j],Frd[i][j],Fld[i][j],Gru[i][j],Grd[i][j],Glu[i][j],Gld[i][j],F_star_L[i][j],F_star_R[i][j],G_star_U[i][j],G_star_D[i][j],F_hll_u[i][j],F_hll_d[i][j],U_star[i][j],S_l_d[i][j],S_r_d[i][j],S_l_u[i][j],S_r_u[i][j],S_d_l[i][j],S_u_l[i][j],S_d_r[i][j],S_u_r[i][j],S_l[i][j],S_r[i][j],S_d[i][j],S_u[i][j],F_star_tmp[i][j]);
                    G_star_2D(Uru[i][j],Ulu[i][j],Fru[i][j],Flu[i][j],Frd[i][j],Fld[i][j],Gld[i][j],Gru[i][j],Glu[i][j],Grd[i][j],F_star_R[i][j],F_star_L[i][j],G_star_D[i][j],G_star_U[i][j],G_hll_r[i][j],G_hll_l[i][j],U_star[i][j],S_l_d[i][j],S_r_d[i][j],S_l_u[i][j],S_r_u[i][j],S_d_l[i][j],S_u_l[i][j],S_d_r[i][j],S_u_r[i][j],S_l[i][j],S_r[i][j],S_d[i][j],S_u[i][j],G_star_tmp[i][j]);
                }
            }
            // both supersonic
            if(S_l[i][j]>=0.0 && S_d[i][j]>=0.0)
            {
                if (S_r[i][j]<=0.0 && S_u[i][j]<=0.0)
                {
                    
                    if(S_l[i][j]>=0.0 && S_d[i][j]>=0.0)
                    {
                        F_star_tmp[i][j] = Fld[i][j];
                        G_star_tmp[i][j] = Gld[i][j];
                    }
                    if(S_l[i][j]>=0.0 && S_u[i][j]<=0.0)
                    {
                        F_star_tmp[i][j] = Flu[i][j];
                        G_star_tmp[i][j] = Glu[i][j];
                    }
                    if(S_r[i][j]<=0.0 && S_d[i][j]>=0.0)
                    {
                        F_star_tmp[i][j] = Frd[i][j];
                        G_star_tmp[i][j] = Grd[i][j];
                    }
                    if(S_r[i][j]<=0.0 && S_u[i][j]<=0.0)
                    {
                        F_star_tmp[i][j] = Fru[i][j];
                        G_star_tmp[i][j] = Gru[i][j];
                    }
                }
            }
            // Supersonic x-direction & subsonic y-direction
            if( S_d[i][j] <= 0.0 && S_u[i][j]>=0.0 )
            {
                if( S_l[i][j] >= 0.0 || S_r[i][j]<=0.0 )
                {
                    F_star_tmp[i][j] = ((S_u[i][j]/(S_u[i][j] - S_d[i][j]))*F_hll_d[i][j]) - ((S_d[i][j]/(S_u[i][j] - S_d[i][j]))*F_hll_u[i][j]);
                    if(S_l[i][j] >= 0.0)
                    {
                        G_star_tmp[i][j] = G_hll_l[i][j];
                    }
                    if(S_r[i][j] <= 0.0)
                    {
                        G_star_tmp[i][j] = G_hll_r[i][j];
                    }
                }
            }
            // Supersonic y-direction & subsonic x-direction
            if( S_l[i][j] <= 0.0 && S_r[i][j]>=0.0 )
            {
                if( S_d[i][j] >= 0.0 || S_u[i][j]<=0.0 )
                {
                    G_star_tmp[i][j] = ((S_r[i][j]/(S_r[i][j] - S_l[i][j]))*G_hll_l[i][j]) - ((S_l[i][j]/(S_r[i][j] - S_l[i][j]))*G_hll_r[i][j]);
                    if(S_d[i][j] >= 0.0)
                    {
                        F_star_tmp[i][j] = F_hll_d[i][j];
                    }
                    if(S_u[i][j] <= 0.0)
                    {
                        F_star_tmp[i][j] = F_hll_u[i][j];
                    }
                }
            }
            
        }
    }
     for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
     {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            F_hll_star[i][j] = F_star_tmp[i][j];
            G_hll_star[i][j] = G_star_tmp[i][j];
        }
     }
    delArr(F_star_tmp,ntot.x);
    delArr(G_star_tmp,ntot.x);
}
void time_advance_hll_flux(double **F_hll_D_1, double **G_hll_L_1,double **F_hll_D_2,double **G_hll_L_2,double **F_hll_D_3,double **G_hll_L_3,double **F_hll_D_4,double **G_hll_L_4,double ** & n_tmp_out, double ** & m_x_tmp_out,double ** & m_y_tmp_out, double ** & e_tmp_out, int my_id, bool notWest, bool notEast,bool notSouth,bool notNorth,int specie)
{    
double aux3, aux4, aux5;
double **varF  = newArr(double, ntot.x,ntot.y);
double  **varG = newArr(double, ntot.x,ntot.y); 
double Fim1x, Gim1y, Gim2y;

/*

//kinetic flux at the left  boundary
if(notWest==false)
{
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
     F_hll_D_1[1][j]  = kflux1_i[j];
     F_hll_D_2[1][j]  = kflux2_i[j];
     F_hll_D_3[1][j]  = kflux3_i[j];
     F_hll_D_4[1][j]  = kflux4_i[j];
    }
}
//kinetic flux at the right boundary
if(notEast==false)
{
    for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
	F_hll_D_1[ntot.x-3][j]  = kflux1_e[j];
        F_hll_D_2[ntot.x-3][j]  = kflux2_e[j];
        F_hll_D_3[ntot.x-3][j]  = kflux3_e[j];
        F_hll_D_4[ntot.x-3][j]  = kflux4_e[j];
     */
     /*
     if(j==nlocal.y/2) 
     {
	if(specie==0)
	cout<<"----------ions----------"<<endl;
	else
	cout<<"----------electrons----------"<<endl;

	cout<<"----------kfluxes_in----------"<<endl; 
        cout<<"kf1_in="<<F_hll_D_1[1][j]<<endl;
        cout<<"kf2_in="<<F_hll_D_2[1][j]<<endl;
        cout<<"kf3_in="<<F_hll_D_3[1][j]<<endl;
        cout<<"kf4_in="<<F_hll_D_4[1][j]<<endl;
        cout<<"----------fluxes_2----------"<<endl;
        cout<<"f1_2="<<F_hll_D_1[2][j]<<endl;
        cout<<"f2_2="<<F_hll_D_1[2][j]<<endl;
        cout<<"f3_2="<<F_hll_D_3[2][j]<<endl;
        cout<<"f4_2="<<F_hll_D_4[2][j]<<endl;
	cout<<"----------kfluxes_N----------"<<endl;
	cout<<"f1_N="<<F_hll_D_1[ntot.x-4][j]<<endl;
        cout<<"f2_N="<<F_hll_D_2[ntot.x-4][j]<<endl;
        cout<<"f3_N="<<F_hll_D_3[ntot.x-4][j]<<endl;
        cout<<"f4_N="<<F_hll_D_4[ntot.x-4][j]<<endl;
	cout<<"----------kfluxes_out----------"<<endl;
        cout<<"kf1_out="<<F_hll_D_1[ntot.x-3][j]<<endl;
        cout<<"kf2_out="<<F_hll_D_2[ntot.x-3][j]<<endl;
        cout<<"kf3_out="<<F_hll_D_3[ntot.x-3][j]<<endl;
        cout<<"kf4_out="<<F_hll_D_4[ntot.x-3][j]<<endl;
	}
     */ 
//   }
//}
   for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            //density
            aux3 = F_hll_D_1[i][j];
            aux4 = (aux3);
            aux3 = F_hll_D_1[i-1][j];
            aux5 = (aux3);
            varF[i][j] = aux4 - aux5;
            aux3 = G_hll_L_1[i][j];
            aux4 = (aux3);
            aux3 = G_hll_L_1[i][j-1];
            aux5 = (aux3);
            varG[i][j] = aux4 - aux5;
            n_tmp_out[i][j]  =   - (1.0/dx*varF[i][j]) - (1.0/dy*varG[i][j]);
            //momentum-x
            aux3 = F_hll_D_2[i][j];
            aux4 = (aux3);
            aux3 = F_hll_D_2[i-1][j];
            aux5 = (aux3);
            varF[i][j] = aux4 - aux5;
            aux3 = G_hll_L_2[i][j];
            aux4 = (aux3);
            aux3 = G_hll_L_2[i][j-1];
            aux5 = (aux3);
            varG[i][j] = aux4 - aux5;
            m_x_tmp_out[i][j]=  - (1.0/dx*varF[i][j]) - (1.0/dy*varG[i][j]);
            //momentum-y
            aux3 = F_hll_D_3[i][j];
            aux4 =  (aux3);
            aux3 = F_hll_D_3[i-1][j];
            aux5 =  (aux3);
            varF[i][j] = aux4 - aux5;
            aux3 = G_hll_L_3[i][j];
            aux4 = (aux3);
            aux3 = G_hll_L_3[i][j-1];
            aux5 = (aux3);
            varG[i][j] = aux4 - aux5;
            m_y_tmp_out[i][j]=   - (1.0/dx*varF[i][j]) - (1.0/dy*varG[i][j]);
            //energy
            aux3 = F_hll_D_4[i][j];
            aux4 = (aux3);
            aux3 = F_hll_D_4[i-1][j];
            aux5 = (aux3);
            varF[i][j] = aux4 - aux5;
            aux3 = G_hll_L_4[i][j];
            aux4 = (aux3);
            aux3 = G_hll_L_4[i][j-1];
            aux5 = (aux3);
            varG[i][j] = aux4 - aux5;
            e_tmp_out[i][j]=  - (1.0/dx*varF[i][j]) - (1.0/dy*varG[i][j]);
        }
    }
    delArr(varF,ntot.x);
    delArr(varG,ntot.x);
}
void time_advance_SIR_flux(double **F_hll_D_1, double **G_hll_L_1,double **F_hll_D_2,double **G_hll_L_2,double **F_hll_D_3,double **G_hll_L_3,double **F_hll_D_4,double **G_hll_L_4,double **F_n_i,double **G_n_i,double **F_mom_x_i,double **G_mom_x_i,double **F_mom_y_i,double **G_mom_y_i,double **F_e_i,double **G_e_i,double ** & n_tmp_out, double ** & m_x_tmp_out,double ** & m_y_tmp_out, double ** & e_tmp_out)
{
    
    double aux1, aux2, aux3, aux4, aux5;
     double **varF  = newArr(double, ntot.x,ntot.y);
    double  **varG = newArr(double, ntot.x,ntot.y);
    
     for (int j=nbuffer;j<ntot.y-nbuffer; j++)
                 {
                 for (int i=nbuffer;i<ntot.x-nbuffer; i++)
                        {
                        aux1 = F_n_i[i][j];
                        aux2 = F_n_i[i][j-1];
                        aux3 = F_hll_D_1[i][j];
                        aux4 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        aux1 = F_n_i[i-1][j];
                        aux2 = F_n_i[i-1][j-1];
                        aux3 = F_hll_D_1[i-1][j];
                        aux5 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        varF[i][j] = aux4 - aux5;
                        aux1 = G_n_i[i-1][j];
                        aux2 = G_n_i[i][j];
                        aux3 = G_hll_L_1[i][j];
                        aux4 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        aux1 = G_n_i[i-1][j-1];
                        aux2 = G_n_i[i][j-1];
                        aux3 = G_hll_L_1[i][j-1];
                        aux5 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        varG[i][j] = aux4 - aux5;
                        n_tmp_out[i][j]  =   - (1.0/dx*varF[i][j]) - (1.0/dy*varG[i][j]);
                        aux1 = F_mom_x_i[i][j];
                        aux2 = F_mom_x_i[i][j-1];
                        aux3 = F_hll_D_2[i][j];
                        aux4 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        aux1 = F_mom_x_i[i-1][j];
                        aux2 = F_mom_x_i[i-1][j-1];
                        aux3 = F_hll_D_2[i-1][j];
                        aux5 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        varF[i][j] = aux4 - aux5;
                        aux1 = G_mom_x_i[i-1][j];
                        aux2 = G_mom_x_i[i][j];
                        aux3 = G_hll_L_2[i][j];
                        aux4 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        aux1 = G_mom_x_i[i-1][j-1];
                        aux2 = G_mom_x_i[i][j-1];
                        aux3 = G_hll_L_2[i][j-1];
                        aux5 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        varG[i][j] = aux4 - aux5;
                        m_x_tmp_out[i][j]=  - (1.0/dx*varF[i][j]) - (1.0/dy*varG[i][j]);
                        aux1 = F_mom_y_i[i][j];
                        aux2 = F_mom_y_i[i][j-1];
                        aux3 = F_hll_D_3[i][j];
                        aux4 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        aux1 = F_mom_y_i[i-1][j];
                        aux2 = F_mom_y_i[i-1][j-1];
                        aux3 = F_hll_D_3[i-1][j];
                        aux5 = ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        varF[i][j] = aux4 - aux5;
                        aux1 = G_mom_y_i[i-1][j];
                        aux2 = G_mom_y_i[i][j];
                        aux3 = G_hll_L_3[i][j];
                        aux4 = ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        aux1 = G_mom_y_i[i-1][j-1];
                        aux2 = G_mom_y_i[i][j-1];
                        aux3 = G_hll_L_3[i][j-1];
                        aux5 = ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        varG[i][j] = aux4 - aux5;
                        m_y_tmp_out[i][j]=   - (1.0/dx*varF[i][j]) - (1.0/dy*varG[i][j]);
                        aux1 = F_e_i[i][j];
                        aux2 = F_e_i[i][j-1];
                        aux3 = F_hll_D_4[i][j];
                        aux4 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        aux1 = F_e_i[i-1][j];
                        aux2 = F_e_i[i-1][j-1];
                        aux3 = F_hll_D_4[i-1][j];
                        aux5 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        varF[i][j] = aux4 - aux5;
                        aux1 = G_e_i[i-1][j];
                        aux2 = G_e_i[i][j];
                        aux3 = G_hll_L_4[i][j];
                        aux4 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        aux1 = G_e_i[i-1][j-1];
                        aux2 = G_e_i[i][j-1];
                        aux3 = G_hll_L_4[i][j-1];
                        aux5 =  ((1.0/6.0)*aux1) + ((1.0/6.0)*aux2) + ((2.0/3.0)*(aux3));
                        varG[i][j] = aux4 - aux5;
                        e_tmp_out[i][j]=  - (1.0/dx*varF[i][j]) - (1.0/dy*varG[i][j]);
                        }
                    }
 
    delArr(varF,ntot.x);
    delArr(varG,ntot.x);
}
void wave_speed_estimate_x(double **n_p_tmp,double **u_xp_tmp,double **u_yp_tmp, double **P_p_tmp, double **n_m_tmp,double **u_xm_tmp, double **u_ym_tmp, double **P_m_tmp, int tmp, double **s_l_d_x, double **s_r_d_x,double **s_l_u_x,double **s_r_u_x,double **s_d_l_x,double **s_u_l_x,double **s_d_r_x, double **s_u_r_x, double **s_m_d_x,double **s_m_u_x,double **s_m_l_x, double **s_m_r_x,double **P_rla_x, double **P_rlb_x, double mtmp)
{
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
        if (tmp==0)
            {
           Wave_speed_linea_stimates_2D(n_p_tmp[i][j],u_xp_tmp[i][j],P_p_tmp[i][j],n_m_tmp[i][j],u_xm_tmp[i][j],P_m_tmp[i][j],s_l_d_x[i][j],s_r_d_x[i][j],s_m_d_x[i][j]);
           Wave_speed_linea_stimates_2D(n_p_tmp[i][j+1],u_xp_tmp[i][j+1],P_p_tmp[i][j+1],n_m_tmp[i][j+1],u_xm_tmp[i][j+1],P_m_tmp[i][j+1],s_l_u_x[i][j],s_r_u_x[i][j],s_m_u_x[i][j]);
            }
            if(tmp==1)
            {
            Wave_speed_pressure_2D(n_p_tmp[i][j],u_xp_tmp[i][j],P_p_tmp[i][j],n_m_tmp[i][j],u_xm_tmp[i][j],P_m_tmp[i][j],s_l_d_x[i][j],s_r_d_x[i][j],s_m_d_x[i][j],mtmp);
           P_rla_x[i][j] = Average_pressure(n_p_tmp[i][j],n_m_tmp[i][j],u_xp_tmp[i][j],u_xm_tmp[i][j],P_p_tmp[i][j],P_m_tmp[i][j],s_l_d_x[i][j],s_r_d_x[i][j],s_m_d_x[i][j],mtmp);
           Wave_speed_pressure_2D(n_p_tmp[i][j+1],u_xp_tmp[i][j+1],P_p_tmp[i][j+1],n_m_tmp[i][j+1],u_xm_tmp[i][j+1],P_m_tmp[i][j+1],s_l_u_x[i][j],s_r_u_x[i][j],s_m_u_x[i][j],mtmp);
           P_rlb_x[i][j] = Average_pressure(n_p_tmp[i][j+1],n_m_tmp[i][j+1],u_xp_tmp[i][j+1],u_xm_tmp[i][j+1],P_p_tmp[i][j+1],P_m_tmp[i][j+1],s_l_u_x[i][j],s_r_u_x[i][j],s_m_u_x[i][j],mtmp);
                }
            }
        }
}
void wave_speed_estimate_y(double **n_p_tmp,double **u_xp_tmp,double **u_yp_tmp, double **P_p_tmp, double **n_m_tmp,double **u_xm_tmp, double **u_ym_tmp, double **P_m_tmp ,int tmp,double **s_l_d_y, double **s_r_d_y, double **s_l_u_y, double **s_r_u_y, double **s_d_l_y, double **s_u_l_y, double **s_d_r_y, double **s_u_r_y, double **s_m_d_y, double **s_m_u_y, double **s_m_l_y, double **s_m_r_y,double **P_rla_y,double  **P_rlb_y, double mtmp)
{
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            if (tmp==0)
            {
                Wave_speed_linea_stimates_2D(n_p_tmp[i][j],u_yp_tmp[i][j],P_p_tmp[i][j],n_m_tmp[i][j],u_ym_tmp[i][j],P_m_tmp[i][j],s_d_l_y[i][j],s_u_l_y[i][j],s_m_l_y[i][j]);
                Wave_speed_linea_stimates_2D(n_p_tmp[i+1][j],u_yp_tmp[i+1][j],P_p_tmp[i+1][j],n_m_tmp[i+1][j],u_ym_tmp[i+1][j],P_m_tmp[i+1][j],s_d_r_y[i][j],s_u_r_y[i][j], s_m_r_y[i][j]);
            }
            if(tmp==1)
            {
                    Wave_speed_pressure_2D(n_p_tmp[i][j],u_yp_tmp[i][j],P_p_tmp[i][j],n_m_tmp[i][j],u_ym_tmp[i][j],P_m_tmp[i][j],s_d_l_y[i][j],s_u_l_y[i][j],s_m_l_y[i][j],mtmp);
                    P_rla_y[i][j] = Average_pressure(n_p_tmp[i][j],n_m_tmp[i][j],u_yp_tmp[i][j],u_ym_tmp[i][j],P_p_tmp[i][j],P_m_tmp[i][j],s_d_l_y[i][j],s_u_l_y[i][j],s_m_l_y[i][j],mtmp);
                    Wave_speed_pressure_2D(n_p_tmp[i+1][j],u_yp_tmp[i+1][j],P_p_tmp[i+1][j],n_m_tmp[i+1][j],u_ym_tmp[i+1][j],P_m_tmp[i+1][j],s_d_r_y[i][j],s_u_r_y[i][j],s_m_r_y[i][j],mtmp);
                    P_rlb_y[i][j] = Average_pressure(n_p_tmp[i+1][j],n_m_tmp[i+1][j],u_yp_tmp[i+1][j],u_ym_tmp[i+1][j],P_p_tmp[i+1][j],P_m_tmp[i+1][j],s_d_r_y[i][j],s_u_r_y[i][j],s_m_r_y[i][j],mtmp);
                }
            }
        }
}

 void fluxes_SIR(double **n_p_tmp,double **u_xp_tmp,double **u_yp_tmp,double **P_p_tmp,double **n_m_tmp,double **u_xm_tmp,double **u_ym_tmp,double **P_m_tmp,double **F_hll_D_1,double **G_star_D_1,double **F_hll_D_2,double **G_star_D_2,double **F_hll_D_3,double **G_star_D_3,double **F_hll_D_4,double **G_star_D_4,double **F_hll_U_1,double **G_star_U_1,double **F_hll_U_2,double **G_star_U_2,double **F_hll_U_3,double **G_star_U_3,double **F_hll_U_4,double **G_star_U_4,double **F_star_L_1,double **G_hll_L_1,double **F_star_L_2,double **G_hll_L_2,double **F_star_L_3,double **G_hll_L_3,double **F_star_L_4,double **G_hll_L_4,double **F_star_R_1,double **G_hll_R_1,double **F_star_R_2,double **G_hll_R_2,double **F_star_R_3,double **G_hll_R_3,double **F_star_R_4,double **G_hll_R_4,double **& F_n_i,double **&  G_n_i,double **& F_mom_x_i,double **& G_mom_x_i,double **& F_mom_y_i,double **& G_mom_y_i,double **& F_e_i,double **& G_e_i,double **s_l_d_x, double  **s_r_d_x, double  **s_l_u_x, double  **s_r_u_x, double **s_d_l_x, double  **s_u_l_x, double  **s_d_r_x, double  **s_u_r_x, double  **s_m_d_x, double  **s_m_u_x, double  **s_m_l_x, double  **s_m_r_x, double  **s_l_d_y, double  **s_r_d_y,double  **s_l_u_y, double  **s_r_u_y,double  **s_d_l_y,double  **s_u_l_y,double   **s_d_r_y,double   **s_u_r_y,double  **s_m_d_y,double  **s_m_u_y,double  **s_m_l_y,double  **s_m_r_y, int my_id,double qtmp, double mtmp)
{
double  **Uld1 = newArr(double, ntot.x,ntot.y);
double  **Uld2 = newArr(double, ntot.x,ntot.y);
double  **Uld3 = newArr(double, ntot.x,ntot.y);
double  **Uld4 = newArr(double, ntot.x,ntot.y);
double  **Ulu1 = newArr(double, ntot.x,ntot.y);
double  **Ulu2 = newArr(double, ntot.x,ntot.y);
double  **Ulu3 = newArr(double, ntot.x,ntot.y);
double  **Ulu4 = newArr(double, ntot.x,ntot.y);
double  **Fld1 = newArr(double, ntot.x,ntot.y);
double  **Fld2 = newArr(double, ntot.x,ntot.y);
double  **Fld3 = newArr(double, ntot.x,ntot.y);
double  **Fld4 = newArr(double, ntot.x,ntot.y);
double  **Flu1 = newArr(double, ntot.x,ntot.y);
double  **Flu2 = newArr(double, ntot.x,ntot.y);
double  **Flu3 = newArr(double, ntot.x,ntot.y);
double  **Flu4 = newArr(double, ntot.x,ntot.y);
double  **Gld1 = newArr(double, ntot.x,ntot.y);
double  **Gld2 = newArr(double, ntot.x,ntot.y);
double  **Gld3 = newArr(double, ntot.x,ntot.y);
double  **Gld4 = newArr(double, ntot.x,ntot.y);
double  **Glu1 = newArr(double, ntot.x,ntot.y);
double  **Glu2 = newArr(double, ntot.x,ntot.y);
double  **Glu3 = newArr(double, ntot.x,ntot.y);
double  **Glu4 = newArr(double, ntot.x,ntot.y);
double  **Urd1 = newArr(double, ntot.x,ntot.y);
double  **Urd2 = newArr(double, ntot.x,ntot.y);
double  **Urd3 = newArr(double, ntot.x,ntot.y);
double  **Urd4 = newArr(double, ntot.x,ntot.y);
double  **Uru1 = newArr(double, ntot.x,ntot.y);
double  **Uru2 = newArr(double, ntot.x,ntot.y);
double  **Uru3 = newArr(double, ntot.x,ntot.y);
double  **Uru4 = newArr(double, ntot.x,ntot.y);
double  **Frd1 = newArr(double, ntot.x,ntot.y);
double  **Frd2 = newArr(double, ntot.x,ntot.y);
double  **Frd3 = newArr(double, ntot.x,ntot.y);
double  **Frd4 = newArr(double, ntot.x,ntot.y);
double  **Fru1 = newArr(double, ntot.x,ntot.y);
double  **Fru2 = newArr(double, ntot.x,ntot.y);
double  **Fru3 = newArr(double, ntot.x,ntot.y);
double  **Fru4 = newArr(double, ntot.x,ntot.y);
double  **Grd1 = newArr(double, ntot.x,ntot.y);
double  **Grd2 = newArr(double, ntot.x,ntot.y);
double  **Grd3 = newArr(double, ntot.x,ntot.y);
double  **Grd4 = newArr(double, ntot.x,ntot.y);
double  **Gru1 = newArr(double, ntot.x,ntot.y);
double  **Gru2 = newArr(double, ntot.x,ntot.y);
double  **Gru3 = newArr(double, ntot.x,ntot.y);
double  **Gru4 = newArr(double, ntot.x,ntot.y);
double  **U_star1 = newArr(double, ntot.x,ntot.y);
double  **U_star2 = newArr(double, ntot.x,ntot.y);
double  **U_star3 = newArr(double, ntot.x,ntot.y);
double  **U_star4 = newArr(double, ntot.x,ntot.y);
double  **s_l = newArr(double, ntot.x,ntot.y);
double  **s_r = newArr(double, ntot.x,ntot.y);
double  **s_d = newArr(double, ntot.x,ntot.y);
double  **s_u = newArr(double, ntot.x,ntot.y);

for (int ii=0;ii<ntot.x;ii++)
{
    for (int jj=0;jj<ntot.y;jj++)
    {
          Uld1[ii][jj] = 0.0;
          Uld2[ii][jj] = 0.0;
          Uld3[ii][jj] = 0.0;
          Uld4[ii][jj] = 0.0;
          Ulu1[ii][jj] = 0.0;
          Ulu2[ii][jj] = 0.0;
          Ulu3[ii][jj] = 0.0;
          Ulu4[ii][jj] = 0.0;
          Fld1[ii][jj] = 0.0;
          Fld2[ii][jj] = 0.0;
          Fld3[ii][jj] = 0.0;
          Fld4[ii][jj] = 0.0;
          Flu1[ii][jj] = 0.0;
          Flu2[ii][jj] = 0.0;
          Flu3[ii][jj] = 0.0;
          Flu4[ii][jj] = 0.0;
          Gld1[ii][jj] = 0.0;
          Gld2[ii][jj] = 0.0;
          Gld3[ii][jj] = 0.0;
          Gld4[ii][jj] = 0.0;
          Glu1[ii][jj] = 0.0;
          Glu2[ii][jj] = 0.0;
          Glu3[ii][jj] = 0.0;
          Glu4[ii][jj] = 0.0;
          Urd1[ii][jj] = 0.0;
          Urd2[ii][jj] = 0.0;
          Urd3[ii][jj] = 0.0;
          Urd4[ii][jj] = 0.0;
          Uru1[ii][jj] = 0.0;
          Uru2[ii][jj] = 0.0;
          Uru3[ii][jj] = 0.0;
          Uru4[ii][jj] = 0.0;
          Frd1[ii][jj] = 0.0;
          Frd2[ii][jj] = 0.0;
          Frd3[ii][jj] = 0.0;
          Frd4[ii][jj] = 0.0;
          Fru1[ii][jj] = 0.0;
          Fru2[ii][jj] = 0.0;
          Fru3[ii][jj] = 0.0;
          Fru4[ii][jj] = 0.0;
          Grd1[ii][jj] = 0.0;
          Grd2[ii][jj] = 0.0;
          Grd3[ii][jj] = 0.0;
          Grd4[ii][jj] = 0.0;
          Gru1[ii][jj] = 0.0;
          Gru2[ii][jj] = 0.0;
          Gru3[ii][jj] = 0.0;
          Gru4[ii][jj] = 0.0;
          U_star1[ii][jj] = 0.0;
          U_star2[ii][jj] = 0.0;
          U_star3[ii][jj] = 0.0;
          U_star4[ii][jj] = 0.0;
          s_l[ii][jj] = 0.0;
          s_r[ii][jj] = 0.0;
          s_d[ii][jj] = 0.0;
          s_u[ii][jj] = 0.0;

    }
}
    //maximal velocities
    Wave_speed_maximal_2D(s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u);
    for(int j=nbuffer-1;j<ntot.y;j++)
    {
        for(int i=nbuffer-1;i<ntot.x;i++)
        {
	 U_F_G_2D(n_p_tmp[i][j],u_xp_tmp[i][j],u_yp_tmp[i][j],P_p_tmp[i][j],n_m_tmp[i][j],u_xm_tmp[i][j],u_ym_tmp[i][j],P_m_tmp[i][j],Uld1[i][j],Uld2[i][j],Uld3[i][j],Uld4[i][j],Urd1[i][j],Urd2[i][j],Urd3[i][j],Urd4[i][j],Fld1[i][j],Fld2[i][j],Fld3[i][j],Fld4[i][j],Frd1[i][j],Frd2[i][j],Frd3[i][j],Frd4[i][j],Gld1[i][j],Gld2[i][j],Gld3[i][j],Gld4[i][j],Grd1[i][j],Grd2[i][j],Grd3[i][j],Grd4[i][j],qtmp, mtmp);  
          U_F_G_2D(n_p_tmp[i][j+1],u_xp_tmp[i][j+1],u_yp_tmp[i][j+1],P_p_tmp[i][j+1],n_m_tmp[i][j+1],u_xm_tmp[i][j+1],u_ym_tmp[i][j+1],P_m_tmp[i][j+1],Ulu1[i][j],Ulu2[i][j],Ulu3[i][j],Ulu4[i][j],Uru1[i][j],Uru2[i][j],Uru3[i][j],Uru4[i][j],Flu1[i][j],Flu2[i][j],Flu3[i][j],Flu4[i][j],Fru1[i][j],Fru2[i][j],Fru3[i][j],Fru4[i][j],Glu1[i][j],Glu2[i][j],Glu3[i][j],Glu4[i][j],Gru1[i][j],Gru2[i][j],Gru3[i][j],Gru4[i][j],qtmp, mtmp);
	}
    }
    //U_star
    U_star_2D(Uru1,Uld1,Urd1,Ulu1,Fru1,Flu1,Frd1,Fld1,Gru1,Grd1,Glu1,Gld1,F_star_R_1,F_star_L_1,G_star_U_1,G_star_D_1,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u,U_star1);
    U_star_2D(Uru2,Uld2,Urd2,Ulu2,Fru2,Flu2,Frd2,Fld2,Gru2,Grd2,Glu2,Gld2,F_star_R_2,F_star_L_2,G_star_U_2,G_star_D_2,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u,U_star2);
    U_star_2D(Uru3,Uld3,Urd3,Ulu3,Fru3,Flu3,Frd3,Fld3,Gru3,Grd3,Glu3,Gld3,F_star_R_3,F_star_L_3,G_star_U_3,G_star_D_3,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u,U_star3);
    U_star_2D(Uru4,Uld4,Urd4,Ulu4,Fru4,Flu4,Frd4,Fld4,Gru4,Grd4,Glu4,Gld4,F_star_R_4,F_star_L_4,G_star_U_4,G_star_D_4,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u,U_star4);
    //F_star
    Flux_HLL_2D(Uru1,Uld1,Urd1,Ulu1,Fru1,Flu1,Frd1,Fld1,Gru1,Grd1,Glu1,Gld1,F_hll_D_1,F_hll_U_1,G_hll_L_1,G_hll_R_1,F_star_L_1,F_star_R_1,G_star_U_1,G_star_D_1,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u,U_star1,F_n_i,G_n_i);
    Flux_HLL_2D(Uru2,Uld2,Urd2,Ulu2,Fru2,Flu2,Frd2,Fld2,Gru2,Grd2,Glu2,Gld2,F_hll_D_2,F_hll_U_2,G_hll_L_2,G_hll_R_2,F_star_L_2,F_star_R_2,G_star_U_2,G_star_D_2,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u,U_star2,F_mom_x_i,G_mom_x_i);
    Flux_HLL_2D(Uru3,Uld3,Urd3,Ulu3,Fru3,Flu3,Frd3,Fld3,Gru3,Grd3,Glu3,Gld3,F_hll_D_3,F_hll_U_3,G_hll_L_3,G_hll_R_3,F_star_L_3,F_star_R_3,G_star_U_3,G_star_D_3,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u,U_star3,F_mom_y_i,G_mom_y_i);
    Flux_HLL_2D(Uru4,Uld4,Urd4,Ulu4,Fru4,Flu4,Frd4,Fld4,Gru4,Grd4,Glu4,Gld4,F_hll_D_4,F_hll_U_4,G_hll_L_4,G_hll_R_4,F_star_L_4,F_star_R_4,G_star_U_4,G_star_D_4,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_l,s_r,s_d,s_u,U_star4,F_e_i,G_e_i);

 delArr(Uld1,ntot.x);
 delArr(Uld2,ntot.x);
 delArr(Uld3,ntot.x);
 delArr(Uld4,ntot.x);
 delArr(Ulu1,ntot.x);
 delArr(Ulu2,ntot.x);
 delArr(Ulu3,ntot.x);
 delArr(Ulu4,ntot.x);
 delArr(Fld1,ntot.x);
 delArr(Fld2,ntot.x);
 delArr(Fld3,ntot.x);
 delArr(Fld4,ntot.x);
 delArr(Flu1,ntot.x);
 delArr(Flu2,ntot.x);
 delArr(Flu3,ntot.x);
 delArr(Flu4,ntot.x);
 delArr(Gld1,ntot.x);
 delArr(Gld2,ntot.x);
 delArr(Gld3,ntot.x);
 delArr(Gld4,ntot.x);
 delArr(Glu1,ntot.x);
 delArr(Glu2,ntot.x);
 delArr(Glu3,ntot.x);
 delArr(Glu4,ntot.x);
 delArr(Urd1,ntot.x);
 delArr(Urd2,ntot.x);
 delArr(Urd3,ntot.x);
 delArr(Urd4,ntot.x);
 delArr(Uru1,ntot.x);
 delArr(Uru2,ntot.x);
 delArr(Uru3,ntot.x);
 delArr(Uru4,ntot.x);
 delArr(Frd1,ntot.x);
 delArr(Frd2,ntot.x);
 delArr(Frd3,ntot.x);
 delArr(Frd4,ntot.x);
 delArr(Fru1,ntot.x);
 delArr(Fru2,ntot.x);
 delArr(Fru3,ntot.x);
 delArr(Fru4,ntot.x);
 delArr(Grd1,ntot.x);
 delArr(Grd2,ntot.x);
 delArr(Grd3,ntot.x);
 delArr(Grd4,ntot.x);
 delArr(Gru1,ntot.x);
 delArr(Gru2,ntot.x);
 delArr(Gru3,ntot.x);
 delArr(Gru4,ntot.x);
 delArr(U_star1,ntot.x);
 delArr(U_star2,ntot.x);
 delArr(U_star3,ntot.x);
 delArr(U_star4,ntot.x);
 delArr(s_l,ntot.x);
 delArr(s_r,ntot.x);
 delArr(s_d,ntot.x);
 delArr(s_u,ntot.x);
}
void Fluxes_Boundary_conditions(double **F_n_i,double **G_n_i,double **F_mom_x_i,double **G_mom_x_i,double **F_mom_y_i,double **G_mom_y_i,double **F_e_i,double **G_e_i, int my_id)
{
if(bc_x==0)
{
                       for (int j=1;j<ntot.y-1;j++)
                       {
                        F_n_i[1][j] = F_n_i[2][j];
                        G_n_i[1][j] = G_n_i[2][j];
                        F_mom_x_i[1][j] = F_mom_x_i[2][j];
                        G_mom_x_i[1][j] = G_mom_x_i[2][j];
                        F_mom_y_i[1][j] = F_mom_y_i[2][j];
                        G_mom_y_i[1][j] = G_mom_y_i[2][j];
                        F_e_i[1][j] = F_e_i[2][j];
                        G_e_i[1][j] = G_e_i[2][j];

                        F_n_i[ntot.x-2][j] = F_n_i[ntot.x-3][j];
                        G_n_i[ntot.x-2][j] = G_n_i[ntot.x-3][j];
                        F_mom_x_i[ntot.x-2][j] = F_mom_x_i[ntot.x-3][j];
                        G_mom_x_i[ntot.x-2][j] = G_mom_x_i[ntot.x-3][j];
                        F_mom_y_i[ntot.x-2][j] = F_mom_y_i[ntot.x-3][j];
                        G_mom_y_i[ntot.x-2][j] = G_mom_y_i[ntot.x-3][j];
                        F_e_i[ntot.x-2][j] = F_e_i[ntot.x-3][j];
                        G_e_i[ntot.x-2][j] = G_e_i[ntot.x-3][j];
                        }
}
if(bc_y==0)
{
                        for (int j=1;j<ntot.x-1;j++)
                        {
                        F_n_i[j][1] = F_n_i[j][2];
                        G_n_i[j][1] = G_n_i[j][2];
                        F_mom_x_i[j][1] = F_mom_x_i[j][2];
                        G_mom_x_i[j][1] = G_mom_x_i[j][2];
                        F_mom_y_i[j][1] = F_mom_y_i[j][2];
                        G_mom_y_i[j][1] = G_mom_y_i[j][2];
                        F_e_i[j][1] = F_e_i[j][2];
                        G_e_i[j][1] = G_e_i[j][2];

                        F_n_i[j][ntot.y-2] = F_n_i[j][ntot.y-3];
                        G_n_i[j][ntot.y-2] = G_n_i[j][ntot.y-3];
                        F_mom_x_i[j][ntot.y-2] = F_mom_x_i[j][ntot.y-3];
                        G_mom_x_i[j][ntot.y-2] = G_mom_x_i[j][ntot.y-3];
                        F_mom_y_i[j][ntot.y-2] = F_mom_y_i[j][ntot.y-3];
                        G_mom_y_i[j][ntot.y-2] = G_mom_y_i[j][ntot.y-3];
                        F_e_i[j][ntot.y-2] = F_e_i[j][ntot.y-3];
                        G_e_i[j][ntot.y-2] = G_e_i[j][ntot.y-3];
                        }
}
/*
if(bc_x==1)
{                       
                        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
                        {
                        F_n_i[1][j] = F_n_i[ntot.x-3][j];
                        G_n_i[1][j] = G_n_i[ntot.x-3][j];
                        F_mom_x_i[1][j] = F_mom_x_i[ntot.x-3][j];
                        G_mom_x_i[1][j] = G_mom_x_i[ntot.x-3][j];
                        F_mom_y_i[1][j] = F_mom_y_i[ntot.x-3][j];
                        G_mom_y_i[1][j] = G_mom_y_i[ntot.x-3][j];
                        F_e_i[1][j] = F_e_i[ntot.x-3][j];
                        G_e_i[1][j] = G_e_i[ntot.x-3][j];
                        }
}
if(bc_y==1)
{                       
                        for (int j=nbuffer;j<ntot.x-nbuffer;j++)
                        {
                        F_n_i[j][1] = F_n_i[j][ntot.y-3];
                        G_n_i[j][1] = G_n_i[j][ntot.y-3];
                        F_mom_x_i[j][1] = F_mom_x_i[j][ntot.y-3];
                        G_mom_x_i[j][1] = G_mom_x_i[j][ntot.y-3];
                        F_mom_y_i[j][1] = F_mom_y_i[j][ntot.y-3];
                        G_mom_y_i[j][1] = G_mom_y_i[j][ntot.y-3];
                        F_e_i[j][1] = F_e_i[j][ntot.y-3];
                        G_e_i[j][1] = G_e_i[j][ntot.y-3];
                        }
        }
*/
}
void fluxes_x(double **n_p_tmp,double **u_xp_tmp,double **u_yp_tmp, double **P_p_tmp, double **n_m_tmp,double **u_xm_tmp, double **u_ym_tmp, double **P_m_tmp ,int tmp,double **F_hll_D_1,double **G_star_D_1,double **F_hll_D_2,double **G_star_D_2,double **F_hll_D_3,double **G_star_D_3,double **F_hll_D_4,double **G_star_D_4,double **F_hll_U_1,double **G_star_U_1,double **F_hll_U_2,double **G_star_U_2,double **F_hll_U_3,double **G_star_U_3,double **F_hll_U_4,double **G_star_U_4, double **s_l_d_x, double  **s_r_d_x, double  **s_l_u_x, double  **s_r_u_x, double **s_d_l_x, double  **s_u_l_x, double  **s_d_r_x, double  **s_u_r_x, double  **s_m_d_x, double  **s_m_u_x, double  **s_m_l_x, double  **s_m_r_x,double qtmp, double mtmp )
{

double  **P_rla_x= new double *[ntot.x]; 
for (int i = 0; i < ntot.x; i++)
P_rla_x [i]= new double[ntot.y]; 

double  **P_rlb_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_rlb_x [i]= new double[ntot.y];

double  **s_l_d_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_l_d_y [i]= new double[ntot.y];
 
double  **s_r_d_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_r_d_y [i]= new double[ntot.y];

double  **s_l_u_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_l_u_y [i]= new double[ntot.y];
 
double  **s_r_u_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_r_u_y [i]= new double[ntot.y];
 
double  **s_d_l_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_d_l_y [i]= new double[ntot.y];

double  **s_u_l_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_u_l_y [i]= new double[ntot.y];

double  **s_d_r_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_d_r_y [i]= new double[ntot.y];

double  **s_u_r_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_u_r_y [i]= new double[ntot.y];

double  **s_m_d_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_d_y [i]= new double[ntot.y];

double  **s_m_u_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_u_y [i]= new double[ntot.y];

double  **s_m_l_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_l_y [i]= new double[ntot.y];
 
double  **s_m_r_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_r_y [i]= new double[ntot.y];

double  **P_rla_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_rla_y [i]= new double[ntot.y];

double  **P_rlb_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_rlb_y [i]= new double[ntot.y];

double  **Uld1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uld1[i]= new double[ntot.y];
 
double  **Uld2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uld2[i]= new double[ntot.y];

double  **Uld3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uld3[i]= new double[ntot.y];

double  **Uld4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uld4[i]= new double[ntot.y];

double  **Ulu1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Ulu1[i]= new double[ntot.y];

double  **Ulu2= new double *[ntot.x];
 for (int i = 0; i < ntot.x; i++)
Ulu2[i]= new double[ntot.y];

double  **Ulu3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Ulu3[i]= new double[ntot.y];

double  **Ulu4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Ulu4[i]= new double[ntot.y];

double  **Fld1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fld1[i]= new double[ntot.y];

double  **Fld2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fld2[i]= new double[ntot.y];

double  **Fld3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fld3[i]= new double[ntot.y];
 
double  **Fld4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fld4[i]= new double[ntot.y];
 
double  **Flu1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Flu1[i]= new double[ntot.y];

double  **Flu2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Flu2[i]= new double[ntot.y];

double  **Flu3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Flu3[i]= new double[ntot.y];

double  **Flu4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Flu4[i]= new double[ntot.y];

double  **Gld1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gld1[i]= new double[ntot.y];

double  **Gld2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gld2[i]= new double[ntot.y];

double  **Gld3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gld3[i]= new double[ntot.y];

double  **Gld4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gld4[i]= new double[ntot.y];

double  **Glu1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Glu1[i]= new double[ntot.y];

double  **Glu2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Glu2[i]= new double[ntot.y];

double  **Glu3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Glu3[i]= new double[ntot.y];

double  **Glu4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Glu4[i]= new double[ntot.y];

double  **Urd1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Urd1[i]= new double[ntot.y];

double  **Urd2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Urd2[i]= new double[ntot.y];

double  **Urd3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Urd3[i]= new double[ntot.y];

double  **Urd4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Urd4[i]= new double[ntot.y];

double  **Uru1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uru1[i]= new double[ntot.y];

double  **Uru2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uru2[i]= new double[ntot.y];

double  **Uru3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uru3[i]= new double[ntot.y];

double  **Uru4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uru4[i]= new double[ntot.y];

double  **Frd1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Frd1[i]= new double[ntot.y];
 
double  **Frd2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Frd2[i]= new double[ntot.y];

double  **Frd3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Frd3[i]= new double[ntot.y];

double  **Frd4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Frd4[i]= new double[ntot.y];

double  **Fru1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fru1[i]= new double[ntot.y];

double  **Fru2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fru2[i]= new double[ntot.y];

double  **Fru3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fru3[i]= new double[ntot.y];

double  **Fru4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fru4[i]= new double[ntot.y];

double  **Grd1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Grd1[i]= new double[ntot.y];

double  **Grd2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Grd2[i]= new double[ntot.y];

double  **Grd3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Grd3[i]= new double[ntot.y];

double  **Grd4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Grd4[i]= new double[ntot.y];
 
double  **Gru1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gru1[i]= new double[ntot.y];

double  **Gru2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gru2[i]= new double[ntot.y];

double  **Gru3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gru3[i]= new double[ntot.y];

double  **Gru4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gru4[i]= new double[ntot.y];

double  **U_star1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star1[i]= new double[ntot.y];

double  **U_star2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star2[i]= new double[ntot.y];

double  **U_star3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star3[i]= new double[ntot.y];

double  **U_star4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star4[i]= new double[ntot.y];



for (int ii=0;ii<ntot.x;ii++)
{
    for (int jj=0;jj<ntot.y;jj++)
    {
           P_rla_x[ii][jj]=0.0;
           P_rlb_x[ii][jj]=0.0;
           s_l_d_y[ii][jj]=0.0;
           s_r_d_y[ii][jj]=0.0;
           s_l_u_y[ii][jj]=0.0;
           s_r_u_y[ii][jj]=0.0;
           s_d_l_y[ii][jj]=0.0;
           s_u_l_y[ii][jj]=0.0;
           s_d_r_y[ii][jj]=0.0;
           s_u_r_y[ii][jj]=0.0;
           s_m_d_y[ii][jj]=0.0;
           s_m_u_y[ii][jj]=0.0;
           s_m_l_y[ii][jj]=0.0;
           s_m_r_y[ii][jj]=0.0;
           P_rla_y[ii][jj]=0.0;
           P_rlb_y[ii][jj]=0.0;
           Uld1[ii][jj]=0.0;
           Uld2[ii][jj]=0.0;
           Uld3[ii][jj]=0.0;
           Uld4[ii][jj]=0.0;
           Ulu1[ii][jj]=0.0;
           Ulu2[ii][jj]=0.0;
           Ulu3[ii][jj]=0.0;
           Ulu4[ii][jj]=0.0;
           Fld1[ii][jj]=0.0;
           Fld2[ii][jj]=0.0;
           Fld3[ii][jj]=0.0;
           Fld4[ii][jj]=0.0;
           Flu1[ii][jj]=0.0;
           Flu2[ii][jj]=0.0;
           Flu3[ii][jj]=0.0;
           Flu4[ii][jj]=0.0;
           Gld1[ii][jj]=0.0;
           Gld2[ii][jj]=0.0;
           Gld3[ii][jj]=0.0;
           Gld4[ii][jj]=0.0;
           Glu1[ii][jj]=0.0;
           Glu2[ii][jj]=0.0;
           Glu3[ii][jj]=0.0;
           Glu4[ii][jj]=0.0;
           Urd1[ii][jj]=0.0;
           Urd2[ii][jj]=0.0;
           Urd3[ii][jj]=0.0;
           Urd4[ii][jj]=0.0;
           Uru1[ii][jj]=0.0;
           Uru2[ii][jj]=0.0;
           Uru3[ii][jj]=0.0;
           Uru4[ii][jj]=0.0;
           Frd1[ii][jj]=0.0;
           Frd2[ii][jj]=0.0;
           Frd3[ii][jj]=0.0;
           Frd4[ii][jj]=0.0;
           Fru1[ii][jj]=0.0;
           Fru2[ii][jj]=0.0;
           Fru3[ii][jj]=0.0;
           Fru4[ii][jj]=0.0;
           Grd1[ii][jj]=0.0;
           Grd2[ii][jj]=0.0;
           Grd3[ii][jj]=0.0;
           Grd4[ii][jj]=0.0;
           Gru1[ii][jj]=0.0;
           Gru2[ii][jj]=0.0;
           Gru3[ii][jj]=0.0;
           Gru4[ii][jj]=0.0;
           U_star1[ii][jj]=0.0;
           U_star2[ii][jj]=0.0;
           U_star3[ii][jj]=0.0;
           U_star4[ii][jj]=0.0;
    }
}

 
    wave_speed_estimate_x(n_p_tmp,u_xp_tmp,u_yp_tmp,P_p_tmp,n_m_tmp,u_xm_tmp,u_ym_tmp,P_m_tmp,tmp,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_x,s_u_l_x,s_d_r_x,s_u_r_x,s_m_d_x,s_m_u_x,s_m_l_x,s_m_r_x,P_rla_x,P_rlb_x,mtmp);
wave_speed_estimate_y(n_p_tmp,u_xp_tmp,u_yp_tmp,P_p_tmp,n_m_tmp,u_xm_tmp,u_ym_tmp,P_m_tmp,tmp,s_l_d_y,s_r_d_y,s_l_u_y,s_r_u_y,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_m_d_y,s_m_u_y,s_m_l_y,s_m_r_y,P_rla_y,P_rlb_y,mtmp);
    
   for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            U_F_G_2D(n_p_tmp[i][j],u_xp_tmp[i][j],u_yp_tmp[i][j],P_p_tmp[i][j],n_m_tmp[i][j],u_xm_tmp[i][j],u_ym_tmp[i][j],P_m_tmp[i][j],Uld1[i][j],Uld2[i][j],Uld3[i][j],Uld4[i][j],Urd1[i][j],Urd2[i][j],Urd3[i][j],Urd4[i][j],Fld1[i][j],Fld2[i][j],Fld3[i][j],Fld4[i][j],Frd1[i][j],Frd2[i][j],Frd3[i][j],Frd4[i][j],Gld1[i][j],Gld2[i][j],Gld3[i][j],Gld4[i][j],Grd1[i][j],Grd2[i][j],Grd3[i][j],Grd4[i][j], qtmp, mtmp);
            U_F_G_2D(n_p_tmp[i][j+1],u_xp_tmp[i][j+1],u_yp_tmp[i][j+1],P_p_tmp[i][j+1],n_m_tmp[i][j+1],u_xm_tmp[i][j+1],u_ym_tmp[i][j+1],P_m_tmp[i][j+1],Ulu1[i][j],Ulu2[i][j],Ulu3[i][j],Ulu4[i][j],Uru1[i][j],Uru2[i][j],Uru3[i][j],Uru4[i][j],Flu1[i][j],Flu2[i][j],Flu3[i][j],Flu4[i][j],Fru1[i][j],Fru2[i][j],Fru3[i][j],Fru4[i][j],Glu1[i][j],Glu2[i][j],Glu3[i][j],Glu4[i][j],Gru1[i][j],Gru2[i][j],Gru3[i][j],Gru4[i][j],qtmp, mtmp);        
}
    }
    if(tmp==0)
    {
        //Down
        Flux_2D_HLL(Uld1,Uld2,Uld3,Uld4,Urd1,Urd2,Urd3,Urd4,Fld1,Fld2,Fld3,Fld4,Frd1,Frd2,Frd3,Frd4,Gld1,Gld2,Gld3,Gld4,Grd1,Grd2,Grd3,Grd4,s_r_d_x,s_l_d_x,s_r_d_y,s_l_d_y,F_hll_D_1,G_star_D_1,F_hll_D_2,G_star_D_2,F_hll_D_3,G_star_D_3,F_hll_D_4,G_star_D_4);
    }
    
    if(tmp==1)
    {
        flux_2D_HLLC(Uld1,Uld2,Uld3,Uld4,Urd1,Urd2,Urd3,Urd4,Fld1,Fld2,Fld3,Fld4,Frd1,Frd2,Frd3,Frd4,Gld1,Gld2,Gld3,Gld4,Grd1,Grd2,Grd3,Grd4,s_r_d_x,s_l_d_x,s_m_d_x,F_hll_D_1,G_star_D_1,F_hll_D_2,G_star_D_2,F_hll_D_3,G_star_D_3,F_hll_D_4,G_star_D_4,P_rla_x,0,mtmp);
    }
    U_HLL(Urd1,Uld1,Frd1,Fld1,s_r_d_x,s_l_d_x,U_star1);
    U_HLL(Urd2,Uld2,Frd2,Fld2,s_r_d_x,s_l_d_x,U_star2);
    U_HLL(Urd3,Uld3,Frd3,Fld3,s_r_d_x,s_l_d_x,U_star3);
    U_HLL(Urd4,Uld4,Frd4,Fld4,s_r_d_x,s_l_d_x,U_star4);
    G_trans_HLL(U_star1,U_star2,U_star3,U_star4,F_hll_D_1,F_hll_D_2,F_hll_D_3,F_hll_D_4,G_star_D_1,0);
    G_trans_HLL(U_star1,U_star2,U_star3,U_star4,F_hll_D_1,F_hll_D_2,F_hll_D_3,F_hll_D_4,G_star_D_2,1);
    G_trans_HLL(U_star1,U_star2,U_star3,U_star4,F_hll_D_1,F_hll_D_2,F_hll_D_3,F_hll_D_4,G_star_D_3,2);
    G_trans_HLL(U_star1,U_star2,U_star3,U_star4,F_hll_D_1,F_hll_D_2,F_hll_D_3,F_hll_D_4,G_star_D_4,3);
    
    //Up
    if(tmp==0)
    {
        Flux_2D_HLL(Ulu1,Ulu2,Ulu3,Ulu4,Uru1,Uru2,Uru3,Uru4,Flu1,Flu2,Flu3,Flu4,Fru1,Fru2,Fru3,Fru4,Glu1,Glu2,Glu3,Glu4,Gru1,Gru2,Gru3,Gru4,s_r_u_x,s_l_u_x,s_r_u_y,s_l_u_y,F_hll_U_1,G_star_U_1,F_hll_U_2,G_star_U_2,F_hll_U_3,G_star_U_3,F_hll_U_4,G_star_U_4);
    }
    if(tmp==1)
    {
        flux_2D_HLLC(Ulu1,Ulu2,Ulu3,Ulu4,Uru1,Uru2,Uru3,Uru4,Flu1,Flu2,Flu3,Flu4,Fru1,Fru2,Fru3,Fru4,Glu1,Glu2,Glu3,Glu4,Gru1,Gru2,Gru3,Gru4,s_r_u_x,s_l_u_x,s_m_u_x,F_hll_U_1,G_star_U_1,F_hll_U_2,G_star_U_2,F_hll_U_3,G_star_U_3,F_hll_U_4,G_star_U_4,P_rlb_x,0,mtmp);
    }
    U_HLL(Uru1,Ulu1,Fru1,Flu1,s_r_u_x,s_l_u_x,U_star1);
    U_HLL(Uru2,Ulu2,Fru2,Flu2,s_r_u_x,s_l_u_x,U_star2);
    U_HLL(Uru3,Ulu3,Fru3,Flu3,s_r_u_x,s_l_u_x,U_star3);
    U_HLL(Uru4,Ulu4,Fru4,Flu4,s_r_u_x,s_l_u_x,U_star4);
    G_trans_HLL(U_star1,U_star2,U_star3,U_star4,F_hll_U_1,F_hll_U_2,F_hll_U_3,F_hll_U_4,G_star_U_1,0);
    G_trans_HLL(U_star1,U_star2,U_star3,U_star4,F_hll_U_1,F_hll_U_2,F_hll_U_3,F_hll_U_4,G_star_U_2,1);
    G_trans_HLL(U_star1,U_star2,U_star3,U_star4,F_hll_U_1,F_hll_U_2,F_hll_U_3,F_hll_U_4,G_star_U_3,2);
    G_trans_HLL(U_star1,U_star2,U_star3,U_star4,F_hll_U_1,F_hll_U_2,F_hll_U_3,F_hll_U_4,G_star_U_4,3);



for (int i = 0; i < ntot.x; i++)
delete[] P_rla_x [i] ; 
delete[] P_rla_x;

for (int i = 0; i < ntot.x; i++)
delete[] P_rlb_x [i] ;
delete[] P_rlb_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_l_d_y [i] ;
delete[] s_l_d_y; 

for (int i = 0; i < ntot.x; i++)
delete[] s_r_d_y [i] ;
delete[] s_r_d_y;

for (int i = 0; i < ntot.x; i++)
delete[] s_l_u_y [i] ;
delete[] s_l_u_y;

for (int i = 0; i < ntot.x; i++)
delete[] s_r_u_y [i] ;
delete[] s_r_u_y;

for (int i = 0; i < ntot.x; i++)
delete[] s_d_l_y [i] ;
delete[] s_d_l_y;

for (int i = 0; i < ntot.x; i++)
delete[] s_u_l_y [i] ;
delete[] s_u_l_y;


for (int i = 0; i < ntot.x; i++)
delete[] s_d_r_y [i] ;
delete[] s_d_r_y;

for (int i = 0; i < ntot.x; i++)
delete[] s_u_r_y [i] ;
delete[] s_u_r_y;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_d_y [i] ;
delete[] s_m_d_y;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_u_y [i] ;
delete[] s_m_u_y;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_l_y [i] ;
delete[] s_m_l_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_r_y [i] ;
delete[] s_m_r_y;

for (int i = 0; i < ntot.x; i++)
delete[] P_rla_y [i] ;
delete[] P_rla_y;

for (int i = 0; i < ntot.x; i++)
delete[] P_rlb_y [i] ;
delete[] P_rlb_y;

for (int i = 0; i < ntot.x; i++)
delete[] Uld1[i] ;
delete[] Uld1 ;


for (int i = 0; i < ntot.x; i++)
delete[] Uld2[i];
delete[] Uld2 ;

for (int i = 0; i < ntot.x; i++)
delete[] Uld3[i] ;
delete[] Uld3;

for (int i = 0; i < ntot.x; i++)
delete[] Uld4[i] ;
delete[] Uld4;
 
for (int i = 0; i < ntot.x; i++)
delete[] Ulu1[i] ;
delete[] Ulu1;
 
for (int i = 0; i < ntot.x; i++)
delete[] Ulu2[i] ;
delete[] Ulu2;

for (int i = 0; i < ntot.x; i++)
delete[] Ulu3[i] ;
delete[] Ulu3;

for (int i = 0; i < ntot.x; i++)
delete[] Ulu4[i] ;
delete[] Ulu4;

for (int i = 0; i < ntot.x; i++)
delete[] Fld1[i] ;
delete[] Fld1;
 
for (int i = 0; i < ntot.x; i++)
delete[] Fld2[i] ;
delete[] Fld2;

for (int i = 0; i < ntot.x; i++)
delete[] Fld3[i] ;
delete[] Fld3;
 
for (int i = 0; i < ntot.x; i++)
delete[] Fld4[i] ;
delete[] Fld4;
 
for (int i = 0; i < ntot.x; i++)
delete[] Flu1[i] ;
delete[] Flu1;

for (int i = 0; i < ntot.x; i++)
delete[] Flu2[i] ;
delete[] Flu2;

for (int i = 0; i < ntot.x; i++)
delete[] Flu3[i] ;
delete[] Flu3;

for (int i = 0; i < ntot.x; i++)
delete[] Flu4[i] ;
delete[] Flu4;

for (int i = 0; i < ntot.x; i++)
delete[] Gld1[i] ;
delete[] Gld1;

for (int i = 0; i < ntot.x; i++)
delete[] Gld2[i] ;
delete[] Gld2;

for (int i = 0; i < ntot.x; i++)
delete[] Gld3[i] ;
delete[] Gld3;

for (int i = 0; i < ntot.x; i++)
delete[] Gld4[i] ;
delete[] Gld4;


for (int i = 0; i < ntot.x; i++)
delete[] Glu1[i] ;
delete[] Glu1;

for (int i = 0; i < ntot.x; i++)
delete[] Glu2[i] ;
delete[] Glu2;

for (int i = 0; i < ntot.x; i++)
delete[] Glu3[i] ;
delete[] Glu3;

for (int i = 0; i < ntot.x; i++)
delete[] Glu4[i] ;
delete[] Glu4;

for (int i = 0; i < ntot.x; i++)
delete[] Urd1[i] ;
delete[] Urd1;

for (int i = 0; i < ntot.x; i++)
delete[] Urd2[i] ;
delete[] Urd2;

for (int i = 0; i < ntot.x; i++)
delete[] Urd3[i] ;
delete[] Urd3;

for (int i = 0; i < ntot.x; i++)
delete[] Urd4[i] ;
delete[] Urd4;

for (int i = 0; i < ntot.x; i++)
delete[] Uru1[i] ;
delete[] Uru1;

for (int i = 0; i < ntot.x; i++)
delete[] Uru2[i] ;
delete[] Uru2;

for (int i = 0; i < ntot.x; i++)
delete[] Uru3[i] ;
delete[] Uru3;

for (int i = 0; i < ntot.x; i++)
delete[] Uru4[i] ;
delete[] Uru4;

for (int i = 0; i < ntot.x; i++)
delete[] Frd1[i] ;
delete[] Frd1;
 
for (int i = 0; i < ntot.x; i++)
delete[] Frd2[i] ;
delete[] Frd2;

for (int i = 0; i < ntot.x; i++)
delete[] Frd3[i] ;
delete[] Frd3;

for (int i = 0; i < ntot.x; i++)
delete[] Frd4[i] ;
delete[] Frd4;

for (int i = 0; i < ntot.x; i++)
delete[] Fru1[i] ;
delete[] Fru1;

for (int i = 0; i < ntot.x; i++)
delete[] Fru2[i] ;
delete[] Fru2;

for (int i = 0; i < ntot.x; i++)
delete[] Fru3[i] ;
delete[] Fru3;

for (int i = 0; i < ntot.x; i++)
delete[] Fru4[i] ;
delete[] Fru4;

for (int i = 0; i < ntot.x; i++)
delete[] Grd1[i] ;
delete[] Grd1;

for (int i = 0; i < ntot.x; i++)
delete[] Grd2[i] ;
delete[] Grd2;

for (int i = 0; i < ntot.x; i++)
delete[] Grd3[i] ;
delete[] Grd3;

for (int i = 0; i < ntot.x; i++)
delete[] Grd4[i] ;
delete[] Grd4;
 
for (int i = 0; i < ntot.x; i++)
delete[] Gru1[i] ;
delete[] Gru1;

for (int i = 0; i < ntot.x; i++)
delete[] Gru2[i] ;
delete[] Gru2;

for (int i = 0; i < ntot.x; i++)
delete[] Gru3[i] ;
delete[] Gru3;

for (int i = 0; i < ntot.x; i++)
delete[] Gru4[i] ;
delete[] Gru4;

for (int i = 0; i < ntot.x; i++)
delete[] U_star1[i] ;
delete[] U_star1;

for (int i = 0; i < ntot.x; i++)
delete[] U_star2[i] ;
delete[] U_star2;

for (int i = 0; i < ntot.x; i++)
delete[] U_star3[i] ;
delete[] U_star3;

for (int i = 0; i < ntot.x; i++)
delete[] U_star4[i] ;
delete[] U_star4;

}
void fluxes_y(double **n_p_tmp,double **u_xp_tmp,double **u_yp_tmp, double **P_p_tmp, double **n_m_tmp,double **u_xm_tmp, double **u_ym_tmp, double **P_m_tmp ,int tmp, double **F_star_L_1,double **G_hll_L_1,double **F_star_L_2,double **G_hll_L_2,double **F_star_L_3,double **G_hll_L_3,double **F_star_L_4,double **G_hll_L_4,double **F_star_R_1,double **G_hll_R_1,double **F_star_R_2,double **G_hll_R_2,double **F_star_R_3,double **G_hll_R_3,double **F_star_R_4,double **G_hll_R_4, double  **s_l_d_y, double  **s_r_d_y,double  **s_l_u_y, double  **s_r_u_y,double  **s_d_l_y,double  **s_u_l_y,double   **s_d_r_y,double   **s_u_r_y,double  **s_m_d_y,double  **s_m_u_y,double  **s_m_l_y,double  **s_m_r_y, double qtmp, double mtmp )
{    

double  **s_l_d_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_l_d_x [i]= new double[ntot.y];

double  **s_r_d_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_r_d_x [i]= new double[ntot.y];

double  **s_l_u_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_l_u_x [i]= new double[ntot.y];

double  **s_r_u_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_r_u_x [i]= new double[ntot.y];

double  **s_d_l_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_d_l_x [i]= new double[ntot.y];

double  **s_u_l_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_u_l_x [i]= new double[ntot.y];

double  **s_d_r_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_d_r_x [i]= new double[ntot.y];

double  **s_u_r_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_u_r_x [i]= new double[ntot.y];

double  **s_m_d_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_d_x [i]= new double[ntot.y];

double  **s_m_u_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_u_x [i]= new double[ntot.y];

double  **s_m_l_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_l_x [i]= new double[ntot.y];

double  **s_m_r_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_r_x [i]= new double[ntot.y];

double  **P_rla_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_rla_x [i]= new double[ntot.y];

double  **P_rlb_x= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_rlb_x [i]= new double[ntot.y];

double  **P_rla_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_rla_y [i]= new double[ntot.y];

double  **P_rlb_y= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_rlb_y [i]= new double[ntot.y];

double  **Uld1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uld1[i]= new double[ntot.y];

double  **Uld2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uld2[i]= new double[ntot.y];

double  **Uld3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uld3[i]= new double[ntot.y];

double  **Uld4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uld4[i]= new double[ntot.y];

double  **Ulu1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Ulu1[i]= new double[ntot.y];

double  **Ulu2= new double *[ntot.x];
 for (int i = 0; i < ntot.x; i++)
Ulu2[i]= new double[ntot.y];

double  **Ulu3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Ulu3[i]= new double[ntot.y];

double  **Ulu4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Ulu4[i]= new double[ntot.y];

double  **Fld1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fld1[i]= new double[ntot.y];

double  **Fld2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fld2[i]= new double[ntot.y];

double  **Fld3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fld3[i]= new double[ntot.y];

double  **Fld4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fld4[i]= new double[ntot.y];

double  **Flu1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Flu1[i]= new double[ntot.y];

double  **Flu2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Flu2[i]= new double[ntot.y];

double  **Flu3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Flu3[i]= new double[ntot.y];

double  **Flu4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Flu4[i]= new double[ntot.y];

double  **Gld1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gld1[i]= new double[ntot.y];

double  **Gld2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gld2[i]= new double[ntot.y];

double  **Gld3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gld3[i]= new double[ntot.y];

double  **Gld4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gld4[i]= new double[ntot.y];

double  **Glu1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Glu1[i]= new double[ntot.y];

double  **Glu2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Glu2[i]= new double[ntot.y];

double  **Glu3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Glu3[i]= new double[ntot.y];

double  **Glu4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Glu4[i]= new double[ntot.y];

double  **Urd1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Urd1[i]= new double[ntot.y];

double  **Urd2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Urd2[i]= new double[ntot.y];

double  **Urd3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Urd3[i]= new double[ntot.y];

double  **Urd4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Urd4[i]= new double[ntot.y];

double  **Uru1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uru1[i]= new double[ntot.y];

double  **Uru2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uru2[i]= new double[ntot.y];

double  **Uru3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uru3[i]= new double[ntot.y];

double  **Uru4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Uru4[i]= new double[ntot.y];

double  **Frd1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Frd1[i]= new double[ntot.y];

double  **Frd2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Frd2[i]= new double[ntot.y];

double  **Frd3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Frd3[i]= new double[ntot.y];

double  **Frd4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Frd4[i]= new double[ntot.y];

double  **Fru1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fru1[i]= new double[ntot.y];

double  **Fru2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fru2[i]= new double[ntot.y];

double  **Fru3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fru3[i]= new double[ntot.y];

double  **Fru4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Fru4[i]= new double[ntot.y];

double  **Grd1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Grd1[i]= new double[ntot.y];

double  **Grd2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Grd2[i]= new double[ntot.y];

double  **Grd3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Grd3[i]= new double[ntot.y];

double  **Grd4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Grd4[i]= new double[ntot.y];

double  **Gru1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gru1[i]= new double[ntot.y];

double  **Gru2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gru2[i]= new double[ntot.y];

double  **Gru3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gru3[i]= new double[ntot.y];

double  **Gru4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Gru4[i]= new double[ntot.y];

double  **U_star1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star1[i]= new double[ntot.y];

double  **U_star2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star2[i]= new double[ntot.y];

double  **U_star3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star3[i]= new double[ntot.y];

double  **U_star4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_star4[i]= new double[ntot.y];

for (int ii=0;ii<ntot.x;ii++)
{   
    for (int jj=0;jj<ntot.y;jj++)
    {    
         s_l_d_x[ii][jj]=0.0;
         s_r_d_x[ii][jj]=0.0;
         s_l_u_x[ii][jj]=0.0;
         s_r_u_x[ii][jj]=0.0;
         s_d_l_x[ii][jj]=0.0;
         s_u_l_x[ii][jj]=0.0;
         s_d_r_x[ii][jj]=0.0;
         s_u_r_x[ii][jj]=0.0;
         s_m_d_x[ii][jj]=0.0;
         s_m_u_x[ii][jj]=0.0;
         s_m_l_x[ii][jj]=0.0;
         s_m_r_x[ii][jj]=0.0;
         P_rla_x[ii][jj]=0.0;
         P_rlb_x[ii][jj]=0.0;
         P_rla_y[ii][jj]=0.0;
         P_rlb_y[ii][jj]=0.0;
         Uld1[ii][jj]=0.0;
         Uld2[ii][jj]=0.0;
         Uld3[ii][jj]=0.0;
         Uld4[ii][jj]=0.0;
         Ulu1[ii][jj]=0.0;
         Ulu2[ii][jj]=0.0;
         Ulu3[ii][jj]=0.0;
         Ulu4[ii][jj]=0.0;
         Fld1[ii][jj]=0.0;
         Fld2[ii][jj]=0.0;
         Fld3[ii][jj]=0.0;
         Fld4[ii][jj]=0.0;
         Flu1[ii][jj]=0.0;
         Flu2[ii][jj]=0.0;
         Flu3[ii][jj]=0.0;
         Flu4[ii][jj]=0.0;
         Gld1[ii][jj]=0.0;
         Gld2[ii][jj]=0.0;
         Gld3[ii][jj]=0.0;
         Gld4[ii][jj]=0.0;
         Glu1[ii][jj]=0.0;
         Glu2[ii][jj]=0.0;
         Glu3[ii][jj]=0.0;
         Glu4[ii][jj]=0.0;
         Urd1[ii][jj]=0.0;
         Urd2[ii][jj]=0.0;
         Urd3[ii][jj]=0.0;
         Urd4[ii][jj]=0.0;
         Uru1[ii][jj]=0.0;
         Uru2[ii][jj]=0.0;
         Uru3[ii][jj]=0.0;
         Uru4[ii][jj]=0.0;
         Frd1[ii][jj]=0.0;
         Frd2[ii][jj]=0.0;
         Frd3[ii][jj]=0.0;
         Frd4[ii][jj]=0.0;
         Fru1[ii][jj]=0.0;
         Fru2[ii][jj]=0.0;
         Fru3[ii][jj]=0.0;
         Fru4[ii][jj]=0.0;
         Grd1[ii][jj]=0.0;
         Grd2[ii][jj]=0.0;
         Grd3[ii][jj]=0.0;
         Grd4[ii][jj]=0.0;
         Gru1[ii][jj]=0.0;
         Gru2[ii][jj]=0.0;
         Gru3[ii][jj]=0.0;
         Gru4[ii][jj]=0.0;
         U_star1[ii][jj]=0.0;
         U_star2[ii][jj]=0.0;
         U_star3[ii][jj]=0.0;
         U_star4[ii][jj]=0.0;
    }
}

    wave_speed_estimate_x(n_p_tmp,u_xp_tmp,u_yp_tmp,P_p_tmp,n_m_tmp,u_xm_tmp,u_ym_tmp,P_m_tmp,tmp,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_x,s_u_l_x,s_d_r_x,s_u_r_x,s_m_d_x,s_m_u_x,s_m_l_x,s_m_r_x,P_rla_x,P_rlb_x,mtmp);
    wave_speed_estimate_y(n_p_tmp,u_xp_tmp,u_yp_tmp,P_p_tmp,n_m_tmp,u_xm_tmp,u_ym_tmp,P_m_tmp,tmp,s_l_d_y,s_r_d_y,s_l_u_y,s_r_u_y,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_m_d_y,s_m_u_y,s_m_l_y,s_m_r_y,P_rla_y,P_rlb_y,mtmp);
    
    for(int j=nbuffer-1;j<=ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer-1;i<=ntot.x-nbuffer;i++)
        {
            
            U_F_G_2D(n_p_tmp[i][j],u_xp_tmp[i][j],u_yp_tmp[i][j],P_p_tmp[i][j],n_m_tmp[i][j],u_xm_tmp[i][j],u_ym_tmp[i][j],P_m_tmp[i][j],Uld1[i][j],Uld2[i][j],Uld3[i][j],Uld4[i][j],Ulu1[i][j],Ulu2[i][j],Ulu3[i][j],Ulu4[i][j],Fld1[i][j],Fld2[i][j],Fld3[i][j],Fld4[i][j],Flu1[i][j],Flu2[i][j],Flu3[i][j],Flu4[i][j],Gld1[i][j],Gld2[i][j],Gld3[i][j],Gld4[i][j],Glu1[i][j],Glu2[i][j],Glu3[i][j],Glu4[i][j],qtmp, mtmp);
            U_F_G_2D(n_p_tmp[i+1][j],u_xp_tmp[i+1][j],u_yp_tmp[i+1][j],P_p_tmp[i+1][j],n_m_tmp[i+1][j],u_xm_tmp[i+1][j],u_ym_tmp[i+1][j],P_m_tmp[i+1][j],Urd1[i][j],Urd2[i][j],Urd3[i][j],Urd4[i][j],Uru1[i][j],Uru2[i][j],Uru3[i][j],Uru4[i][j],Frd1[i][j],Frd2[i][j],Frd3[i][j],Frd4[i][j],Fru1[i][j],Fru2[i][j],Fru3[i][j],Fru4[i][j],Grd1[i][j],Grd2[i][j],Grd3[i][j],Grd4[i][j],Gru1[i][j],Gru2[i][j],Gru3[i][j],Gru4[i][j],qtmp, mtmp);
        }
    }
    if(tmp==0)
    {
        Flux_2D_HLL(Uld1,Uld2,Uld3,Uld4,Ulu1,Ulu2,Ulu3,Ulu4,Fld1,Fld2,Fld3,Fld4,Flu1,Flu2,Flu3,Flu4,Gld1,Gld2,Gld3,Gld4,Glu1,Glu2,Glu3,Glu4,s_u_l_x,s_d_l_x,s_u_l_y,s_d_l_y,F_star_L_1,G_hll_L_1,F_star_L_2,G_hll_L_2,F_star_L_3,G_hll_L_3,F_star_L_4,G_hll_L_4);
    }
    else if(tmp==1)
    {
        flux_2D_HLLC(Uld1,Uld2,Uld3,Uld4,Ulu1,Ulu2,Ulu3,Ulu4,Fld1,Fld2,Fld3,Fld4,Flu1,Flu2,Flu3,Flu4,Gld1,Gld2,Gld3,Gld4,Glu1,Glu2,Glu3,Glu4,s_u_l_y,s_d_l_y,s_m_l_y,F_star_L_1,G_hll_L_1,F_star_L_2,G_hll_L_2,F_star_L_3,G_hll_L_3,F_star_L_4,G_hll_L_4,P_rla_y,1,mtmp);    
    }
    else
    {
        cout <<"NA\n";
    }
    U_HLL(Ulu1,Uld1,Glu1,Gld1,s_u_l_y,s_d_l_y,U_star1);
    U_HLL(Ulu2,Uld2,Glu2,Gld2,s_u_l_y,s_d_l_y,U_star2);
    U_HLL(Ulu3,Uld3,Glu3,Gld3,s_u_l_y,s_d_l_y,U_star3);
    U_HLL(Ulu4,Uld4,Glu4,Gld4,s_u_l_y,s_d_l_y,U_star4);
    F_trans_HLL(U_star1,U_star2,U_star3,U_star4,G_hll_L_1,G_hll_L_2,G_hll_L_3,G_hll_L_4,F_star_L_1,0);
    F_trans_HLL(U_star1,U_star2,U_star3,U_star4,G_hll_L_1,G_hll_L_2,G_hll_L_3,G_hll_L_4,F_star_L_2,1);
    F_trans_HLL(U_star1,U_star2,U_star3,U_star4,G_hll_L_1,G_hll_L_2,G_hll_L_3,G_hll_L_4,F_star_L_3,2);
    F_trans_HLL(U_star1,U_star2,U_star3,U_star4,G_hll_L_1,G_hll_L_2,G_hll_L_3,G_hll_L_4,F_star_L_4,3);
  
    if(tmp==0)
    {
        //Right
        Flux_2D_HLL(Urd1,Urd2,Urd3,Urd4,Uru1,Uru2,Uru3,Uru4,Frd1,Frd2,Frd3,Frd4,Fru1,Fru2,Fru3,Fru4,Grd1,Grd2,Grd3,Grd4,Gru1,Gru2,Gru3,Gru4,s_u_r_x,s_d_r_x,s_u_r_y,s_d_r_y,F_star_R_1,G_hll_R_1,F_star_R_2,G_hll_R_2,F_star_R_3,G_hll_R_3,F_star_R_4,G_hll_R_4);
    }
    else if(tmp==1)
    {
        flux_2D_HLLC(Urd1,Urd2,Urd3,Urd4,Uru1,Uru2,Uru3,Uru4,Frd1,Frd2,Frd3,Frd4,Fru1,Fru2,Fru3,Fru4,Grd1,Grd2,Grd3,Grd4,Gru1,Gru2,Gru3,Gru4,s_u_r_y,s_d_r_y,s_m_r_y,F_star_R_1,G_hll_R_1,F_star_R_2,G_hll_R_2,F_star_R_3,G_hll_R_3,F_star_R_4,G_hll_R_4,P_rlb_y,1,mtmp);
    }
    else{
        cout <<"NA\n";
    }
    U_HLL(Uru1,Urd1,Gru1,Grd1,s_u_r_y,s_d_r_y,U_star1);
    U_HLL(Uru2,Urd2,Gru2,Grd2,s_u_r_y,s_d_r_y,U_star2);
    U_HLL(Uru3,Urd3,Gru3,Grd3,s_u_r_y,s_d_r_y,U_star3);
    U_HLL(Uru4,Urd4,Gru4,Grd4,s_u_r_y,s_d_r_y,U_star4);
    F_trans_HLL(U_star1,U_star2,U_star3,U_star4,G_hll_R_1,G_hll_R_2,G_hll_R_3,G_hll_R_4,F_star_R_1,0);
    F_trans_HLL(U_star1,U_star2,U_star3,U_star4,G_hll_R_1,G_hll_R_2,G_hll_R_3,G_hll_R_4,F_star_R_2,1);
    F_trans_HLL(U_star1,U_star2,U_star3,U_star4,G_hll_R_1,G_hll_R_2,G_hll_R_3,G_hll_R_4,F_star_R_3,2);
    F_trans_HLL(U_star1,U_star2,U_star3,U_star4,G_hll_R_1,G_hll_R_2,G_hll_R_3,G_hll_R_4,F_star_R_4,3);


for (int i = 0; i < ntot.x; i++)
delete[] s_l_d_x [i] ;
delete[] s_l_d_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_r_d_x [i] ;
delete[] s_r_d_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_l_u_x [i] ;
delete[] s_l_u_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_r_u_x [i] ;
delete[] s_r_u_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_d_l_x [i] ;
delete[] s_d_l_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_u_l_x [i] ;
delete[] s_u_l_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_d_r_x [i] ;
delete[] s_d_r_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_u_r_x [i] ;
delete[] s_u_r_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_d_x [i] ;
delete[] s_m_d_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_u_x [i] ;
delete[] s_m_u_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_l_x [i] ;
delete[] s_m_l_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_r_x [i] ;
delete[] s_m_r_x;

for (int i = 0; i < ntot.x; i++)
delete[] P_rla_x [i] ;
delete[] P_rla_x;

for (int i = 0; i < ntot.x; i++)
delete[] P_rlb_x [i] ;
delete[] P_rlb_x;

for (int i = 0; i < ntot.x; i++)
delete[] P_rla_y [i] ;
delete[] P_rla_y;

for (int i = 0; i < ntot.x; i++)
delete[] P_rlb_y [i] ;
delete[] P_rlb_y;

for (int i = 0; i < ntot.x; i++)
delete[] Uld1[i] ;
delete[] Uld1 ;


for (int i = 0; i < ntot.x; i++)
delete[] Uld2[i];
delete[] Uld2 ;

for (int i = 0; i < ntot.x; i++)
delete[] Uld3[i] ;
delete[] Uld3;

for (int i = 0; i < ntot.x; i++)
delete[] Uld4[i] ;
delete[] Uld4;

for (int i = 0; i < ntot.x; i++)
delete[] Ulu1[i] ;
delete[] Ulu1;

for (int i = 0; i < ntot.x; i++)
delete[] Ulu2[i] ;
delete[] Ulu2;

for (int i = 0; i < ntot.x; i++)
delete[] Ulu3[i] ;
delete[] Ulu3;

for (int i = 0; i < ntot.x; i++)
delete[] Ulu4[i] ;
delete[] Ulu4;

for (int i = 0; i < ntot.x; i++)
delete[] Fld1[i] ;
delete[] Fld1;

for (int i = 0; i < ntot.x; i++)
delete[] Fld2[i] ;
delete[] Fld2;

for (int i = 0; i < ntot.x; i++)
delete[] Fld3[i] ;
delete[] Fld3;

for (int i = 0; i < ntot.x; i++)
delete[] Fld4[i] ;
delete[] Fld4;

for (int i = 0; i < ntot.x; i++)
delete[] Flu1[i] ;
delete[] Flu1;

for (int i = 0; i < ntot.x; i++)
delete[] Flu2[i] ;
delete[] Flu2;

for (int i = 0; i < ntot.x; i++)
delete[] Flu3[i] ;
delete[] Flu3;

for (int i = 0; i < ntot.x; i++)
delete[] Flu4[i] ;
delete[] Flu4;

for (int i = 0; i < ntot.x; i++)
delete[] Gld1[i] ;
delete[] Gld1;

for (int i = 0; i < ntot.x; i++)
delete[] Gld2[i] ;
delete[] Gld2;

for (int i = 0; i < ntot.x; i++)
delete[] Gld3[i] ;
delete[] Gld3;

for (int i = 0; i < ntot.x; i++)
delete[] Gld4[i] ;
delete[] Gld4;


for (int i = 0; i < ntot.x; i++)
delete[] Glu1[i] ;
delete[] Glu1;

for (int i = 0; i < ntot.x; i++)
delete[] Glu2[i] ;
delete[] Glu2;

for (int i = 0; i < ntot.x; i++)
delete[] Glu3[i] ;
delete[] Glu3;

for (int i = 0; i < ntot.x; i++)
delete[] Glu4[i] ;
delete[] Glu4;

for (int i = 0; i < ntot.x; i++)
delete[] Urd1[i] ;
delete[] Urd1;

for (int i = 0; i < ntot.x; i++)
delete[] Urd2[i] ;
delete[] Urd2;

for (int i = 0; i < ntot.x; i++)
delete[] Urd3[i] ;
delete[] Urd3;

for (int i = 0; i < ntot.x; i++)
delete[] Urd4[i] ;
delete[] Urd4;

for (int i = 0; i < ntot.x; i++)
delete[] Uru1[i] ;
delete[] Uru1;

for (int i = 0; i < ntot.x; i++)
delete[] Uru2[i] ;
delete[] Uru2;

for (int i = 0; i < ntot.x; i++)
delete[] Uru3[i] ;
delete[] Uru3;

for (int i = 0; i < ntot.x; i++)
delete[] Uru4[i] ;
delete[] Uru4;

for (int i = 0; i < ntot.x; i++)
delete[] Frd1[i] ;
delete[] Frd1;

for (int i = 0; i < ntot.x; i++)
delete[] Frd2[i] ;
delete[] Frd2;

for (int i = 0; i < ntot.x; i++)
delete[] Frd3[i] ;
delete[] Frd3;

for (int i = 0; i < ntot.x; i++)
delete[] Frd4[i] ;
delete[] Frd4;

for (int i = 0; i < ntot.x; i++)
delete[] Fru1[i] ;
delete[] Fru1;

for (int i = 0; i < ntot.x; i++)
delete[] Fru2[i] ;
delete[] Fru2;

for (int i = 0; i < ntot.x; i++)
delete[] Fru3[i] ;
delete[] Fru3;

for (int i = 0; i < ntot.x; i++)
delete[] Fru4[i] ;
delete[] Fru4;

for (int i = 0; i < ntot.x; i++)
delete[] Grd1[i] ;
delete[] Grd1;

for (int i = 0; i < ntot.x; i++)
delete[] Grd2[i] ;
delete[] Grd2;

for (int i = 0; i < ntot.x; i++)
delete[] Grd3[i] ;
delete[] Grd3;

for (int i = 0; i < ntot.x; i++)
delete[] Grd4[i] ;
delete[] Grd4;

for (int i = 0; i < ntot.x; i++)
delete[] Gru1[i] ;
delete[] Gru1;

for (int i = 0; i < ntot.x; i++)
delete[] Gru2[i] ;
delete[] Gru2;

for (int i = 0; i < ntot.x; i++)
delete[] Gru3[i] ;
delete[] Gru3;

for (int i = 0; i < ntot.x; i++)
delete[] Gru4[i] ;
delete[] Gru4;

for (int i = 0; i < ntot.x; i++)
delete[] U_star1[i] ;
delete[] U_star1;

for (int i = 0; i < ntot.x; i++)
delete[] U_star2[i] ;
delete[] U_star2;

for (int i = 0; i < ntot.x; i++)
delete[] U_star3[i] ;
delete[] U_star3;

for (int i = 0; i < ntot.x; i++)
delete[] U_star4[i] ;
delete[] U_star4;

}
void fluxes(double **ni_tmp, double **ux_tmp,double **uy_tmp,double **pi_tmp,double ** U_rk_tmp_1_0,double ** U_rk_tmp_2_0 ,double ** U_rk_tmp_3_0,double ** U_rk_tmp_4_0, int bc_x, int bc_y,int my_id, double qtmp, double mtmp,bool notWest, bool notEast,bool notSouth,bool notNorth, int specie, double * flux_cath)
{

double *kf1_end,*kf2_end,*kf3_end,*kf4_end;
double *kf1_in,*kf2_in,*kf3_in,*kf4_in;
double kf_tmp;
double *diffn;
kf1_end =  (double*)malloc(ntot.y*sizeof(double)); 
kf2_end =  (double*)malloc(ntot.y*sizeof(double));
kf3_end =  (double*)malloc(ntot.y*sizeof(double));
kf4_end =  (double*)malloc(ntot.y*sizeof(double));

kf1_in =  (double*)malloc(ntot.y*sizeof(double));
kf2_in =  (double*)malloc(ntot.y*sizeof(double));
kf3_in =  (double*)malloc(ntot.y*sizeof(double));
kf4_in =  (double*)malloc(ntot.y*sizeof(double));


diffn = (double*)malloc(ntot.y*sizeof(double));

double  **n_p = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
n_p[i] = new double[ntot.y];

double  **n_m = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
n_m[i] = new double[ntot.y];

double  **u_xp = new double *[ntot.x]; 
for (int i = 0; i < ntot.x; i++)
u_xp[i] = new double[ntot.y];

double  **u_xm = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
u_xm[i] = new double[ntot.y];

double  **u_yp = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
u_yp[i] = new double[ntot.y];

double  ** u_ym = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
u_ym[i] = new double[ntot.y];

double  **P_p = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_p[i] = new double[ntot.y];

double  **P_m = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_m[i] = new double[ntot.y];

double  **F_hll_D_1 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_hll_D_1[i] = new double[ntot.y];

double  **G_star_D_1 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_star_D_1[i] = new double[ntot.y];

double  **F_hll_D_2 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_hll_D_2[i] = new double[ntot.y];

double  **G_star_D_2 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_star_D_2[i] = new double[ntot.y];

double  **F_hll_D_3 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_hll_D_3[i] = new double[ntot.y];

double  **G_star_D_3 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_star_D_3[i] = new double[ntot.y];

double  **F_hll_D_4 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_hll_D_4[i] = new double[ntot.y];

double  **G_star_D_4 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_star_D_4[i] = new double[ntot.y];

double  **F_hll_U_1 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_hll_U_1[i] = new double[ntot.y];

double  **G_star_U_1 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_star_U_1[i] = new double[ntot.y];

double  **F_hll_U_2 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_hll_U_2[i] = new double[ntot.y];

double  **G_star_U_2 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_star_U_2[i] = new double[ntot.y];

double  **F_hll_U_3 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_hll_U_3[i] = new double[ntot.y];

double  **G_star_U_3 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_star_U_3[i] = new double[ntot.y];

double  **F_hll_U_4 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_hll_U_4[i] = new double[ntot.y];

double  **G_star_U_4 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_star_U_4[i] = new double[ntot.y];

double  **F_star_R_1 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_star_R_1[i] = new double[ntot.y];

double  **G_hll_R_1 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_hll_R_1[i] = new double[ntot.y];

double  **F_star_R_2 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_star_R_2[i] = new double[ntot.y];

double  **G_hll_R_2 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_hll_R_2[i] = new double[ntot.y];

double  **F_star_R_3 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_star_R_3[i] = new double[ntot.y];

double  **G_hll_R_3 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_hll_R_3[i] = new double[ntot.y];

double  **F_star_R_4 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_star_R_4[i] = new double[ntot.y];

double  **G_hll_R_4 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_hll_R_4[i] = new double[ntot.y];

double  **G_hll_L_1 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_hll_L_1[i] = new double[ntot.y];

double  **F_star_L_1 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_star_L_1[i] = new double[ntot.y];

double  **G_hll_L_2 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_hll_L_2[i] = new double[ntot.y];

double  **F_star_L_2 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_star_L_2[i] = new double[ntot.y];

double  **G_hll_L_3 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_hll_L_3[i] = new double[ntot.y];

double  **F_star_L_3 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_star_L_3[i] = new double[ntot.y];

double  **G_hll_L_4 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_hll_L_4[i] = new double[ntot.y];

double  ** F_star_L_4 = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_star_L_4[i] = new double[ntot.y];

double  **F_n_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_n_i[i] = new double[ntot.y];

double  **G_n_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_n_i[i] = new double[ntot.y];

double  **F_mom_x_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_mom_x_i[i] = new double[ntot.y];

double  **G_mom_x_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_mom_x_i[i] = new double[ntot.y];

double  **F_mom_y_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_mom_y_i[i] = new double[ntot.y];

double  **G_mom_y_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_mom_y_i[i] = new double[ntot.y];

double  **F_e_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
F_e_i[i] = new double[ntot.y];

double  **G_e_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
G_e_i[i] = new double[ntot.y];

double  **s_l_d_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_l_d_y[i] = new double[ntot.y];

double  **s_r_d_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_r_d_y[i] = new double[ntot.y];

double  **s_l_u_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_l_u_y[i] = new double[ntot.y];

double  **s_r_u_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_r_u_y[i] = new double[ntot.y];

double  **s_d_l_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_d_l_y[i] = new double[ntot.y];

double  **s_u_l_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_u_l_y[i] = new double[ntot.y];

double   **s_d_r_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_d_r_y[i] = new double[ntot.y];

double   **s_u_r_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_u_r_y[i] = new double[ntot.y];

double  **s_m_d_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_d_y[i] = new double[ntot.y];

double  **s_m_u_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_u_y[i] = new double[ntot.y];

double  **s_m_l_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_l_y[i] = new double[ntot.y];

double  **s_m_r_y = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_r_y[i] = new double[ntot.y];

double **s_l_d_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_l_d_x[i] = new double[ntot.y];

double  **s_r_d_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_r_d_x[i] = new double[ntot.y];

double  **s_l_u_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_l_u_x[i] = new double[ntot.y];

double  **s_r_u_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_r_u_x[i] = new double[ntot.y];

double **s_d_l_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_d_l_x[i] = new double[ntot.y];

double  **s_u_l_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_u_l_x[i] = new double[ntot.y];

double  **s_d_r_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_d_r_x[i] = new double[ntot.y];

double  **s_u_r_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_u_r_x[i] = new double[ntot.y];

double  **s_m_d_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_d_x[i] = new double[ntot.y];

double  **s_m_u_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_u_x[i] = new double[ntot.y];

double  **s_m_l_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_l_x[i] = new double[ntot.y];

double  **s_m_r_x = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
s_m_r_x[i] = new double[ntot.y];

for (int ii=0;ii<ntot.x;ii++)
{
    for (int jj=0;jj<ntot.y;jj++)
    {
         n_p[ii][jj]=0.0;
         n_m[ii][jj]=0.0;
         u_xp[ii][jj]=0.0;
         u_xm[ii][jj]=0.0;
         u_yp[ii][jj]=0.0;
         u_ym[ii][jj]=0.0;
         P_p[ii][jj]=0.0;
         P_m[ii][jj]=0.0;
         F_hll_D_1[ii][jj]=0.0;
         G_star_D_1[ii][jj]=0.0;
         F_hll_D_2[ii][jj]=0.0;
         G_star_D_2[ii][jj]=0.0;
         F_hll_D_3[ii][jj]=0.0;
         G_star_D_3[ii][jj]=0.0;
         F_hll_D_4[ii][jj]=0.0;
         G_star_D_4[ii][jj]=0.0;
         F_hll_U_1[ii][jj]=0.0;
         G_star_U_1[ii][jj]=0.0;
         F_hll_U_2[ii][jj]=0.0;
         G_star_U_2[ii][jj]=0.0;
         F_hll_U_3[ii][jj]=0.0;
         G_star_U_3[ii][jj]=0.0;
         F_hll_U_4[ii][jj]=0.0;
         G_star_U_4[ii][jj]=0.0;
         F_star_R_1[ii][jj]=0.0;
         G_hll_R_1[ii][jj]=0.0;
         F_star_R_2[ii][jj]=0.0;
         G_hll_R_2[ii][jj]=0.0;
         F_star_R_3[ii][jj]=0.0;
         G_hll_R_3[ii][jj]=0.0;
         F_star_R_4[ii][jj]=0.0;
         G_hll_R_4[ii][jj]=0.0;
         G_hll_L_1[ii][jj]=0.0;
         F_star_L_1[ii][jj]=0.0;
         G_hll_L_2[ii][jj]=0.0;
         F_star_L_2[ii][jj]=0.0;
         G_hll_L_3[ii][jj]=0.0;
         F_star_L_3[ii][jj]=0.0;
         G_hll_L_4[ii][jj]=0.0;
         F_star_L_4[ii][jj]=0.0;
         F_n_i[ii][jj]=0.0;
         G_n_i[ii][jj]=0.0;
         F_mom_x_i[ii][jj]=0.0;
         G_mom_x_i[ii][jj]=0.0;
         F_mom_y_i[ii][jj]=0.0;
         G_mom_y_i[ii][jj]=0.0;
         F_e_i[ii][jj]=0.0;
         G_e_i[ii][jj]=0.0;
       	 s_l_d_y[ii][jj]=0.0;
         s_r_d_y[ii][jj]=0.0;
         s_l_u_y[ii][jj]=0.0;
         s_r_u_y[ii][jj]=0.0;
         s_d_l_y[ii][jj]=0.0;
         s_u_l_y[ii][jj]=0.0;
         s_d_r_y[ii][jj]=0.0;
         s_u_r_y[ii][jj]=0.0;
         s_m_d_y[ii][jj]=0.0;
         s_m_u_y[ii][jj]=0.0;
         s_m_l_y[ii][jj]=0.0;
         s_m_r_y[ii][jj]=0.0;
         s_l_d_x[ii][jj]=0.0;
         s_r_d_x[ii][jj]=0.0;
         s_l_u_x[ii][jj]=0.0;
         s_r_u_x[ii][jj]=0.0;
         s_d_l_x[ii][jj]=0.0;
         s_u_l_x[ii][jj]=0.0;
         s_d_r_x[ii][jj]=0.0;
         s_u_r_x[ii][jj]=0.0;
         s_m_d_x[ii][jj]=0.0;
         s_m_u_x[ii][jj]=0.0;
         s_m_l_x[ii][jj]=0.0;
         s_m_r_x[ii][jj]=0.0;
	}
}
if(limiter==0){n_m = ni_tmp; n_p = ni_tmp; u_xp=ux_tmp; u_xm=ux_tmp; u_yp=uy_tmp; u_ym=uy_tmp; P_p=pi_tmp; P_m=pi_tmp;}
    else if(limiter==1)
    {
        limiter_2D(ni_tmp,ux_tmp,uy_tmp,pi_tmp,n_p,n_m,u_xp,u_xm,u_yp,u_ym,P_p,P_m,variabletype,limitertype,1,bc_x, bc_y, my_id,notWest,notEast,notSouth,notNorth,specie,flux_cath);
    }
        fluxes_y(n_p,u_xp,u_yp,P_p,n_m,u_xm,u_ym,P_m ,solver,F_star_L_1,G_hll_L_1,F_star_L_2,G_hll_L_2,F_star_L_3,G_hll_L_3,F_star_L_4,G_hll_L_4,F_star_R_1,G_hll_R_1,F_star_R_2,G_hll_R_2,F_star_R_3,G_hll_R_3,F_star_R_4,G_hll_R_4,s_l_d_y,s_r_d_y,s_l_u_y,s_r_u_y,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_m_d_y,s_m_u_y,s_m_l_y,s_m_r_y,qtmp, mtmp);

        if(limiter!=0)
        {
        limiter_2D(ni_tmp,ux_tmp,uy_tmp,pi_tmp,n_p,n_m,u_xp,u_xm,u_yp,u_ym,P_p,P_m,variabletype,limitertype,0,bc_x, bc_y,my_id,notWest,notEast,notSouth,notNorth,specie,flux_cath);
        }
	fluxes_x(n_p,u_xp,u_yp,P_p,n_m,u_xm,u_ym,P_m,solver,F_hll_D_1,G_star_D_1,F_hll_D_2,G_star_D_2,F_hll_D_3,G_star_D_3,F_hll_D_4,G_star_D_4,F_hll_U_1,G_star_U_1,F_hll_U_2,G_star_U_2,F_hll_U_3,G_star_U_3,F_hll_U_4,G_star_U_4,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_x,s_u_l_x,s_d_r_x,s_u_r_x,s_m_d_x,s_m_u_x,s_m_l_x,s_m_r_x, qtmp,mtmp);
	if(notWest==false)
        {
/*
                if(specie==0)
                {
                 for (int j=0;j<ntot.y;j++)
                        {
                        flux_cath[j] = 0.0;
                        diffn[j] = 0.0;
			}
                        flux_cathode_tot = 0.0;
                 for (int j=0;j<ntot.y;j++)
                        {
                        kf_tmp = n_m[2][j]*u_xm[2][j];
                        flux_cath[j] = kf_tmp;
                        diffn[j] = ni_tmp[ntot.x-3][j];
                        }
                        if(my_id==0)
                        cout<<"Jm="<<ech*kf_tmp<<endl;
                }
                if(specie==1)
                {
                        for (int j=0;j<ntot.y;j++)
                        {
                        kf_tmp = n_m[2][j]*u_xm[2][j];
                        flux_cath[j] = (kf_tmp - flux_cath[j]);
                        flux_cathode_tot += flux_cath[j]/((double)(Ny));
                        diffn[j] = ni_tmp[ntot.x-3][j] - diffn[j];
                        }
                        if(my_id==0)
                        cout<<"flux_anode="<<flux_cathode_tot<<endl;
                }
*/
        }
	if(notWest==false)
        {
/*
         for (int j=0;j<ntot.y;j++)
           {    
                kf1_in[j] = flux_bc(0.0,n_p[2][j],0.0,u_xp[2][j],0.0,u_yp[2][j],0.0,P_p[2][j],0,mtmp,1);
                kf2_in[j] = flux_bc(0.0,n_p[2][j],0.0,u_xp[2][j],0.0,u_yp[2][j],0.0,P_p[2][j],1,mtmp,1);
                kf3_in[j] = flux_bc(0.0,n_p[2][j],0.0,u_xp[2][j],0.0,u_yp[2][j],0.0,P_p[2][j],2,mtmp,1);
                kf4_in[j] = flux_bc(0.0,n_p[2][j],0.0,u_xp[2][j],0.0,u_yp[2][j],0.0,P_p[2][j],3,mtmp,1);
	    }
	}
	if(notEast==false)
	{
	 for (int j=0;j<ntot.y;j++)
           {
                kf1_end[j] = flux_bc(n_m[ntot.x-3][j],0.0,u_xm[ntot.x-3][j],0.0,u_ym[ntot.x-3][j],0.0,P_m[ntot.x-3][j],0.0,0,mtmp,0);
            	kf2_end[j] = flux_bc(n_m[ntot.x-3][j],0.0,u_xm[ntot.x-3][j],0.0,u_ym[ntot.x-3][j],0.0,P_m[ntot.x-3][j],0.0,1,mtmp,0);
 		kf3_end[j] = flux_bc(n_m[ntot.x-3][j],0.0,u_xm[ntot.x-3][j],0.0,u_ym[ntot.x-3][j],0.0,P_m[ntot.x-3][j],0.0,2,mtmp,0);
		kf4_end[j] = flux_bc(n_m[ntot.x-3][j],0.0,u_xm[ntot.x-3][j],0.0,u_ym[ntot.x-3][j],0.0,P_m[ntot.x-3][j],0.0,3,mtmp,0);  
	   }
*/	
	}

        if(SIR==0)
        {
        //compute the fluxes at the vertex
        time_advance_hll_flux(F_hll_D_1,G_hll_L_1,F_hll_D_2,G_hll_L_2,F_hll_D_3,G_hll_L_3,F_hll_D_4,G_hll_L_4,U_rk_tmp_1_0,U_rk_tmp_2_0,U_rk_tmp_3_0,U_rk_tmp_4_0,my_id,notWest,notEast,notSouth,notNorth,specie);
	}
	 else
        {
        fluxes_SIR(n_p,u_xp,u_yp,P_p,n_m,u_xm,u_ym,P_m,F_hll_D_1,G_star_D_1,F_hll_D_2,G_star_D_2,F_hll_D_3,G_star_D_3,F_hll_D_4,G_star_D_4,F_hll_U_1,G_star_U_1,F_hll_U_2,G_star_U_2,F_hll_U_3,G_star_U_3,F_hll_U_4,G_star_U_4,F_star_L_1,G_hll_L_1,F_star_L_2,G_hll_L_2,F_star_L_3,G_hll_L_3,F_star_L_4,G_hll_L_4,F_star_R_1,G_hll_R_1,F_star_R_2,G_hll_R_2,F_star_R_3,G_hll_R_3,F_star_R_4,G_hll_R_4, F_n_i, G_n_i, F_mom_x_i, G_mom_x_i, F_mom_y_i, G_mom_y_i, F_e_i, G_e_i,s_l_d_x,s_r_d_x,s_l_u_x,s_r_u_x,s_d_l_x,s_u_l_x,s_d_r_x,s_u_r_x,s_m_d_x,s_m_u_x,s_m_l_x,s_m_r_x,s_l_d_y,s_r_d_y,s_l_u_y,s_r_u_y,s_d_l_y,s_u_l_y,s_d_r_y,s_u_r_y,s_m_d_y,s_m_u_y,s_m_l_y,s_m_r_y,my_id,qtmp,mtmp);

        Fluxes_Boundary_conditions(F_n_i,G_n_i,F_mom_x_i,G_mom_x_i,F_mom_y_i,G_mom_y_i,F_e_i,G_e_i,my_id);
        time_advance_SIR_flux(F_hll_D_1,G_hll_L_1,F_hll_D_2,G_hll_L_2,F_hll_D_3,G_hll_L_3,F_hll_D_4,G_hll_L_4,F_n_i,G_n_i,F_mom_x_i,G_mom_x_i,F_mom_y_i,G_mom_y_i,F_e_i,G_e_i,U_rk_tmp_1_0,U_rk_tmp_2_0,U_rk_tmp_3_0,U_rk_tmp_4_0);
      }

free(kf1_end);
free(kf2_end);
free(kf3_end);
free(kf4_end);        
free(kf1_in);
free(kf2_in);
free(kf3_in);
free(kf4_in);
free(diffn);

for (int i = 0; i < ntot.x; i++)
delete[] n_p[i]   ;
delete [] n_p ;

for (int i = 0; i < ntot.x; i++)
delete[] n_m[i]  ; 
delete[] n_m ;

for (int i = 0; i < ntot.x; i++)
delete[]u_xp[i]  ;
delete[]u_xp ;

for (int i = 0; i < ntot.x; i++)
delete[] u_xm[i]  ; 
delete[] u_xm ;

for (int i = 0; i < ntot.x; i++)
delete[] u_yp[i]   ;
delete[] u_yp ;

for (int i = 0; i < ntot.x; i++)
delete[] u_ym[i]  ; 
delete[] u_ym ;

for (int i = 0; i < ntot.x; i++)
delete[] P_p[i] ;  
delete[] P_p ;

for (int i = 0; i < ntot.x; i++)
delete[] P_m[i]  ; 
delete[] P_m ;

for (int i = 0; i < ntot.x; i++)
delete[] F_hll_D_1[i]  ; 
delete[] F_hll_D_1 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_star_D_1[i] ;  
delete[] G_star_D_1 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_hll_D_2[i]  ; 
delete[] F_hll_D_2 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_star_D_2[i]   ;
delete[] G_star_D_2 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_hll_D_3[i]   ;
delete[] F_hll_D_3 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_star_D_3[i] ;  
delete[] G_star_D_3 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_hll_D_4[i]  ; 
delete[] F_hll_D_4 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_star_D_4[i]   ;
delete[] G_star_D_4 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_hll_U_1[i] ;  
delete[] F_hll_U_1 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_star_U_1[i]  ; 
delete[] G_star_U_1 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_hll_U_2[i] ;  
delete[] F_hll_U_2 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_star_U_2[i] ;  
delete[] G_star_U_2 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_hll_U_3[i]  ; 
delete[] F_hll_U_3 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_star_U_3[i]  ; 
delete[] G_star_U_3 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_hll_U_4[i] ;  
delete[] F_hll_U_4 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_star_U_4[i]   ;
delete[] G_star_U_4 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_star_R_1[i]  ; 
delete[] F_star_R_1 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_hll_R_1[i]  ; 
delete[] G_hll_R_1 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_star_R_2[i] ;  
delete[] F_star_R_2 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_hll_R_2[i]  ; 
delete[] G_hll_R_2 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_star_R_3[i] ;  
delete[] F_star_R_3 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_hll_R_3[i]  ; 
delete[] G_hll_R_3 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_star_R_4[i] ;  
delete[] F_star_R_4 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_hll_R_4[i]  ; 
delete[] G_hll_R_4 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_hll_L_1[i] ;  
delete[] G_hll_L_1 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_star_L_1[i] ;  
delete[] F_star_L_1 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_hll_L_2[i] ;  
delete[] G_hll_L_2 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_star_L_2[i] ;  
delete[] F_star_L_2 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_hll_L_3[i]   ;
delete[] G_hll_L_3 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_star_L_3[i] ; 
delete[] F_star_L_3 ;

for (int i = 0; i < ntot.x; i++)
delete[] G_hll_L_4[i] ;  
delete[] G_hll_L_4 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_star_L_4[i]  ; 
delete[] F_star_L_4 ;

for (int i = 0; i < ntot.x; i++)
delete[] F_n_i[i]  ; 
delete[] F_n_i ;

for (int i = 0; i < ntot.x; i++)
delete[] G_n_i[i]  ; 
delete[] G_n_i ;

for (int i = 0; i < ntot.x; i++)
delete[] F_mom_x_i[i]  ; 
delete[] F_mom_x_i ;

for (int i = 0; i < ntot.x; i++)
delete[] G_mom_x_i[i] ;  
delete[] G_mom_x_i ;

for (int i = 0; i < ntot.x; i++)
delete[] F_mom_y_i[i]   ;
delete[] F_mom_y_i ;

for (int i = 0; i < ntot.x; i++)
delete[] G_mom_y_i[i] ;  
delete[] G_mom_y_i ;

for (int i = 0; i < ntot.x; i++)
delete[] F_e_i[i]   ;
delete[] F_e_i ;

for (int i = 0; i < ntot.x; i++)
delete[] G_e_i[i] ;  
delete[] G_e_i ;

for (int i = 0; i < ntot.x; i++)
delete[] s_l_d_y[i]   ;
delete[] s_l_d_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_r_d_y[i] ;  
delete[] s_r_d_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_l_u_y[i]  ; 
delete[] s_l_u_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_r_u_y[i]  ; 
delete[] s_r_u_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_d_l_y[i] ;  
delete[] s_d_l_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_u_l_y[i]  ; 
delete[] s_u_l_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_d_r_y[i];  
delete[] s_d_r_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_u_r_y[i] ;  
delete[] s_u_r_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_d_y[i] ;  
delete[] s_m_d_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_u_y[i]   ;
delete[] s_m_u_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_l_y[i] ;
delete[] s_m_l_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_r_y[i] ;  
delete[] s_m_r_y ;

for (int i = 0; i < ntot.x; i++)
delete[] s_l_d_x[i] ;  
delete[] s_l_d_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_r_d_x[i]  ; 
delete[] s_r_d_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_l_u_x[i]  ; 
delete[] s_l_u_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_r_u_x[i]  ; 
delete[] s_r_u_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_d_l_x[i] ;  
delete[] s_d_l_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_u_l_x[i] ; 
delete[] s_u_l_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_d_r_x[i] ;  
delete[] s_d_r_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_u_r_x[i];  
delete[] s_u_r_x ;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_d_x[i];
delete[] s_m_d_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_u_x[i];
delete[] s_m_u_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_l_x[i]; 
delete[] s_m_l_x;

for (int i = 0; i < ntot.x; i++)
delete[] s_m_r_x[i];
delete[] s_m_r_x;  

}
void initialize(double **n_i,double **u_ix,double **u_iy,double **P_i,double **m_x_i,double **m_y_i,double **e_i)
{
    //initialization
    for (int i=0;i<ntot.x;i++)
    {
        for (int j=0;j<ntot.y;j++)
        {
            n_i[i][j]=0.0;
            u_ix[i][j]=0.0;
            u_iy[i][j]=0.0;
            m_x_i[i][j]=0.0;
            m_y_i[i][j]=0.0;
            P_i[i][j]=0.0;
            e_i[i][j]=0.0;
        }
        
    }
    
}
void primitive_variable_BC(double **n_i, double **m_x_i, double **m_y_i, double **e_i, int bc_x, int bc_y, int my_id)
{
    	if(bc_x==0)
    	{
        for (int j=0;j<ntot.y;j++)
        {
            n_i[1][j] = n_i[nbuffer][j];
            m_x_i[1][j] = m_x_i[nbuffer][j];
            m_y_i[1][j] = m_y_i[nbuffer][j];
            e_i[1][j] = e_i[nbuffer][j];
            n_i[0][j] = n_i[1][j];
            m_x_i[0][j] = m_x_i[1][j];
            m_y_i[0][j] = m_y_i[1][j];
            e_i[0][j] = e_i[1][j];
        }
    	}
        if(bc_x==0)
        {
        for (int j=0;j<ntot.y;j++)
        {   
            n_i[ntot.x-2][j] = n_i[ntot.x-3][j];
            e_i[ntot.x-2][j] = e_i[ntot.x-3][j];
            m_x_i[ntot.x-2][j] = m_x_i[ntot.x-3][j];
            m_y_i[ntot.x-2][j] = m_y_i[ntot.x-3][j];
            n_i[ntot.x-1][j] = n_i[ntot.x-3][j];
            e_i[ntot.x-1][j] = e_i[ntot.x-3][j];
            m_x_i[ntot.x-1][j] = m_x_i[ntot.x-3][j];
            m_y_i[ntot.x-1][j] = m_y_i[ntot.x-3][j];
        }
        }
          if(bc_y==0)
           {
           for (int j=0;j<ntot.x;j++)
           {
            n_i[j][1] = n_i[j][nbuffer];
            m_x_i[j][1] = m_x_i[j][nbuffer];
            m_y_i[j][1] = m_y_i[j][nbuffer];
            e_i[j][1] = e_i[j][nbuffer];
            n_i[j][0] = n_i[j][1];
            m_x_i[j][0] = m_x_i[j][1];
            m_y_i[j][0] = m_y_i[j][1];
            e_i[j][0] = e_i[j][1];
           }
          } 
          if(bc_y==0)
           {
           for (int j=0;j<ntot.x;j++)
           {
            n_i[j][ntot.y-2] = n_i[j][ntot.y-3];
            e_i[j][ntot.y-2] = e_i[j][ntot.y-3];
            m_x_i[j][ntot.y-2] = m_x_i[j][ntot.y-3];
            m_y_i[j][ntot.y-2] = m_y_i[j][ntot.y-3];
            n_i[j][ntot.y-1] = n_i[j][ntot.y-3];
            e_i[j][ntot.y-1] = e_i[j][ntot.y-3];
            m_x_i[j][ntot.y-1] = m_x_i[j][ntot.y-3];
            m_y_i[j][ntot.y-1] = m_y_i[j][ntot.y-3];
           }
          }	
       	//}

}
void InitialConditions (int my_id, int p,  GridFluid *** &gf,double mtmp, int specie,bool notWest, bool notEast,bool notSouth,bool notNorth){

   int ntot_x, ntot_vx;
   int nx, nvx;
   double P0,qp;
   double dv,vx,f,f_random,As,Ar;
   double rho1,rho2,u01,v0,u02,w0,delta_s,T_alpha;
   delta_s = 0.035;
   w0 = 0.01;
   static int flag = 0;
   int  emission_line;
   double b0 = 60.0e-4;
   double bd = 10.0e-4;
   double bmax = B_0;
   double atmp, btmp, stmp;
   double *x_nodes, *y_nodes;
   int *ind_x_nodes, *ind_y_nodes;
   int np_emiss_tmp;
    x_nodes =  (double*)malloc(ntot.x*sizeof(double));
    y_nodes =  (double*)malloc(ntot.y*sizeof(double));
    ind_x_nodes = (int*)malloc(ntot.x*sizeof(int));
    ind_y_nodes = (int*)malloc(ntot.y*sizeof(int));
    double xtmp,ytmp;
    double ind_xtmp, ind_ytmp;


   printf("ntot.x=%d, ntot.y=%d\n",ntot.x,ntot.y);

//KH initial condition 
   for(int i=0;i<ntot.x;i++)
   {
     for(int j=0;j<ntot.y;j++)
     {
        ind_xtmp = (nstart.x+(dx/2)) + i;
        ind_ytmp = (nstart.y+(dy/2)) + j;
        xtmp = xmin + dx*ind_xtmp;
        ytmp = ymin + dy*ind_ytmp;
        
        if (fabs(ytmp-0.5) < 0.25)
        {
            rho2 = 2.0; 
            P0   = 2.5;
            u02   = 0.5;
            v0   = 0.0;
        }
        else if (fabs(ytmp-0.5) >= 0.25)
        {
            rho1 = 1.0;
            P0   = 2.5;
            u01   = -0.5;
            v0   = 0.0;
         } 
         gf[i][j][0].n_i = (0.5*(rho2-rho1))*( tanh((ytmp-(1.0/4.0))/delta_s) - tanh((ytmp-(3.0/4.0))/delta_s)) + rho1 ;
         gf[i][j][0].u_ix =  (0.5*(u02-u01))*( tanh((ytmp-(1.0/4.0))/delta_s) - tanh((ytmp-(3.0/4.0))/delta_s) - 1 );
         gf[i][j][0].u_iy = w0*sin(4*M_PI*xtmp);
         gf[i][j][0].P_i = P0;
     	 gf[i][j][0].m_x_i = gf[i][j][0].n_i*gf[i][j][0].u_ix;
         gf[i][j][0].m_y_i = gf[i][j][0].n_i*gf[i][j][0].u_iy;
         gf[i][j][0].e_i = P0/(gam-1.0) + 0.5*mtmp*gf[i][j][0].n_i*((gf[i][j][0].u_ix * gf[i][j][0].u_ix)+(gf[i][j][0].u_iy*gf[i][j][0].u_iy));
	}
   }
 free(ind_x_nodes);
 free(ind_y_nodes);
 free(x_nodes);
 free(y_nodes);
}
void Copy_array(GridFluid ** gfc, double ** &n_tp,double ** &ux_tp, double ** &uy_tp,double ** &p_tp)
{
 for (int ii=0;ii<ntot.x;ii++)
    {
        for (int jj=0;jj<ntot.y;jj++)
        {
	   n_tp[ii][jj] = gfc[ii][jj].n_i;
	   ux_tp[ii][jj] = gfc[ii][jj].u_ix;
	   uy_tp[ii][jj] = gfc[ii][jj].u_iy;
       p_tp[ii][jj] = gfc[ii][jj].P_i;
	}
     }
}

void BC_x (double ** vdf,int bc_x, int bc_y)
{
   for(int itmp=nbuffer;itmp<ntot.y-nbuffer;itmp++)
   {
     if(bc_x==0)
     {
         vdf[0][itmp] = vdf[nbuffer][itmp];
         vdf[1][itmp] = vdf[nbuffer][itmp];
         vdf[ntot.x-nbuffer][itmp] = vdf[ntot.x-nbuffer-1][itmp];
         vdf[ntot.x-1][itmp] = vdf[ntot.x-nbuffer-1][itmp];
     }
      else if(bc_x==1)
      {
          vdf[0][itmp] = vdf[ntot.x-nbuffer-2][itmp];
          vdf[1][itmp] = vdf[ntot.x-nbuffer-1][itmp];
          vdf[ntot.x-nbuffer][itmp] = vdf[nbuffer][itmp];
          vdf[ntot.x-1][itmp] = vdf[nbuffer+1][itmp];
      }
     else

             cout <<"Something wrong" <<endl;
   }
}
void BC_y (double ** vdf,int bc_x, int bc_y)
{
   for(int itmp=nbuffer;itmp<ntot.x-nbuffer;itmp++)
   {
     if(bc_y==0)
     {
         vdf[itmp][0] = vdf[itmp][nbuffer];
         vdf[itmp][1] = vdf[itmp][nbuffer];
         vdf[itmp][ntot.y-nbuffer] = vdf[itmp][ntot.y-nbuffer-1];
         vdf[itmp][ntot.y-1] = vdf[itmp][ntot.y-nbuffer-1];
         
     }
      else if(bc_y==1)
      {
          vdf[itmp][0] = vdf[itmp][ntot.y-nbuffer-2];
          vdf[itmp][1] = vdf[itmp][ntot.y-nbuffer-1];
          vdf[itmp][ntot.y-nbuffer]  = vdf[itmp][nbuffer];
          vdf[itmp][ntot.y-1] = vdf[itmp][nbuffer+1];
      }
     else

             cout <<"Something wrong" <<endl;
   }
    
    //corners
    if(bc_y==1)
    {
        //top-left corner
        vdf[0][0] = vdf[ntot.x-4][0];
        vdf[1][0] = vdf[ntot.x-3][0];
        vdf[0][1] = vdf[ntot.x-4][1];
        vdf[1][1] = vdf[ntot.x-3][1];
        
        //top-right corner
        vdf[ntot.x-nbuffer][0] = vdf[2][0];
        vdf[ntot.x-1][0] = vdf[3][0];
        vdf[ntot.x-nbuffer][1] = vdf[2][1];
        vdf[ntot.x-1][1] = vdf[3][1];
        
        //bottom-left corner
        vdf[0][ntot.y-2] = vdf[ntot.x-4][ntot.y-2];
        vdf[1][ntot.y-2] = vdf[ntot.x-3][ntot.y-2];
        vdf[0][ntot.y-1] = vdf[ntot.x-4][ntot.y-1];
        vdf[1][ntot.y-1] = vdf[ntot.x-3][ntot.y-1];
        
        //bottom-right corner
        vdf[ntot.x-nbuffer][ntot.y-2] = vdf[2][ntot.y-2];
        vdf[ntot.x-1][ntot.y-2] = vdf[3][ntot.y-2];
        vdf[ntot.x-nbuffer][ntot.y-1] = vdf[2][ntot.y-1];
        vdf[ntot.x-1][ntot.y-1] = vdf[3][ntot.y-1];
    }
    if(bc_y==0)
    {
        //top-left corner
        vdf[0][0] = vdf[0][2];
        vdf[1][0] = vdf[1][3];
        vdf[0][1] = vdf[0][2];
        vdf[1][1] = vdf[1][3];
        
        //top-right corner
        vdf[ntot.x-nbuffer][0] = vdf[ntot.x-nbuffer][2];
        vdf[ntot.x-1][0] = vdf[ntot.x-1][2];
        vdf[ntot.x-nbuffer][1] = vdf[ntot.x-nbuffer][3];
        vdf[ntot.x-1][1] = vdf[ntot.x-1][3];
        
        //bottom-left corner
        vdf[0][ntot.y-2] = vdf[0][ntot.y-4];
        vdf[1][ntot.y-2] = vdf[1][ntot.y-4];
        vdf[0][ntot.y-1] = vdf[0][ntot.y-3];
        vdf[1][ntot.y-1] = vdf[1][ntot.y-3];
        
        //bottom-right corner
        vdf[ntot.x-nbuffer][ntot.y-2] = vdf[ntot.x-nbuffer][ntot.y-4];
        vdf[ntot.x-1][ntot.y-2] = vdf[ntot.x-1][ntot.y-4];
        vdf[ntot.x-nbuffer][ntot.y-1] = vdf[ntot.x-nbuffer][ntot.y-3];
        vdf[ntot.x-1][ntot.y-1] = vdf[ntot.x-1][ntot.y-3];
    }
}



void time_advance_twofluid(int my_id,GridFluid ***gfi,GridFluid ***gfe, int bc_x, int bc_y,bool notWest,bool notEast,bool notSouth,bool notNorth,int time_solver)
{

double ** n_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
n_i[i]= new double[ntot.y]; 

double ** u_ix = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
u_ix[i]= new double[ntot.y]; 

double ** u_iy = new double *[ntot.x]; 
for (int i = 0; i < ntot.x; i++)
u_iy[i]= new double[ntot.y]; 

double ** m_x_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
m_x_i[i]= new double[ntot.y]; 

double ** m_y_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
m_y_i[i]= new double[ntot.y]; 

double ** P_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
P_i[i]= new double[ntot.y]; 

double ** e_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
e_i[i]= new double[ntot.y]; 

double ** E_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
E_i[i]= new double[ntot.y]; 

double ** c_i = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
c_i[i]= new double[ntot.y]; 

double ** n_e = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
n_e[i]= new double[ntot.y]; 

double ** u_ex = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
u_ex[i]= new double[ntot.y]; 

double ** u_ey = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
u_ey[i]= new double[ntot.y]; 

double ** m_x_e = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
m_x_e[i]= new double[ntot.y]; 

double ** m_y_e = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
m_y_e[i]= new double[ntot.y]; 

double ** P_e = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
 P_e[i]= new double[ntot.y]; 

double ** e_e = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
e_e[i]= new double[ntot.y]; 

double ** E_e = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
E_e[i]= new double[ntot.y]; 

double ** c_e = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
c_e[i]= new double[ntot.y]; 

double ** rho_tmp = new double *[ntot_nodes.x];
for (int i = 0; i < ntot_nodes.x; i++)
rho_tmp[i]= new double[ntot_nodes.y]; 

double ** phi_tmp = new double *[ntot_nodes.x];
for (int i = 0; i < ntot_nodes.x; i++)
phi_tmp[i]= new double[ntot_nodes.y]; 

double ** rho = new double *[ntot.x+1];
for (int i = 0; i < ntot.x+1; i++)
rho[i]= new double[ntot.y+1]; 

double ** phi = new double *[ntot.x+1];
for (int i = 0; i < ntot.x+1; i++)
phi[i]= new double[ntot.y+1]; 

double **Ex_tmp = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Ex_tmp[i]= new double[ntot.y]; 

double **Ey_tmp = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Ey_tmp[i]= new double[ntot.y]; 

double **Br_tmp = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Br_tmp[i]= new double[ntot.y]; 

double **S_tmp = new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
S_tmp[i]= new double[ntot.y]; 

double **U_rk_1_0_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_1_0_i[i]= new double[ntot.y]; 

double **U_rk_2_0_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_2_0_i[i]= new double[ntot.y]; 

double **U_rk_3_0_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_3_0_i[i]= new double[ntot.y]; 

double **U_rk_4_0_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_4_0_i[i]= new double[ntot.y]; 

double **U_rk_1_1_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_1_1_i[i]= new double[ntot.y]; 

double **U_rk_2_1_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_2_1_i[i]= new double[ntot.y]; 

double **U_rk_3_1_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_3_1_i[i]= new double[ntot.y]; 

double **U_rk_4_1_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_4_1_i[i]= new double[ntot.y]; 

double **U_rk_1_2_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_1_2_i[i]= new double[ntot.y]; 

double **U_rk_2_2_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_2_2_i[i]= new double[ntot.y]; 

double **U_rk_3_2_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_3_2_i[i]= new double[ntot.y]; 

double **U_rk_4_2_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_4_2_i[i]= new double[ntot.y]; 

double **U_rk_1_3_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_1_3_i[i]= new double[ntot.y]; 

double **U_rk_2_3_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_2_3_i[i]= new double[ntot.y]; 

double **U_rk_3_3_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_3_3_i[i]= new double[ntot.y]; 

double **U_rk_4_3_i= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_4_3_i[i]= new double[ntot.y]; 

double **ni_tmp= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
ni_tmp[i]= new double[ntot.y]; 

double **ux_i_tmp= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
ux_i_tmp[i]= new double[ntot.y]; 

double **uy_i_tmp= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
uy_i_tmp[i]= new double[ntot.y]; 

double **pi_tmp= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
pi_tmp[i]= new double[ntot.y]; 

double **U_rk_1_0_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_1_0_e[i]= new double[ntot.y]; 

double **U_rk_2_0_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_2_0_e[i]= new double[ntot.y]; 

double **U_rk_3_0_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_3_0_e[i]= new double[ntot.y]; 

double **U_rk_4_0_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_4_0_e[i]= new double[ntot.y]; 

double **U_rk_1_1_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_1_1_e[i]= new double[ntot.y]; 

double **U_rk_2_1_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_2_1_e[i]= new double[ntot.y]; 

double **U_rk_3_1_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_3_1_e[i]= new double[ntot.y]; 

double **U_rk_4_1_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_4_1_e[i]= new double[ntot.y]; 

double **U_rk_1_2_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_1_2_e[i]= new double[ntot.y]; 

double **U_rk_2_2_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_2_2_e[i]= new double[ntot.y]; 

double **U_rk_3_2_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_3_2_e[i]= new double[ntot.y]; 

double **U_rk_4_2_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)

U_rk_4_2_e[i]= new double[ntot.y]; 
double **U_rk_1_3_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_1_3_e[i]= new double[ntot.y]; 

double **U_rk_2_3_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_2_3_e[i]= new double[ntot.y]; 

double **U_rk_3_3_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_3_3_e[i]= new double[ntot.y]; 

double **U_rk_4_3_e= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
U_rk_4_3_e[i]= new double[ntot.y]; 

double **ne_tmp= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
ne_tmp[i]= new double[ntot.y]; 

double **ux_e_tmp= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
ux_e_tmp[i]= new double[ntot.y]; 

double **uy_e_tmp= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
uy_e_tmp[i]= new double[ntot.y]; 

double **pe_tmp= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
pe_tmp[i]= new double[ntot.y]; 

double **Source1= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Source1[i]= new double[ntot.y]; 

double **Source2= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Source2[i]= new double[ntot.y]; 

double **Source3= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Source3[i]= new double[ntot.y]; 

double **Source4= new double *[ntot.x];
for (int i = 0; i < ntot.x; i++)
Source4[i]= new double[ntot.y]; 


 for (int ii=0;ii<ntot.x;ii++)
    {
        for (int jj=0;jj<ntot.y;jj++)
        {
            n_i[ii][jj]=0.0;    
            u_ix[ii][jj]=0.0;
            u_iy[ii][jj]=0.0;
            m_x_i[ii][jj]=0.0;
            m_y_i[ii][jj]=0.0;
            P_i[ii][jj]=0.0;
            e_i[ii][jj]=0.0;
            E_i[ii][jj]=0.0;
            c_i[ii][jj]=0.0;
            rho_tmp[ii][jj]=0.0;
            phi_tmp[ii][jj]=0.0;
            Ex_tmp[ii][jj]=0.0;
            Ey_tmp[ii][jj]=0.0;
            Br_tmp[ii][jj]=0.0; 
	}
    }
    for (int ii=0;ii<ntot.x;ii++)
    {
        for (int jj=0;jj<ntot.y;jj++)
        {
            n_e[ii][jj]=0.0;	
            u_ex[ii][jj]=0.0;
            u_ey[ii][jj]=0.0;
            m_x_e[ii][jj]=0.0;
            m_y_e[ii][jj]=0.0;
            P_e[ii][jj]=0.0;
            e_e[ii][jj]=0.0;
            E_e[ii][jj]=0.0;
            c_e[ii][jj]=0.0;
        }
    }
                               
 for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
	        //ions
            n_i[i][j]  = gfi[i][j][0].n_i  ;
            m_x_i[i][j]  = gfi[i][j][0].m_x_i ;
            m_y_i[i][j]  = gfi[i][j][0].m_y_i  ;
            e_i[i][j]  = gfi[i][j][0].e_i ;
            u_ix[i][j] = gfi[i][j][0].u_ix;
            u_iy[i][j] = gfi[i][j][0].u_iy;
            P_i[i][j]  = gfi[i][j][0].P_i;
	        //electrons
            n_e[i][j]  = gfe[i][j][0].n_i  ;
            m_x_e[i][j]  = gfe[i][j][0].m_x_i ;
            m_y_e[i][j]  = gfe[i][j][0].m_y_i  ;
            e_e[i][j]  = gfe[i][j][0].e_i ;
            u_ex[i][j] = gfe[i][j][0].u_ix;
            u_ey[i][j] = gfe[i][j][0].u_iy;
            P_e[i][j]  = gfe[i][j][0].P_i;
	
        }
    }

	        BC_x(n_i,bc_x,bc_y);
            BC_x(m_x_i,bc_x,bc_y);
            BC_x(m_y_i,bc_x,bc_y);
            BC_x(P_i,bc_x,bc_y);
    
            BC_y(n_i,bc_x,bc_y);
            BC_y(m_x_i,bc_x,bc_y);
            BC_y(m_y_i,bc_x,bc_y);
            BC_y(P_i,bc_x,bc_y);
 
            BC_x(n_e,bc_x,bc_y);
            BC_x(m_x_e,bc_x,bc_y);
            BC_x(m_y_e,bc_x,bc_y);
            BC_x(P_e,bc_x,bc_y);
    
            BC_y(n_e,bc_x,bc_y);
            BC_y(m_x_e,bc_x,bc_y);
            BC_y(m_y_e,bc_x,bc_y);
            BC_y(P_e,bc_x,bc_y);

	 for (int ii=0;ii<ntot.x;ii++)
	 {   
    		for (int jj=0;jj<ntot.y;jj++)
    		{    
         	U_rk_1_0_i[ii][jj]=0.0;
         	U_rk_2_0_i[ii][jj]=0.0;
         	U_rk_3_0_i[ii][jj]=0.0;
         	U_rk_4_0_i[ii][jj]=0.0;
         	U_rk_1_1_i[ii][jj]=0.0;
         	U_rk_2_1_i[ii][jj]=0.0;
         	U_rk_3_1_i[ii][jj]=0.0;
         	U_rk_4_1_i[ii][jj]=0.0;
         	U_rk_1_2_i[ii][jj]=0.0;
         	U_rk_2_2_i[ii][jj]=0.0;
         	U_rk_3_2_i[ii][jj]=0.0;
         	U_rk_4_2_i[ii][jj]=0.0;
         	U_rk_1_3_i[ii][jj]=0.0;
         	U_rk_2_3_i[ii][jj]=0.0;
         	U_rk_3_3_i[ii][jj]=0.0;
         	U_rk_4_3_i[ii][jj]=0.0;
         	ni_tmp[ii][jj]=0.0;
         	ux_i_tmp[ii][jj]=0.0;
         	uy_i_tmp[ii][jj]=0.0;
         	pi_tmp[ii][jj]=0.0;
         	U_rk_1_0_e[ii][jj]=0.0;
         	U_rk_2_0_e[ii][jj]=0.0;
         	U_rk_3_0_e[ii][jj]=0.0;
         	U_rk_4_0_e[ii][jj]=0.0;
         	U_rk_1_1_e[ii][jj]=0.0;
         	U_rk_2_1_e[ii][jj]=0.0;
         	U_rk_3_1_e[ii][jj]=0.0;
         	U_rk_4_1_e[ii][jj]=0.0;
         	U_rk_1_2_e[ii][jj]=0.0;
         	U_rk_2_2_e[ii][jj]=0.0;
         	U_rk_3_2_e[ii][jj]=0.0;
         	U_rk_4_2_e[ii][jj]=0.0;
         	U_rk_1_3_e[ii][jj]=0.0;
         	U_rk_2_3_e[ii][jj]=0.0;
         	U_rk_3_3_e[ii][jj]=0.0;
         	U_rk_4_3_e[ii][jj]=0.0;
         	ne_tmp[ii][jj]=0.0;
         	ux_e_tmp[ii][jj]=0.0;
         	uy_e_tmp[ii][jj]=0.0;
         	pe_tmp[ii][jj]=0.0;
   		}
	} 

  
 for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            U_rk_1_0_i[i][j]=n_i[i][j];
            U_rk_2_0_i[i][j]=m_x_i[i][j];
            U_rk_3_0_i[i][j]=m_y_i[i][j];
            U_rk_4_0_i[i][j]=e_i[i][j];
        }
    }
     for (int i=nbuffer;i<ntot.x-nbuffer; i++)
    {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
            if (U_rk_1_0_i[i][j]<=epsil)
            {
                ni_tmp[i][j]=epsil;
            }
            else
            {
                ni_tmp[i][j]=U_rk_1_0_i[i][j];
            }
            ux_i_tmp[i][j]=U_rk_2_0_i[i][j]/(ni_tmp[i][j]);
            uy_i_tmp[i][j]=U_rk_3_0_i[i][j]/(ni_tmp[i][j]);
            pi_tmp[i][j]=(gam-1.0)*(U_rk_4_0_i[i][j]-0.5*mass_i*ni_tmp[i][j]*(pow(ux_i_tmp[i][j],2) + pow(uy_i_tmp[i][j],2)));
            if(pi_tmp[i][j]<=0)
	    {
	    pi_tmp[i][j]=max(epsil,pi_tmp[i][j]);
	    }
 	    if(((ux_i_tmp[i][j]*dt/dx) + (uy_i_tmp[i][j]*dt/dy)) >= 1.0)
            {
            cout<<"my_id"<<my_id<<endl;    
	    cout<<"wrong CFL_i="<< ((ux_i_tmp[i][j]*dt/dx) + (uy_i_tmp[i][j]*dt/dy)) <<"at i="<<i<<",j="<<j<<endl;	
       	    }
	}
    }
    
               BC_x(ni_tmp,bc_x,bc_y);
               BC_x(ux_i_tmp,bc_x,bc_y);
               BC_x(uy_i_tmp,bc_x,bc_y);
               BC_x(pi_tmp,bc_x,bc_y);
    
               BC_y(ni_tmp,bc_x,bc_y);
               BC_y(ux_i_tmp,bc_x,bc_y);
               BC_y(uy_i_tmp,bc_x,bc_y);
               BC_y(pi_tmp,bc_x,bc_y);
  fluxes(ni_tmp,ux_i_tmp,uy_i_tmp,pi_tmp,U_rk_1_1_i,U_rk_2_1_i,U_rk_3_1_i,U_rk_4_1_i,bc_x,bc_y,my_id,ech,mass_i,notWest,notEast,notSouth,notNorth,0,flux_cathode);

     for (int i=nbuffer;i<ntot.x-nbuffer; i++)
    {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
            U_rk_1_1_i[i][j] = U_rk_1_0_i[i][j] + (dt*U_rk_1_1_i[i][j]);// + (dt*S_tmp[i][j]);

            U_rk_2_1_i[i][j] = U_rk_2_0_i[i][j] + (dt*U_rk_2_1_i[i][j]);// + ((((ech/mass_i)*dt*ni_tmp[i][j])*(Ex_tmp[i][j] + (uy_i_tmp[i][j]*Br_tmp[i][j]))));

            U_rk_3_1_i[i][j] = U_rk_3_0_i[i][j] + (dt*U_rk_3_1_i[i][j]);// + ((((ech/mass_i)*dt*ni_tmp[i][j])*(Ey_tmp[i][j] - (ux_i_tmp[i][j]*Br_tmp[i][j]))));

            U_rk_4_1_i[i][j] = U_rk_4_0_i[i][j] + (dt*U_rk_4_1_i[i][j]);// + (((ech*dt*ni_tmp[i][j])*((ux_i_tmp[i][j]*Ex_tmp[i][j]) + (uy_i_tmp[i][j]*Ey_tmp[i][j]))));
       }
    }
    //primitive_variable_BC(U_rk_1_1_i,U_rk_2_1_i,U_rk_3_1_i,U_rk_4_1_i,bc_x,bc_y,my_id);
        

	for (int i=nbuffer;i<ntot.x-nbuffer; i++)
    	{   
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        { 
            if (U_rk_1_1_i[i][j]<=epsil)
            {
            ni_tmp[i][j]=epsil;
            }
            else
            {
            ni_tmp[i][j]=U_rk_1_1_i[i][j];
            }
            ux_i_tmp[i][j]=U_rk_2_1_i[i][j]/(ni_tmp[i][j]);
            uy_i_tmp[i][j]=U_rk_3_1_i[i][j]/(ni_tmp[i][j]);
            pi_tmp[i][j]=((gam-1.0)*U_rk_4_1_i[i][j]) - (0.5*mass_i*ni_tmp[i][j]*(pow(ux_i_tmp[i][j],2) + pow(uy_i_tmp[i][j],2)));
	    if(pi_tmp[i][j]<=0)
            {
            pi_tmp[i][j]=max(epsil,pi_tmp[i][j]);
            }
	   if(((ux_i_tmp[i][j]*dt/dx) + (uy_i_tmp[i][j]*dt/dy)) >= 1.0)
            {
            cout<<"my_id"<<my_id<<endl;   
	    cout<<"wrong CFL_i="<< ((ux_i_tmp[i][j]*dt/dx) + (uy_i_tmp[i][j]*dt/dy)) <<"at x="<<i<<",y="<<j<<endl; 
	    }
       }
    }
   //........Runge-Kutta stage 1...electron...........................................
  //
   for(int j=nbuffer;j<ntot.y-nbuffer;j++)
    {
        for(int i=nbuffer;i<ntot.x-nbuffer;i++)
        {
            U_rk_1_0_e[i][j]=n_e[i][j];
            U_rk_2_0_e[i][j]=m_x_e[i][j];
            U_rk_3_0_e[i][j]=m_y_e[i][j];
            U_rk_4_0_e[i][j]=e_e[i][j];
        }
    }
            

 for (int i=nbuffer;i<ntot.x-nbuffer; i++)
        {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
	   if (U_rk_1_0_e[i][j]<=epsil)
            {
                ne_tmp[i][j]=epsil;
            }
            else
            {
                ne_tmp[i][j]=U_rk_1_0_e[i][j];
            }
            ux_e_tmp[i][j]=U_rk_2_0_e[i][j]/(ne_tmp[i][j]);
            uy_e_tmp[i][j]=U_rk_3_0_e[i][j]/(ne_tmp[i][j]);
            pe_tmp[i][j]=((gam-1.0)*U_rk_4_0_e[i][j]) - (0.5*mass_e*ne_tmp[i][j]*(pow(ux_e_tmp[i][j],2) + pow(uy_e_tmp[i][j],2)));
	    if(pe_tmp[i][j]<=0)
            {
            pe_tmp[i][j]=max(epsil,pe_tmp[i][j]);
            }
	    if(((ux_e_tmp[i][j]*dt/dx) + (uy_e_tmp[i][j]*dt/dy)) >= 1.0)
            {
                cout<<"my_id"<<my_id<<endl;    
		cout<<"wrong CFL_e="<< ((ux_e_tmp[i][j]*dt/dx) + (uy_e_tmp[i][j]*dt/dy)) <<"at x="<<i<<",y="<<j<<endl;
            }
	}
    }
    
               BC_x(ne_tmp,bc_x,bc_y);
               BC_x(ux_e_tmp,bc_x,bc_y);
               BC_x(uy_e_tmp,bc_x,bc_y);
               BC_x(pe_tmp,bc_x,bc_y);
    
               BC_y(ne_tmp,bc_x,bc_y);
               BC_y(ux_e_tmp,bc_x,bc_y);
               BC_y(uy_e_tmp,bc_x,bc_y);
               BC_y(pe_tmp,bc_x,bc_y);
    fluxes(ne_tmp,ux_e_tmp,uy_e_tmp,pe_tmp,U_rk_1_1_e,U_rk_2_1_e,U_rk_3_1_e,U_rk_4_1_e,bc_x,bc_y,my_id,ech,mass_e,notWest,notEast,notSouth,notNorth,1,flux_cathode);
     for (int i=nbuffer;i<ntot.x-nbuffer; i++)
        {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
            U_rk_1_1_e[i][j] = U_rk_1_0_e[i][j] + (dt*U_rk_1_1_e[i][j]);
            U_rk_2_1_e[i][j] = U_rk_2_0_e[i][j] + (dt*U_rk_2_1_e[i][j]);// + Source2[i][j];
            U_rk_3_1_e[i][j] = U_rk_3_0_e[i][j] + (dt*U_rk_3_1_e[i][j]) ;//+ Source3[i][j];
            U_rk_4_1_e[i][j] = U_rk_4_0_e[i][j] + (dt*U_rk_4_1_e[i][j]);// + Source4[i][j];
        }
    }
    //primitive_variable_BC(U_rk_1_1_e,U_rk_2_1_e,U_rk_3_1_e,U_rk_4_1_e,bc_x,bc_y,my_id);
            
 for (int i=nbuffer;i<ntot.x-nbuffer; i++)
        {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
	    if (U_rk_1_1_e[i][j]<=epsil)
            {
            ne_tmp[i][j]=epsil;
            }
            else
            {
            ne_tmp[i][j]=U_rk_1_1_e[i][j];
            }
            ux_e_tmp[i][j]=U_rk_2_1_e[i][j]/(ne_tmp[i][j]);
            uy_e_tmp[i][j]=U_rk_3_1_e[i][j]/(ne_tmp[i][j]);
            pe_tmp[i][j]=((gam-1.0)*U_rk_4_1_e[i][j]) - (0.5*mass_e*ne_tmp[i][j]*(pow(ux_e_tmp[i][j],2) + pow(uy_e_tmp[i][j],2)));
	    if(pe_tmp[i][j]<=0)
            {
            pe_tmp[i][j]=max(epsil,pe_tmp[i][j]);
            }
	    if(((ux_e_tmp[i][j]*dt/dx) + (uy_e_tmp[i][j]*dt/dy)) >= 1.0)
             {
                cout<<"my_id"<<my_id<<endl;   
		cout<<"wrong CFL_e="<< ((ux_e_tmp[i][j]*dt/dx) + (uy_e_tmp[i][j]*dt/dy)) <<"at x="<<i<<",y="<<j<<endl;
       	     }
 	}
    }
 /////////////////////////////////RK2/////////////////////////////////////
               BC_x(ni_tmp,bc_x,bc_y);
               BC_x(ux_i_tmp,bc_x,bc_y);
               BC_x(uy_i_tmp,bc_x,bc_y);
               BC_x(pi_tmp,bc_x,bc_y);
    
               BC_y(ni_tmp,bc_x,bc_y);
               BC_y(ux_i_tmp,bc_x,bc_y);
               BC_y(uy_i_tmp,bc_x,bc_y);
               BC_y(pi_tmp,bc_x,bc_y);
    
    fluxes(ni_tmp,ux_i_tmp,uy_i_tmp,pi_tmp,U_rk_1_2_i,U_rk_2_2_i,U_rk_3_2_i,U_rk_4_2_i,bc_x,bc_y,my_id,ech,mass_i,notWest,notEast,notSouth,notNorth,0,flux_cathode);

     for (int i=nbuffer;i<ntot.x-nbuffer; i++)
    {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
            U_rk_1_2_i[i][j] = ((3.0/4.0)*U_rk_1_0_i[i][j]) + ((1.0/4.0)*U_rk_1_1_i[i][j]) + ((dt/4.0)*U_rk_1_2_i[i][j]);// + ((dt/4.0)*S_tmp[i][j]);

            U_rk_2_2_i[i][j] = ((3.0/4.0)*U_rk_2_0_i[i][j]) + ((1.0/4.0)*U_rk_2_1_i[i][j]) + ((dt/4.0)*U_rk_2_2_i[i][j]);//  + (((ech/mass_i)*(dt/4.0)*ni_tmp[i][j])*(Ex_tmp[i][j]  + (uy_i_tmp[i][j]*Br_tmp[i][j])));

            U_rk_3_2_i[i][j] = ((3.0/4.0)*U_rk_3_0_i[i][j]) + ((1.0/4.0)*U_rk_3_1_i[i][j]) + ((dt/4.0)*U_rk_3_2_i[i][j]);// + (((ech/mass_i)*(dt/4.0)*ni_tmp[i][j])*(Ey_tmp[i][j]  - (ux_i_tmp[i][j]*Br_tmp[i][j])));

            U_rk_4_2_i[i][j] = ((3.0/4.0)*U_rk_4_0_i[i][j]) + ((1.0/4.0)*U_rk_4_1_i[i][j]) + ((dt/4.0)*U_rk_4_2_i[i][j]);// + ((dt/4.0)*(ech*ni_tmp[i][j])*((ux_i_tmp[i][j]*Ex_tmp[i][j]) + (uy_i_tmp[i][j]*Ey_tmp[i][j])));
        }
    }
    //primitive_variable_BC(U_rk_1_2_i,U_rk_2_2_i,U_rk_3_2_i,U_rk_4_2_i,bc_x,bc_y,my_id);
        
 for (int i=nbuffer;i<ntot.x-nbuffer; i++)
        {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
	    if (U_rk_1_2_i[i][j]<=epsil)
            {
                ni_tmp[i][j]=epsil;
            }
            else
            {
                ni_tmp[i][j]=U_rk_1_2_i[i][j];
            }
            ux_i_tmp[i][j]=U_rk_2_2_i[i][j]/(ni_tmp[i][j]);
            uy_i_tmp[i][j]=U_rk_3_2_i[i][j]/(ni_tmp[i][j]);
            pi_tmp[i][j]=((gam-1.0)*U_rk_4_2_i[i][j]) - (0.5*mass_i*ni_tmp[i][j]*(pow(ux_i_tmp[i][j],2) + pow(uy_i_tmp[i][j],2)));
	    if(pi_tmp[i][j]<=0)
            {
            pi_tmp[i][j]=max(epsil,pi_tmp[i][j]);
            }
	    if(((ux_i_tmp[i][j]*dt/dx) + (uy_i_tmp[i][j]*dt/dy)) >= 1.0)
	    {
		cout<<"my_id"<<my_id<<endl;	
                cout<<"wrong CFL_i="<< ((ux_i_tmp[i][j]*dt/dx) + (uy_i_tmp[i][j]*dt/dy)) <<"at x="<<i<<",y="<<j<<endl;
            }
	}
    }
//........Runge-Kutta stage 2....electrons..........................................
               BC_x(ne_tmp,bc_x,bc_y);
               BC_x(ux_e_tmp,bc_x,bc_y);
               BC_x(uy_e_tmp,bc_x,bc_y);
               BC_x(pe_tmp,bc_x,bc_y);
    
               BC_y(ne_tmp,bc_x,bc_y);
               BC_y(ux_e_tmp,bc_x,bc_y);
               BC_y(uy_e_tmp,bc_x,bc_y);
               BC_y(pe_tmp,bc_x,bc_y);
    fluxes(ne_tmp,ux_e_tmp,uy_e_tmp,pe_tmp,U_rk_1_2_e,U_rk_2_2_e,U_rk_3_2_e,U_rk_4_2_e,bc_x,bc_y,my_id,ech,mass_e,notWest,notEast,notSouth,notNorth,1,flux_cathode);
     

    for (int i=nbuffer;i<ntot.x-nbuffer; i++)
    {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
            U_rk_1_2_e[i][j] = ((3.0/4.0)*U_rk_1_0_e[i][j]) + ((1.0/4.0)*U_rk_1_1_e[i][j]) + ((dt/4.0)*U_rk_1_2_e[i][j]);
            U_rk_2_2_e[i][j] = ((3.0/4.0)*U_rk_2_0_e[i][j]) + ((1.0/4.0)*U_rk_2_1_e[i][j]) + ((dt/4.0)*U_rk_2_2_e[i][j]);
            U_rk_3_2_e[i][j] = ((3.0/4.0)*U_rk_3_0_e[i][j]) + ((1.0/4.0)*U_rk_3_1_e[i][j]) + ((dt/4.0)*U_rk_3_2_e[i][j]);
            U_rk_4_2_e[i][j] = ((3.0/4.0)*U_rk_4_0_e[i][j]) + ((1.0/4.0)*U_rk_4_1_e[i][j]) + ((dt/4.0)*U_rk_4_2_e[i][j]);
        }
    }
   // primitive_variable_BC(U_rk_1_2_e,U_rk_2_2_e,U_rk_3_2_e,U_rk_4_2_e,bc_x,bc_y,my_id);
        

 for (int i=nbuffer;i<ntot.x-nbuffer; i++)
        {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
	  if (U_rk_1_2_e[i][j]<=epsil)
            {
                ne_tmp[i][j]=epsil;
            }
            else
            {
                ne_tmp[i][j]=U_rk_1_2_e[i][j];
            }
            ux_e_tmp[i][j]=U_rk_2_2_e[i][j]/(ne_tmp[i][j]);
            uy_e_tmp[i][j]=U_rk_3_2_e[i][j]/(ne_tmp[i][j]);
            pe_tmp[i][j]=((gam-1.0)*U_rk_4_2_e[i][j]) - (0.5*mass_e*ne_tmp[i][j]*(pow(ux_e_tmp[i][j],2) + pow(uy_e_tmp[i][j],2)));
	    if(pe_tmp[i][j]<=0)
            {
            pe_tmp[i][j]=max(epsil,pe_tmp[i][j]);
	    }
	    if(((ux_e_tmp[i][j]*dt/dx) + (uy_e_tmp[i][j]*dt/dy)) >= 1.0)
            {
            cout<<"my_id"<<my_id<<endl;    
	    cout<<"wrong CFL_e="<< ((ux_e_tmp[i][j]*dt/dx) + (uy_e_tmp[i][j]*dt/dy)) <<"at x="<<i<<",y="<<j<<endl;
            }
	}
    }
//////////////////////////////////////RK3/////////////////////////////////////
               BC_x(ni_tmp,bc_x,bc_y);
               BC_x(ux_i_tmp,bc_x,bc_y);
               BC_x(uy_i_tmp,bc_x,bc_y);
               BC_x(pi_tmp,bc_x,bc_y);
    
               BC_y(ni_tmp,bc_x,bc_y);
               BC_y(ux_i_tmp,bc_x,bc_y);
               BC_y(uy_i_tmp,bc_x,bc_y);
               BC_y(pi_tmp,bc_x,bc_y);
    fluxes(ni_tmp,ux_i_tmp,uy_i_tmp,pi_tmp,U_rk_1_3_i,U_rk_2_3_i,U_rk_3_3_i,U_rk_4_3_i,bc_x,bc_y,my_id,ech,mass_i,notWest,notEast,notSouth,notNorth,0,flux_cathode);
    for (int i=nbuffer;i<ntot.x-nbuffer; i++)
    {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {

            n_i[i][j] =  ((1.0/3.0)*U_rk_1_0_i[i][j]) + ((2.0/3.0)*U_rk_1_2_i[i][j]) + ((2.0*dt/3.0)*U_rk_1_3_i[i][j]);
            m_x_i[i][j] = ((1.0/3.0)*U_rk_2_0_i[i][j]) + ((2.0/3.0)*U_rk_2_2_i[i][j]) + ((2.0*dt/3.0)*U_rk_2_3_i[i][j]);
            m_y_i[i][j] = ((1.0/3.0)*U_rk_3_0_i[i][j]) + ((2.0/3.0)*U_rk_3_2_i[i][j])  + ((2.0*dt/3.0)*U_rk_3_3_i[i][j]);
            e_i[i][j] = ((1.0/3.0)*U_rk_4_0_i[i][j]) + ((2.0/3.0)*U_rk_4_2_i[i][j]) + ((2.0*dt/3.0)*U_rk_4_3_i[i][j]);
    	}
    //primitive_variable_BC(n_i,m_x_i,m_y_i,e_i,bc_x,bc_y,my_id);
    }        
 for (int i=nbuffer;i<ntot.x-nbuffer; i++)
        {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
		if (n_i[i][j]<=epsil)
                {
                n_i[i][j]=epsil;
                }
                else
                {
                n_i[i][j]=n_i[i][j];
                }
                u_ix[i][j]=m_x_i[i][j]/(n_i[i][j]);
                u_iy[i][j]=m_y_i[i][j]/(n_i[i][j]);
                P_i[i][j]=((gam-1.0)*e_i[i][j]) - (0.5*mass_i*n_i[i][j]*(pow(u_ix[i][j],2) + pow(u_iy[i][j],2)));
		if(P_i[i][j]<=0)
            	{
            	P_i[i][j]=max(epsil,P_i[i][j]);
            	}
		if(((u_ix[i][j]*dt/dx) + (u_iy[i][j]*dt/dy)) >= 1.0)
                {
                cout<<"my_id"<<my_id<<endl;
		cout<<"wrong CFL_i="<< ((u_ix[i][j]*dt/dx) + (u_iy[i][j]*dt/dy)) <<"at x="<<i<<",y="<<j<<endl;
	        }
	}
        } 
//........Runge-Kutta stage 3...electrons...........................................
    
    BC_x(ne_tmp,bc_x,bc_y);
    BC_x(ux_e_tmp,bc_x,bc_y);
    BC_x(uy_e_tmp,bc_x,bc_y);
    BC_x(pe_tmp,bc_x,bc_y);
    
    BC_y(ne_tmp,bc_x,bc_y);
    BC_y(ux_e_tmp,bc_x,bc_y);
    BC_y(uy_e_tmp,bc_x,bc_y);
    BC_y(pe_tmp,bc_x,bc_y);
    fluxes(ne_tmp,ux_e_tmp,uy_e_tmp,pe_tmp,U_rk_1_3_e,U_rk_2_3_e,U_rk_3_3_e,U_rk_4_3_e,bc_x,bc_y,my_id,ech,mass_e,notWest,notEast,notSouth,notNorth,1,flux_cathode);
        for (int i=nbuffer;i<ntot.x-nbuffer; i++)
        {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
            n_e[i][j] =  ((1.0/3.0)*U_rk_1_0_e[i][j]) + ((2.0/3.0)*U_rk_1_2_e[i][j]) + ((2.0*dt/3.0)*U_rk_1_3_e[i][j]);//  + Source1[i][j] ;

            m_x_e[i][j] = ((1.0/3.0)*U_rk_2_0_e[i][j]) + ((2.0/3.0)*U_rk_2_2_e[i][j]) + ((2.0*dt/3.0)*U_rk_2_3_e[i][j]);// + Source2[i][j];

            m_y_e[i][j] = ((1.0/3.0)*U_rk_3_0_e[i][j]) + ((2.0/3.0)*U_rk_3_2_e[i][j]) + ((2.0*dt/3.0)*U_rk_3_3_e[i][j]);// + Source3[i][j];

            e_e[i][j] = ((1.0/3.0)*U_rk_4_0_e[i][j]) + ((2.0/3.0)*U_rk_4_2_e[i][j]) + ((2.0*dt/3.0)*U_rk_4_3_e[i][j]);// + Source4[i][j];
        }
    }
    //primitive_variable_BC(n_e,m_x_e,m_y_e,e_e,bc_x,bc_y,my_id);
	 
	for (int i=nbuffer;i<ntot.x-nbuffer; i++)
        {
        for (int j=nbuffer;j<ntot.y-nbuffer;j++)
        {
                if (n_e[i][j]<=epsil)
                {
                n_e[i][j]=epsil;
                }
                else
                {
                n_e[i][j]=n_e[i][j];
                }
                u_ex[i][j]=m_x_e[i][j]/(n_e[i][j]);
                u_ey[i][j]=m_y_e[i][j]/(n_e[i][j]);
                P_e[i][j]=((gam-1.0)*e_e[i][j]) - (0.5*mass_e*n_e[i][j]*(pow(u_ex[i][j],2) + pow(u_ey[i][j],2)));
                if(P_e[i][j]<=0)
            	{
            	P_e[i][j]=max(epsil,P_e[i][j]);
		}
	        if(((u_ex[i][j]*dt/dx) + (u_ey[i][j]*dt/dy)) >= 0.3)
                {
                cout<<"my_id"<<my_id<<endl;
		cout<<"wrong CFL_e="<< ((u_ex[i][j]*dt/dx) + (u_ey[i][j]*dt/dy)) <<"at x="<<i<<",y="<<j<<endl;
		}
	 }
         }
   
      for (int ii=nbuffer;ii<ntot.x-nbuffer;ii++)
      {
        for (int jj=nbuffer;jj<ntot.y-nbuffer;jj++)
        {
	  gfi[ii][jj][0].n_i = n_i[ii][jj];
	  gfi[ii][jj][0].u_ix = u_ix[ii][jj];
	  gfi[ii][jj][0].u_iy = u_iy[ii][jj];
          gfi[ii][jj][0].P_i = P_i[ii][jj];
	  gfi[ii][jj][0].m_x_i = m_x_i[ii][jj];
	  gfi[ii][jj][0].m_y_i= m_y_i[ii][jj];
          gfi[ii][jj][0].e_i = e_i[ii][jj];
	}
      }
      for (int ii=nbuffer;ii<ntot.x-nbuffer;ii++)
      {
        for (int jj=nbuffer;jj<ntot.y-nbuffer;jj++)
        {
	  gfe[ii][jj][0].n_i = n_e[ii][jj];
	  gfe[ii][jj][0].u_ix = u_ex[ii][jj];
	  gfe[ii][jj][0].u_iy = u_ey[ii][jj];
          gfe[ii][jj][0].P_i = P_e[ii][jj];
	  gfe[ii][jj][0].m_x_i = m_x_e[ii][jj];
	  gfe[ii][jj][0].m_y_i= m_y_e[ii][jj];
          gfe[ii][jj][0].e_i = e_e[ii][jj];
	}
      }

for (int i = 0; i < ntot.x; i++)
   delete[] n_i[i] ;
delete[] n_i;
   
for (int i = 0; i < ntot.x; i++)
   delete[] u_ix[i] ;    
delete[] u_ix;

for (int i = 0; i < ntot.x; i++)   
   delete[] u_iy[i] ;
delete[] u_iy;
 
for (int i = 0; i < ntot.x; i++)  
   delete[] m_x_i[i] ;
delete[] m_x_i;

for (int i = 0; i < ntot.x; i++)
   delete[] m_y_i[i] ;
delete[] m_y_i;

for (int i = 0; i < ntot.x; i++)
   delete[] P_i[i] ;
delete[] P_i;

for (int i = 0; i < ntot.x; i++)
   delete[] e_i[i] ;
delete[] e_i;

for (int i = 0; i < ntot.x; i++)
   delete[] E_i[i] ;
delete[] E_i;

for (int i = 0; i < ntot.x; i++)
   delete[] c_i[i] ;
delete[] c_i;

for (int i = 0; i < ntot.x; i++)
   delete[] n_e[i] ;
delete[] n_e;

for (int i = 0; i < ntot.x; i++)
   delete[] u_ex[i] ;
 delete[] u_ex;

for (int i = 0; i < ntot.x; i++)
   delete[] u_ey[i] ;
 delete[] u_ey;

for (int i = 0; i < ntot.x; i++)
   delete[] m_x_e[i] ;
delete[] m_x_e;

for (int i = 0; i < ntot.x; i++)
   delete[] m_y_e[i] ;
 delete[] m_y_e;

for (int i = 0; i < ntot.x; i++)
   delete[] P_e[i] ;
delete[] P_e;

for (int i = 0; i < ntot.x; i++)
   delete[] e_e[i] ;
delete[] e_e;

for (int i = 0; i < ntot.x; i++)
   delete[] E_e[i] ;
 delete[] E_e;

for (int i = 0; i < ntot.x; i++)
   delete[] c_e[i] ; 
delete[] c_e; 

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_1_0_i[i] ;
delete[] U_rk_1_0_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_2_0_i[i] ;
delete[] U_rk_2_0_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_3_0_i[i] ;
delete[] U_rk_3_0_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_4_0_i[i] ;
delete[] U_rk_4_0_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_1_1_i[i] ;
delete[] U_rk_1_1_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_2_1_i[i] ;
delete[] U_rk_2_1_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_3_1_i[i] ;
delete[] U_rk_3_1_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_4_1_i[i] ;
delete[] U_rk_4_1_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_1_2_i[i] ;
delete[] U_rk_1_2_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_2_2_i[i] ;
delete[] U_rk_2_2_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_3_2_i[i] ;
delete[] U_rk_3_2_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_4_2_i[i] ;
delete[] U_rk_4_2_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_1_3_i[i] ;
delete[] U_rk_1_3_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_2_3_i[i] ;
delete[] U_rk_2_3_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_3_3_i[i] ;
delete[] U_rk_3_3_i;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_4_3_i[i] ;
delete[] U_rk_4_3_i;

for (int i = 0; i < ntot.x; i++)
   delete[] ni_tmp[i] ;
delete[] ni_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] ux_i_tmp[i] ;
delete[] ux_i_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] uy_i_tmp[i] ;
delete[] uy_i_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] pi_tmp[i] ;
delete[] pi_tmp;


for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_1_0_e[i] ;
delete[] U_rk_1_0_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_2_0_e[i] ;
delete[] U_rk_2_0_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_3_0_e[i] ;
delete[] U_rk_3_0_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_4_0_e[i] ;
 delete[] U_rk_4_0_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_1_1_e[i] ;
delete[] U_rk_1_1_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_2_1_e[i] ;
delete[] U_rk_2_1_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_3_1_e[i] ;
delete[] U_rk_3_1_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_4_1_e[i] ;
delete[] U_rk_4_1_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_1_2_e[i] ;
delete[] U_rk_1_2_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_2_2_e[i] ;
delete[] U_rk_2_2_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_3_2_e[i] ;
delete[] U_rk_3_2_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_4_2_e[i] ;
delete[] U_rk_4_2_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_1_3_e[i] ;
delete[] U_rk_1_3_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_2_3_e[i] ;
 delete[] U_rk_2_3_e;


for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_3_3_e[i] ;
delete[] U_rk_3_3_e;

for (int i = 0; i < ntot.x; i++)
   delete[] U_rk_4_3_e[i] ;
delete[] U_rk_4_3_e;

for (int i = 0; i < ntot.x; i++)
   delete[] ne_tmp[i] ;
delete[] ne_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] ux_e_tmp[i] ;
delete[] ux_e_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] uy_e_tmp[i] ;
 delete[] uy_e_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] pe_tmp[i] ;
 delete[] pe_tmp;

for (int i = 0; i < ntot_nodes.x; i++)
   delete[] rho_tmp[i] ;
delete[] rho_tmp;

for (int i = 0; i < ntot_nodes.x; i++)
   delete[] phi_tmp[i] ;
 delete[] phi_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] Ex_tmp[i] ;
delete[] Ex_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] Ey_tmp[i] ;
delete[] Ey_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] Br_tmp[i] ;
delete[] Br_tmp;

for (int i = 0; i < ntot.x+1; i++)
   delete[] phi[i] ;
delete[] phi;

for (int i = 0; i < ntot.x+1; i++)
   delete[] rho[i] ;
delete[] rho;

for (int i = 0; i < ntot.x; i++)
   delete[] S_tmp[i] ; 
delete[] S_tmp;

for (int i = 0; i < ntot.x; i++)
   delete[] Source1[i] ;
 delete[] Source1;

for (int i = 0; i < ntot.x; i++)
   delete[] Source2[i] ;
delete[] Source2;

for (int i = 0; i < ntot.x; i++)
   delete[] Source3[i] ;
 delete[] Source3;

for (int i = 0; i < ntot.x; i++)
   delete[] Source4[i] ;
delete[] Source4;
}
