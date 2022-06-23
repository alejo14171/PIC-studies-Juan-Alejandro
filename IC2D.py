import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import math

import clases

def BorisA(Ex,Ey,Ez,Bx,By,Bz,velx,vely,velz,q,m,dt):
  #print(Ex,Ey,Ez,Bx,By,Bz,velx,vely,velz)

  #Hallar la velocidad final 

  #paso 1 ecuacion (3) del documento (Inicicializar velocidades)
  ux0 = velx
  uy0 = vely
  uz0 = velz

  #Velocidad un más adelante "(Ex*q/(2.0*m))""
  u_menosx= ux0 + (Ex*q/(2.0*m))
  u_menosy= uy0 + (Ey*q/(2.0*m))
  u_menosz= uz0 + (Ez*q/(2.0*m))

  #Ángulo de fase de la rotación 

  #Paso 2 ecuacion (6) del documento
  tx=q*dt*Bx/(2.0*m)
  ty=q*dt*By/(2.0*m)
  tz=q*dt*Bz/(2.0*m)

  #paso 3 ecuacion (8) del documento. 
  u_primax=u_menosx + ((u_menosy*tz)-(u_menosz*ty))
  u_primay=u_menosy + ((u_menosz*tx)-(u_menosx*tz))
  u_primaz=u_menosz + ((u_menosx*ty)-(u_menosy*tx))

  # Agrupando terminos:
  sx=2*tx/(1.0+(tx*tx))
  sy=2*ty/(1.0+(ty*ty))
  sz=2*tz/(1.0+(tz*tz))
  # Paso 4 ecuacion(9)
  u_masx = u_menosx + ((u_primay*sz) - (u_primaz*sy))
  u_masy = u_menosy + ((u_primaz*sx) - (u_primax*sz))
  u_masz = u_menosz + ((u_primax*sy) - (u_primay*sx))
  #Paso 5
  u_finalx= u_masx + (Ex*q/(2*m))
  u_finaly= u_masy + (Ey*q/(2*m))
  u_finalz= u_masz + (Ez*q/(2*m))
  
  return u_finalx,u_finaly,u_finalz

#Método de leapfrog para hallar las componentes de la siguiente posición

def leapfrog(posx,posy,posz,ux,uy,uz,dt):
  pos_finalx=posx+ux*dt
  pos_finaly=posy+uy*dt
  pos_finalz=posz+uz*dt
  return pos_finalx,pos_finaly,pos_finalz
  

def interpolation(posx,posy,posz):
    #x-direction
    ix =  int((posx-xmin)/dx) #index x
    xleft = posx    #position of the left node
    xright = xleft+dx  #position of the right node
    ileft  = ix    #index of the left node
    iright = ix+1  # #index of the right node
    hxleft  = (posx-xx[ix])/dx  #particle fractional x-distance from the nearest node
    #y-direction
    jy = int((posy-ymin)/dy) #index y
    ydown = posy #position of the down node
    yup = ydown+dy #position of the up node
    jdown  = jy  #index of the down node
    jup = jy+1   #index of up node
    hydown  = (posy-yy[jy])/dy; #particle fractional y-distance from the nearest node
    #print(ileft,iright,hxleft,jdown,jup,hydown)
    return ileft,iright,hxleft,jdown,jup,hydown
  
def get_fields_nodes(ex,ey,ez,bx,by,bz,ileft,iright,hxleft,jdown,jup,hydown):
    #weight functions
    w1 = (1.0-hxleft)*(1.0-hydown)
    w2 = hxleft*(1.0-hydown)
    w3 = (1.0-hxleft)*hydown
    w4 = hxleft*hydown
    # linear interpolations of fields on the particle position
    Ex = w1*ex[ileft,jdown] +  w2*ex[iright,jdown] + w3*ex[ileft,jup] + w4*ex[iright,jup]
    Ey = w1*ey[ileft,jdown] +  w2*ey[iright,jdown] + w3*ey[ileft,jup] + w4*ey[iright,jup]
    Ez = w1*ez[ileft,jdown] +  w2*ez[iright,jdown] + w3*ez[ileft,jup] + w4*ez[iright,jup]
    Bx = w1*bx[ileft,jdown] +  w2*bx[iright,jdown] + w3*bx[ileft,jup] + w4*bx[iright,jup]
    By = w1*by[ileft,jdown] +  w2*by[iright,jdown] + w3*by[ileft,jup] + w4*by[iright,jup]
    Bz = w1*bz[ileft,jdown] +  w2*bz[iright,jdown] + w3*bz[ileft,jup] + w4*bz[iright,jup]
    return Ex,Ey,Ez,Bx,By,Bz
    
def zigzag(xi,yi,zi,xf,yf,zf,dx,dy,dt):

    #xgc lo dejamos en standby

    ii=np.zeros(2)
    jj=np.zeros(2)
    Fx=np.zeros(2)
    Fy=np.zeros(2)
    Wx=np.zeros(2)
    Wy=np.zeros(2)
    J1=np.zeros([2,2])
    #JY1=np.zeros([2,2])
    J2=np.zeros([2,2])
    #JY2=np.zeros([2,2])
    


    ii[0]=math.floor(xi,dx)
    jj[0]=math.floor(yi,dy)


    ii[1]=math.floor(xf/dx)  
    jj[1]=math.floor(yf/dy)
        
        
    val1x=min(ii[0]*dx,ii[1]*dx)+dx
    val2x=max(max(ii[0]*dx,ii[1]*dx),(xi+xf)/2.0)
    xr=min(val1x,val2x)

    val1y=min(jj[0]*dy,jj[1]*dy)+dy
    val2y=max(max(jj[0]*dy,jj[1]*dy),(yi+yf)/2.0)
    yr=min(val1y,val2y)

                 
    Fx[0]=(xr-xi)/dt   
    Fy[0]=(yr-yi)/dt    

    Fx[1]=(xf-xr)/dt  
    Fy[1]=(yf-yr)/dt   
        

    Wx[0]=((xi+xr)/(2*dx))-ii[0]
    Wy[0]=((yi+yr)/(2*dy))-jj[0]
    

    Wx[1]=((xr+xf)/(2*dx))-ii[1]
    Wy[1]=((yr+yf)/(2*dy))-jj[1]

    J1[1,0]=(1/(dx*dy))*Fx[0]*(1-Wy[0])
    J1[1,2]=(1/(dx*dy))*Fx[0]*Wy[0]
    J1[0,1]=(1/(dx*dy))*Fy[0]*(1-Wx[0])
    J1[2,1]=(1/(dx*dy))*Fy[0]*Wx[0]
    
    J2[1,0]=(1/(dx*dy))*Fx[1]*(1-Wy[1])
    J2[1,2]=(1/(dx*dy))*Fx[1]*Wy[1]
    J2[0,1]=(1/(dx*dy))*Fy[1]*(1-Wx[1])
    J2[2,1]=(1/(dx*dy))*Fy[1]*Wx[1]

def densidad_corriente_zigzag(xgc, pos_antigua, pos_nueva):
    #num_save = 3
    # Atemp_local
    # Atemp_total
    # ntotal = Nx*Ny*Nz*num_save
    # Atemp_local = ???
    # Atemp_total = ???
    dV = dx*dy

    #Atemp_local = np.zeros(ntotal)
    #Atemp_total = np.zeros(ntotal)

    xgc.jx = np.zeros([Nx, Ny])
    xgc.jy = np.zeros([Nx, Ny])

    
    


    


#definicion de malla computacional
xmin = -1.0*np.pi
xmax = 1.0*np.pi
ymin = -1.0*np.pi
ymax = 1.0*np.pi
Nx = 64
Ny = 64
Lx =  xmax-xmin
Ly =  ymax-ymin
dx = Lx/Nx
dy = Ly/Ny
dt = np.pi/10
tfinal = 1000.0*np.pi
npart = 1  #number of particles
nout =  int(tfinal/dt)#numero de iteraciones+
#particles = [npart,{x,y,z,vx,vy,vz},nout]
particulas = np.zeros([npart,6,nout])
energia_cinetica = np.zeros(nout)
epsilon_r=np.zeros(nout)
q=1.0
m=1.0

#crear malla y campos
xx = np.zeros(xmin,xmax,Nx+1)
yy = np.zeros(ymin,ymax,Ny+1)
x,y = np.meshgrid(xx,yy)



#distribucion inicial de particulas (t=0) para 1 particula
for i in range(npart):
  particulas[i,0,0] = 2.0 #posicion x[numero de particulas,posx,posy,posz,vx,vy,vz]
  particulas[i,1,0] = 0.0 #posicion y
  particulas[i,2,0] = 0.0 #posicion z
  particulas[i,3,0] = 0.75 #posicion vx
  particulas[i,4,0] = 0.0 #posicion vy
  particulas[i,5,0] = 0.0 #posicion vz


#campo electromagnetico (B & E on the nodes)
bx= np.zeros([Nx+1,Ny+1])
by= np.zeros([Nx+1,Ny+1])
bz= np.zeros([Nx+1,Ny+1])
ex= np.zeros([Nx+1,Ny+1])
ey= np.zeros([Nx+1,Ny+1])
ez= np.zeros([Nx+1,Ny+1])

for iy in range(0,Ny+1):
    for ix in range(0,Nx+1):
        tmp_far = x[ix,iy]**2 + y[ix,iy]**2
        #print(tmp_far)
        if(tmp_far==0):
            bz[ix,iy] = 0.0
            ex[ix,iy] = 0.0
            ey[ix,iy] = 0.0
        else:
            bz[ix,iy] = np.sqrt(tmp_far)
            ex[ix,iy] = (0.01*x[ix,iy])/pow(tmp_far,1.5)
            ey[ix,iy] = (0.01*y[ix,iy])/pow(tmp_far,1.5)





#version (2) leap-frog: posicion a tiempos enteros (n) y velocidad a tiempo intermedio (n+1/2)
Bx=0;By=0;Bz=0;Ex=0.0;Ey=0.0;Ez=0.0
for i in range(nout):
    #cada partícula
    for j in range(npart):
        #step 1-a v^n ----- v^&{n+1/2} (solo en t=0)
        if(i==0):
            x0 = particulas[j,0,i]
            y0 = particulas[j,1,i]
            z0 = particulas[j,2,i]
            ux0 = particulas[j,3,i]
            uy0 = particulas[j,4,i]
            uz0=particulas[j,5,i]
            #found the closest node indexes according to particle position
            ileft,iright,hxleft,jdown,jup,hydown = interpolation(x0,y0,z0)
            #linear interpolation of fields (nodes) to particle position
            Ex,Ey,Ez,Bx,By,Bz  = get_fields_nodes(ex,ey,ez,bx,by,bz,ileft,iright,hxleft,jdown,jup,hydown)
            #Advance particle velocity from n to n+1/2
            [u_finalx,u_finaly,u_finalz] = BorisA(Ex,Ey,Ez,Bx,By,Bz,ux0,uy0,uz0,q,m,0.5*dt)
        else:
            #step 1-b v^n+1/2 ----- v^&{n+3/2} (when t!=0)
            x0 = particulas[j,0,i-1]
            y0 = particulas[j,1,i-1]
            z0 = particulas[j,2,i-1]
            ux0 = particulas[j,3,i-1]
            uy0 = particulas[j,4,i-1]
            uz0=particulas[j,5,i-1]
            #found the closest node indexes according to particle position
            ileft,iright,hxleft,jdown,jup,hydown = interpolation(x0,y0,z0)
            #linear interpolation of fields (nodes) to particle position
            Ex,Ey,Ez,Bx,By,Bz  = get_fields_nodes(ex,ey,ez,bx,by,bz,ileft,iright,hxleft,jdown,jup,hydown)
            #print(Ex,Ey,Ez,Bx,By,Bz)
            #step 2 v^n+1/2 ----- v^&{n+3/2}
            [u_finalx,u_finaly,u_finalz] = BorisA(Ex,Ey,Ez,Bx,By,Bz,ux0,uy0,uz0,q,m,dt)
        #x^{n} ----- x^&{n+1}
        
        [pos_finalx,pos_finaly,pos_finalz] = leapfrog(x0,y0,z0,u_finalx,u_finaly,u_finalz,dt)
        #guardar nuevas posiciones y velocidades
        particulas[j,0,i] = pos_finalx
        particulas[j,1,i] = pos_finaly
        particulas[j,2,i] = pos_finalz
        particulas[j,3,i] = u_finalx
        particulas[j,4,i] = u_finaly
        particulas[j,5,i] = u_finalz

        #zigzag(pos_finalx, pos_finaly, pos_finalz, dx, dy)
        #Energía

        
        energia_cinetica [i] = 0.5*m*(u_finalx**2+u_finaly**2+u_finalz**2)
        epsilon_r [i] = (energia_cinetica[i]- energia_cinetica[0] )/energia_cinetica[0]









#Gráficas


fig, (ax1, ax2, ax3) = plt.subplots(3,1)

ax1.plot(energia_cinetica)
ax1.set_title('Energía cinética a través del tiempo ')
ax1.set_xlabel('t')
ax1.set_ylabel('Ke')

datax = []
datay = []

for ii in range(0,i):
    datax.append(particulas[0,0,ii])
    datay.append(particulas[0,1,ii])

dataxNP = np.array(datax)
datayNP = np.array(datay)

ax2.plot(dataxNP, datayNP,'tab:green')
ax2.set_title('Trayectoria de la particula')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')

ax3.plot(epsilon_r)
ax3.set_title('Error relativo')
ax3.set_xlabel('t')
ax3.set_ylabel('Ke')

plt.show()