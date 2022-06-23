import numpy as np
import matplotlib.pyplot as plt

def BorisA(Ex,Ey,Ez,Bx,By,Bz,velx,vely,velz,q,m,dt):
  #print(Ex,Ey,Ez,Bx,By,Bz,velx,vely,velz)
  #paso 1 ecuacion (3) del documento
  ux0 = velx
  uy0 = vely
  uz0 = velz
  u_menosx= ux0 + (Ex*q/(2.0*m))
  u_menosy= uy0 + (Ey*q/(2.0*m))
  u_menosz= uz0 + (Ez*q/(2.0*m))
  #Paso 2 ecuacion (6) del documento
  tx=q*dt*Bx/(2.0*m)
  ty=q*dt*By/(2.0*m)
  tz=q*dt*Bz/(2.0*m)
  #paso 3 ecuacion (8) del documento. (aqui esta pendiente de acabar)
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
    return Ex,Ey,Ez,Bz,By,Bz
    

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
nout =  int(tfinal/dt)#numero de iteraciones
#particles = [npart,{x,y,z,vx,vy,vz},nout]
particulas= np.zeros([npart,6,nout])
q=1.0
m=1.0

#crear malla y campos
xx = np.linspace(xmin,xmax,Nx+1)
yy = np.linspace(ymin,ymax,Ny+1)
x,y = np.meshgrid(xx,yy)

#distribucion inicial de particulas (t=0)
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
    #print(i)
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
            Ex,Ey,Ez,Bz,By,Bz  = get_fields_nodes(ex,ey,ez,bx,by,bz,ileft,iright,hxleft,jdown,jup,hydown)
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
            Ex,Ey,Ez,Bz,By,Bz  = get_fields_nodes(ex,ey,ez,bx,by,bz,ileft,iright,hxleft,jdown,jup,hydown)
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


plt.figure(10);plt.clf();plt.imshow(bz,extent=[-1.0*np.pi,1.0*np.pi,-1.0*np.pi,1.0*np.pi],aspect='auto')

for ii in range(0,i):
    plt.plot(particulas[0,0,ii],particulas[0,1,ii],'ok')
