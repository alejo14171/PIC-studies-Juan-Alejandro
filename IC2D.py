import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import math

#import clases, funciones
from funciones import *

# definicion de malla computacional
xmin = -1.0*np.pi
xmax = 1.0*np.pi
ymin = -1.0*np.pi
ymax = 1.0*np.pi
Nx = 64
Ny = 64
Lx = xmax-xmin
Ly = ymax-ymin
dx = Lx/Nx
dy = Ly/Ny
dt = np.pi/10
tfinal = 1000.0*np.pi
npart = 1  # number of particles
nout = int(tfinal/dt)  # numero de iteraciones+
#particles = [npart,{x,y,z,vx,vy,vz},nout]
particulas = np.zeros([npart, 6, nout])
energia_cinetica = np.zeros(nout)
epsilon_r = np.zeros(nout)
q = 1.0
m = 1.0

# crear malla y campos
xx = np.linspace(xmin, xmax, Nx+1)
yy = np.linspace(ymin, ymax, Ny+1)
x, y = np.meshgrid(xx, yy)


# distribucion inicial de particulas (t=0) para 1 particula
for i in range(npart):
    # posicion x[numero de particulas,posx,posy,posz,vx,vy,vz]
    particulas[i, 0, 0] = 2.0
    particulas[i, 1, 0] = 0.0  # posicion y
    particulas[i, 2, 0] = 0.0  # posicion z
    particulas[i, 3, 0] = 0.75  # posicion vx
    particulas[i, 4, 0] = 0.0  # posicion vy
    particulas[i, 5, 0] = 0.0  # posicion vz


# campo electromagnetico (B & E on the nodes)
bx = np.zeros([Nx+1, Ny+1])
by = np.zeros([Nx+1, Ny+1])
bz = np.zeros([Nx+1, Ny+1])
ex = np.zeros([Nx+1, Ny+1])
ey = np.zeros([Nx+1, Ny+1])
ez = np.zeros([Nx+1, Ny+1])

for iy in range(0, Ny+1):
    for ix in range(0, Nx+1):
        tmp_far = x[ix, iy]**2 + y[ix, iy]**2
        # print(tmp_far)
        if(tmp_far == 0):
            bz[ix, iy] = 0.0
            ex[ix, iy] = 0.0
            ey[ix, iy] = 0.0
        else:
            bz[ix, iy] = np.sqrt(tmp_far)
            ex[ix, iy] = (0.01*x[ix, iy])/pow(tmp_far, 1.5)
            ey[ix, iy] = (0.01*y[ix, iy])/pow(tmp_far, 1.5)


# version (2) leap-frog: posicion a tiempos enteros (n) y velocidad a tiempo intermedio (n+1/2)
Bx = 0
By = 0
Bz = 0
Ex = 0.0
Ey = 0.0
Ez = 0.0
for i in range(nout):
    # cada partícula
    for j in range(npart):
        # step 1-a v^n ----- v^&{n+1/2} (solo en t=0)
        if(i == 0):
            x0 = particulas[j, 0, i]
            y0 = particulas[j, 1, i]
            z0 = particulas[j, 2, i]
            ux0 = particulas[j, 3, i]
            uy0 = particulas[j, 4, i]
            uz0 = particulas[j, 5, i]
            # found the closest node indexes according to particle position
            ileft, iright, hxleft, jdown, jup, hydown = interpolation(
                x0, y0, z0)
            # linear interpolation of fields (nodes) to particle position
            Ex, Ey, Ez, Bx, By, Bz = get_fields_nodes(
                ex, ey, ez, bx, by, bz, ileft, iright, hxleft, jdown, jup, hydown)
            # Advance particle velocity from n to n+1/2
            [u_finalx, u_finaly, u_finalz] = BorisA(
                Ex, Ey, Ez, Bx, By, Bz, ux0, uy0, uz0, q, m, 0.5*dt)
        else:
            # step 1-b v^n+1/2 ----- v^&{n+3/2} (when t!=0)
            x0 = particulas[j, 0, i-1]
            y0 = particulas[j, 1, i-1]
            z0 = particulas[j, 2, i-1]
            ux0 = particulas[j, 3, i-1]
            uy0 = particulas[j, 4, i-1]
            uz0 = particulas[j, 5, i-1]
            # found the closest node indexes according to particle position
            ileft, iright, hxleft, jdown, jup, hydown = interpolation(
                x0, y0, z0)
            # linear interpolation of fields (nodes) to particle position
            Ex, Ey, Ez, Bx, By, Bz = get_fields_nodes(
                ex, ey, ez, bx, by, bz, ileft, iright, hxleft, jdown, jup, hydown)
            # print(Ex,Ey,Ez,Bx,By,Bz)
            # step 2 v^n+1/2 ----- v^&{n+3/2}
            [u_finalx, u_finaly, u_finalz] = BorisA(
                Ex, Ey, Ez, Bx, By, Bz, ux0, uy0, uz0, q, m, dt)
        # x^{n} ----- x^&{n+1}

        [pos_finalx, pos_finaly, pos_finalz] = leapfrog(
            x0, y0, z0, u_finalx, u_finaly, u_finalz, dt)
        [J1x1, J1x2, J1y1, J1y2, J2x1, J2x2, J2y1, J2y2] = zig_zag(
            x0, y0, pos_finalx, pos_finaly)
        # guardar nuevas posiciones y velocidades
        particulas[j, 0, i] = pos_finalx
        particulas[j, 1, i] = pos_finaly
        particulas[j, 2, i] = pos_finalz
        particulas[j, 3, i] = u_finalx
        particulas[j, 4, i] = u_finaly
        particulas[j, 5, i] = u_finalz

        #zigzag(pos_finalx, pos_finaly, pos_finalz, dx, dy)
        # Energía
        energia_cinetica[i] = 0.5*m*(u_finalx**2+u_finaly**2+u_finalz**2)
        epsilon_r[i] = (energia_cinetica[i] -
                        energia_cinetica[0])/energia_cinetica[0]


# Gráficas


fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

ax1.plot(energia_cinetica)
ax1.set_title('Energía cinética a través del tiempo ')
ax1.set_xlabel('t')
ax1.set_ylabel('Ke')

datax = []
datay = []

for ii in range(0, i):
    datax.append(particulas[0, 0, ii])
    datay.append(particulas[0, 1, ii])

dataxNP = np.array(datax)
datayNP = np.array(datay)

ax2.plot(dataxNP, datayNP, 'tab:green')
ax2.set_title('Trayectoria de la particula')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')

ax3.plot(epsilon_r)
ax3.set_title('Error relativo')
ax3.set_xlabel('t')
ax3.set_ylabel('Ke')

plt.show()
