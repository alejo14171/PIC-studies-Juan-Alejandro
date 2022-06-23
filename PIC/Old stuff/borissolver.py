from re import M
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from numpy import array, transpose, sum, append, delete, shape, sqrt, linspace, sin



# Caso 1:  Campo Electrico uniforme
Ex,Ey,Ez=[1],[0],[0]
Bx,By,Bz=[1],[0],[0]
x, y, z  = [1], [0], [0]
Vx, Vy, Vz= [1], [0] , [0]
M=1
Q=1
Δt=1*10**-9
N=1
T=20

#paso de 1/2 de la velocidad


for i in range(0,N):
            
            QM = Q/M
            
            Vx[i] += QM* (Ex[i] + Bz[i]* Vy[i] - By[i]* Vz[i])* Δt /2
            Vy[i] += QM* (Ey[i] + Bx[i]* Vz[i] - Bz[i]* Vx[i])* Δt /2
            Vz[i] += QM* (Ez[i] + By[i]* Vx[i] - Bx[i]* Vy[i])* Δt /2

# Step - 2 : Boris Algorithm
t  = 0
tan = [0]*N
hx, hy, hz = [0]* N, [0]* N, [0]* N
Sx, Sy, Sz = [0]* N, [0]* N, [0]* N
Ux, Uy, Uz = [0]* N, [0]* N, [0]* N
Ux_d, Uy_d, Uz_d = [0]* N, [0]* N, [0]*N

        
tan = (Q*Δt)/ (2* M)
        
        
while (t <= T):
            
    for i in range(0, N):
            
        hx[i]   = tan* Bx[i]
        hy[i]   = tan* By[i]
        hz[i]   = tan* Bz[i]
                
        Sx[i]   = 2* hx[i]/ (1 + ( hx[i]**2 + hy[i]**2 + hz[i]**2 ) )
        Sy[i]   = 2* hy[i]/ (1 + ( hx[i]**2 + hy[i]**2 + hz[i]**2 ) )
        Sz[i]   = 2* hz[i]/ (1 + ( hx[i]**2 + hy[i]**2 + hz[i]**2 ) )
                
        Ux[i]   = Vx[i] + tan* Ex[i]
        Uy[i]   = Vy[i] + tan* Ey[i]
        Uz[i]   = Vz[i] + tan* Ez[i]
                
        Ux_d[i] = Ux[i] + Sz[i]* ( Uy[i] + Uz[i]*hx[i] - hz[i]*Ux[i] ) - Sy[i]* ( Uz[i] + Ux[i]*hy[i] - hx[i]*Uy[i] )
        Uy_d[i] = Uy[i] + Sx[i]* ( Uz[i] + Ux[i]*hy[i] - hx[i]*Uy[i] ) - Sz[i]* ( Ux[i] + Uy[i]*hz[i] - hy[i]*Uz[i] )
        Uz_d[i] = Uz[i] + Sy[i]* ( Ux[i] + Uy[i]*hz[i] - hy[i]*Uz[i] ) - Sx[i]* ( Uy[i] + Uz[i]*hx[i] - hz[i]*Ux[i] )
                
        Vx[i] = Ux_d[i] + tan* Ex[i]
        Vy[i] = Uy_d[i] + tan* Ey[i]
        Vz[i] = Uz_d[i] + tan* Ez[i]
                
        x[i] += Vx[i]* Δt
        y[i] += Vy[i]* Δt
        z[i] += Vz[i]* Δt


        fig=plt.figure()
        ax= plt.axes(projection="3d")

        ax.plot3D(x[i],y[i],z[i],'red')

        

    t += Δt
    





