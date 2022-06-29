import matplotlib as plt
import numpy as np
import random as rnd
import math

# Constantes
q = 1

# Números aleatorios pero peques Δ
Δt = 0.000001 
Δx = 0.1
Δy = 0.1 

vx = rnd.random() #velocidad en x en t + Δt/2
vy = rnd.random() #velocidad en y en t + Δt/2

x1 = 0  # posición x en t
y1 = 0  # posición y en t

x2 = x1 + vx*Δt
y2 = y1 + vy*Δt

# largest integer values not greater than x/Δx and y/Δy respectively
i1 = math.floor(x1/Δx)
i2 = math.floor(x2/Δx)
j1 = math.floor(y1/Δy)
j2 = math.floor(y2/Δy)

"""WHEN A PARTICLE REMAINS IN THE SAME CELL square during its movement, a charge
 flux of a particle can be computed from the start point (x1, y1) and end point
 (x2, y2) of the particle movement. Charge conservation of the assigned
current densities is realized with the procedure described by Eastwood """

# First order shape factor
Wx = ((x1+x2)/2*Δx) - i1
Wy = ((y1+y2)/2*Δx) - j1

# Charge flux
Fx = q((x2-x1)/Δt)
Fy = q((y2-y1)/Δt)

# Constant
c = 1/Δx*Δy

# Current density
Jx1_Eastwood = c*Fx*(1-Wy)
Jx2_Eastwood = c*Fx*(Wy)
Jy1_Eastwood = c*Fy*(1-Wx)
Jy1_Eastwood = c*Fy*(Wy)


"""WHEN A PARTICLE MOVES ACROSS CELL MESHES, we decompose the particle movement with a special 
assignment pattern as shown in Figs. 3 and 4. We have to think of the following four cases 
depending on positions of the particle at time t and t + Δt: (a) i1 = i2 and j1 ≠ j2, 
(b) i1 ≠ i2 and j1 = j2, (c) i1 = i2 and j1 ≠ j2, and (d) i1 = i2 and j1 = j2."""

# Practicamente, para resumir estos casos a dos lineas (sin usar if) se usa la siguiente 
# expresión teniendo en cuenta que el movimiento de la partícula se descompone en dos:
# desde (x1,y1) hasta (xr, yr) y desde (xr, yr) hasta (x2, y2). Las coordenadas (xr, yr)
# son:

xr = min(min(i1*Δx, i2*Δx) + Δx, max(max(i1*Δx, i2*Δx), (x1+x2)/2))
yr = min(min(j1*Δy, j2*Δy) + Δy, max(max(j1*Δy, j2*Δy), (y1+y2)/2))

# Así mismo como se descompuse la trayectoria, también se descompone el flujo así:
# (F1, F2) = ([Fx1, Fy1], [Fx2, Fy2]):

Fx1 = q*(xr-x1)/Δt
Fy1 = q*(yr-y1)/Δt
Fx2 = q*(x2-xr)/Δt
Fy2 = q*(y2-yr)/Δt

# Antes de hayar el flujo de carga, se calculan los pesos respectivos a cada trayectoria
# Por tanto, dos para F1 y dos para F2:

Wx1 = ((x1+xr)/2*Δx) - i1
Wy1 = ((y1+yr)/2*Δy) - j1
Wx2 = ((xr+x2)/2*Δx) - i2
Wy2 = ((yr+y2)/2*Δx) - j2

# Por último, se agregan los segmentos de carga a 8 puntos en la malla, los cuales son:

J1x1 = c*Fx1*(1-Wy1) #Jx(i1+1/2, j1)
J1x2 = c*Fx1*(Wy1)   #Jx(i1+1/2, j1+1)
J1y1 = c*Fy1*(1-Wx1) #Jy(i1, j1+1/2)
J1y2 = c*Fy1*(Wx1)   #Jy(i1+1, j1+1/2)

J2x1 = c*Fx2*(1-Wy2) #Jx(i2+1/2, j2)
J2x2 = c*Fx2*(Wy2)   #Jx(i2+1/2, j2+1)
J2y1 = c*Fy2*(1-Wx2) #Jy(i2, j2+1/2)
J2y2 = c*Fy2*(Wx2)   #Jy(i2+1, j2+1/2)

# Finally we superpose the charge fluxes contributed by each particle to obtain the total
#  current densities.


def zig_zag(x1, y1, x2, y2, Δx, Δy, Δt, q): # Δx, Δy, Δt, q pueden ser globales

    i1 = math.floor(x1/Δx)
    i2 = math.floor(x2/Δx)
    j1 = math.floor(y1/Δy)
    j2 = math.floor(y2/Δy)

    xr = min(min(i1*Δx, i2*Δx) + Δx, max(max(i1*Δx, i2*Δx), (x1+x2)/2))
    yr = min(min(j1*Δy, j2*Δy) + Δy, max(max(j1*Δy, j2*Δy), (y1+y2)/2))

    Fx1 = q*(xr-x1)/Δt
    Fy1 = q*(yr-y1)/Δt
    Fx2 = q*(x2-xr)/Δt
    Fy2 = q*(y2-yr)/Δt

    Wx1 = ((x1+xr)/2*Δx) - i1
    Wy1 = ((y1+yr)/2*Δy) - j1
    Wx2 = ((xr+x2)/2*Δx) - i2
    Wy2 = ((yr+y2)/2*Δx) - j2

    J1x1 = c*Fx1*(1-Wy1) #Jx(i1+1/2, j1)
    J1x2 = c*Fx1*(Wy1)   #Jx(i1+1/2, j1+1)
    J1y1 = c*Fy1*(1-Wx1) #Jy(i1, j1+1/2)
    J1y2 = c*Fy1*(Wx1)   #Jy(i1+1, j1+1/2)

    J2x1 = c*Fx2*(1-Wy2) #Jx(i2+1/2, j2)
    J2x2 = c*Fx2*(Wy2)   #Jx(i2+1/2, j2+1)
    J2y1 = c*Fy2*(1-Wx2) #Jy(i2, j2+1/2)
    J2y2 = c*Fy2*(Wx2)   #Jy(i2+1, j2+1/2)

    return(J1x1, J1x2, J1y1, J1y2, J2x1, J2x2, J2y1, J2y2)