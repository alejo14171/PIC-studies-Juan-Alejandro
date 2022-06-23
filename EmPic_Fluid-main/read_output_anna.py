import numpy as np
import struct
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import glob
import sys


def read_parallel_info():
    # read parallel information ----------------
    file_prl = open('./parallel_info.dat','rb')

    skip = struct.unpack("f",file_prl.read(4))

    npe = int((struct.unpack("f",file_prl.read(4)))[0])
    nvar = int((struct.unpack("f",file_prl.read(4)))[0])

    skip = struct.unpack("f",file_prl.read(4))

    file_prl.close()
    return npe, nvar


def read_grid():
    # read nx, ny, nz and grid-------------------------
    file_grid = open('./grid.dat', 'rb')

    skip = struct.unpack("f",file_grid.read(4))

    nx = int((struct.unpack("f",file_grid.read(4)))[0])
    ny = int((struct.unpack("f",file_grid.read(4)))[0])

    skip = struct.unpack("f",file_grid.read(4))

    xgrid = np.zeros(nx)
    ygrid = np.zeros(ny)

    skip = struct.unpack("f",file_grid.read(4))

    for i in range(nx):
        xgrid[i] = (struct.unpack("f",file_grid.read(4)))[0]
    for i in range(ny):
        ygrid[i] = (struct.unpack("f",file_grid.read(4)))[0]

    skip = struct.unpack("f",file_grid.read(4))

    file_grid.close()
    return xgrid,ygrid

def read_EBM():
    #read the EBM_info.dat
    #which includes time,radius,Ur

    file_EBM = np.array(np.loadtxt('./EBM_info.dat'))

    if len(file_EBM.shape)==1:
        file_EBM = np.reshape(file_EBM,(int(len(file_EBM)/3),3))

    t_EBM = file_EBM[:,0]
    radius = file_EBM[:,1]
    Ur_EBM = file_EBM[:,2]

    return t_EBM,radius,Ur_EBM
def read_uu(filename,nx,ny,nvar):
    file_uu = open(filename, 'rb')
    uu = np.zeros([nx,ny,nvar])

    skip = (struct.unpack("f",file_uu.read(4)))[0]
    t = (struct.unpack("f",file_uu.read(4)))[0] 
    skip = (struct.unpack("f",file_uu.read(4)))[0]

    for ivar in range(nvar):
            for iy in range(ny):
                for ix in range(nx):
                    uu[ix,iy,ivar] = (struct.unpack(\
                        "d",file_uu.read(8)))[0] 

    return t, uu


xgrid, ygrid = read_grid()
nx,ny= len(xgrid), len(ygrid)
print('nx, ny =', nx, ny)


npe,nvar = read_parallel_info()                                 
print('npe, nvar = ', npe, nvar)

t_EBM, radius, Ur = read_EBM()

files = sorted(glob.glob('./out*dat'))
nout = len(files)

time = np.zeros(nout)
rho = np.zeros([nx,ny,nout])
ux = np.zeros([nx,ny,nout])
uy = np.zeros([nx,ny,nout])
uz = np.zeros([nx,ny,nout])
Bx = np.zeros([nx,ny,nout])
By = np.zeros([nx,ny,nout])
Bz = np.zeros([nx,ny,nout])
p = np.zeros([nx,ny,nout])

for nt in range(nout):
    t, uu = read_uu(files[nt],nx,ny,nvar)

    print('t = {:.3f}'.format(t))

    time[nt] = t 

    rho[:,:,nt] = uu[:,:,0]

    ux[:,:,nt] = uu[:,:,1]
    uy[:,:,nt] = uu[:,:,2]
    uz[:,:,nt] = uu[:,:,3]

    Bx[:,:,nt] = uu[:,:,4]
    By[:,:,nt] = uu[:,:,5]
    Bz[:,:,nt] = uu[:,:,6]

    p[:,:,nt] = uu[:,:,7]

avg_rho = np.zeros(nout)
avg_bx = np.zeros(nout)
avg_bz = np.zeros(nout)
avg_uz = np.zeros(nout)
avg_p = np.zeros(nout)
avg_by = np.zeros(nout)
var_by = np.zeros(nout)
var_bz = np.zeros(nout)
var_bx = np.zeros(nout)
var_rho = np.zeros(nout)
var_p = np.zeros(nout)
var_uz = np.zeros(nout)

for k in range(nout):
    avg_rho[k] = np.sum(rho[:,:,k])
    avg_bx[k] = np.sum(Bx[:,:,k])
    avg_bz[k] = np.sum(Bz[:,:,k])
    avg_by[k] = np.sum(By[:,:,k])
    avg_uz[k] = np.sum(uz[:,:,k])
    avg_p[k] = np.sum(p[:,:,k])

avg_rho= avg_rho/(nx*ny)
avg_bx= avg_bx/(nx*ny)
avg_bz= avg_bz/(nx*ny)
avg_by= avg_by/(nx*ny)
avg_uz= avg_uz/(nx*ny)
avg_p= avg_p/(nx*ny)

for k in range(nout):
    var_rho[k] = np.sum((rho[:,:,k]-avg_rho[k])**2)
    var_bx[k] = np.sum((Bx[:,:,k]-avg_bx[k])**2)
    var_bz[k] = np.sum((By[:,:,k]-avg_by[k])**2)
    var_by[k] = np.sum((Bz[:,:,k]-avg_bz[k])**2)
    var_uz[k] = np.sum((uz[:,:,k]-avg_uz[k])**2)
    var_p[k] = np.sum((p[:,:,k]-avg_p[k])**2)

var_bz = var_bz/(nx*ny)
var_by = var_by/(nx*ny)
var_bx = var_bx/(nx*ny)
var_rho = var_rho/(nx*ny)
var_uz = var_uz/(nx*ny)
var_p = var_p/(nx*ny)

bzl=np.zeros([nx,ny,nout],dtype=complex)
byl=np.zeros([nx,ny,nout],dtype=complex)
bxl=np.zeros([nx,ny,nout],dtype=complex)
rhol=np.zeros([nx,ny,nout],dtype=complex)

for i in range(nout):
    byl[:,:,i] = (np.fft.fft2(By[:,:,i]))/(nx*ny)
    bxl[:,:,i] = (np.fft.fft2(Bx[:,:,i]-1))/(nx*ny)
    bzl[:,:,i] = (np.fft.fft2(Bz[:,:,i]))/(nx*ny)
    rhol[:,:,i] = (np.fft.fft2(rho[:,:,i]-avg_rho[k]))/(nx*ny)
