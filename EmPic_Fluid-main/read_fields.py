import numpy as np                                                  
import struct                                                       
import matplotlib.pyplot as plt                                     
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import glob
import sys
 
#read particles from processor 0 
files = sorted(glob.glob('bx_c*'))
nout = len(files)
Nx=64
Ny=64
bx_c = np.zeros([Nx,Ny,nout])

for i in range(nout):                                               
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      bx_c[:,:,i] = tmp


files = sorted(glob.glob('by_c*'))
nout = len(files)
by_c = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      by_c[:,:,i] = tmp

files = sorted(glob.glob('bz_c*'))
nout = len(files)
bz_c = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      bz_c[:,:,i] = tmp


files = sorted(glob.glob('ex_c*'))
nout = len(files)
ex_c = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      ex_c[:,:,i] = tmp


files = sorted(glob.glob('ey_c*'))
nout = len(files)
ey_c = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      ey_c[:,:,i] = tmp

files = sorted(glob.glob('ez_c*'))
nout = len(files)
ez_c = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      ez_c[:,:,i] = tmp



#read particles from processor 0
files = sorted(glob.glob('bx_n*'))
nout = len(files)

bx_n = np.zeros([Nx+1,Ny+1,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      bx_n[:,:,i] = tmp


files = sorted(glob.glob('by_n*'))
nout = len(files)
by_n = np.zeros([Nx+1,Ny+1,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      by_n[:,:,i] = tmp

files = sorted(glob.glob('bz_n*'))
nout = len(files)
bz_n = np.zeros([Nx+1,Ny+1,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      bz_n[:,:,i] = tmp


files = sorted(glob.glob('ex_n*'))
nout = len(files)
ex_n = np.zeros([Nx+1,Ny+1,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      ex_n[:,:,i] = tmp


files = sorted(glob.glob('ey_n*'))
nout = len(files)
ey_n = np.zeros([Nx+1,Ny+1,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      ey_n[:,:,i] = tmp

files = sorted(glob.glob('ez_n*'))
nout = len(files)
ez_n = np.zeros([Nx+1,Ny+1,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      ez_n[:,:,i] = tmp
