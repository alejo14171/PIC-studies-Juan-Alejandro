import numpy as np                                                  
import struct                                                       
import matplotlib.pyplot as plt                                     
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import glob
import sys
 
#read particles from processor 0 
files = sorted(glob.glob('Density*'))
nout = len(files)
Nx=512
Ny=512
rho = np.zeros([Nx,Ny,nout])

for i in range(nout):                                               
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      rho[:,:,i] = tmp


files = sorted(glob.glob('Ux*'))
nout = len(files)
Ux = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      Ux[:,:,i] = tmp

files = sorted(glob.glob('Uy*'))
nout = len(files)
Uy = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      Uy[:,:,i] = tmp


files = sorted(glob.glob('Uz*'))
nout = len(files)
Uz = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      Uz[:,:,i] = tmp


files = sorted(glob.glob('Fvpape*'))
nout = len(files)
fvpape = np.zeros([Nx,Ny,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx,Ny],order = 'C')
      fvpape[:,:,i] = tmp      
