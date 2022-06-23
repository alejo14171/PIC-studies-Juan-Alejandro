import numpy as np                                                  
import struct                                                       
import matplotlib.pyplot as plt                                     
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import glob
import sys
 
#read particles from processor 0 
files = sorted(glob.glob('./output/ion/*_00.dat'))
nout = len(files)  

#uu[#particles,#files,#frames]
# #particles in each processor
#files =[x_p,y_p,vx_p,vy_p,vz_p,Ex_p,Ey_p,Ez_p,Bx_p,By_p,Bz_p,Vpar_p,Vperp_p

len_array = np.shape(np.loadtxt(files[0]) )
uu = np.zeros([len_array[0],len_array[1],nout]) 
for i in range(nout-1):
      print(files[i])                                                 
      uu[:,:,i] = np.loadtxt(files[i]) 
