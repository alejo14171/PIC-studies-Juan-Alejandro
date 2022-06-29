from cmath import sqrt
import numpy as np
class PuntosDeMalla:
    def __init__(self): #Por convención se le pone self
        self.x;
        self.y;
        self.z;
        self.phi;
        self.den_i;
        self.ux_i;
        self.uy_i;
        self.uz_i;
        self.pxx_i;
        self.pxy_i;
        self.pyy_i;
        self.pxz_i;
        self.pyz_i;
        self.pzz_i;
        self.ex;
        self.ey;
        self.ez;
        self.bx;
        self.by;
        self.bz;
        self.jx;
        self.jy;
        self.jz;
        self.fvpape;
        self.fxvx;

class PuntosDeCampo:
    def __init__(self):
        self.x;
        self.y;
        self.z;
        self.ex;
        self.ey;
        self.ez;
        self.bx; 
        self.by;
        self.bz;
        self.jx;
        self.jy;
        self.jz;

class TresDVector:
    def __init__(self, *args):
        self.coordenadas = np.zeros(3)
    
    def llenar(self, x, y, z):
        self.coordenadas[0] = x
        self.coordenadas[1] = y
        self.coordenadas[2] = z

    def setX(self, x):
        self.coordenadas[0] = x

    def setY(self, y):
        self.coordenadas[0] = y

    def setX(self, z):
        self.coordenadas[0] = z

    def magnitud(self):
        x = self.coordenadas[0]
        y = self.coordenadas[1]
        z = self.coordenadas[2]
        return sqrt(x^2 + y^2 + z^2)

    def printV(self):
        print('x', self.coordenadas[0])
        print('y', self.coordenadas[1])
        print('z', self.coordenadas[2])