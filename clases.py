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
        self.coordenadas = np.array(args)
        self.coordenadas = np.zeros(3)
    
