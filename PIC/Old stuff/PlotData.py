

from numpy import array, transpose, sum, append, delete, shape, sqrt, linspace, sin
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


class Energy:

    '''
        Esta clase va a calcular la Energía total de la particula
    '''

    def __init__(self, Q, M, N, x, y, z, Vx, Vy, Vz):

        self.Q = Q
        self.M,  self.N  = M,  N
        self.x,  self.y,  self.z  = x,  y,  z
        self.Vx, self.Vy, self.Vz = Vx, Vy, Vz


    def Pot(self, x, y, z):
        
        W = 0

        for i in range(0, self.N):
        
            # Case - 1 : External Force Potential
            # W -= x[i]* 1 + y[i]* 0 + z[i]* 0
            
            # Case - 2 : Nonuniform E field
            W -= 0.5* sin(2* x[i])
            
            # Case - 3 : Interaction Potential
            # for j in range(0, self.N):
            #     if (i != j):

            #         r = ( (x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2 )**0.5

            #         W += 0.5* self.Q[i]* self.Q[j]/ r

        return W


    def EnergyCalculations(self):

        # kinetic energy
        KE = 0
        for i in range(0, len(self.M)):
            KE += 0.5* self.M[i]* ( self.Vx[:,i]**2 + self.Vy[:,i]**2 + self.Vz[:,i]**2 )

        # electrostatic potential energy
        p, _ = shape( self.x )
        PE = array([])
        for i in range(0, p):
            PE_dt = self.Pot( self.x[i,:], self.y[i,:], self.z[i,:] )
            PE = append(PE, PE_dt)

        # Total energy
        E = KE + PE

        # Energy per particle
        E  = E / self.N
        KE = KE/ self.N
        PE = PE/ self.N

        return E, KE, PE

#%%

class PlotData:

    def __init__(self):

        # Read Data Files
        M    = pd.read_csv(r'Data/M.txt',  header = None, sep = ' ')
        Q    = pd.read_csv(r'Data/Q.txt',  header = None, sep = ' ')
        x    = pd.read_csv(r'Data/x.txt',  header = None, sep = ' ')
        y    = pd.read_csv(r'Data/y.txt',  header = None, sep = ' ')
        z    = pd.read_csv(r'Data/z.txt',  header = None, sep = ' ')
        Vx   = pd.read_csv(r'Data/vx.txt', header = None, sep = ' ')
        Vy   = pd.read_csv(r'Data/vy.txt', header = None, sep = ' ')
        Vz   = pd.read_csv(r'Data/vz.txt', header = None, sep = ' ')
        data = pd.read_csv('inputs.txt', header = None, sep = '=')
        data = data[1].tolist()

        self.M  = M.values.tolist()[0]
        self.Q  = Q.values.tolist()[0]
        self.N  = int(data[0])
        self.Δt = float(data[2])
        self.x  = array(x)
        self.y  = array(y)
        self.z  = array(z)
        self.Vx = array(Vx)
        self.Vy = array(Vy)
        self.Vz = array(Vz)
        
        self.Velocity()
        
        p, _ = shape(self.x)
        self.time = linspace(0, (p-1)* data[2], p)

        self.E, self.KE, self.PE = Energy(self.Q, self.M, self.N, self.x, self.y, self.z, self.Vx, self.Vy, self.Vz).EnergyCalculations()


    def Velocity(self):
        
        Vel = lambda Pos, Δt : [(Pos[i+2] - Pos[i])/ (2* Δt) for i in range(0, len(Pos) - 2)]
        
        Dmy_Vx, Dmy_Vy, Dmy_Vz = [0]* self.N, [0]* self.N, [0]* self.N
        
        for i in range(0, self.N):
            
            Dmy_Vx[i] = [self.Vx[0][i]] + Vel(self.x[:,i], self.Δt)
            Dmy_Vy[i] = [self.Vy[0][i]] + Vel(self.y[:,i], self.Δt)
            Dmy_Vz[i] = [self.Vz[0][i]] + Vel(self.z[:,i], self.Δt)
            
        self.Vx = transpose( array(Dmy_Vx) )
        self.Vy = transpose( array(Dmy_Vy) )
        self.Vz = transpose( array(Dmy_Vz) )
            
        self.x = self.x[0:int(len(self.x[:,0])-1),:]
        self.y = self.y[0:int(len(self.y[:,0])-1),:]
        self.z = self.z[0:int(len(self.z[:,0])-1),:]





























    def AnimamteParticles(self, Skip = 20):

        gs = GridSpec(1, 2)

        fig = plt.figure( figsize = (80, 60) )
        ax1 = fig.add_subplot( gs[0, 0], projection = '3d' )
        ax2 = fig.add_subplot( gs[0, 1] )

        for frames, t in enumerate( self.time[::Skip] ):

            fS = frames* Skip

            # PLOT - 1: Particles
            ax1.clear()

            ax1.scatter( self.x[fS], self.y[fS], self.z[fS], marker = 'o' )

            Vmag = ( self.Vx[fS]**2 + self.Vy[fS]**2 + self.Vz[fS]**2 )**0.5

            ax1.quiver(  self.x[fS],  self.y[fS],  self.z[fS], \
                        self.Vx[fS], self.Vy[fS], self.Vz[fS], \
                        normalize = True)

            ax1.set_xlabel('x', fontsize = 16)
            ax1.set_ylabel('y', fontsize = 16)
            ax1.set_zlabel('z', fontsize = 16)
            ax1.set_xlim([-5, 5])
            ax1.set_ylim([-5, 5])
            ax1.set_zlim([-5, 5])

            # PLOT - 2: Energy
            if (frames == 0):

                e  = self.E[fS]
                ke = self.KE[fS]
                pe = self.PE[fS]
                to = t

                ax2.scatter(to,  e, label = r'$E_{sys}$')
                ax2.scatter(to, ke, label = r'KE')
                ax2.scatter(to, pe, label = r'PE')

                ax2.legend()

            else:

                ax2.plot([to,t], [ e, self.E[fS]  ], color = 'blue')
                ax2.plot([to,t], [ke, self.KE[fS] ], color = 'orange')
                ax2.plot([to,t], [pe, self.PE[fS] ], color = 'green')

                e  = self.E[fS]
                ke = self.KE[fS]
                pe = self.PE[fS]
                to = t

            plt.pause(1e-11)


    def PlotEnergy(self):

        fig = plt.figure('Energy conservation')
        plt.plot( self.time, self.E,  label = r'$E_{sys}$')
        plt.plot( self.time, self.KE, label = r'KE')
        plt.plot( self.time, self.PE, label = r'PE')
        plt.legend()
        plt.grid()
        plt.show()




if __name__ == '__main__':

    DoMyWork = PlotData()
    # DoMyWork.AnimamteParticles(Skip = 75)
    DoMyWork.PlotEnergy()