def BorisA(Ex,Ey,Ez,Bx,By,Bz,velx,vely,velz,q,m,dt):
  #print(Ex,Ey,Ez,Bx,By,Bz,velx,vely,velz)

  #Hallar la velocidad final 

  #paso 1 ecuacion (3) del documento (Inicicializar velocidades)
  ux0 = velx
  uy0 = vely
  uz0 = velz

  #Velocidad un más adelante "(Ex*q/(2.0*m))""
  u_menosx= ux0 + (Ex*q/(2.0*m))
  u_menosy= uy0 + (Ey*q/(2.0*m))
  u_menosz= uz0 + (Ez*q/(2.0*m))

  #Ángulo de fase de la rotación 

  #Paso 2 ecuacion (6) del documento
  tx=q*dt*Bx/(2.0*m)
  ty=q*dt*By/(2.0*m)
  tz=q*dt*Bz/(2.0*m)

  #paso 3 ecuacion (8) del documento. 
  u_primax=u_menosx + ((u_menosy*tz)-(u_menosz*ty))
  u_primay=u_menosy + ((u_menosz*tx)-(u_menosx*tz))
  u_primaz=u_menosz + ((u_menosx*ty)-(u_menosy*tx))

  # Agrupando terminos:
  sx=2*tx/(1.0+(tx*tx))
  sy=2*ty/(1.0+(ty*ty))
  sz=2*tz/(1.0+(tz*tz))
  # Paso 4 ecuacion(9)
  u_masx = u_menosx + ((u_primay*sz) - (u_primaz*sy))
  u_masy = u_menosy + ((u_primaz*sx) - (u_primax*sz))
  u_masz = u_menosz + ((u_primax*sy) - (u_primay*sx))
  #Paso 5
  u_finalx= u_masx + (Ex*q/(2*m))
  u_finaly= u_masy + (Ey*q/(2*m))
  u_finalz= u_masz + (Ez*q/(2*m))
  
  return u_finalx,u_finaly,u_finalz

#Método de leapfrog para hallar las componentes de la siguiente posición

def leapfrog(posx,posy,posz,ux,uy,uz,dt):
  pos_finalx=posx+ux*dt
  pos_finaly=posy+uy*dt
  pos_finalz=posz+uz*dt
  return pos_finalx,pos_finaly,pos_finalz
  

def interpolation(posx,posy,posz):
    #x-direction
    ix =  int((posx-xmin)/dx) #index x
    xleft = posx    #position of the left node
    xright = xleft+dx  #position of the right node
    ileft  = ix    #index of the left node
    iright = ix+1  # #index of the right node
    hxleft  = (posx-xx[ix])/dx  #particle fractional x-distance from the nearest node
    #y-direction
    jy = int((posy-ymin)/dy) #index y
    ydown = posy #position of the down node
    yup = ydown+dy #position of the up node
    jdown  = jy  #index of the down node
    jup = jy+1   #index of up node
    hydown  = (posy-yy[jy])/dy; #particle fractional y-distance from the nearest node
    #print(ileft,iright,hxleft,jdown,jup,hydown)
    return ileft,iright,hxleft,jdown,jup,hydown
  
def get_fields_nodes(ex,ey,ez,bx,by,bz,ileft,iright,hxleft,jdown,jup,hydown):
    #weight functions
    w1 = (1.0-hxleft)*(1.0-hydown)
    w2 = hxleft*(1.0-hydown)
    w3 = (1.0-hxleft)*hydown
    w4 = hxleft*hydown
    # linear interpolations of fields on the particle position
    Ex = w1*ex[ileft,jdown] +  w2*ex[iright,jdown] + w3*ex[ileft,jup] + w4*ex[iright,jup]
    Ey = w1*ey[ileft,jdown] +  w2*ey[iright,jdown] + w3*ey[ileft,jup] + w4*ey[iright,jup]
    Ez = w1*ez[ileft,jdown] +  w2*ez[iright,jdown] + w3*ez[ileft,jup] + w4*ez[iright,jup]
    Bx = w1*bx[ileft,jdown] +  w2*bx[iright,jdown] + w3*bx[ileft,jup] + w4*bx[iright,jup]
    By = w1*by[ileft,jdown] +  w2*by[iright,jdown] + w3*by[ileft,jup] + w4*by[iright,jup]
    Bz = w1*bz[ileft,jdown] +  w2*bz[iright,jdown] + w3*bz[ileft,jup] + w4*bz[iright,jup]
    return Ex,Ey,Ez,Bx,By,Bz
    
def zigzag(xi,yi,zi,xf,yf,zf,dx,dy,dt):

    #xgc lo dejamos en standby

    ii=np.zeros(2)
    jj=np.zeros(2)
    Fx=np.zeros(2)
    Fy=np.zeros(2)
    Wx=np.zeros(2)
    Wy=np.zeros(2)
    J1=np.zeros([2,2])
    #JY1=np.zeros([2,2])
    J2=np.zeros([2,2])
    #JY2=np.zeros([2,2])
    


    ii[0]=math.floor(xi,dx)
    jj[0]=math.floor(yi,dy)


    ii[1]=math.floor(xf/dx)  
    jj[1]=math.floor(yf/dy)
        
        
    val1x=min(ii[0]*dx,ii[1]*dx)+dx
    val2x=max(max(ii[0]*dx,ii[1]*dx),(xi+xf)/2.0)
    xr=min(val1x,val2x)

    val1y=min(jj[0]*dy,jj[1]*dy)+dy
    val2y=max(max(jj[0]*dy,jj[1]*dy),(yi+yf)/2.0)
    yr=min(val1y,val2y)

                 
    Fx[0]=(xr-xi)/dt   
    Fy[0]=(yr-yi)/dt    

    Fx[1]=(xf-xr)/dt  
    Fy[1]=(yf-yr)/dt   
        

    Wx[0]=((xi+xr)/(2*dx))-ii[0]
    Wy[0]=((yi+yr)/(2*dy))-jj[0]
    

    Wx[1]=((xr+xf)/(2*dx))-ii[1]
    Wy[1]=((yr+yf)/(2*dy))-jj[1]

    J1[1,0]=(1/(dx*dy))*Fx[0]*(1-Wy[0])
    J1[1,2]=(1/(dx*dy))*Fx[0]*Wy[0]
    J1[0,1]=(1/(dx*dy))*Fy[0]*(1-Wx[0])
    J1[2,1]=(1/(dx*dy))*Fy[0]*Wx[0]
    
    J2[1,0]=(1/(dx*dy))*Fx[1]*(1-Wy[1])
    J2[1,2]=(1/(dx*dy))*Fx[1]*Wy[1]
    J2[0,1]=(1/(dx*dy))*Fy[1]*(1-Wx[1])
    J2[2,1]=(1/(dx*dy))*Fy[1]*Wx[1]

def densidad_corriente_zigzag(puntosDeMalla, pos_inicial, pos_final):
    #num_save = 3
    # Atemp_local
    # Atemp_total
    # ntotal = Nx*Ny*Nz*num_save
    # Atemp_local = ???
    # Atemp_total = ???
    dV = dx*dy

    #Atemp_local = np.zeros(ntotal)
    #Atemp_total = np.zeros(ntotal)

    puntosDeMalla = np.zeros([Nx, Ny]) #xgc
    