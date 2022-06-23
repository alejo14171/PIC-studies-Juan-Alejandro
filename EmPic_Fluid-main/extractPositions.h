//==========================================
//	CELLS
//==========================================
// For non-periodic or periodic boundary condition
//   - in x-direction
//   iflag = 0 : non-periodic (wall) 
//   iflag = 1 : periodic
//==========================================
void extractPositions_x_cell
(GridPoints ***gc,
 double xa, 
// double &hxright, 
 double &hxleft,
 int &iright,    int &ileft,
 int iflag)
{
       int ix; 
       double xleft, xright;
       if(iflag<0 || iflag>1) cout << "error in iflag extractPositions_x" <<endl;

       if(Nx!=1)
       {
          if(xa-xmin>=0.5*dx) 
                  ix = (int)((xa-xmin)/dx-0.5); //index
          else 
                  ix = -1;

          if(ix> Nx) cout << " cell: ix above Nx" << ix << " >" << Nx << endl;
          if(ix< -1) cout << " cell:ix below -1" << ix << " <0" << endl;

          if(ix<0)
          {
             xright = gc[0][0][0].x;
             xleft  = xright-dx;

             iright = 0;

             if(iflag==0)      ileft  = 0;
             else if(iflag==1) ileft  = Nx-1;
          }
          else if(ix>=Nx-1)
          {
             xleft = gc[Nx-1][0][0].x;
             xright  = xleft+dx;

             ileft = Nx-1;
             if(iflag==0)      iright  = Nx-1;
             else if(iflag==1) iright  = 0;
          }
          else
          {
             xleft = gc[ix][0][0].x;
             xright = xleft+dx;

             ileft  = ix;
             iright = ix+1;
          }

          hxleft  = (xa-xleft)/dx;
//          hxright = (xright-xa)/dx;
       }
       else
       {
          iright = 0;
          ileft  = 0;
          hxleft  = 1.0;
//          hxright = 0.0;
       }
}

//==========================================
// For periodic boundary condition
//   - in y-direction
//==========================================
void extractPositions_y_cell
(GridPoints ***gc,
  double ya, 
//  double &hyright,
  double &hyleft,
  int &iright,    int &ileft,
 int iflag)
{
       int iy; 
       double yleft, yright;

       if(Ny!=1)
       {
          if(ya-ymin>=0.5*dy) 
                  iy = (int)((ya-ymin)/dy-0.5); //indey
          else 
                  iy = -1;

          if(iy<0)
          {
             yright = gc[0][0][0].y;
             yleft  = yright-dy;

             iright = 0;
             ileft  = Ny-1;
          }
          else if(iy>=Ny-1)
          {
             yleft = gc[0][Ny-1][0].y;
             yright  = yleft+dy;

             ileft  = Ny-1;
             iright = 0;
          }
          else
          {
             yleft = gc[0][iy][0].y;
             yright = yleft+dy;

             ileft  = iy;
             iright = iy+1;
          }

          hyleft  = (ya-yleft)/dy;
          //hyright = (yright-ya)/dy;

	  // should not happen
          if(hyleft<0.0 || hyleft>1.0) cout << " cel shout hyleft " << hyleft <<endl;  
//          if(hyright<0.0 || hyright>1.0) cout << " cel shout hyright " << hyright <<endl;  
       }
       else // if only 1 cell
       {
          iright = 0;
          ileft  = 0;
          hyleft  = 1.0;
//          hyright = 0.0;
       }
}


//==========================================
//	NODES
//==========================================
// periodic boundary condition in x-direction
//==========================================
void extractPositions_x_node
(FieldPoints ***field,
 double xa, 
 double &hxleft,
 int &iright,    int &ileft)
{
       int ix; 
       double xleft, xright;

       // If information in X
       if(Nx!=1)
       {
          if(xa-xmin>=0.0)
                  ix = (int)((xa-xmin)/dx); //index
          else 
                  ix = -1;

          // Sometimes xa = xmax gives ix = Nx
          if(ix>= Nx)
	  {
             cout << " node: ix above Nx" << ix << ">-" << Nx << " xa:" << xa <<  endl;
	     ix = Nx-1;
	  }

          // this should not happen but just incase
          if(ix< 0) 
          {
             cout << " node: ix below 0"    << ix << "<0 xa:" << xa   <<endl;
             ix = 0;
          }


          // 0<=ix<Nx
          xleft = field[ix][0][0].x;
          xright = xleft+dx;

          ileft  = ix;   // >=0
          iright = ix+1; // <=Nx

          hxleft  = (xa-xleft)/dx;

          // sometimes roundoff errors 
          if(hxleft>1.0)// || hxright<0.0)
          {
              hxleft  = 1.0;
          }

          if(hxleft<0.0)// || hxright>1.0)
          {
              hxleft  = 0.0;
          }
       }
       else // if only 1 cell
       {
          iright = 0;
          ileft  = 0;
          hxleft  = 1.0;
       }
}

//==========================================
// For periodic boundary condition
//   - in y-direction
//==========================================
void extractPositions_y_node
(FieldPoints ***field,
  double ya, 
  double &hyleft,
  int &iright,    int &ileft)
{
       int iy; 
       double yleft, yright;

       if(Ny!=1)
       {
          if(ya-ymin>=0.0) 
                  iy = (int)((ya-ymin)/dy); //indey
          else 
                  iy = -1;

          // Sometimes ya = ymax gives iy = Ny
          if(iy>= Ny)
	  {
             cout << " node: iy above Ny" << iy << ">-" << Ny << " ya:" << ya <<  endl;
	     iy = Ny-1;
	  }

          // this should not happen but just incase
          if(iy< 0) 
          {
             cout << " node: iy below 0"    << iy << "<0 ya:" << ya   <<endl;
             iy = 0;
          }

          // 0<=iy<Ny
          yleft = field[0][iy][0].y;
          yright = yleft+dy;

          ileft  = iy;   // >= 0 
          iright = iy+1; // <= Ny

          hyleft  = (ya-yleft)/dy;

          // sometimes roundoff errors 
          if(hyleft>1.0)// || hyright<0.0)
          {
              hyleft  = 1.0;
          }

          if(hyleft<0.0)// || hyright>1.0)
          {
              hyleft  = 0.0;
          }
       }
       else // if only 1 cell
       {
          iright = 0;
          ileft  = 0;
          hyleft  = 1.0;
       }
}


