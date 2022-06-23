void setProperties()
{
    double datom;
    char a[50];

    ifstream infile;
    infile.open("parameters.in");
    infile>>a>>tmax;
    infile>>a>>neMax;
    infile>>a>>dt;
    infile>>a>>xmax;
    infile>>a>>xmin;
    infile>>a>>ymax;
    infile>>a>>ymin;
    infile>>a>>zmax;
    infile>>a>>zmin;
    infile>>a;
    infile>>a>>Nx;
    infile>>a>>Ny;
    infile>>a>>Nz;
    infile>>a>>Npart_cell;
    infile>>a;
    infile>>a>> mi_me;
    infile>>a>> vth_i;
    infile>>a>> vth_e;
    infile>>a;
    infile>>a>>solver_Maxwell;
    infile>>a;
    infile>>a>>solver_Boris;
    infile>>a;
    infile>>a>>solver;
    infile>>a>>gam;
    infile>>a>>bc_x;
    infile>>a>>bc_y;
    infile>>a>>SIR;
    infile>>a>>limiter;
    infile>>a>>variabletype;
    infile>>a>>limitertype;
    infile>>a>>epsil;    
    infile.close();

}

