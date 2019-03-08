#ifndef BURGERS_H
#define BURGERS_H



#define F77NAME(x) x##_
extern "C" {
    /* Level 1 functions */
    double F77NAME(ddot)(
                         const int& n,
                         const double *x, const int& incx,
                         const double *y, const int& incy
                         );
    
    double F77NAME(daxpy)(
                          const int& n, const double& alpha,
                          const double *x,const int& incx,
                          const double *y, const int& incy
                          );
};



class burgers
{
public:
    burgers(Model* m);
    //   ~burgers();
    
    // Declaring functions
    void Initial_velocity();
    void Integrate_velocity();
    void Energy();
    
private:
    Model* m;
    // Numerics
    double*x;
    double*y;
    
    
    double*u;
    double*v;
    
    
    double x0;
    double y0;
    double Lx;
    double Ly;
    double T;
    int  Nx;
    int  Ny;
    int Nt;
    double dx;
    double dy;
    double dt;
    
    //    //double dt = T / Nt; // taken out of initilizing velocity
    //    double uArray[Nx*Ny];  // u velocity (space x and y)
    //    //double unArray[Nx*Ny]; // previous velocity u
    //
    //    double vArray[Nx*Ny]; // v velocity (space x and y)
    //    //double vnArray[Nx*Ny]; // previous velocity v
    //
    //    double rhsxArray[Nx*Ny];
    //    double rhsyArray[Nx*Ny];
    //    double xArray[Nx*Ny];
    
    
    // Physics
    double ax;
    double ay;
    double b;
    double c;
    
    //
    
};



burgers::burgers(Model* m_)
{
    m = m_;
    Nx = m->GetNx();
    Ny = m->GetNy();
    
    
    
    ax = m->GetAx();
    ay = m->GetAy();
    b = m->GetB();
    c = m->GetC();
    Nt = m->GetNt();
    
    dx = m->GetDx();
    dy = m->GetDy();
    dt = m->GetDt();
    
    T = m->GetT();
    x0 = m->GetX0();
    y0 = m->GetY0();
    
    Lx = m->GetLx(); // domain of x
    Ly = m->GetLy(); // domain of y

    
    
    
}


void burgers::Initial_velocity()
{
    // Initialisation
    
    double L=Lx; // setting the value of L
    
    
    
    
    // defining step size in x and y.
    double dx = L / (Nx -1);
    double dy = L / (Ny -1);
    
    
    
    // Initilizing the arrays for the space and velocity
    x=new double[Nx];
    y=new double[Ny];
    u=new double[Nx*Ny];
    v=new double[Nx*Ny];
    
    //x[0] = -L / 2;
    //y[0] = -L / 2;
    
    
    // initilizing space matrices
    // x axis
    for(int temp = 0; temp<Nx; temp++) {
        x[temp] = -L/2 + temp*dx;
    }
    // y axis
    for(int temp = 0; temp<Ny; temp++) {
        y[temp] = -L/2 + dy*temp;
    }
    
    
    // Initial Conditions for velocity at t=0;
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            double r = pow(pow(x[i], 2) + pow(y[j], 2), 1/2);
            if(r <= 1) {
                u[i*Ny+j] = 2 * pow(1 - r, 4) * (4 * r + 1);
                v[i*Ny+j] = u[i*Ny+j];
            } else {
                u[i*Ny+j] = 0;
                v[i*Ny+j] = 0;
            }
        }
    }
    
    // Boundary Conditions
    for(int temp = 0; temp < Ny; temp++) {
        u[0*Ny+temp] = 0.0;
        u[Nx*Ny+temp] = 0.0;
        v[0*Ny+temp] = 0;
        v[Nx*Ny+temp] = 0;
    }
    
    for(int temp = 0; temp < Nx; temp++) {
        u[temp*Ny+1] = 0.0;
        u[temp*Ny+Ny] = 0.0;
        v[temp*Ny+1] = 0;
        v[temp*Ny+Ny] = 0;
    }
    
    
    
    //        for(int temp = 0; temp < Nx; temp++) {
    //            for(int temp2=0; temp2 < Ny; temp2++){
    //                cout << u[temp][temp2] << endl;
    //            }
    //        }
    
}


void burgers::Integrate_velocity()
{
    
    double L=Lx; // setting the value of L
    
    // defining step size in x and y.
    double dx = L / (Nx -1);
    double dy = L / (Ny -1);
    double dt=T/Nt;
    
    
    // Initilizing the matrices to renwew v and u
    //double un[Nx*Ny];
    //double vn[Nx*Ny];
    
    double u_n_plus;
    double v_n_plus;

    
    
    const double C_i_minus_j = c / dx / dx + ax / dx;                               // i-1,j
    const double C_i_plus_j = c / dx / dx;                                          // i+1,j
    const double C_i_j  = -2.0*c*(1/dx/dx + 1/dy/dy) - ax/dx - ay/dy + 1/dt;        // i,j
    const double C_i_j_minus = c / dy / dy + ay / dy;                               // i,j-1
    const double C_i_j_plus = c / dy / dy;                                          // i,j+1
    
    
    
    for (int i=1;i<Nx-1;i++){
        for (int j=1;j<Ny-1;j++){
            for (int t_couter=1;t_couter<Nt;t_couter++){
                
                // linear terms
                u_n_plus=C_i_minus_j*u[(i-1)*Ny+j]+C_i_plus_j*u[(i+1)*Ny+j]+C_i_j*u[i*Ny+j]+C_i_j_minus*u[i*Ny+(j-1)]+C_i_j_plus*u[i*Ny+(j+1)];
                
                // non linear terms
                u_n_plus+=b*u[i*Ny+j]*(-u[i*Ny+j]+u[(i-1)*Ny+j])/dx+b*v[i*Ny+j]*(u[i*Ny+(j-1)]-u[i*Ny+j])/dy;
                
                
                v_n_plus=C_i_minus_j*v[(i-1)*Ny+j]+C_i_plus_j*v[(i+1)*Ny+j]+C_i_j*v[i*Ny+j]+C_i_j_minus*v[i*Ny+(j-1)]+C_i_j_plus*v[i*Ny+(j+1)];
                // non linear terms
                v_n_plus+=b*v[i*Ny+j]*(-v[i*Ny+j]+v[(i-1)*Ny+j])/dx+b*u[i*Ny+j]*(v[i*Ny+(j-1)]-v[i*Ny+j])/dy;
                
                u[i*Ny+j]+=dt*u_n_plus;
                v[i*Ny+j]+=dt*v_n_plus;
        
        
            }
        }
    }
    
}

void burgers::Energy(){
    double Energy;
    
    
    

  F77NAME(ddot)(Nx*Ny,u,1,u,1);
  F77NAME(ddot)(Nx*Ny,v,1,v,1);
    
    F77NAME(daxpy)(Nx*Ny,1.0,u,1,v,1);
    
    cout<<*v<<endl;
    //Energy=dx*dy*v/2;
    
}


#endif // BURGERS_H
