#include <chrono>
#include <cmath>


#include "Model.h"
#include "Burgers.h"

using namespace std;


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
    
    // Initilizing the arrays for the space and velocity
    x=new double[Nx];
    y=new double[Ny];
    u=new double[Nx*Ny];
    v=new double[Nx*Ny];
    
    
    // initilizing space matrices
    // x axis
    for(int temp = 0; temp<Nx; temp++) {
        x[temp] = -Lx/2 + temp*dx;
    }
    // y axis
    for(int temp = 0; temp<Ny; temp++) {
        y[temp] = -Ly/2 + dy*temp;
    }
    
    double r;
    
    // Initial Conditions for velocity at t=0;
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            r = sqrt((x[i]*x[i] + y[j]*y[j]));
            
            
            if(r <= 1) {
                u[i*Ny+j] = 2.0 *(1.0 - r)*(1.0 - r)*(1.0 - r)*(1.0 - r) * (4.0 * r + 1.0);
                v[i*Ny+j] = 2.0 *(1.0 - r)*(1.0 - r)*(1.0 - r)*(1.0 - r) * (4.0 * r + 1.0);
                
                //v[i*Ny+j] = u[i*Ny+j];
            } else {
                u[i*Ny+j] = 0;
                v[i*Ny+j] = 0;
            }
        }
    }
    
    // Boundary Conditions
    for(int temp = 0; temp < Ny; temp++) {
        u[0*Ny+temp] = 0.0;
        u[(Nx-1)*Ny+temp] = 0.0;
        v[0*Ny+temp] = 0.0;
        v[(Nx-1)*Ny+temp] = 0.0;
    }
    
    for(int temp = 0; temp < Nx; temp++) {
        u[temp*Ny+0] = 0.0;
        u[temp*Ny+(Ny-1)] = 0.0;
        v[temp*Ny+0] = 0.0;
        v[temp*Ny+(Ny-1)] = 0.0;
    }
    
    
    
    
    
    //            for(int temp = 0; temp < Nx; temp++) {
    //                for(int temp2=0; temp2 < Ny; temp2++){
    //                    cout << u[temp*Ny+temp2] << endl;
    //                }
    //            }
//    for (int i=0;i<Nx;i++){
//        for (int j=0;j<Ny;j++){
//            cout.precision(4);
//            cout<<fixed<<u[i*Ny+j];
//        }
//        cout<<" ";
//    }
//    
    
    
    //
}


void burgers::Integrate_velocity()
{
    
    
    
    // double* temp;
    
    
    // Calculating the new values of each space vector
    double u_n_plus;
    double v_n_plus;
    
    
    // Initilizing the matrices to renwew v and u
    auto* u_new = new double[Nx*Ny];
    auto* v_new = new double[Nx*Ny];
    
    
    const double C_i_minus_j = c/dx/dx + ax/dx;                               // i-1,j
    const double C_i_plus_j = c / dx / dx;                                          // i+1,j
    const double C_i_j  = -2.0*c*(1/dx/dx + 1/dy/dy) - ax/dx - ay/dy + 1.0/dt;        // i,j
    const double C_i_j_minus = ay/dy + c/dy/dy ;                               // i,j-1
    const double C_i_j_plus = c/dy/dy;                                          // i,j+1
    
    //    cout<<C_i_minus_j<<endl;
    //    cout<<C_i_plus_j<<endl;
    //    cout<<C_i_j<<endl;
    //    cout<<C_i_j_minus<<endl;
    //    cout<<C_i_j_plus<<endl;
    
    
    cout<<"\tCurrent Energy: ";
    Energy();
    cout<<endl;
    
    for (int t_couter=1;t_couter<Nt;t_couter++){
        for (int i=1;i<Nx-1;i++){
            for (int j=1;j<Ny-1;j++){
                
                
                
                // linear terms
                u_n_plus=C_i_minus_j*u[(i-1)*Ny+j]
                +C_i_plus_j*u[(i+1)*Ny+j]
                +C_i_j*u[i*Ny+j]
                +C_i_j_minus*u[i*Ny+(j-1)]
                +C_i_j_plus*u[i*Ny+(j+1)];
                
                
                // non linear terms
                u_n_plus+=b*u[i*Ny+j]*(u[(i-1)*Ny+j]-u[i*Ny+j])/dx
                +b*v[i*Ny+j]*(u[i*Ny+(j-1)]-u[i*Ny+j])/dy;
                
                
                v_n_plus=C_i_minus_j*v[(i-1)*Ny+j]
                +C_i_plus_j*v[(i+1)*Ny+j]
                +C_i_j*v[i*Ny+j]
                +C_i_j_minus*v[i*Ny+(j-1)]
                +C_i_j_plus*v[i*Ny+(j+1)];
                
                
                
                // non linear terms
                v_n_plus+=b*u[i*Ny+j]*(v[i*Ny+(j-1)]-v[i*Ny+j])/dx
                +b*v[i*Ny+j]*(v[(i-1)*Ny+j]-v[i*Ny+j])/dy;
                
                
                u_new[i*Ny+j]=dt*u_n_plus;
                v_new[i*Ny+j]=dt*v_n_plus;
                
                
                
                
            }
        }
        
        
        //    u=u_new;
        //    v=v_new;
        
        for (int i=1;i<Nx-1;i++){
            for (int j=1;j<Ny-1;j++){
                //            temp=&u_new[i*Ny+j];
                //            u_new[i*Ny+j]=u[i*Ny+j];
                //            u[i*Ny+j]=*temp;
                //
                //
                //            temp=&v_new[i*Ny+j];
                //            v_new[i*Ny+j]=v[i*Ny+j];
                //            v[i*Ny+j]=*temp;
                
                u[i*Ny+j]=u_new[i*Ny+j];
                v[i*Ny+j]=v_new[i*Ny+j];
                
                //cout<<v[i*Ny+j]<<endl;
                
            }
        }
        
        //
        //    if (t_couter % 10 ==0){
        //
        //        cout<<"Time Iterations Passed: "<<t_couter<<"\tCurrent Energy: ";
        //        Energy();
        //        cout<<"\tVelocity[57]: "<<u[190*10+50]<<endl;
        //    }
        
        if (t_couter==1){
        for (int i=0;i<Nx;i++){
            for (int j=0;j<Ny;j++){
                cout.precision(4);
                cout<<fixed<<u[i*Ny+j];
            }
            cout<<" ";
        }
        
        }
        
    }
    
    
}

void burgers::Energy(){
    double Energy=0.0;
    
    
    for (int i=0;i<Nx*Ny;i++){
        
        Energy+=u[i]*u[i]+v[i]*v[i];
        
    }
    Energy*=0.5*dx*dy;
    
    cout<<Energy;
}
