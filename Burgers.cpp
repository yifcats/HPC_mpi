#include <chrono>
#include <cmath>
#include <mpi.h>
#include "Burgers.h"
#include "Model.h"

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

    Px = m->GetPx(); // domain of x
    Py = m->GetPy(); // domain of y

}

void burgers::Initial_velocity()
{

  
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // gets the current rank of a proccessor
    // my rank is an interger starts from 0
    
    // describes the size of the sub arrays mean!
    Nx_sub = Nx / Px;
    Ny_sub = Nx / Py;

    int rem_x = Nx % Px;
    int rem_y = Nx % Py;

    position_x = myrank % Px + 1; //s gives the location of the first matrix value globally
    position_y = myrank / Px + 1;

    int Position_global_x; // finding the actual x and y values in global frame
    int Position_global_y;

    if(rem_x == 0 && rem_y == 0) {
        Position_global_x = position_x * Nx_sub;
        Position_global_y = position_y * Ny_sub;
    } else {
        if(position_x <= rem_x) {
            Position_global_x = (position_x - 1) * Nx_sub + (position_x - 1);
        } else {
            Position_global_x = (position_x - 1) * Nx_sub + (rem_x); // to go to next line
        }

        if(position_y <= rem_y) {
            Position_global_y = (position_y - 1) * Ny_sub + (position_y - 1);
        } else {
            Position_global_y = (position_y - 1) * Ny_sub + (rem_y);
        }
    }
//
//    cout << myrank << endl;
//    cout << Position_global_x << endl;
//    cout << Position_global_y << endl;

    // getting correct dimensions of the sub matrix
    if(position_x <= rem_x) {
        Nx_sub++;
    }

    if(position_y <= rem_y) {
        Ny_sub++;
    }

    // Initilizing the arrays for the space and velocity
    //    x=new double[Nx];
    //    y=new double[Ny];
    //    u=new double[Nx*Ny];
    //    v=new double[Nx*Ny];

    x = new double[Nx_sub];
    y = new double[Ny_sub];

    u = new double[Nx_sub * Ny_sub];
    v = new double[Nx_sub * Ny_sub];

    // initilizing space matrices
    // x axis - sub arrays for mpi
    for(int temp = 0; temp < Nx_sub; temp++) {
        x[temp] = -Lx / 2 + (temp + Position_global_x) * dx;
    }
    // y axis
    for(int temp = 0; temp < Ny_sub; temp++) {
        y[temp] = -Ly / 2 + dy * (Position_global_y + temp);
    }

    double r;

    // Initial Conditions for velocity at t=0;
    for(int i = 0; i < Nx_sub; i++) {
        for(int j = 0; j < Ny_sub; j++) {
            r = sqrt((x[i] * x[i] + y[j] * y[j]));

            if(r <= 1) {
                u[i * Ny_sub + j] = 2.0 * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (4.0 * r + 1.0);
                v[i * Ny_sub + j] = 2.0 * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (4.0 * r + 1.0);

                // v[i*Ny+j] = u[i*Ny+j];
            } else {
                u[i * Ny_sub + j] = 0;
                v[i * Ny_sub + j] = 0;
            }
        }
    }

    // Boundary Conditions
    //    for(int temp = 0; temp < Ny_sub; temp++) {
    //        u[0 * Ny_sub + temp] = 0.0;
    //        u[(Nx_sub - 1) * Ny_sub + temp] = 0.0;
    //        v[0 * Ny_sub + temp] = 0.0;
    //        v[(Nx_sub - 1) * Ny_sub + temp] = 0.0;
    //    }
    //
    //    for(int temp = 0; temp < Nx_sub; temp++) {
    //        u[temp * Ny_sub + 0] = 0.0;
    //        u[temp * Ny_sub + (Ny_sub - 1)] = 0.0;
    //        v[temp * Ny_sub + 0] = 0.0;
    //        v[temp * Ny_sub + (Ny_sub - 1)] = 0.0;
    //    }

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
    
    cout << "Initialised velocity for: " << myrank << endl;
}

void burgers::Integrate_velocity()
{
    
    cout << "Running integration for: " << myrank << endl;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // gets the current rank of a proccessor
    // my rank is an interger starts from 0

    // describes the size of the sub arrays mean!

    int rem_x = Nx % Px;
    int rem_y = Nx % Py;

    int position_x = myrank % Px + 1; // gives the location of the first matrix value globally
    int position_y = myrank / Px + 1;

    int Position_global_x; // finding the actual x and y values in global frame
    int Position_global_y;

    if(rem_x == 0 && rem_y == 0) {
        Position_global_x = position_x * Nx_sub;
        Position_global_y = position_y * Ny_sub;
    } else {
        if(position_x <= rem_x) {
            Position_global_x = (position_x - 1) * Nx_sub + (position_x - 1);
        } else {
            Position_global_x = (position_x - 1) * Nx_sub + (rem_x); // to go to next line
        }

        if(position_y <= rem_y) {
            Position_global_y = (position_y - 1) * Ny_sub + (position_y - 1);
        } else {
            Position_global_y = (position_y - 1) * Ny_sub + (rem_y);
        }
    }

//    cout << myrank << endl;
//    cout << Position_global_x << endl;
//    cout << Position_global_y << endl;

    // getting correct dimensions of the sub matrix
    if(position_x <= rem_x) {
        Nx_sub++;
    }

    if(position_y <= rem_y) {
        Ny_sub++;
    }

    // redifined check later how to pass
    // Nx_sub = Nx / Px;
    // Ny_sub = Nx / Py;

    // cout<<"hindjknfdkndnfkdnfkdnkdnkdnkdnk"<<Nx_sub<<endl;
    // double* temp;

    // Calculating the new values of each space vector
    double u_n_plus;
    double v_n_plus;

    // Initilizing the matrices to renwew v and u
    auto* u_new = new double[Nx_sub * Ny_sub];
    auto* v_new = new double[Nx_sub * Ny_sub];

    const double C_i_minus_j = c / dx / dx + ax / dx;                                           // i-1,j
    const double C_i_plus_j = c / dx / dx;                                                      // i+1,j
    const double C_i_j = -2.0 * c * (1 / dx / dx + 1 / dy / dy) - ax / dx - ay / dy + 1.0 / dt; // i,j
    const double C_i_j_minus = ay / dy + c / dy / dy;                                           // i,j-1
    const double C_i_j_plus = c / dy / dy;                                                      // i,j+1

    //    cout<<C_i_minus_j<<endl;
    //    cout<<C_i_plus_j<<endl;
    //    cout<<C_i_j<<endl;
    //    cout<<C_i_j_minus<<endl;
    //    cout<<C_i_j_plus<<endl;

    //cout << "\tCurrent Energy: ";
    //Energy();
    //cout << endl;

    int Lower_i;
    int Upper_i;
    int Lower_j;
    int Upper_j;
    
    
    
    // Left top corner
    if(position_x == 1 && position_y == 1) {
        Lower_i = 1;
        Upper_i = Nx_sub;

        Lower_j = 1;
        Upper_j = Ny_sub;
    }

    // Left bottom corner
    else if(position_x == Px && position_y == 1) {
        Lower_i = 0;
        Upper_i = Nx_sub - 1;

        Lower_j = 0;
        Upper_j = Ny_sub - 1;

    }

    // Left bottom corner
    else if(position_x == 1 && position_y == Py) {
        Lower_i = 1;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub - 1;

    }

    // Left bottom corner
    else if(position_x == Px && position_y == Py) {
        Lower_i = 0;
        Upper_i = Nx_sub - 1;
        Lower_j = 0;
        Upper_j = Ny_sub - 1;
    }

    else if(position_x != 1 && position_x != Px && position_y == 1) {
        Lower_i = 0;
        Upper_i = Nx_sub;
        Lower_j = 1;
        Upper_j = Ny_sub;
    }

    else if(position_x != 1 && position_x != Px && position_y == Py) {
        Lower_i = 0;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub - 1;
    }

    else if(position_x == 1 && position_y != 1 && position_y != Py) {
        Lower_i = 1;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub;
    }

    else if(position_x == Px && position_y != 1 && position_y != Py) {
        Lower_i = 0;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub - 1;
    }

    else {
        Lower_i = 0;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub;
    }

    u_left = new double[Ny_sub];
    v_left = new double[Ny_sub];

    u_right = new double[Ny_sub];
    v_right = new double[Ny_sub];

    u_up = new double[Nx_sub];
    v_up = new double[Nx_sub];

    u_down = new double[Nx_sub];
    v_down = new double[Nx_sub];
    
    
    u_left_send = new double[Ny_sub];
    v_left_send = new double[Ny_sub];

    u_right_send = new double[Ny_sub];
    v_right_send = new double[Ny_sub];

    u_up_send = new double[Nx_sub];
    v_up_send = new double[Nx_sub];

    u_down_send = new double[Nx_sub];
    v_down_send = new double[Nx_sub];


    
    double u_i_minus_j;
    double u_i_plus_j;
      
    double u_i_j_minus;
    double u_i_j_plus;
    
    double v_i_minus_j;
    double v_i_plus_j;
      
    double v_i_j_minus;
    double v_i_j_plus;

    for(int t_couter = 1; t_couter < 3; t_couter++) {
        // MPI_Send(&Lower_i, 1, MPI_INT, &myrank,23, MPI_COMM_WOLRD);
                cout << "Creating boundaries for: " << myrank << endl;
                Create_boundaris();
                cout << "Created boundaries for: " << myrank << endl;
                MPI_Barrier(MPI_COMM_WORLD);
                //MPI_Bcast(&Lower_i, 1, MPI_INT, 0, MPI_COMM_WORLD);
                Communication();
                
                cout << "Communication for: " << myrank << endl;
                
                

        for(int i = Lower_i; i < Upper_i; i++) {
            for(int j = Lower_j; j < Upper_j; j++) {

               
                
                // setting u and v for the case with full boundaries
                
                
                
                // close boundary conditions
                
                
                if (i==0){
                    // LEFT
                    u_i_minus_j=u_left[j]; // singular value
                    v_i_minus_j=v_left[j];
                } else if (i==Nx_sub){
                    // RIGHT
                    u_i_plus_j=u_right[j];
                    v_i_plus_j=v_right[j];
                } else{
                    u_i_minus_j=u[(i-1)*Ny+j];
                    u_i_plus_j=u[(i+1)*Ny+j];
                
                    u_i_j_minus=u[i*Ny+(j-1)];
                    u_i_j_plus=u[i*Ny+(j+1)];
                
                
                    v_i_minus_j=v[(i-1)*Ny+j];
                    v_i_plus_j=v[(i+1)*Ny+j];
                
                    v_i_j_minus=v[i*Ny+(j-1)];
                    v_i_j_plus=v[i*Ny+(j+1)]; 
                }
                

                
                
                if (j==0){
                    // UP
                    u_i_j_minus=u_up[i];
                    v_i_j_minus=v_up[i];
                }  else if (j==Ny_sub){
                    // DOWN
                    u_i_j_plus=u_down[i];
                    v_i_j_plus=v_down[i];
                } else {
                    u_i_minus_j=u[(i-1)*Ny+j];
                    u_i_plus_j=u[(i+1)*Ny+j];
                
                    u_i_j_minus=u[i*Ny+(j-1)];
                    u_i_j_plus=u[i*Ny+(j+1)];
                
                
                    v_i_minus_j=v[(i-1)*Ny+j];
                    v_i_plus_j=v[(i+1)*Ny+j];
                
                    v_i_j_minus=v[i*Ny+(j-1)];
                    v_i_j_plus=v[i*Ny+(j+1)]; 
                    
                }
                
               
                
                
                // Liner Terms
                
                u_n_plus=C_i_minus_j*u_i_minus_j+
                         C_i_plus_j*u_i_plus_j+
                         C_i_j*u[i*Ny+j]+
                         C_i_j_minus*u_i_j_minus+
                         C_i_j_plus*u_i_j_plus;
                         
                         
                v_n_plus=C_i_minus_j*v_i_minus_j+
                         C_i_plus_j*v_i_plus_j+
                         C_i_j*v[i*Ny+j]+
                         C_i_j_minus*v_i_j_minus+
                         C_i_j_plus*v_i_j_plus;        
                
                
                
                // non liner terms
                
                u_n_plus+=b*u[i*Ny+j]*(u_i_minus_j-u[i*Ny+j])/dx
                                +b*v[i*Ny+j]*(u_i_j_minus-u[i*Ny+j])/dy;
                
                v_n_plus+=b*u[i*Ny+j]*(v_i_j_plus-v[i*Ny+j])/dx
                                +b*v[i*Ny+j]*(v_i_minus_j-v[i*Ny+j])/dy;
                
                
               
                //                // linear terms
//                                u_n_plus=C_i_minus_j*u[(i-1)*Ny+j]
//                                +C_i_plus_j*u[(i+1)*Ny+j]
//                                +C_i_j*u[i*Ny+j]
//                                +C_i_j_minus*u[i*Ny+(j-1)]
//                                +C_i_j_plus*u[i*Ny+(j+1)];
            
                //
                //                // non linear terms
//                                u_n_plus+=b*u[i*Ny+j]*(u[(i-1)*Ny+j]-u[i*Ny+j])/dx
//                                +b*v[i*Ny+j]*(u[i*Ny+(j-1)]-u[i*Ny+j])/dy;
                //
                //
//                                v_n_plus=C_i_minus_j*v[(i-1)*Ny+j]
//                                +C_i_plus_j*v[(i+1)*Ny+j]
//                                +C_i_j*v[i*Ny+j]
//                                +C_i_j_minus*v[i*Ny+(j-1)]
//                                +C_i_j_plus*v[i*Ny+(j+1)];
                
                
                
                //                // non linear terms
//                                v_n_plus+=b*u[i*Ny+j]*(v[i*Ny+(j-1)]-v[i*Ny+j])/dx
//                                +b*v[i*Ny+j]*(v[(i-1)*Ny+j]-v[i*Ny+j])/dy;
                //
                //
                            u_new[i*Ny_sub+j]=dt*u_n_plus;
                               v_new[i*Ny_sub+j]=dt*v_n_plus;
            }
        }


        for(int i = Lower_i; i < Upper_i; i++) {
            for(int j = Lower_j; j < Upper_j; j++) {
                //            temp=&u_new[i*Ny+j];
                //            u_new[i*Ny+j]=u[i*Ny+j];
                //            u[i*Ny+j]=*temp;
                //
                //
                //            temp=&v_new[i*Ny+j];
                //            v_new[i*Ny+j]=v[i*Ny+j];
                //            v[i*Ny+j]=*temp;

                u[i * Ny_sub + j] = u_new[i * Ny_sub + j];
                v[i * Ny_sub + j] = v_new[i * Ny_sub + j];

                // cout<<v[i*Ny+j]<<endl;
            }
        }

        //
        //    if (t_couter % 10 ==0){
        //
        //        cout<<"Time Iterations Passed: "<<t_couter<<"\tCurrent Energy: ";
        //        Energy();
        //        cout<<"\tVelocity[57]: "<<u[190*10+50]<<endl;
        //    }

        
            for(int i = 0; i < Nx_sub; i++) {
                for(int j = 0; j < Ny_sub; j++) {
                    cout.precision(4);
                    //cout << fixed << u[i * Ny_sub + j];
                }
                //cout << " ";
            }
        cout << "Time is: "<< t_couter << endl;
    }
}

void burgers::Create_boundaris(){

  
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // gets the current rank of a proccessor
    // my rank is an interger starts from 0

    // describes the size of the sub arrays mean!

    int rem_x = Nx % Px;
    int rem_y = Nx % Py;

    int Position_global_x; // finding the actual x and y values in global frame
    int Position_global_y;

    if(rem_x == 0 && rem_y == 0) {
        Position_global_x = position_x * Nx_sub;
        Position_global_y = position_y * Ny_sub;
    } else {
        if(position_x <= rem_x) {
            Position_global_x = (position_x - 1) * Nx_sub + (position_x - 1);
        } else {
            Position_global_x = (position_x - 1) * Nx_sub + (rem_x); // to go to next line
        }

        if(position_y <= rem_y) {
            Position_global_y = (position_y - 1) * Ny_sub + (position_y - 1);
        } else {
            Position_global_y = (position_y - 1) * Ny_sub + (rem_y);
        }
    }

    //cout << myrank << endl;
    //cout << Position_global_x << endl;
    //cout << Position_global_y << endl;

    // getting correct dimensions of the sub matrix
    if(position_x <= rem_x) {
        Nx_sub++;
    }

    if(position_y <= rem_y) {
        Ny_sub++;
    }

    if(position_x - 1 > 0) {
        rank_left = myrank - 1;
        
       // MPI_Recv(u_left, Ny_sub, MPI_DOUBLE, rank_left, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //MPI_Recv(v_left, Ny_sub, MPI_DOUBLE, rank_left, 2, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(int j = 0; j < Ny_sub; j++) {
//            u_left[j] = u[Nx_sub * Ny_sub + j];
//            v_left[j] = v[Nx_sub * Ny_sub + j];
//            
            
            u_left_send[j]=u[0*Ny_sub+j];
            v_left_send[j]=u[0*Ny_sub+j];
        }
    }

    // u_1=u[Nx_sub*Ny_sub+Ny_sub];

    if(position_x + 1 <= Px) {
        rank_right = myrank + 1;

        for(int j = 0; j < Ny_sub; j++) {
//            u_right[j] = u[0 * Ny_sub + j];
//            v_right[j] = v[0 * Ny_sub + j];
            
            u_right_send[j]=u[Nx_sub*Ny_sub+j];
            v_right_send[j]=u[Nx_sub*Ny_sub+j];
//            u_right_send[j]=(double)j;
//            v_right_send[j]=(double)j;
        }
    }

    if(position_y + 1 <= Py) {
        rank_down = myrank + Px;
        for(int i = 0; i < Ny_sub; i++) {
//            u_down[i] = u[i * Ny_sub + 0];
//            v_down[i] = v[i * Ny_sub + 0];
        
        
            u_down_send[i]=u[i*Ny_sub+Ny_sub];
            v_down_send[i]=u[i*Ny_sub+Ny_sub];
        }
    }

    if(position_y - 1 > 0) {
        rank_up = myrank - Px;
        for(int i = 0; i < Ny_sub; i++) {
//            u_up[i] = u[i * Ny_sub + Ny_sub];
//            v_up[i] = v[i * Ny_sub + Ny_sub];
//            
            u_up_send[i]=u[i*Ny_sub+0];
            v_up_send[i]=u[i*Ny_sub+0];
            
//            u_down[i]=u_up_send[i]/
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}


void burgers::Communication(){
                

                        
                if(position_x - 1 > 0) {
                    cout<<"my rank: "<<myrank<<" receiving left from "<<rank_left<<endl;
                MPI_Recv(u_left, Ny_sub, MPI_DOUBLE, rank_left, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(v_left, Ny_sub, MPI_DOUBLE, rank_left, 2, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                cout<<"my rank: "<<myrank<<" received left"<<endl;
                for (int j=0;j<Ny_sub;j++){
                    cout<<u_left[j]<<endl;
                }
                    cout<<"my rank: "<<myrank<<" sending left to "<<rank_left<<endl;
                MPI_Send(u_left_send, Ny_sub, MPI_DOUBLE, rank_left, 3, MPI_COMM_WORLD);
                MPI_Send(v_left_send, Ny_sub, MPI_DOUBLE, rank_left, 4, MPI_COMM_WORLD);
                cout<<"my rank: "<<myrank<<" sent left"<<endl;
                }
                
              
                    
                if(position_x +1<= Px) {
                    cout<<"my rank: "<<myrank<<" sending right to "<<rank_right<<endl;
                MPI_Send(u_right_send, Ny_sub, MPI_DOUBLE, rank_right, 1,MPI_COMM_WORLD);
                MPI_Send(v_right_send, Ny_sub, MPI_DOUBLE, rank_right, 2,MPI_COMM_WORLD);
                cout<<"my rank: "<<myrank<<" sent right"<<endl;
                    cout<<"my rank: "<<myrank<<" receiving right from "<<rank_right<<endl;
                
                MPI_Recv(u_right, Ny_sub, MPI_DOUBLE,rank_right, 3,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(v_right, Ny_sub, MPI_DOUBLE, rank_right,4,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                cout<<"my rank: "<<myrank<<" received right"<<endl;
        
                
                }
                cout<<"my rank: "<<myrank<<"Got here"<<endl;
                if(position_y + 1 <= Py) {
                MPI_Send(u_down_send, Nx_sub, MPI_DOUBLE, rank_down, 13,MPI_COMM_WORLD);
                MPI_Send(v_down_send, Nx_sub, MPI_DOUBLE, rank_down, 14,MPI_COMM_WORLD);
                
                MPI_Recv(u_down, Nx_sub, MPI_DOUBLE, rank_down,11,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(v_down, Nx_sub, MPI_DOUBLE,  rank_down,12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                
                }

                if(position_y - 1 > 0) {
                MPI_Recv(u_up, Nx_sub, MPI_DOUBLE,  rank_up, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(v_up, Nx_sub, MPI_DOUBLE,  rank_up, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                MPI_Send(u_up_send, Nx_sub, MPI_DOUBLE, rank_up, 11,MPI_COMM_WORLD);
                MPI_Send(v_up_send, Nx_sub, MPI_DOUBLE, rank_up, 12,MPI_COMM_WORLD);
                }
                
                //MPI_Barrier(MPI_COMM_WORLD);
}





void burgers::Energy()
{
    cout << "Running energy: "<<myrank<< endl;
    
    double Energy;
    cout << "Barrier " << endl;

    for(int i = 0; i < Nx_sub * Ny_sub; i++) {
    
        //MPI_Recv(&Energy,1,MPI_DOUBLE,myrank,50,MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receives energy
        
        Energy += u[i] * u[i] + v[i] * v[i];
    }
    MPI_Send(&Energy,1,MPI_DOUBLE, 0, 50, MPI_COMM_WORLD); // passes energy
    MPI_Barrier(MPI_COMM_WORLD);
    Energy *= 0.5 * dx * dy;

    cout << "Energy is: " << Energy << endl;;
    
    
}
