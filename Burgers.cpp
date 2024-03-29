#include "Burgers.h"
#include "Model.h"
#include <chrono>
#include <cmath>
#include <mpi.h>

#include <fstream>
#include <iostream>

//#define F77NAME(x) x##_
//extern "C" {
//// Copies a vector to another vector (double-precision). Y = X
//double F77NAME(dcopy)(const int& n, const double* x, const int& incx, const double* y, const int& incy);
//}

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
    


    rank_error_check = MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // zero-based

    if(rank_error_check == MPI_ERR_COMM) {
        cout << "Invalid communicator" << endl;
    }
    
}

void burgers::Initial_velocity()
{
    // cout<<"go"<<endl;

    //MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // gets the current rank of a proccessor
    // my rank is an interger starts from 0

    // describes the size of the sub arrays mean!
    Nx_sub = Nx / Px;
    Ny_sub = Nx / Py;

    rem_x = Nx % Px;
    rem_y = Nx % Py;

    position_x = myrank % Px; // s gives the location of the first matrix value globally
    position_y = myrank / Px;

    // finding the actual x and y values in global frame
    if(rem_x == 0 && rem_y == 0) {
        Position_global_x = (position_x) * Nx_sub;
        Position_global_y = (position_y) * Ny_sub;
    } else {
        if(position_x < rem_x) {
            Position_global_x = (position_x) * Nx_sub + (position_x );
        } else {
            Position_global_x = (position_x) * Nx_sub + (rem_x); // to go to next line
        }

        if(position_y < rem_y) {
            Position_global_y = (position_y) * Ny_sub + (position_y );
        } else {
            Position_global_y = (position_y) * Ny_sub + (rem_y);
        }
    }

    // cout<<"global Position"<<Position_global_x<<endl;

    // getting correct dimensions of the sub matrix which is dependend on the 
    // remineder for a grid which is not divisable by px and py. 
    // the reminder is added from the first prossor to until the if condition is unsatisfied
    if(position_x < rem_x) {
        Nx_sub++;
    }

    if(position_y < rem_y) {
        Ny_sub++;
    }


    // Initilizing the size of the arrays
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
        y[temp] = Ly / 2 - dy * (Position_global_y + temp);
    }


    double r;

    int i_j;
    // Initial Conditions for velocity at t=0;
    for(int i = 0; i < Nx_sub; i++) {
        for(int j = 0; j < Ny_sub; j++) {

            i_j = i * Ny_sub + j;
            r = sqrt((x[i] * x[i] + y[j] * y[j]));

            if(r <= 1) {
                u[i_j] = 2.0 * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (4.0 * r + 1.0);
                v[i_j] = u[i_j];
            } else {
                u[i_j] = 0.0;
                v[i_j] = 0.0;
            }
        }
    }
    
    
    // creating matrices to be used in the integrations phase
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
    
    
    
    
    
    delete[] x; 
    delete[] y;

}

void burgers::Integrate_velocity()
{   
    

    // cout << "Running integration for: " << myrank << endl;
    //MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // gets the current rank of a proccessor
    // my rank is an interger starts from 0

    // Calculating the new values of each space vector
    double u_n_plus = 0.0;
    double v_n_plus = 0.0;

    // Initilizing the matrices to renwew v and u
    auto* u_new = new double[Nx_sub * Ny_sub];
    auto* v_new = new double[Nx_sub * Ny_sub];

    const double C_i_minus_j = dt * (c / dx / dx + ax / dx);                                           // i-1,j
    const double C_i_plus_j = dt * (c / dx / dx);                                                      // i+1,j
    const double C_i_j = dt * (-2.0 * c * (1 / dx / dx + 1 / dy / dy) - ax / dx - ay / dy + 1.0 / dt); // i,j
    const double C_i_j_plus = dt * (ay / dy + c / dy / dy);                                            // i,j+1
    const double C_i_j_minus = dt * (c / dy / dy);                                                     // i,j-1

   

    
    //

    
    
    
    Set_integration_boundaries();

    for(int t_couter = 1; t_couter <= Nt; t_couter++) {

        // creating boundaries for the sub matrices
        Create_sending_boundaris();
        // swapping processor own boundaries with surrounding boundaries
        Communication();

        for(int i = Lower_i; i < Upper_i; i++) {
            for(int j = Lower_j; j < Upper_j; j++) {

                // Indicies
                i_j = i * Ny_sub + j; // i,j

                i_plus_j = i_j + Ny_sub;  // i+1,j
                i_minus_j = i_j - Ny_sub; // i-1,j

                i_j_plus = i_j + 1;  // i,j+1
                i_j_minus = i_j - 1; // i,j-1

                Current_velocity(i,j);

                // Liner Terms

                u_n_plus = C_i_minus_j * u_i_minus_j + C_i_plus_j * u_i_plus_j + C_i_j * u[i_j] +
                    C_i_j_minus * u_i_j_minus + C_i_j_plus * u_i_j_plus;

                v_n_plus = C_i_minus_j * v_i_minus_j + C_i_plus_j * v_i_plus_j + C_i_j * v[i_j] +
                    C_i_j_minus * v_i_j_minus + C_i_j_plus * v_i_j_plus;

                // non liner terms

                u_n_plus += dt * (b * u[i_j] * (u_i_minus_j - u[i_j]) / dx + b * v[i_j] * (u_i_j_plus - u[i_j]) / dy);

                v_n_plus += dt * (b * u[i_j] * (v_i_j_plus - v[i_j]) / dx + b * v[i_j] * (v_i_minus_j - v[i_j]) / dy);

                // Filling in a matrix of the new values of u for the new time step
                u_new[i_j] = u_n_plus;
                v_new[i_j] = v_n_plus;
            }
        }

        // MPI_Barrier(MPI_COMM_WORLD);
        for(int i = Lower_i; i < Upper_i; i++) {
            for(int j = Lower_j; j < Upper_j; j++) {

                // new indicies counter
                i_j_2 = i * Ny_sub + j;

                // over writing the values of u for the n+1 time step
                u[i_j_2] = u_new[i_j_2];
                v[i_j_2] = v_new[i_j_2];

                // F77NAME(dcopy)(Nx_sub*Ny_sub,u_new, 1,u, 1);
            }
        }

        if(t_couter % 10 == 0 && myrank == 0) {
            cout << "Time Iterations Passed: " << t_couter << endl;
        }
    }

    
    Assembling_Submatices();
    Print_matrix();
    
    
    delete[] u_new;
    delete[] v_new;
    
    
    
    delete[] u_left;
    delete[] v_left;

    delete[] u_right;
    delete[] v_right;

    delete[] u_up; 
    delete[] v_up;

    delete[] u_down;
    delete[] v_down;

    delete[] u_left_send; 
    delete[] v_left_send;

    delete[] u_right_send;
    delete[] v_right_send;

    delete[] u_up_send; 
    delete[] v_up_send; 

    delete[] u_down_send;
    delete[] v_down_send;
    
    //delete[] U_assembled;
    
    
}



void burgers::Set_integration_boundaries(){
    // This function outputs the lower and upper bounds of i and j to be used in integration for loop
    // It makes sure that external boundaries are uneffected to make sure that boundary condition given 
    // is satisfied. And to insure that there will be no segmentation errors.
    
    
    // Left top corner
    if(position_x == 0 && position_y == 0) {
        Lower_i = 1;
        Upper_i = Nx_sub;

        Lower_j = 1;
        Upper_j = Ny_sub;
    }

    // Right TOP corner
    else if(position_x == Px-1 && position_y == 0) {
        Lower_i = 0;
        Upper_i = Nx_sub - 1;

        Lower_j = 1;
        Upper_j = Ny_sub;

    }

    // Left bottom corner
    else if(position_x == 0 && position_y == Py-1) {
        Lower_i = 1;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub - 1;

    }

    // Right bottom corner
    else if(position_x == Px-1 && position_y == Py-1) {
        Lower_i = 0;
        Upper_i = Nx_sub - 1;
        Lower_j = 0;
        Upper_j = Ny_sub - 1;
    }

    // Top middle (anything between first and last processor
    else if(position_x != 0 && position_x != Px-1 && position_y == 0) {
        Lower_i = 0;
        Upper_i = Nx_sub;
        Lower_j = 1;
        Upper_j = Ny_sub;
    }
    // bottom middle (anything between first and last processor)
    else if(position_x != 0 && position_x != Px-1 && position_y == Py-1) {
        Lower_i = 0;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub - 1;
    }
    // Left middle (anything between first and last processor)
    else if(position_x == 0 && position_y != 0 && position_y != Py-1) {
        Lower_i = 1;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub;
    }
    // right middle (anything between first and last processor)
    else if(position_x == Px-1 && position_y != 0 && position_y != Py-1) {
        Lower_i = 0;
        Upper_i = Nx_sub - 1;
        Lower_j = 0;
        Upper_j = Ny_sub;
    }

    else {
        Lower_i = 0;
        Upper_i = Nx_sub;
        Lower_j = 0;
        Upper_j = Ny_sub;
    }
    
    
}



void burgers::Current_velocity(int i, int j)
{
    // This function reduces effiecny but increases legability which was felt nessessary

    // This function recives the velocities form surrounding processors and re defines its values accordingly
    // this function checks for the possible combinations of where the processor is situated.

    if(i == 0 || i == Nx_sub - 1 || j == 0 || j == Ny_sub - 1) {
        if(i == 0 && j == 0) {
            // LEFT i-1 incicies
            u_i_minus_j = u_left[j];
            v_i_minus_j = v_left[j];

            // UP j-1
            u_i_j_minus = u_up[i];
            v_i_j_minus = v_up[i];

            u_i_plus_j = u[i_plus_j];
            v_i_plus_j = v[i_plus_j];

            u_i_j_plus = u[i_j_plus];
            v_i_j_plus = v[i_j_plus];

        } else if(i == Nx_sub - 1 && j == Ny_sub - 1) {
            // RIGHT i+1 indices
            u_i_plus_j = u_right[j];
            v_i_plus_j = v_right[j];

            // DOWN
            u_i_j_plus = u_down[i];
            v_i_j_plus = v_down[i];

            u_i_minus_j = u[i_minus_j];
            v_i_minus_j = v[i_minus_j];

            u_i_j_minus = u[i_j_minus];
            v_i_j_minus = u[i_j_minus];
        } else if(i == Nx_sub - 1 && j == 0) {
            // RIGHT i+1 indices
            u_i_plus_j = u_right[j];
            v_i_plus_j = v_right[j];

            // UP j-1
            u_i_j_minus = u_up[i];
            v_i_j_minus = v_up[i];

            u_i_minus_j = u[i_minus_j];
            v_i_minus_j = v[i_minus_j];

            u_i_j_plus = u[i_j_plus];
            v_i_j_plus = v[i_j_plus];

        } else if(i == 0 && j == Ny_sub - 1) {
            // DOWN
            u_i_j_plus = u_down[i];
            v_i_j_plus = v_down[i];

            // LEFT i-1 incicies
            u_i_minus_j = u_left[j];
            v_i_minus_j = v_left[j];

            u_i_plus_j = u[i_plus_j];
            v_i_plus_j = v[i_plus_j];

            u_i_j_minus = u[i_j_minus];
            v_i_j_minus = u[i_j_minus];
        } else if(i == 0 && j != 0 && j != Ny_sub - 1) {

            // LEFT i-1 incicies
            u_i_minus_j = u_left[j];
            v_i_minus_j = v_left[j];

            u_i_j_plus = u[i_j_plus];
            v_i_j_plus = v[i_j_plus];

            u_i_plus_j = u[i_plus_j];
            v_i_plus_j = v[i_plus_j];

            u_i_j_minus = u[i_j_minus];
            v_i_j_minus = v[i_j_minus];

        } else if(j == 0 && i != 0 && i != Nx_sub - 1) {

            // UP j-1
            u_i_j_minus = u_up[i];
            v_i_j_minus = v_up[i];

            v_i_j_plus = v[i_j_plus];
            u_i_j_plus = u[i_j_plus];

            u_i_minus_j = u[i_minus_j];
            v_i_minus_j = v[i_minus_j];

            u_i_plus_j = u[i_plus_j];
            v_i_plus_j = v[i_plus_j];

        } else if(j == Ny_sub - 1 && i != 0 && i != Ny_sub - 1) {
            // DOWN
            u_i_j_plus = u_down[i];
            v_i_j_plus = v_down[i];

            u_i_minus_j = u[i_minus_j];
            v_i_minus_j = v[i_minus_j];

            u_i_plus_j = u[i_plus_j];
            v_i_plus_j = v[i_plus_j];

            u_i_j_minus = u[i_j_minus];
            v_i_j_minus = v[i_j_minus];

        } else if(i == Nx_sub - 1 && j != 0 && j != Ny_sub) {
            // RIGHT i+1 indices
            u_i_plus_j = u_right[j];
            v_i_plus_j = v_right[j];

            v_i_j_plus = v[i_j_plus];
            u_i_j_plus = u[i_j_plus];

            u_i_minus_j = u[i_minus_j];
            v_i_minus_j = v[i_minus_j];

            u_i_j_minus = u[i_j_minus];
            v_i_j_minus = v[i_j_minus];
        }

    } else {

        v_i_j_plus = v[i_j_plus];
        u_i_j_plus = u[i_j_plus];

        u_i_minus_j = u[i_minus_j];
        v_i_minus_j = v[i_minus_j];

        u_i_plus_j = u[i_plus_j];
        v_i_plus_j = v[i_plus_j];

        u_i_j_minus = u[i_j_minus];
        v_i_j_minus = v[i_j_minus];
    }
}






// Assembling By rows
void burgers::Assembling_Submatices(){
    
    if (myrank==0){
    
    array_global_y=new int[Px*Py];
    array_global_x=new int[Px*Py];
           
    array_size_Nx_sub=new int[Px*Py];
    array_size_Ny_sub=new int[Px*Py];
        
        
    for (int rank_now = 1; rank_now < Px*Py; rank_now++){
    MPI_Recv(&array_global_x[rank_now],Px,MPI_DOUBLE,rank_now,888,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&array_size_Nx_sub[rank_now],Px,MPI_DOUBLE,rank_now,777,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    //cout<<"rank now: "<<rank_now<<endl;
    }
    
    array_global_x[0]=0;
    array_size_Nx_sub[0]=Nx_sub;
    
    for (int rank_now = 1; rank_now < Px*Py; rank_now++){
    MPI_Recv(&array_global_y[rank_now],Py,MPI_DOUBLE,rank_now,666,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&array_size_Ny_sub[rank_now],Py,MPI_DOUBLE,rank_now,555,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    array_global_y[0]=0;
    array_size_Ny_sub[0]=Ny_sub;
    
    
    
    }else{
     
        MPI_Send(&Position_global_x,1,MPI_DOUBLE,0,888,MPI_COMM_WORLD);   
        MPI_Send(&Nx_sub,1,MPI_DOUBLE,0,777,MPI_COMM_WORLD); 

        MPI_Send(&Position_global_y,1,MPI_DOUBLE,0,666,MPI_COMM_WORLD);   
        MPI_Send(&Ny_sub,1,MPI_DOUBLE,0,555,MPI_COMM_WORLD); 
        
        //cout<<"myrank "<<myrank<<endl;
 
    }
    
    
    if (myrank==0){
        U_assembled=new double[Nx*Ny];
        V_assembled=new double[Nx*Ny];

    for(int rank_now = 1; rank_now < Px * Py; rank_now++) {
        for (int i=0;i<array_size_Nx_sub[rank_now];i++){
            
        MPI_Recv(&U_assembled[(array_global_x[rank_now]+i)*Ny+array_global_y[rank_now]],array_size_Ny_sub[rank_now],MPI_DOUBLE,rank_now,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        MPI_Recv(&V_assembled[(array_global_x[rank_now]+i)*Ny+array_global_y[rank_now]],array_size_Ny_sub[rank_now],MPI_DOUBLE,rank_now,1999,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        }
    }
    
    
    for (int i=0;i<Nx_sub;i++){
        for (int j=0;j<Ny_sub;j++){
        U_assembled[(Position_global_x+i)*Ny+Position_global_y+j]=u[i*Ny_sub+j];
        
        }
        
    }
    
    
        
    }else {
        
        //double* Assemble_send=new double[Nx_sub]
        for (int i=0;i<Nx_sub;i++){
        //Assemble_send[i]=u[]
        MPI_Send(&u[i*Ny_sub],Ny_sub,MPI_DOUBLE,0,999,MPI_COMM_WORLD);
        MPI_Send(&v[i*Ny_sub],Ny_sub,MPI_DOUBLE,0,1999,MPI_COMM_WORLD);
        }
        
    }
    
    
    
    
    delete[] array_global_y;
    delete[] array_global_x;

           
    delete[] array_size_Nx_sub;
    delete[] array_size_Ny_sub;
 
    
}














void burgers::Print_matrix()
{
    

    //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    if (myrank==0){
    // Open file
    ofstream f_out("Velocity_Feild_Diff.txt");
    
    if(!f_out.good()) { // check if the file was open seccesfully
        cout << "Error: unable to open output file: sine.txt" << endl;
    } else {
      
        f_out.precision(5);
        for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            

                f_out <<fixed<<U_assembled[i * Ny + j] << " ";
                
            }
            f_out << endl;
            
        }

        f_out.precision(5);
        for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            

                f_out <<fixed<<V_assembled[i * Ny + j] << " ";
                
            }
            f_out << endl;
            
        }

    }

    // Close file
    f_out.close();

}

}









void burgers::Create_sending_boundaris(){

    // The following function computes the boundaries of the current processor and saves this arrays 
    // to be sent using the next function 
    // the amount of boundaries nessessary depend on the processors surroundings 
    //MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // gets the current rank of a proccessor

    for(int j = 0; j < Ny_sub; j++) {

        u_left_send[j] = u[0 * Ny_sub + j];
        v_left_send[j] = v[0 * Ny_sub + j];

    }
    //}

    // u_1=u[Nx_sub*Ny_sub+Ny_sub];

    // if(position_x + 1 <= Px) {
    // rank_right = myrank + 1;
    for(int j = 0; j < Ny_sub; j++) {
        //            u_right[j] = u[0 * Ny_sub + j];
        //            v_right[j] = v[0 * Ny_sub + j];

        u_right_send[j] = u[(Nx_sub - 1) * Ny_sub + j]; // Change here -1
        v_right_send[j] = v[(Nx_sub - 1) * Ny_sub + j]; // Change here -1
                                                        //            u_right_send[j]=(double)j;
                                                        //            v_right_send[j]=(double)j;
        // u[(Nx_sub-1)*Ny_sub+j]=u_right[j];
        // v[(Nx_sub-1)*Ny_sub+j]=v_right[j];
    }
    //}

    // if(position_y + 1 <= Py) {
    // rank_down = myrank + Px;
    for(int i = 0; i < Nx_sub; i++) { // Change here
                                      //            u_down[i] = u[i * Ny_sub + 0];
                                      //            v_down[i] = v[i * Ny_sub + 0];

        u_down_send[i] = u[i * Ny_sub + Ny_sub - 1]; // Change here -1
        v_down_send[i] = v[i * Ny_sub + Ny_sub - 1]; // Change here -1

        // u[i*Ny_sub+Ny_sub-1]=u_down[i];
        // v[i*Ny_sub+Ny_sub-1]=v_down[i];
    }
    //}

    // if(position_y - 1 > 0) {
    // rank_up = myrank - Px;
    for(int i = 0; i < Nx_sub; i++) { // Change here
                                      //            u_up[i] = u[i * Ny_sub + Ny_sub];
                                      //            v_up[i] = v[i * Ny_sub + Ny_sub];
                                      //
        u_up_send[i] = u[i * Ny_sub + 0];
        v_up_send[i] = v[i * Ny_sub + 0];

        // u[i*Ny_sub+0]=u_up[i];
        // v[i*Ny_sub+0]=v_up[i];
        //            u_down[i]=u_up_send[i]/
    }
}




void burgers::Communication()
{
    
    

    if(position_x  > 0) {
        rank_left = myrank - 1;
        // cout<<"my rank: "<<myrank<<" receiving left from "<<rank_left<<endl;
        MPI_Recv(u_left, Ny_sub, MPI_DOUBLE, rank_left, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(v_left, Ny_sub, MPI_DOUBLE, rank_left, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // cout<<"my rank: "<<myrank<<" received left"<<endl;
        // cout<<"my rank: "<<myrank<<" sending left to "<<rank_left<<endl;
        MPI_Send(u_left_send, Ny_sub, MPI_DOUBLE, rank_left, 3, MPI_COMM_WORLD);
        MPI_Send(v_left_send, Ny_sub, MPI_DOUBLE, rank_left, 4, MPI_COMM_WORLD);
        // cout<<"my rank: "<<myrank<<" sent left"<<endl;
    }

    if(position_x + 1 < Px) {
        rank_right = myrank + 1;
        // cout<<"my rank: "<<myrank<<" sending right to "<<rank_right<<endl;
        MPI_Send(u_right_send, Ny_sub, MPI_DOUBLE, rank_right, 1, MPI_COMM_WORLD);
        MPI_Send(v_right_send, Ny_sub, MPI_DOUBLE, rank_right, 2, MPI_COMM_WORLD);
        // cout<<"my rank: "<<myrank<<" sent right"<<endl;
        // cout<<"my rank: "<<myrank<<" receiving right from "<<rank_right<<endl;

        MPI_Recv(u_right, Ny_sub, MPI_DOUBLE, rank_right, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(v_right, Ny_sub, MPI_DOUBLE, rank_right, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // cout<<"my rank: "<<myrank<<" received right"<<endl;
    }
    // cout<<"my rank: "<<myrank<<"Got here"<<endl;
    if(position_y + 1 < Py) {
        rank_down = myrank + Px;
        MPI_Send(u_down_send, Nx_sub, MPI_DOUBLE, rank_down, 11, MPI_COMM_WORLD);
        MPI_Send(v_down_send, Nx_sub, MPI_DOUBLE, rank_down, 12, MPI_COMM_WORLD);

        MPI_Recv(u_down, Nx_sub, MPI_DOUBLE, rank_down, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(v_down, Nx_sub, MPI_DOUBLE, rank_down, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(position_y > 0) {
        rank_up = myrank - Px;
        MPI_Recv(u_up, Nx_sub, MPI_DOUBLE, rank_up, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(v_up, Nx_sub, MPI_DOUBLE, rank_up, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(u_up_send, Nx_sub, MPI_DOUBLE, rank_up, 13, MPI_COMM_WORLD);
        MPI_Send(v_up_send, Nx_sub, MPI_DOUBLE, rank_up, 14, MPI_COMM_WORLD);
    }

    // MPI_Barrier(MPI_COMM_WORLD);
}





void burgers::Energy()
{
    // cout << "Running energy: "<<myrank<< endl;
    //MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // gets the current rank of a proccessor
    if (myrank==0){
    double Energy = 0.0;
    //double Energy_new = 0.0;
    // cout << "Barrier " << endl;
    // (myrank==0){
    //for(int i = 0; i < Nx * Ny; i++) {
        for(int i = 0; i < Nx * Ny; i++) {

        // MPI_Recv(&Energy,1,MPI_DOUBLE,myrank,50,MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receives energy

            Energy += U_assembled[i] * U_assembled[i] + V_assembled[i] * V_assembled[i];
        
        //Energy += U_assembled[i]*U_assembled[i]+V_assembled[i]*V_assembled[i];
        // Energy *= 0.5 * dx * dy;
        }
        cout.precision(10);
        Energy *= 0.5 * dx * dy;
        cout << "Energy is: " << Energy << endl;
        
    }

//            cout.precision(10);
//        Energy *= 0.5 * dx * dy;
//        cout << "Energy is: " << Energy << endl;
//    }
    //cout << "rank: " << myrank << "\t Energy: " << 0.5 * Energy * dx * dy << endl;
//
//    if(myrank != 0) {
//        MPI_Send(&Energy, 1, MPI_DOUBLE, 0, 555, MPI_COMM_WORLD); // passes energy
//    } else{
//        for(int i = 1; i < Px * Py; i++) {
//            MPI_Recv(&Energy_new, 1, MPI_DOUBLE, i, 555, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            Energy += Energy_new;
//        }
//        cout.precision(10);
//        Energy *= 0.5 * dx * dy;
//        cout << "Energy is: " << Energy << endl;
//        
//    }

    

    delete[] u;
    delete[] v;
}
