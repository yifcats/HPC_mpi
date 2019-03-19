#ifndef BURGERS_H
#define BURGERS_H


#include "Model.h"


class burgers{

public:
    burgers(Model* m);
    //~burgers();
    
    // Declaring functions
    void Initial_velocity();
    void Integrate_velocity();
    void Energy();
    void Create_sending_boundaris();
    void Communication();
    void Set_integration_boundaries();
    void Assembling_Submatices();
    
    void Current_velocity(int i, int j);
    
    
    void Print_matrix(); // matrix printing for part 1b)
//    int GetRank();
    
private:
    Model* m;
    // Numerics
    double*x;
    double*y;
    
    // pointer to u and v velocity arrays
    double*u;
    double*v;
    
    // Creating the boundary Arrays to be recived
    
    double*u_left;
    double*v_left;
    
    double*u_right;
    double*v_right;
    
    double*u_up;
    double*v_up;
    
    
    double*u_down;
    double*v_down;
    
    
    
    double*u_left_send;
    double*v_left_send;
    
    double*u_right_send;
    double*v_right_send;
    
    double*u_up_send;
    double*v_up_send;
    
    
    double*u_down_send;
    double*v_down_send;
    
    
    
    //int myrank;
    int myrank;
    
    int rank_left;
    int rank_right;
    int rank_up;
    int rank_down;
    
    /////////////////////////
    int Nx_sub;
    int Ny_sub;
    
    int position_x;
    int position_y;
    
    int Position_global_x=0;
    int Position_global_y=0;
    
    int rem_x;
    int rem_y;
    
    // Assembling Marix
    
    double* U_assembled;
    double* V_assembled;
    
    int* array_global_y;
    int* array_global_x;
           
    int* array_size_Nx_sub;
    int* array_size_Ny_sub;
    
    
    
    // Initilizing boundaries for integrations baised on locations of the processors
    int Lower_i = 0;
    int Upper_i = 0;
    int Lower_j = 0;
    int Upper_j = 0;
    
    
    
    
    
    
    
    
    
    
    
    
    // current values for u and v within the loop using boundary conditions and processors vals sending/reciving
    
    double u_i_minus_j=0.0;
    double u_i_plus_j=0.0;
      
    double u_i_j_minus=0.0;
    double u_i_j_plus=0.0;
    
    double v_i_minus_j=0.0;
    double v_i_plus_j=0.0;
      
    double v_i_j_minus=0.0;
    double v_i_j_plus=0.0;
    
    // indicies
    int i_plus_j=0;
    int i_minus_j=0;
    int i_j=0;
    int i_j_plus=0;
    int i_j_minus=0;
    
    
    int i_j_2=0;
    
    
    
    
    
    
    
    
    
    
    
    ////////////////////////
    
    
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
    
    // Physics
    double ax;
    double ay;
    double b;
    double c;
    
    // position of processor
    

    
    
    
    int Px;
    int Py;
//    
//    int Nx_sub;
//    int Ny_sub;
    
};



#endif // BURGERS_H
