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
    void Create_boundaris();
    void Communication();
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
