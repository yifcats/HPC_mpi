#ifndef BURGERS_H
#define BURGERS_H

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
    
    // pointer to u and v velocity arrays
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
    
    // Physics
    double ax;
    double ay;
    double b;
    double c;
    
    //
    
};



#endif // BURGERS_H
