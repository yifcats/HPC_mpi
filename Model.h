#ifndef CLASS_MODEL
#define CLASS_MODEL



#include <iostream>
#include <cctype>
//#include <stdlib.h>
//#include <stdio.h>
// #include <stdexcept>
#include <mpi.h>


using namespace std;

class Model {
public:
    
    Model(int argc, char* argv[]);
    void PrintParameters();
    //void ValidateParameters(bool valid);
    void IsValid(int argc, char* argv[]);
    // ~Model();
    
        

    
    // Getters
    bool   IsVerbose() const { return verbose; }
    bool   IsHelp()    const { return help; }
    double GetX0()     const { return x0; }
    double GetY0()     const { return y0; }
    double GetLx()     const { return Lx; }
    double GetLy()     const { return Ly; }
    double GetT()      const { return T; }
    int    GetNx()     const { return Nx; }
    int    GetNy()     const { return Ny; }
    int    GetNt()     const { return Nt; }
    double GetDx()     const { return dx; }
    double GetDy()     const { return dy; }
    double GetDt()     const { return dt; }
    double GetAx()     const { return ax; }
    double GetAy()     const { return ay; }
    double GetB()      const { return b; }
    double GetC()      const { return c; }
    
    // Add any other getters here...
    int GetPx()      const { return Px; }
    int GetPy()      const { return Py; }
private:
    
    void ParseParameters(int argc, char* argv[]);
    void ValidateParameters(int argc, char* argv[]);
    void Print_help();
    
    bool verbose=false;
    bool help=false; // input a description of the parameters
    // Numerics
    double x0;
    double y0;
    double Lx=10;
    double Ly=10;
    double T=1;
    int    Nx;
    int    Ny;
    int    Nt;
    double dx;
    double dy;
    double dt;
    
    // Physics
    double ax;
    double ay;
    double b;
    double c;
    
    // Add any additional parameters here...
    int Px;
    int Py;
    
};

#endif

