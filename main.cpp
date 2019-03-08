#ifndef CLASS_MODEL
#define CLASS_MODEL



#include <iostream>

//#include <mpi.h>

using namespace std;

class Model {
public:
    
    Model(int argc, char* argv[]){
        
        for (int i=0;i<argc;i++){
            if (*argv[i]=='h'){ // if( string(argv[i]=="-h")
                Print_help();
                help=true;
            }
            else if (*argv[i]=='v'){
                verbose=true;
            }
        }
        if (help!=true){
            ParseParameters(argc,argv);
            Lx=10;
            Ly=10;
            T=0;
            
            if (verbose==true){
                PrintParameters();
            }
        }
    }
    
    // ~Model();
    
    void PrintParameters(){
        cout<<"hi"<<endl; // print all the prameters
    }
    
    bool IsValid(){ // return value of valid if all
        
        bool valid;
        return valid;
    }
    
    
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
    
private:
    
    
    // Inputing Parameters
    void ParseParameters(int argc, char* argv[]){
        if (argc!=6 && argc!=5){ // checking if the correct number of arguments have been passed by the user
            cout<<"Invialid Parameter Inputs"<<endl;
            Print_help();
            return;
        }
        
        ax=stod(argv[1]);
        ay=stod(argv[2]);
        b=stod(argv[3]);
        c=stod(argv[4]);
        
        
        
        cout<<ax<<ay<<b<<c<<endl;
    }
    
    void ValidateParameters(int argc, char* argv[]){ // checking each entery is vaild
        
    }
    
    void Print_help(){ // just prints out the help
        cerr<<"\nHelp Documentation:"
        <<"\nInput: ax, ay, b, c, User"
        <<"\n\nUser: <Options>"
        <<"\n\t h \t help"
        <<"\n\t v \t verbose\n"<<endl;
    }
    
    
    bool verbose;
    bool help; // input a description of the parameters
    
    // Numerics
    double x0;
    double y0;
    double Lx;
    double Ly;
    double T;
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
};

#endif

