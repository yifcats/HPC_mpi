#ifndef CLASS_MODEL
#define CLASS_MODEL



#include <iostream>
#include <cctype>
// #include <stdexcept>

#include <mpi.h>


using namespace std;

class Model {
public:
    
    Model(int argc, char* argv[]){
        ValidateParameters(argc,argv); // print status of mission
        if (IsValid(argc,argv)){
            
            for (int i=0;i<argc;i++){
                if (*argv[i]=='h'){ // if( string(argv[i]=="-h")
                    Print_help();
                    help=true;
                }
                else if (*argv[i]=='v'){
                    verbose=true;
                } else {
                    help=false;
                    verbose=false;
                }
            }
            
            
            if (!help){ // print help if entery is inputed
                ParseParameters(argc,argv);
                Lx=10;
                Ly=10;
                T=0;
            }
            if (verbose){ // print paramteters if verbose is true.
                PrintParameters();
            }
            
        } else { // print help if invallid parameters are entered
            Print_help();
            
        }
    }
    
    // ~Model();
    
    void PrintParameters(){ // print all the prameters
        cout<<"\n Numerics:\n"<<endl;
        cout<<"x0:\t"<<x0<<endl;
        cout<<"y0:\t"<<y0<<endl;
        cout<<"Lx:\t"<<Lx<<endl;
        cout<<"Ly:\t"<<Ly<<endl;
        cout<<"T:\t"<<T<<endl;
        cout<<"Nx:\t"<<Nx<<endl;
        cout<<"Ny:\t"<<Ny<<endl;
        cout<<"Nt:\t"<<Nt<<endl;
        cout<<"dx:\t"<<dx<<endl;
        cout<<"dy:\t"<<dy<<endl;
        cout<<"dt:\t"<<dt<<endl;
        
        cout<<"\n Physics:\n"<<endl;
        cout<<"ax:\t"<<ax<<endl;
        cout<<"ay:\t"<<ay<<endl;
        cout<<"b:\t"<<b<<endl;
        cout<<"c:\t"<<c<<endl;
    }
    
    bool IsValid(int argc, char* argv[]){ // return value of valid if all of bellow is true which will alow to run the functions.
        bool valid=true;
        bool value=true;
        int size;
        
        if (argc!=6 && argc!=5){ // checking if the correct number of arguments have been passed by the user
            valid=false;
            
        } else {
            if (argc==6){
                size=1;
                if (*argv[5]!='h' && *argv[5]!='v'){ // isalpha(argv[5]) checking that character is inputed is within the once diffined for user help or v
                    valid=false;
                    cout<<value<<endl;
                }
                
            }else { // else size ==5
                size=0;
            }
            
            for (int i=1;i<argc-size;i++) {
                
                if (isdigit(*argv[i])){ // check if all values are digits (2-5)
                    value=true;
                    
                } else {
                    value=false;
                }
                
                valid*=value;
                //cout << "Testing Parameter ["<< i << "]: "<< value << endl; // say when it fails
            }
           //cout << "Testing Success: "<< valid << endl; // enter if valid or not
        }
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
    void ParseParameters(int argc, char* argv[]){ // parsing the parameters
        
        try{
            ax=stod(argv[1]); // convertes strings to doubles
            ay=stod(argv[2]);
            b=stod(argv[3]);
            c=stod(argv[4]);
            
        }
        
        catch (const invalid_argument){ // checking for invalid agument errors
            cout << "Last Check Found An Error"<<endl;
            Print_help();
            return;
        }
        
        cout<<ax<<ay<<endl;
        
//        if (argc!=6 && argc!=5){ // checking if the correct number of arguments have been passed by the user
//          cout<<"Invialid Parameter Inputs"<<endl;
//         Print_help();
//        return;
//         }
//        ax=stod(argv[1]); // convertes strings to doubles
//        ay=stod(argv[2]);
//        b=stod(argv[3]);
//        c=stod(argv[4]);
    }
    
    void ValidateParameters(int argc, char* argv[]){ // checking each entery is vaild
        cout << "Testing Success: "<<IsValid(argc, argv)<<endl; // enter if valid or not
    }
    
    void Print_help(){ // just prints out the help
        cerr<<"\nHelp Documentation:"
        <<"\n\n Input: {ax, ay, b, c, Optional(User)}"
        <<"\n\n Data Type:"
        <<"\t ax: 'double'"
        <<"\t ay: 'double'"
        <<"\n\n\t\t b: 'double'"
        <<"\t c: 'double'"
        <<"\n\n\t\t User: 'char'"
        <<"\n\n\n User: <Options>"
        <<"\n\n\t h \t help"
        <<"\n\t v \t verbose\n"<<endl;
    }
    
    bool verbose;
    bool help; // input a description of the parameters
    
    // Numerics
    double x0;
    double y0;
    double Lx=10;
    double Ly=10;
    double T;
    int    Nx=20;
    int    Ny=20;
    int    Nt=30;
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

