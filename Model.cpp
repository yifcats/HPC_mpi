//
//  model.cpp
//  HPC_CW
//
//  Created by Yiftach Cats on 10/03/2019.
//  Copyright © 2019 Yiftach Cats. All rights reserved.
//

#include "Model.h"


using namespace std;

Model::Model(int argc, char* argv[]){
    
    
    IsValid(argc, argv);
    ParseParameters(argc, argv);
    
    
    
    Print_help();
    PrintParameters();
    
    
    
    //ValidateParameters(argc,argv); // print status of mission
    


}



void Model::PrintParameters(){ // print all the prameters
    if (verbose){
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
        
}


void Model::IsValid(int argc, char* argv[]){ // return value of valid if all of bellow is true which will alow to run the functions.
    
    int mpi_size;
    
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    
    int size_mesh=stoi(argv[8],NULL,10)*stoi(argv[9],NULL,10);
    
    if (mpi_size != size_mesh){ // the px and py must be a factor of mpi_size
        cout<<"\n\n\n\nInvalid: The Number of Processors must equal the product of Px and Py\n\n\n\n" << endl;
        exit(EXIT_SUCCESS);
    } 
    
    if (argc!=11 && argc!=10){ // checking if the correct number of arguments have been passed by the user 7 or 7 + argument (v or h)
        cout<<"\n\n\n\nInvalid: Number of Inputs Exceeds the Maximum allowed\n\n\n\n" << endl;
        exit(EXIT_SUCCESS);
        
    }
    
    
    if (argc==11){
        if (*argv[10]=='v'){
        verbose=true;
    
    } else if (*argv[10]=='h'){
        help=true;
    
    } else {
        cout<<"\n\n\n\nInvalid User Inputs must be single character h: help v: verbose\n\n\n\n" << endl;
        help=false;
        verbose=false;
        exit(EXIT_SUCCESS);
    }
    
    }
    

}


// Inputing Parameters
void Model::ParseParameters(int argc, char* argv[]){ // parsing the parameters
    
    try{
        ax=stod(argv[1]); // convertes strings to doubles
        ay=stod(argv[2]);
        b=stod(argv[3]);
        c=stod(argv[4]);
        
        Nt=stoi(argv[5],NULL,10);
        Nx=stoi(argv[6],NULL,10);
        Ny=stoi(argv[7],NULL,10);
        
        Px=stoi(argv[8],NULL,10);
        Py=stoi(argv[9],NULL,10);
        
        
    }
    
    catch (const invalid_argument){ // checking for invalid agument errors
        cout << "\n\n\n\nInvalid Entery for ax, ay, b and c. Data Type double\n\n\n\n"<<endl;
        Print_help();
        exit(EXIT_SUCCESS);
        return;
    }
    
    // Initilising the step sizes
    dx=Lx/(Nx-1);
    dy=Ly/(Ny-1);
    dt=T/Nt;
    
}




//void Model::ValidateParameters(int argc, char* argv[]){ // checking each entery is vaild
//   // cout << "Testing Success: "<<IsValid(argc, argv)<<endl; // enter if valid or not
//}


void Model::Print_help(){ // just prints out the help

    if (help){
    cerr<<"\nHelp Documentation:"
    <<"\n\n Input: {ax, ay, b, c, Px, Py, Optional(User)}"
    <<"\n\n Data Type:"
    <<"\t ax: 'double'"
    <<"\t ay: 'double'"
    <<"\n\n\t\t b: 'double'"
    <<"\t c: 'double'"
    <<"\n\n\t\t Px: 'int'"
    <<"\t Py: 'int'"
    <<"\n\n\t\t User: 'char'"
    <<"\n\n\n User: <Options>"
    <<"\n\n\t h \t help"
    <<"\n\t v \t verbose\n"<<endl;
    }
}




