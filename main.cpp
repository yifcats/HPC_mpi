#include <chrono>
#include <cmath>


#include "Model.h"
#include "Burgers.h"


#include <mpi.h>


int main(int argc, char* argv[]){
   
   int mpi_init_err=MPI_Init(&argc,&argv);

    if(mpi_init_err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI" << endl;
        return -1;
    }
   
   
    
    
    Model m(argc, argv);
    burgers b(& m);
    
    // Call code to initialise the problem here
    b.Initial_velocity();
    
    
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
    b.Integrate_velocity();
    //Call code to perform time integration here

    hrc::time_point end = hrc::now();
    //MPI_Barrier(MPI_COMM_WORLD); // waits until other process have done the same
    //Calculate final energy and write output
    b.Energy();
    
    
   
    auto diff = end - start;
    cout <<"\n\nexexution time: "<< chrono::duration <double, milli> (diff).count() << " ms" << endl;
    
    MPI_Finalize();
    return 0;
}

