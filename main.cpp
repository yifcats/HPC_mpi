#include <chrono>
#include <cmath>


#include "Model.h"
#include "Burgers.h"



int main(int argc, char* argv[]){
    
    
    Model m(argc, argv);
    burgers b(&m);
    
    // Call code to initialise the problem here
    b.Initial_velocity();
    b.Integrate_velocity();
    
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
    
    //Call code to perform time integration here
    
    hrc::time_point end = hrc::now();
    
    //Calculate final energy and write output
    b.Energy();
    
    
   
    
    auto diff = end - start;
    cout <<"\n\nexexution time: "<< chrono::duration <double, milli> (diff).count() << " ms" << endl;
    return 0;
}

