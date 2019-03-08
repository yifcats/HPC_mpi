#include <stdio.h>
#include <cmath>

int main(int argc, char **argv)
{
    
    
    // MATLAB TO HPC
    
    // Test Parameters
    
    int const nx=20;
    int const ny=20;
    int const L=10;
    int const T=1;
    int const nt=30;
    
    double ax=1;
    double ay=1;
    double b=0.1;
    double c=0.02;
    
    double uArray[nx][ny];
    double unArray[nx][ny];
    double vArray[nx][ny];
    double vnArray[nx][ny];
    double rhsxArray[nx][ny];
    double rhsyArray[nx][ny];
    double xArray[nx][ny];
    
    double dx = L / (nx -1);
    double dy = L / (ny -1);
    double dt = T / nt;
    
    double x[nx];
    double y[ny];
    double u[nx][ny];
    double v[nx][ny];
    double un[nx][ny];
    double vn[nx][ny];
    
    
    x[0]=-L/2;
    y[0]=-L/2;
    
    for (int temp=1; temp++; x[nx]>0){
        x[temp]=x[temp-1]+dx;
    }
    
    for (int temp=1; temp++; y[ny]>0){
        y[temp]=y[temp-1]+dy;
    }
    // Velocity profile
    // Initial Conditions
    
    for  (int i = 0; i < nx; i++){
        for  (int j = 0; j < ny; j++){
            double r = pow(pow(x[i],2)+pow(y[j],2),2);
            if (r<=1){
                u[i][j]= 2 * pow(1-r,4) * (4 *r +1);
                v[i][j]= u[i][j];
            }
            else {
                u[i][j]=0;
                v[i][j]=0;
            }
        }
    }
    
    // Boundary Conditions
    
    for (int temp=0;  temp++; temp < ny){
        u[0][temp]=0.0;
        u[nx][temp]=0.0;
        v[0][temp]=0;
        v[nx][temp]=0;
    }
    
    for (int temp=0;  temp++; temp < nx){
        u[temp][1]=0.0;
        u[temp][ny]=0.0;
        v[temp][1]=0;
        v[temp][ny]=0;
    }
    
    
    // integration
    int mArray[nx-2];
    for (int m =1; m < nx - 1; m++ ) {
        mArray[m]=m;
    }
    
    int nArray[ny-2];
    for (int n =1; n < ny - 1; n++ ) {
        nArray[n]=n;
    }
    
    for (int it=0; it < nt+1; it++){
        
        
        for (int tempx=0;  tempx++; tempx < nx){
            for (int tempy=0;  tempy++; tempy < ny){
                un[tempx][tempy] = u[tempx][tempy];
                vn[tempx][tempy] = v[tempx][tempy];
            }
        }
        
        for (int i =1; i < nx-1; i++ ) {
            
            for (int j =1; j < ny-1; j++ ) {
                
                u[i][j] = un[i][j] - (dt*(un[i][j]-un[i-1][j])*(ax+b*un[i][j])/dx) - (dt*(un[i][j]-un[i][j-1])*(ay+b*vn[i][j])/dy) + (c*dt*(un[i+1][j]-2*un[i][j]+un[i-1][j])/(dx*dx)) + (c*dt*(un[i][j-1]-2*un[i][j]+un[i][j+1])/(dy*dy));
                
                v[i][j] = vn[i][j] - (dt*(vn[i][j]-vn[i-1][j])*(ax+b*un[i][j])/dx) - (dt*(vn[i][j]-vn[i][j-1])*(ay+b*vn[i][j])/dy) + (c*dt*(vn[i+1][j]-2*vn[i][j]+vn[i-1][j])/(dx*dx)) + (c*dt*(vn[i][j-1]-2*vn[i][j]+vn[i][j+1])/(dy*dy));
            }
        }
        
        for (int temp=0;  temp++; temp < ny){
            u[0][temp]=0.0;
            u[nx][temp]=0.0;
            v[0][temp]=0;
            v[nx][temp]=0;
        }
        
        for (int temp=0;  temp++; temp < nx){
            u[temp][1]=0.0;
            u[temp][ny]=0.0;
            v[temp][1]=0;
            v[temp][ny]=0;
        }
        
    }
    
    return 0;
}
