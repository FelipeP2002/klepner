  
#include <iostream>
#include <cmath>
#include <vector>
#include "integration.h"

const double E = 0.6;
const double T0 = 0.0;
const double TF = 35;
const double H = 0.0005;
const int NSTEPS = (TF)/H;
const int DIM = 2;

template<typename state_t>
    void integrate(state_t & data);

typedef std::vector<double> state_t;
/* y[0]=q_1
 * y[1]=q_2
 * y[2]=q_1'=p_1
 * y[3]=q_2'=p_2
 */
int main(void)
{
    state_t y(2*DIM);
    y = {1 - E , 0.0 , 0.0 , std::sqrt( (1 + E)/(1 - E) ) };  // initial conditions

    integrate(y);

    return 0;
}

template<typename state_t>
    void integrate(state_t & data )
{
    
    std::cout.precision(7);
    const int nsteps = (TF-T0)/H;

    Velvet_initial_conditions(data, H);
    
    for (int ii = 0; ii <= nsteps; ++ii){
        double t = T0 + ii*H;
        std::cout<< t << "  ";
        for(double val : data ){std::cout<< val << "  ";}
        std::cout<<"\n";
        
        Verlet(data,H);

    }
}
