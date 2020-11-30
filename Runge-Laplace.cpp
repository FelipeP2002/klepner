  
#include <iostream>
#include <cmath>
#include <vector>
#include "integration.h"

const double E = 0.6;
const double T0 = 0.0;
const double TF = 20*M_PI;
const double H = 0.0005;
const int NSTEPS = (TF)/H;
const int DIM = 4;

template<typename system_t, typename state_t>
    void integrate(system_t deriv, state_t & data, double tmin, double tmax, double dt);

typedef std::vector<double> state_t;
/* y[0]=q_1
 * y[1]=q_2
 * y[2]=q_1'=p_1
 * y[3]=q_2'=p_2
 */

int main(void)
{
    state_t y(DIM);
    y = {1 - E , 0.0 , 0.0 , std::sqrt( (1 + E)/(1 - E) ) };  // initial conditions

    integrate(fderiv, y, T0, TF, H);

    return 0;
}

template<typename system_t, typename state_t>
    void integrate(system_t deriv, state_t & data, double tmin, double tmax, double dt)
{
    
    std::cout.precision(7);
    const int nsteps = (tmax-tmin)/dt;

    state_t A(2);
    
    for (int ii = 0; ii <= nsteps; ++ii){
        double t = tmin + ii*dt;

        R_L_R(data, A, t);

        std::cout<< t << "  ";
        for(double val : A){std::cout<< val << "  ";}
        std::cout<<"\n";
        
        rk4(deriv, data, t, dt);

    }
}
