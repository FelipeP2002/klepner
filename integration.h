#include <iostream>

#ifndef __INTEGRATION_H_
#define __INTEGRATION_H_

typedef std::vector<double> state_t;
void fderiv(const state_t & y, state_t & dydt, double t)
{
    dydt[0] = y[2];                                         // q_1=p_1
    dydt[1] = y[3];                                         // q_2'=p_2
    dydt[2] = -y[0]/(std::hypot( y[0] , y[1]) * std::hypot( y[0] , y[1]) * std::hypot( y[0] , y[1])); // q_1''
    dydt[3] = -y[1]/(std::hypot( y[0] , y[1]) * std::hypot( y[0] , y[1]) * std::hypot( y[0] , y[1])); // q_2''
}

template<typename systema_t, typename state_t>
void euler_central(systema_t deriv, state_t & data, double t, double h)
{
    state_t dydt(data.size());
    state_t y_1(data.size());
    state_t y_2(data.size());
    deriv(data, dydt, t);

    for(int ii = 0; ii < data.size(); ++ii) {
        y_1[ii] = data[ii] + h*dydt[ii];
    }
    
    for(int ii = 0; ii < data.size(); ++ii) {
        y_2[ii] = (data[ii] + y_1[ii])/2;
    }
    
    deriv(y_2, dydt, t);
    
    for(int ii = 0; ii < data.size(); ++ii) {
        data[ii] = data[ii] + h*dydt[ii];
    }
}

template<typename system_t, typename state_t>
void euler(system_t deriv, state_t & data, double t, double h)
{
    state_t dydt(data.size());
    deriv(data, dydt, t);

    for(int ii = 0; ii < data.size(); ++ii) {
        data[ii] = data[ii] + h*dydt[ii];
    }
}


template<typename system_t, typename state_t>
    void rk4(system_t deriv, state_t & data, double t, double h)
{
    state_t dydt(data.size());
    state_t k1(data.size()), k2(data.size()),
        k3(data.size()), k4(data.size()), aux(data.size());

    // k1
    deriv(data, dydt, t);
    for(int ii = 0; ii < data.size(); ++ii) {
        k1[ii] = h*dydt[ii];
    }
    // k2 aux
    for(int ii = 0; ii < data.size(); ++ii) {
        aux[ii] = data[ii] + k1[ii]/2;
    }
    // k2
    deriv(aux, dydt, t + h/2);
    for(int ii = 0; ii < data.size(); ++ii) {
        k2[ii] = h*dydt[ii];
    }
    // k3 aux
    for(int ii = 0; ii < data.size(); ++ii) {
        aux[ii] = data[ii] + k2[ii]/2;
    }
    // k3
    deriv(aux, dydt, t + h/2);
    for(int ii = 0; ii < data.size(); ++ii) {
        k3[ii] = h*dydt[ii];
    }
    // k4 aux
    for(int ii = 0; ii < data.size(); ++ii) {
        aux[ii] = data[ii] + k3[ii];
    }
    // k4
    deriv(aux, dydt, t + h);
    for(int ii = 0; ii < data.size(); ++ii) {
        k4[ii] = h*dydt[ii];
    }

    // write new data
    for(int ii = 0; ii < data.size(); ++ii) {
        data[ii] = data[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    }

}

template<typename systema_t, typename state_t>
void Velvet(systema_t deriv, state_t & data, double t, double h)
{
    state_t dydt(data.size());

    deriv(data, dydt, t);

    for(int ii = 0; ii < 2 ; ++ii ){
        data[ii + 2] += h*dydt[ii + 2];};

    for(int jj = 0; jj < 2 ; ++jj){
        data[jj] += h*data[jj + 2];};
}

template<typename systema_t, typename state_t>
void Velvet_initial_conditions(systema_t deriv, state_t & data, double t, double h)
{
    state_t dydt(data.size());

    deriv(data, dydt, t);

    for(int ii = 0; ii < 2 ; ++ii ){
        data[ii + 2] -= h*dydt[ii + 2];};

}

template< typename S, typename state_t>
S Hamilton(S t, state_t  data)
{
    double H{};

    H = (data[2]*data[2] + data[3]*data[3])/2 - 1/std::hypot(data[0] , data[1]);
    
    return H;
}

template< typename S, typename state_t>
S Momentum(S t, state_t data)
{
    double L{};

    L = data[0]*data[3] - data[1]*data[2];

    return L;

}


#endif // __INTEGRATION_H_
