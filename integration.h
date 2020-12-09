#include <iostream>
#include <cmath>
#include <vector>

#ifndef __INTEGRATION_H_
#define __INTEGRATION_H_

typedef std::vector <double> state_t;
void fderiv(const state_t & y, state_t & dydt, double t)
{
    dydt[0] = y[2];                                         // q_1'=p_1
    dydt[1] = y[3];                                         // q_2'=p_2
    dydt[2] = -y[0]/(std::hypot(y[0] , y[1]) * std::hypot(y[0] , y[1]) * std::hypot(y[0] , y[1])); // q_1''
    dydt[3] = -y[1]/(std::hypot(y[0] , y[1]) * std::hypot(y[0] , y[1]) * std::hypot(y[0] , y[1])); // q_2''
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

    H = (data[2]*data[2] + data[3]*data[3])/2 - 1/std::hypot(data[0], data[1]);
    
    return H;
}

template< typename S,  typename state_t>
S Momentum(S t, state_t data)
{
    double L{};

    L = data[0]*data[3] - data[1]*data[2];

    return L;

}


template< typename state_t>
state_t R_L_R(state_t data, state_t & A, double t)
{

    A[0] = data[3] * Momentum(t, data) - data[0]/std::hypot(data[0] , data[1]);
    A[1] = -data[2] * Momentum( t , data) - data[1]/std::hypot(data[0] , data[1]);

    return A;
}

template<typename system_t, typename state_t>
void int_4(system_t deriv  ,state_t & data, double t, double h){

state_t dydt(data.size());
state_t map(8);

double K = std::cbrt(2.0);

map[0] = 1.0/(2.0*(2.0-K));//c_1
map[2] = (1.0-K)/(2.0*(2.0-K));//c_2
map[4] = map[2];//c_3
map[6] = map[0];//c_4

map[1] = 1.0/(2.0-K);//d_1
map[3] = -K/(2.0-K);//d_2
map[5] = map[1];//d_3
map[7] = 0;//d_4


for(int kk = 0; kk < 4; ++kk)
{

  deriv(data, dydt, t);

    for(int ii = 0; ii < 2 ; ++ii){
      data[ii] = data[ii] + h*map[0 + 2*kk]*dydt[ii];
    }

  deriv(data, dydt, t);

    for(int jj = 2; jj < 4 ; ++jj){
     data[jj] = data[jj] + h*map[1 + 2*kk]*dydt[jj];
    }

  }

}

template<typename system_t, typename state_t>
void int_6(system_t deriv  ,state_t & data, double t, double h){

state_t dydt(data.size());
state_t map(16);


map[0] = 0.392256805238780;//c_1
map[2] = 0.510043411918458;//c_2
map[4] = -0.471053385409758;//c_3
map[6] = 0.0687531682525198;//c_4
map[8] = map[6];//c_5
map[10] = map[4];//c_6
map[12]= map[2];//c_7
map[14] = map[0];//c_8

map[1] = 0.784513610477560;//d_1
map[3] = 0.235573213359357;//d_2
map[5] = -1.17767998417887;//d_3
map[7] = 1.31518632068391;//d_4
map[9] = map[5];//d_5
map[11] = map[3];//d_6
map[13] = map[1];//d_7
map[15] = 0;//d_8

for(int kk = 0; kk < 8; ++kk)
{

  deriv(data, dydt, t);

    for(int ii = 0; ii < 2 ; ++ii){
      data[ii] = data[ii] + h*map[2*kk]*dydt[ii];
    }

  deriv(data, dydt, t);

    for(int jj = 2; jj < 4 ; ++jj){
     data[jj] = data[jj] + h*map[1 + 2*kk]*dydt[jj];
    }

  }

}
#endif // __INTEGRATION_H_
