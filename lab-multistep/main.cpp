#include<iostream>
#include "multistep-funcs.h"
#include<fstream>

int main() {
    // std::cout << "Hello, world!" << std::endl;

    // Oscillator_Multi::StateAndArg init;
    // init.state = Eigen::Vector<double, Oscillator_Multi::dim>(0.0, 1.0);
    // init.arg = 0.0;

    IntegrationParameters params;
    params.epsilon = pow(10, -6);
    params.maxIter = 10;

    // Oscillator_Multi osc;
    // Linear_Multi lin; 
    // Linear_Multi::StateAndArg init_lin;
    // init_lin.state = Eigen::Vector<double, Linear_Multi::dim>(0.0);
    // init_lin.arg = 0.0;

    // std::cout << integrate<BDF4, Oscillator_Multi, RK4Table>(init, 5, params, osc)[100].state[0] << std::endl;
    // std::cout << integrate<BDF4, Linear_Multi, RK4Table>(init_lin, 5, params, lin)[100].state[0] << std::endl;
    // std::cout << integrate<BDF4, Oscillator_Multi, RK4Table>(init, 5, params, osc)[100].arg << std::endl;

    std::ofstream out("multi-data.txt");
    double endTime = 5.;

    Oscillator_Multi::StateAndArg init;
    init.state = Eigen::Vector<double, Oscillator_Multi::dim>(0.0, 1.0);
    init.arg = 0.0;
    Oscillator_Multi osc; 

    // std::cout << integrate<BDF4, Oscillator_Multi, RK4Table>(init, endTime, params, osc)[5/0.124-1].state[0] << std::endl;

    Linear_Multi::StateAndArg init_linear;
    init_linear.state = Eigen::Vector<double, Linear::dim>(0.0);
    init_linear.arg = 0.0;
    Linear_Multi linear;

    for(int i = 5; i <= 1000; i++) {
        params.step = 5./i;
        auto res_osc = integrate<BDF4, Oscillator_Multi, RK4Table>(init, endTime, params, osc);
        auto res_lin = integrate<BDF4, Linear_Multi, RK4Table>(init_linear, endTime, params, linear);

        double max_osc = 0;
        double max_lin = 0;
        for(int j = 0; j < i-2; j++) {
            if(max_osc < fabs(res_osc[j].state[0] - sin(res_osc[j].arg))) {
                max_osc = fabs(res_osc[j].state[0] - sin(res_osc[j].arg));
            }
            if(max_lin < fabs(res_lin[j].state[0] - pow(res_lin[j].arg, 4)/4)) {
                max_lin = fabs(res_lin[j].state[0] - pow(res_lin[j].arg, 4)/4);
            }
        }
        out << endTime/i << '\t' << max_osc << '\t' << max_lin << std::endl;
    }
}