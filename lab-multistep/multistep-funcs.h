#include<iostream>
#include<array>
#include<vector>
#include<Eigen/Dense>
#include "onestep-rk.h"

/* Это коэффициенты для метода на 4 шагах*/
struct BDF4{
    static constexpr unsigned int size = 4;
    static constexpr std::array<double, size+1> alpha = {-1./4, 4./3, -3, 4, 25./12};
};

class Oscillator_Multi {
    
public:
    
    static constexpr unsigned int dim = 2;  // размерность задачи
    
    using Argument = double;  // тип аргумента, тип t
    
    using State = Eigen::Vector<double, dim>;  // состояние
    
    struct StateAndArg{
        State state;
        Argument arg;
    };

    /*** Вычисляет разницу двух состояний (решения нелиненого уравнения) ***/
    double calcDif(const State& first, const State& second) const {
        double max = 0;
        for(int i = 0; i < dim; i++) {
            if(fabs(first[i] - second[i]) > max) {
                max = fabs(first[i] - second[i]);
            }
        }
        return max;
    }
    
    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.state(1), -stateAndArg.state(0)};
    } 
};

class Linear_Multi {

    public:

    static constexpr unsigned int dim = 1;

    using Argument = double;

    using State = Eigen::Vector<double, dim>;

    struct StateAndArg{
        State state;
        Argument arg;
    };

    /*** Вычисляет разницу двух состояний (решения нелиненого уравнения) ***/
    double calcDif(const State& first, const State& second) const {
        double max = 0;
        for(int i = 0; i < dim; i++) {
            if(fabs(first[i] - second[i]) > max) {
                max = fabs(first[i] - second[i]);
            }
        }
        return max;
    }

    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{pow(stateAndArg.arg, 3)};
    } 
};

struct IntegrationParameters{
    double step;  // шаг интегрирования
    double epsilon;  // точность решения нелинейного уравнения
    double maxIter;  // максимальное количество итераций для решения нелинейного уравнения
};

/***
BDF - структура с коэффициентами метода
RHS - правая часть Д.У.
RKTable -  таблица Бутчера метода Рунге-Кутты для разгона.
***/

template<typename BDF, typename RHS, typename RKTable>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
    const typename RHS::StateAndArg& initialState,
    const typename RHS::Argument& endTime,
    const IntegrationParameters& parameters,
    const RHS& rhs
)
{
    std::vector<typename RHS::StateAndArg> res;

    typename RHS::StateAndArg init = initialState;
    typename RHS::Argument end = init.arg + parameters.step*4;

    std::vector<typename RHS::StateAndArg> y2 = integrateRK<RK4Table, RHS>(init, end, parameters.step, rhs);
    res.push_back(y2[0]);
    res.push_back(y2[1]);
    res.push_back(y2[2]);
    res.push_back(y2[3]);
    y2[4] = y2[3];
    
    for(typename RHS::Argument t = init.arg + parameters.step*4; t < endTime; t += parameters.step) {
        double dif = 5*parameters.epsilon;
        typename RHS::StateAndArg y;
        y.arg = end;
        y.state = y2[4].state;
        for(int j = 0; j < parameters.maxIter; j++) {
            if(dif > parameters.epsilon) {
                for(int m = 0; m < rhs.dim; m++) {
                    y.state[m] = parameters.step*rhs.calc(y2[4])[m]/BDF4::alpha[4]; 
                    for(unsigned int s = 0; s < BDF4::size; s++) {
                        y.state[m] += BDF4::alpha[s]*y2[s].state[m]/BDF4::alpha[4];
                    }
                }
            }
            else {
                break;
            }
            dif = rhs.calcDif(y.state, y2[4].state);
            y2[4] = y;
        }
        res.push_back(y);

        y2[0] = y2[1];
        y2[1] = y2[2];
        y2[2] = y2[3];
        y2[3] = y2[4];

        init.state = y.state;
        init.arg = end;
        y2[4].arg = end;
        end = init.arg + parameters.step;
    }
    return res;
}
