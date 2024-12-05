#include<iostream>
#include<array>
#include<vector>
#include<Eigen/Dense>
#include <cmath>

// Таблица Бутчера для метода Рунге-Кутты 4 порядка
struct RK4Table{
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = { {{0, 0, 0, 0},
     {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}} };
    static constexpr std::array<double, stages> cColumn = {0, 0.5, 0.5, 1};
    static constexpr std::array<double, stages> bString = {1.0/6, 1.0/3, 1.0/3, 1.0/6};
};

//Реализация класса правой части дифференциального уравнения. То есть класс f(t, y) для y' = f(t,y)
class Oscillator {
    
public:
    
    static constexpr unsigned int dim = 2;  // размерность задачи
    
    using Argument = double;  // тип аргумента, тип t
    
    using State = Eigen::Vector<double, dim>;  // состояние
    
    struct StateAndArg{
        State state;
        Argument arg;
    };
    
    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.state(1), -stateAndArg.state(0)};
    } 
};

class Linear {

    public:

    static constexpr unsigned int dim = 1;

    using Argument = double;

    using State = Eigen::Vector<double, dim>;

    struct StateAndArg{
        State state;
        Argument arg;
    };

    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        // std::cout << pow(stateAndArg.arg, 3) << " " << stateAndArg.arg << std::endl;
        return Eigen::Vector<double, dim>{pow(stateAndArg.arg, 3)};
    } 
};

//Сигнатура для метода интегрирования:
template<typename Table, typename RHS>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrateRK(
    const typename RHS::StateAndArg& initialState,
    const typename RHS::Argument& endTime,
    double step,
    const RHS& rhs
) {
    const unsigned int dim = RHS::dim;
    const unsigned int stages = Table::stages;
    std::array<Eigen::Vector<double, dim>, stages> k;

    typename RHS::StateAndArg currState = initialState;
    std::vector<typename RHS::StateAndArg> res;
    res.push_back(initialState);

    for(double t = initialState.arg; t < endTime; t += step) {
        Eigen::Vector<double, dim> w;
        for(int m = 0; m < dim; m++) {
                w[m] = 0;
        }

        for(int j = 0; j < stages; j++) {
            Eigen::Vector<double, dim> summ;
            for(int m = 0; m < dim; m++) {
                summ[m] = 0;
                for(int l = 0; l < j; l++) {
                    summ[m] += Table::table[j][l]*k[l][m];
                }
            }
            
            typename RHS::StateAndArg tempState = currState;
            tempState.arg = t + Table::cColumn[j]*step;
            Eigen::Vector<double, dim> tmp;
            for(int m = 0; m < dim; m++) {
                tmp[m] = currState.state[m] + step*summ[m];
            }
            tempState.state = tmp;
            k[j] = rhs.calc(tempState);

            for(int m = 0; m < dim; m++) {
                w[m] += step * k[j][m] * Table::bString[j];
            }
        }
        typename RHS::StateAndArg nextState;
        for(int m = 0; m < dim; m++) {
            nextState.state[m] = currState.state[m] + w[m];
        }
        nextState.arg = t + step;
        // std::cout << "CHECK " << t << " " << nextState.arg << " " << currState.state[0] << " " << nextState.state[0] <<
        // " " <<  w[0] << std::endl; 
        currState = nextState;
        res.push_back(currState);
        // std::cout << "CHECK 2 " << res.back().arg << " " << res.back().state[0] << std::endl;
    }
    return res;
}