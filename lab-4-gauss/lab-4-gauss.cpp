#include <iostream>
#include <array>
#include <type_traits>
#include <cmath>
#include <fstream>

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

/* Функция производит интегрирование на одном отрезке */
template<typename Callable, std::size_t N>
decltype(auto) integrate(   
    const Callable& func,  // Интегрируемая функция
    const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
    const typename ArgumentGetter<Callable>::Argument& end  // конец отрезка
                        )
    {
        // Узлы и веса Гаусса (два варианта на выбор)

        std::array<typename ArgumentGetter<Callable>::Argument, N> nodes;
        std::array<typename ArgumentGetter<Callable>::Argument, N> weights;

        if(N == 4)
        {
            nodes[0] = 0.86114;
            nodes[1] = 0.33998;
            nodes[2] = -0.33998;
            nodes[3] = -0.86114;

            weights[0] = 0.65216;
            weights[1] = 0.34786;
            weights[2] = 0.34786;
            weights[3] = 0.65216;
        }
        else if(N == 6)
        {
            nodes[0] = -0.23862;
            nodes[1] = 0.23862;
            nodes[2] = -0.66121;
            nodes[3] = 0.66121;
            nodes[4] = -0.93247;
            nodes[5] = 0.93247;

            weights[0] = 0.46791;
            weights[1] = 0.46791;
            weights[2] = 0.36076;
            weights[3] = 0.36076;
            weights[4] = 0.17132;
            weights[5] = 0.17132;
        }
        
        // Подсчет по квадратурной формуле (b-a)*sum(wk*f(xk))
        Dif<typename ArgumentGetter<Callable>::Argument> res = 0.0;
        for(int k = 0; k < N; k++)
        {
            // Правильный x из nodes[k]
            Dif<typename ArgumentGetter<Callable>::Argument> xk = 
                                        (start + end)/2 + ((end - start)/2)*nodes[k];

            // Прибавляем к res
            res += weights[k] * func(xk);
        }
        res *= (end - start)/2;

        return res;
    }

/* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
template<typename Callable, std::size_t N>
decltype(auto) integrate(   
    const Callable& func,  // Интегрируемая функция
    const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
    const typename ArgumentGetter<Callable>::Argument& end,  // конец отрезка
    const Dif<typename ArgumentGetter<Callable>::Argument>& dx  // Длина подотрезка
                        )
    {
        // Узлы и веса Гаусса (два варианта на выбор)

        std::array<typename ArgumentGetter<Callable>::Argument, N> nodes;
        std::array<typename ArgumentGetter<Callable>::Argument, N> weights;

        if(N == 4)
        {
            nodes[0] = 0.8611363;
            nodes[1] = 0.3399810;
            nodes[2] = -0.3399810;
            nodes[3] = -0.8611363;

            weights[0] = 0.6521451;
            weights[1] = 0.3478548;
            weights[2] = 0.3478548;
            weights[3] = 0.6521451;
        }
        else if(N == 6)
        {
            nodes[0] = -0.2386142;
            nodes[1] = 0.2386142;
            nodes[2] = -0.6612094;
            nodes[3] = 0.6612094;
            nodes[4] = -0.9324700;
            nodes[5] = 0.9324700;

            weights[0] = 0.4679140;
            weights[1] = 0.4679140;
            weights[2] = 0.3607616;
            weights[3] = 0.3607616;
            weights[4] = 0.1713245;
            weights[5] = 0.1713245;
        }
        
        // Подсчет по квадратурной формуле (b-a)*sum(wk*f(xk))
        Dif<typename ArgumentGetter<Callable>::Argument> res = 0;

        for(double j = start; j <= end-dx; j+=dx)
        {
            Dif<typename ArgumentGetter<Callable>::Argument> start_j = j;
            Dif<typename ArgumentGetter<Callable>::Argument> end_j = j+dx;
            
            Dif<typename ArgumentGetter<Callable>::Argument> res_j = 0;

            for(int k = 0; k < N; k++)
            {
                // Правильный x из nodes[k]
                Dif<typename ArgumentGetter<Callable>::Argument> xk = 
                                            (start_j + end_j)/2 + ((end_j - start_j)/2)*nodes[k];

                // Прибавляем к res
                res_j += weights[k] * func(xk);
            }
            res_j *= (end_j - start_j)/2;
            res += res_j;
        }
        // std::cout << res << std::endl;
        return res;
    }

int main()
{
    // Подготовка файла
    std::ofstream out("lab-4-data-6.txt");

    for(int i = 0; i < 15; i++)
    {
        double h = 10/(pow(2, i));

        // std::cout << h << std::endl;

        out << h << "\t" << fabs(1.83907152908 - integrate<double(double), 6>(sin, 0., 10., h)) << std::endl; 
    }

}