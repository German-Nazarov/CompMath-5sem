#include<iostream>
#include<cmath>

/**
    Решает уравнение Кеплера методом Ньютона
    * ecc - эксцентриситет, принадлежит (0, 1)
    * meanAnomaly - средняя аномалия, М (радианы)
    * maxIter - максимальное количество итераций
    * tol - точность, с которой нужно отыскать решение

    Рекомендуемое поведение. Если решение не нашлось за maxIter итераций - выбрасывать исключение.
    Если приближения к решению между итерациями меняются не более, чем на tol, то решение достигнуто.

    E - e*sin(E) = M -> Метод Ньютона: f(E) = E - e*sin(E) - M = 0
    f'(E) = 1 - e*cos(E)
    E_i+1 = E_i - f(E_i)/(f'(E_i)) 
**/
double keplerSolver(double ecc, double meanAnomaly, unsigned int maxIter, double tol) {
    double E = meanAnomaly; // Эксцентрическая аномалия
    for(unsigned int i = 0; i < maxIter; i++) {
        double E_next = E - (E - ecc*sin(E) - meanAnomaly)/(1 - ecc*cos(E));
        if(fabs(E_next - E) <= tol) {
            return E;
        }
        else {
            E = E_next;
        }
    }
    return E;
}

/**
    Метод простой итерации 
    f(x) = 0  <=>  x_i+1 = x_i - tau*f(x_i)
    f(x) = x^2 + tg^2(x) - 1
    0 <= f'(x) <= 2(1 + 1/cos^2(1)) для x из [0, 1]
    tau_opt = 2/(M+m) = cos^2(1)/(1 + cos^2(1))
**/
template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename Callable, typename RealType>
decltype(auto) solve(   
    const Callable& func,                                             // функция F
    const RealType& tau,                                              // шаг тау
    const typename ArgumentGetter<Callable>::Argument& initialGuess,  // начальное приближение
    const unsigned int nIteration                                     // количество итераций
                    ) 
{
    typename ArgumentGetter<Callable>::Argument x = initialGuess;
    for(unsigned int i = 0; i < nIteration; i++) {
        x = x - tau*func(x);
    }
    return x;
}