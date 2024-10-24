#include <iostream>
#include <fstream>
#include <array>
#include <Eigen/Dense>
#include <cmath>

template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points, const int derOrder){
    //Матрица коэффициентов для A1, A2, A3, ... (если следовать терминологии семинара) для составления лин. системы
    Eigen::MatrixXd syst(N+1,N+1);
    Eigen::MatrixXd rightSide(N+1, 1); //Правый столбец для лин. системы

    //Инициализация правого столбца и первого столбца матрицы
    for(int i = 0; i < N+1; i++){
        if(i == derOrder){
            rightSide(i, 0) = 1;
            syst(i, 0) = 0;
        }
        else if(i == 0){
            rightSide(i, 0) = 0;
            syst(i, 0) = 1;
        }
        else{
            rightSide(i, 0) = 0;
            syst(i, 0) = 0;
        }
    } 

    //Проверка
    // std::cout << rightSide << std::endl; 

    //Заполнение матрицы коэфф-ов для лин. системы
    //Заполнение матрицы коэфф-ов для лин. системы
    double f = 1;
    for(int j = 0; j < N+1; j++){
        if(j == 0)
        {
            f = 1;
        }
        else
        {
            f *= j;
        }

        for(int i = 0; i < N; i++){
            syst(j, i+1) = (1/f) * pow(points[i], j);
        }
    } 

    //Решение линейной системы
    Eigen::MatrixXd solution = syst.colPivHouseholderQr().solve(rightSide);

    //Проверка
    std::cout << solution << std::endl; 
    
    //Перевод из матрицы в нужные типы переменных 
    DerivativeCoef<RealType, N> df;
    df.centralCoef = solution(0,0);
    for(int i = 0; i < N; i++){
        df.otherCoefs[i] = solution(i+1, 0);
    }
    
    return df;
}


int main()
{
    //Начальные данные
    const unsigned int N = 3;
    const int derOrder = 2;
    std::array<double, N> points{{-1, 1, 2}};
    
    std::array<double, 16> hs;
    for(int i = 0; i < hs.size(); i++){
        hs[i] = pow(10, -i);
    }

    //Расчет коэффициентов
    DerivativeCoef<double, N> df = calcDerivativeCoef<double, N>(points, derOrder);
    double x0 = 1; double D = 0;
    
    //Подготовка файла
    std::ofstream out("points_3_extra.txt");

    for(int j = 0; j < hs.size(); j++){
        D = exp(x0)*df.centralCoef;
        for(int i = 0; i < N; i++){
            D += df.otherCoefs[i]*exp(x0 + hs[j]*points[i]);
        }

        //Вывод в файл (степень 10, значение производной, разницу м/у посчитанным и теорией)
        out << j << " " << D/pow(hs[j], derOrder) << " " << exp(1) - D/pow(hs[j], derOrder) << std::endl;
    }
    return 0;
}
