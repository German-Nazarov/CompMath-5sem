#include <array>
#include <iostream>
#include <iterator>
#include <cmath>
#include <fstream>

template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& o, const std::array<T, N>& arr)
{
    copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
    return o;
}

/**
* xType - тип аргумента x.
* yType - тип значения функции y
* N - количество точек для интерполяции
*
* Рекомедую обратить внимание. Разность (xType - xType) не обязана быть типом xType
*/
template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {

    std::array<xType, N> points;
    std::array<yType, N> values;

    public:
    NewtonInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values)
    {   
        this->points = points;
        this->values = values;
    }
    
    yType interpolate(const xType& x) const
    {
        // Инициализация и заполнение массива для разделенных разностей, который в начале совпадает с values
        std::array<yType, N> div_diff = values;

        for(int j = 0; j < N; j++)
        {
            // std::cout << "divided differences (step " << j+1 << "): " << div_diff << std::endl;
            yType old_div = div_diff[j];
            for(int i = j+1; i < N; i++)
            {   
                yType temp = div_diff[i];
                div_diff[i] = (div_diff[i] - old_div)/(points[i] - points[i-(j+1)]);
                old_div = temp;
            }
        }
        // std::cout << std::endl;
        // std::cout << "divided differences (step " << N << "): " << div_diff << std::endl; // проверка разделенных разностей

        // Вычисление значения полинома Ньютона в точке x при помощи схемы Горнера
        yType D = div_diff[0]; // инициализация

        for(int k = 1; k < N; k++)
        {
            // Отдельный подсчет произведения П(x-xj)
            auto product = x - points[0];
            for(int j = 1; j < k; j++)
            {
                product *= x - points[j];
            }

            D += div_diff[k] * product;
        }

        return D;
    }


};

//Функции для формирования равномерных узлов или узлов Чебышева из заданного отрезка

template<typename xType, unsigned int N>
std::array<xType, N> uniform_points(const xType& lead, const xType& end) // равномерное расположение
{
    std::array<xType, N> p;
    p[0] = lead;
    for(int i = 1; i < N; i++)
    {
        p[i] = p[i-1] + (end - lead)/(N-1);
    }

    return p;
};

template<typename xType, unsigned int N>
std::array<xType, N> chebyshev_points(const xType& lead, const xType& end) // чебышевское расположение
{
    std::array<xType, N> p;
    for(int i = 0; i < N; i++)
    {
        p[i] = (end + lead)/2 + ((end - lead)/2)*cos((M_PI*(2*i + 1))/(2*(N-1) + 2));
    }

    return p;
};

int main()
{   
    // Выбор N
    const unsigned int N = 5;

    // Подготовка файла
    std::ofstream out("lab-2-interp_data-n5-u.txt");

    // Равномерное расположение узлов
    for(int i = 0; i <= 5; i++)
    {   
        double lead = 0.0;
        double end = 2.0/(pow(2, i));

        std::array<double, N> points = uniform_points<double, N>(lead, end);
        
        std::array<double, N> values;
        for(int j = 0; j < N; j++)
        {   
            values[j] = exp(points[j]);
        }

        NewtonInterpolator<double, double, N> NewtonInterp(points, values);

        // Ошибка интерполяции
        double err = 0.0;

        for(int k = 0; k < 1000; k++)
        {
            double x = (end - lead) / (k+1);
            double difference = std::fabs(exp(x) - NewtonInterp.interpolate(x));
            if(err <= difference)
            {
                err = difference;
            }
        }

        out << end-lead << "  " << err << std::endl; 
    }

    // Подготовка файла
    std::ofstream out_ch("lab-2-interp_data-n5-ch.txt");

    // Чебышевское расположение узлов
    for(int i = 0; i <= 5; i++)
    {   
        double lead = 0.0;
        double end = 2.0/(pow(2, i));

        std::array<double, N> points = chebyshev_points<double, N>(lead, end);
        
        std::array<double, N> values;
        for(int j = 0; j < N; j++)
        {   
            values[j] = exp(points[j]);
        }

        NewtonInterpolator<double, double, N> NewtonInterp(points, values);

        // Ошибка интерполяции
        double err = 0.0;

        for(int k = 0; k < 1000; k++)
        {
            double x = (end - lead) / (k+1);
            double difference = std::fabs(exp(x) - NewtonInterp.interpolate(x));
            if(err <= difference)
            {
                err = difference;
            }
        }

        out_ch << end-lead << "  " << err << std::endl; 
    }
}