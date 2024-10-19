#include <iostream>
#include <array>
#include <vector>
#include <type_traits>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>

template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& o, const std::array<T, N>& arr)
{
    copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
    return o;
}

template <typename S>
std::ostream& operator<<(std::ostream& os,
                    const std::vector<S>& vector) {
  
    // Printing all the elements using <<
    for (auto i : vector) 
        os << i << " " << std::endl;
    return os;
}

/** класс для работы с трехдиагональной матрицей **/
template<typename Type>
class ThreeDiagonalMatrix {
    /*** Здесь какие-то поля и методы ***/
    private:
    std::vector<std::array<Type, 3>> data;

    public:
    void fill_line(std::array<Type, 3> &arr)
    {
        data.push_back(arr);
    }

    std::vector<std::array<Type, 3>> get_data() const
    {
        return data;
    }
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

/** Функция для решения методом прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve( const ThreeDiagonalMatrix<mType>& matrix,
                                            const std::vector<cType>& column)
{
    std::vector<std::array<mType, 3>> matrix_data = matrix.get_data();

    std::size_t M = column.size();

    // Прямой ход
    std::vector<DivisType<cType, mType>> p;
    std::vector<DivisType<cType, mType>> q;

    p.push_back(-1 * (matrix_data[0][2] / matrix_data[0][1]));
    q.push_back(column[0]/matrix_data[0][1]);

    for(std::size_t i = 1; i < M-1; ++i)
    {
        p.push_back((-1*matrix_data[i][2])/(matrix_data[i][0] * p[i-1]  + matrix_data[i][1]));
        q.push_back((column[i] - matrix_data[i][0]*q[i-1])/(matrix_data[i][0] * p[i-1]  + matrix_data[i][1]));
    }

    // Обратный ход -> определение коэфф-ов ci для интерполянта: si(x) = ai + bi(x-xi) + ci/2(x-xi)^2 _ di/6(x-xi)^3
    std::vector<DivisType<cType, mType>> c;

    c.push_back((column[M-1] - matrix_data[M-1][0]*q[M-2])/(matrix_data[M-1][0] * p[M-2]  + matrix_data[M-1][1]));
    for(long long int i = M-2; i >= 0; i--)
    {
        c.push_back(p[i] * c[M-i-2] + q[i]);
    }

    // c.push_back(p[0] * c[M-2] + q[0]);
    c.push_back(0);

    std::reverse(c.begin(), c.end());

    c.push_back(0);

    return (c);
}


/**
* xType - тип аргумента x.
* yType - тип значения функции y
*/
template<typename xType, typename yType>
class CubicSpline {
    // Трёхдиагональная матрица и столбец
    ThreeDiagonalMatrix<xType> matrix;
    std::vector<DivisType<xType, yType>> column;
    std::vector<xType> points;

    // Коэфф-ы интерполянта
    std::vector<DivisType<xType, yType>> a;
    std::vector<DivisType<xType, yType>> b;
    std::vector<DivisType<xType, yType>> c;
    std::vector<DivisType<xType, yType>> d;

    public:
    CubicSpline(const std::vector<xType> &points, const std::vector<yType>& values)
    {
        this->points = points;

        // Сразу же сформируем трёхдиагональную матрицу и столбец
        for(std::size_t i = 1; i < points.size()-1; ++i)
        {   
            xType h1 = points[i] - points[i-1];
            xType h2 = points[i+1] - points[i];

            // Разделенная разность для столбца
            DivisType<xType, yType> div_diff = 
                    ((values[i+1] - values[i])/h2 - (values[i] - values[i-1])/h1)/(h2+h1);
            this->column.push_back(6*div_diff);

            // Ненулевые элементы строчки трёхдиагональной матрицы
            if(i == 1)
            {
                std::array<xType, 3> line = {-500, 2, h2/(h1+h2)};
                matrix.fill_line(line);
            }
            else if(i == points.size()-2)
            {
                std::array<xType, 3> line = {h1/(h1+h2), 2, -500};
                matrix.fill_line(line);
            }
            else
            {
                std::array<xType, 3> line = {h1/(h1+h2), 2, h2/(h1+h2)};
                matrix.fill_line(line);
            }
        }

        

        // Коэфф-ы ai, bi  и di для i-го интерполянта si
        a = values;
        c = solve(matrix, column);

        b.push_back(0);
        b.push_back((c[1] * (points[1] - points[0]))/3 + (values[1] - values[0])/(points[1] - points[0]));
        for(std::size_t k = 1; k < c.size()-1; k++)
        {   
            b.push_back((c[k+1] * (points[k+1] - points[k]))/3 +
                 (values[k+1] - values[k])/(points[k+1] - points[k]) + (c[k]/6)*(points[k+1] - points[k]));
        }

        d.push_back(0);
        d.push_back(c[1]/(points[1] - points[0]));
        for(std::size_t k = 1; k < c.size()-1; k ++)
        {
            d.push_back((c[k+1] - c[k])/(points[k+1] - points[k]));
        }
    }
                        
    yType interpolate(const xType& x) const
    {
        std::size_t index = 0;

        if(x > points[0])
        {
            for(std::size_t i = 0; i < a.size()-1; i++)
            {
                if(x > points[i] && x <= points[i+1])
                {
                    index = i+1;
                }
            }
        }
        else if(x == points[0])
        {
            index = 0;
        }
        // std::cout << a[index] << std::endl;
        yType s = a[index] + b[index]*(x - points[index]) + 0.5*c[index]*pow((x - points[index]), 2) +
                                (1/6)*d[index]*pow((x - points[index]), 3);

        return s;
    }
};

// Формирование равномерных узлов
template<typename xType>
std::vector<xType> uniform_points(const xType& lead, const xType& end, unsigned int N) // равномерное расположение
{
    std::vector<xType> p;
    p.push_back(lead);
    for(int i = 1; i < N; i++)
    {
        p.push_back(p[i-1] + (end - lead)/(N-1));
    }

    return p;
};

int main()
{
    std::ofstream out("lab-3-spline_data.txt");

    std::array<unsigned int, 6> Ns = {5, 10, 20, 40, 80, 160};

    for(int j = 0; j < 6; j++)
    {
        std::vector<double> points = uniform_points<double>(0, 10, Ns[j]);
        std::vector<double> values;
        for(int l = 0; l < Ns[j]; l++)
        {   
            values.push_back(exp(points[l]));
        }

        CubicSpline<double, double> spline(points, values);

        double err = 0;
        double x = 2 * 10./1000;
        for(int l = 2; l < 900; l++)
        {
            x = x + 10./1000;
            // std::cout << "XXXXX: " << x << std::endl;
            double res = spline.interpolate(x);
            // std::cout << Ns[j] << " " << res << " " << exp(x) << std::endl;
            if(err <= std::fabs(res - exp(x)))
            {
                err = std::fabs(res - exp(x));
            }
        }
        out << Ns[j] << "\t" << err << std::endl; 
    }
}