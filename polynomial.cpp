#include "polynomial.h"
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
polynomial::polynomial() = default;
polynomial::polynomial(size_t degree, const coeff_type& value): coeffs_(degree + 1, value)
{trim();}
polynomial::polynomial(const container& coeffs): coeffs_(coeffs)
{trim();}
void polynomial::trim()
{
    while (!coeffs_.empty() && coeffs_.back() == 0)
    {coeffs_.pop_back();}
}
size_t polynomial::degree() const
{
    if (coeffs_.empty()){return 0;}
    size_t i = coeffs_.size();
    while (i > 0 && coeffs_[i - 1] == 0){--i;}
    return i == 0 ? 0 : i - 1;
}
const big_int& polynomial::operator[](size_t i) const
{return coeffs_.at(i);}//ат для границ возвр конст ссылку на коэф с индексом i
big_int& polynomial::operator[](size_t i)
{
    if (i >= coeffs_.size())
    {coeffs_.resize(i + 1, 0);}
    return coeffs_[i];
}
polynomial polynomial::operator+(const polynomial& other) const
{
    const size_t size = std::max(coeffs_.size(), other.coeffs_.size());
    container result(size, 0);
    for (size_t i = 0; i < coeffs_.size(); ++i)
    {result[i] += coeffs_[i];}
    for (size_t i = 0; i < other.coeffs_.size(); ++i)
    {result[i] += other.coeffs_[i];}
    return polynomial(result);
}
polynomial polynomial::operator-(const polynomial& other) const
{
    const size_t size = std::max(coeffs_.size(), other.coeffs_.size());
    container result(size, 0);
    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        result[i] += coeffs_[i];
    }
    for (size_t i = 0; i < other.coeffs_.size(); ++i)
    {
        result[i] -= other.coeffs_[i];//коэф второго полинома вычитаются из реза
    }
    return polynomial(result);
}
polynomial polynomial::operator*(const polynomial& other) const
{
    if (coeffs_.empty() || other.coeffs_.empty())
    {
        return polynomial();
    }
    container result(
        coeffs_.size() + other.coeffs_.size() - 1,
        0);
    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        for (size_t j = 0; j < other.coeffs_.size(); ++j)
        {
            result[i + j] +=
                coeffs_[i] * other.coeffs_[j];
        }
    }
    return polynomial(result);
}
polynomial polynomial::operator*(const coeff_type& scalar) const
{
    if (scalar == 0 || coeffs_.empty())
    {
        return polynomial();
    }
    container result(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        result[i] = coeffs_[i] * scalar;
    }
    return polynomial(result);
}
polynomial polynomial::mod(const coeff_type& mod) const
{
    if (mod <= 0)
    {
        throw std::invalid_argument(
            "polynomial::mod: invalid modulus");
    }
    container result(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        result[i] = coeffs_[i] % mod;
        if (result[i] < 0)
        {
            result[i] += mod;
        }
    }
    return polynomial(result);
}
polynomial polynomial::center_lift(const coeff_type& mod) const
{
    container result(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        result[i] = ::center_lift(coeffs_[i], mod);//для каждого коэф глобальная функция
    }
    return polynomial(result);
}
polynomial polynomial::mul_mod(
    const polynomial& other,
    const coeff_type& mod,
    size_t N) const
{
    if (N == 0)
    {
        throw std::invalid_argument(
            "mul_mod: N must be > 0");
    }
    if (mod <= 0)
    {
        throw std::invalid_argument(
            "mul_mod: modulus must be positive");
    }
    container result(N, 0);//вектор коэф произведения
    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        for (size_t j = 0; j < other.coeffs_.size(); ++j)
        {
            const size_t index = (i + j) % N;//индекс в рез полиноме
            big_int term =
                (coeffs_[i] * other.coeffs_[j]) % mod;
            if (term < 0)
            {
                term += mod;
            }
            result[index] = (result[index] + term) % mod;//Добавляет term к текущему значению в result[index] и снова берёт по модулю mod.
        }
    }
    return polynomial(result);
}
bool polynomial::operator==(const polynomial& other) const
{
    container a = coeffs_;
    container b = other.coeffs_;
    while (!a.empty() && a.back() == 0)
    {
        a.pop_back();
    }
    while (!b.empty() && b.back() == 0)
    {
        b.pop_back();
    }
    return a == b;
}
bool polynomial::operator!=(const polynomial& other) const
{
    return !(*this == other);
}
std::string polynomial::to_string() const
{
    std::ostringstream stream;
    if (coeffs_.empty())
    {
        return "0";
    }
    for (size_t i = 0; i < coeffs_.size(); ++i)
    {
        if (i > 0)
        {
            stream << " + ";
        }
        stream << coeffs_[i] << "x^" << i;
    }
    return stream.str();
}
polynomial polynomial::random(
    size_t degree,
    const std::vector<coeff_type>& coeff_set)
{
    if (coeff_set.empty())
    {
        throw std::invalid_argument("empty coeff set");
    }
    container coeffs(degree + 1);
    for (size_t i = 0; i < coeffs.size(); ++i)
    {
        const size_t index =
            random_bigint(0, coeff_set.size() - 1)
            .convert_to<size_t>();

        coeffs[i] = coeff_set[index];
    }
    return polynomial(coeffs);//возвр полином созданный из вектора коэф
}
void polynomial::resize(size_t N, const coeff_type& fill)
{
    coeffs_.resize(N, fill);
}
static std::pair<polynomial, polynomial> divmod_poly(//возвр частное остаток
    const polynomial& a,
    const polynomial& b,
    const big_int& mod)
{
    if (b.get_coeffs().empty())
    {
        throw std::invalid_argument(
            "Division by zero polynomial");
    }
    polynomial A = a.mod(mod);
    polynomial B = b.mod(mod);
    if (A.get_coeffs().empty() || A.degree() < B.degree())
    {
        return { polynomial(), A };//частное 0 остаток А
    }
    polynomial quotient;
    polynomial remainder = A;
    const big_int inverse_leading =
        mod_inverse(B.get_coeffs().back(), mod);//вычисляет обр элемент старшего коэф делителя B по модулю mod
    while (!remainder.get_coeffs().empty() &&
        remainder.degree() >= B.degree())//пока остаток не ноль и степень B не меньше степени делителя
    {
        const size_t shift =
            remainder.degree() - B.degree();//сдвиг разность степенией и делителя
        big_int coefficient =
            (remainder.get_coeffs().back() * inverse_leading) % mod;//вычисляет коэф частного: старший коэф остатка,умножить на обратный старшего коэф делителя по модулю
        if (coefficient < 0)
        {
            coefficient += mod;
        }
        polynomial::container term_coeffs(shift + 1, 0);
        term_coeffs[shift] = coefficient;
        quotient =
            (quotient + polynomial(term_coeffs)).mod(mod);//создвкт полином,добавляет к частному м по модулю
        polynomial subtraction = B * coefficient;//умножает делитель В на коэф,затем вставляет шифт нулей в начало вектора коэф, чтобы сдвинуть полином на шифт степений вверх
        polynomial::container subtraction_coeffs =
            subtraction.get_coeffs();
        subtraction_coeffs.insert(
            subtraction_coeffs.begin(),
            shift,
            big_int(0));
        subtraction.set_coeffs(subtraction_coeffs);
        remainder = (remainder - subtraction).mod(mod);
    }
    return { quotient, remainder };
}
polynomial polynomial::inverse_mod(//нахождение обратного в кольце
    const coeff_type& mod,
    size_t N) const
{
    if (mod <= 1)
    {
        throw std::invalid_argument(
            "inverse_mod: invalid modulus");
    }
    if (N == 0)
    {
        throw std::invalid_argument(
            "inverse_mod: N must be > 0");
    }
    std::vector<std::vector<big_int>> matrix(
        N,
        std::vector<big_int>(N + 1, 0));//Создаёт матрицу размера N x (N+1) для решения системы линейных уравнений методом Гаусса. Матрица будет использоваться для поиска обратного полинома: мы решаем систему a * x ≡ 1 (mod x^N - 1).
    for (size_t row = 0; row < N; ++row)
    {
        for (size_t col = 0; col < N; ++col)
        {
            const size_t index = (row + N - col) % N;//индекс коэф полинома а который соотв матричному элементу
            if (index < coeffs_.size())
            {
                matrix[row][col] = coeffs_[index] % mod;
                if (matrix[row][col] < 0)
                {
                    matrix[row][col] += mod;
                }
            }
        }
        matrix[row][N] = (row == 0) ? 1 : 0;//в последний столб правая часть для первой строки 1(коэф при х^0 для остальных 0)
    }
    size_t pivot_row = 0;//тек строка где ищем ведущий
    for (size_t col = 0;
        col < N && pivot_row < N;
        ++col)
    {
        size_t pivot = pivot_row;//начинаем поиск ведущ эл в тек столбце col начиная с тек строки
        while (pivot < N && matrix[pivot][col] == 0)//идем вниз по столбцу пока не найдем строку где элемент не равен нулю
        {
            ++pivot;
        }
        if (pivot == N)
        {
            continue;
        }
        std::swap(matrix[pivot], matrix[pivot_row]);//Меняем местами строку с найденным ведущим элементом (pivot) и текущую строку (pivot_row), чтобы переместить её наверх.
        const big_int inverse =
            mod_inverse(matrix[pivot_row][col], mod);//вычисляем мультипл обрат для ведущ знач (чтоыб вед элемент равным 1)
        for (size_t j = col; j <= N; ++j)
        {
            matrix[pivot_row][j] =
                (matrix[pivot_row][j] * inverse) % mod;//делаем ведущий равным 1 приводим к единичной матрице
            if (matrix[pivot_row][j] < 0)
            {
                matrix[pivot_row][j] += mod;
            }
        }
        for (size_t row = 0; row < N; ++row)
        {
            if (row == pivot_row || matrix[row][col] == 0)
            {
                continue;
            }
            const big_int factor = matrix[row][col];//вычисляем коэф и вычитаем из строки ведущую умноженную на фактор, чтобы обнулить элемент в стобце col  
            for (size_t j = col; j <= N; ++j)
            {
                matrix[row][j] =
                    (matrix[row][j] -
                        factor * matrix[pivot_row][j]) % mod;
                if (matrix[row][j] < 0)
                {
                    matrix[row][j] += mod;
                }
            }
        }
        ++pivot_row;
    }
    if (pivot_row < N)//не удалось найти ведущий обратного не сущ
    {
        throw std::runtime_error("Polynomial not invertible");
    }
    container result(N, 0);//содержит коэф обрат полинома
    for (size_t i = 0; i < N; ++i)
    {
        result[i] = matrix[i][N];//извлекаем послед столбец матрицы, после приведения последний содержит решение, т.е коэф иск полинома
    }
    return polynomial(result);
}
