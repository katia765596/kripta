#include "polynomial.h"
#include <algorithm>
#include <stdexcept>
polynomial::polynomial()
{
}
polynomial::polynomial(uint64_t value)
{
    if (value != 0)
    {
        data.push_back(value);
    }
}
polynomial::polynomial(const std::vector<uint64_t>& value)
    : data(value)
{
    normalize();
}
polynomial polynomial::from_binary(const std::string& value)
{
    polynomial result;
    if (value.empty())
    {
        return result;
    }
    for (size_t i = 0; i < value.size(); ++i)
    {
        char current = value[value.size() - 1 - i];//посл симвл соотв младшему биту
        if (current == '1')
        {
            result.set_bit(i, true);
        }
        else if (current != '0')
        {
            throw std::invalid_argument("invalid binary polynomial");
        }
    }
    return result;
}
polynomial polynomial::one()
{
    return polynomial(1);
}
polynomial polynomial::x()
{
    return polynomial(2);//x^1 - бит 1 установлен, соотв знач 2
}
bool polynomial::is_zero() const
{
    return data.empty();
}
bool polynomial::is_one() const
{
    return data.size() == 1 && data[0] == 1;
}
size_t polynomial::degree() const
{
    if (data.empty())
    {
        return 0;
    }
    uint64_t value = data.back();//последнее слово в векторе, которое содержит старш бит полинома
    size_t bits = 0;
    while (value != 0)
    {
        value >>= 1;//цикл сдвигает битс вправо пока не станет нуль и подсчитывает кол-во бит в слове
        ++bits;
    }
    return (data.size() - 1) * 64 + bits - 1;//к кол-ву бит в пред словах прибавить битс-1 (индекс последнего установленного бита в посл слове)
}
bool polynomial::get_bit(size_t index) const
{
    size_t word = index / 64;
    size_t bit = index % 64;
    if (word >= data.size())
    {
        return false;
    }
    return ((data[word] >> bit) & 1ULL) != 0;//извл бит с помощью сдвига вправо и маски возвр тру если бит равен 1
}
void polynomial::set_bit(size_t index, bool value)
{
    size_t word = index / 64;
    size_t bit = index % 64;
    if (word >= data.size())
    {
        data.resize(word + 1, 0);
    }
    uint64_t mask = 1ULL << bit;
    if (value)
    {
        data[word] |= mask;//установить бит
    }
    else
    {
        data[word] &= ~mask;//сбросить бит
    }
    normalize();
}
void polynomial::normalize()
{
    while (!data.empty() && data.back() == 0)
    {
        data.pop_back();
    }
}
uint64_t polynomial::to_uint64() const
{
    if (data.size() > 1)
    {
        throw std::overflow_error("polynomial does not fit uint64");
    }
    return data.empty() ? 0 : data[0];
}
std::string polynomial::to_binary() const
{
    if (is_zero())
    {
        return "0";
    }
    std::string result;
    size_t current_degree = degree();
    for (size_t i = current_degree + 1; i > 0; --i)
    {
        result.push_back(get_bit(i - 1) ? '1' : '0');
    }
    return result;
}
bool polynomial::operator==(const polynomial& other) const
{
    return data == other.data;
}
bool polynomial::operator!=(const polynomial& other) const
{
    return !(*this == other);
}
polynomial polynomial::add(const polynomial& other) const
{
    size_t size = std::max(data.size(), other.data.size());
    polynomial result;
    result.data.resize(size, 0);
    for (size_t i = 0; i < size; ++i)
    {
        uint64_t first = i < data.size() ? data[i] : 0;
        uint64_t second = i < other.data.size() ? other.data[i] : 0;
        result.data[i] = first ^ second;
    }
    result.normalize();
    return result;
}
void polynomial::xor_shifted(const polynomial& other, size_t shift)
{//this ^= other * x^shift.
    if (other.is_zero())
    {
        return;
    }
    size_t word_shift = shift / 64;
    size_t bit_shift = shift % 64;
    size_t required_size = word_shift + other.data.size();//мин размер вектора для хранения рез
    if (bit_shift != 0)
    {
        ++required_size;
    }
    if (data.size() < required_size)
    {
        data.resize(required_size, 0);
    }
    for (size_t i = 0; i < other.data.size(); ++i)
    {
        data[word_shift + i] ^= other.data[i] << bit_shift;//для каждого i выполняем xor тек слова со знач офердата и сдвиг байтов внутри влево на бит_шифт бит

        if (bit_shift != 0)
        {
            data[word_shift + i + 1] ^= other.data[i] >> (64 - bit_shift);//часть битов выходящая за пределы 64-битного слова,перенесена в след слово (сдвиг с переносом)
        }
    }
    normalize();
}
polynomial polynomial::multiply(const polynomial& other) const
{
    polynomial result;
    if (is_zero() || other.is_zero())
    {
        return result;
    }
    const polynomial* first = this;
    const polynomial* second = &other;
    if (first->degree() < second->degree())
    {
        std::swap(first, second);//полином с большей степенью -первый
    }
    for (size_t word = 0; word < second->data.size(); ++word)
    {//внешний цикл по словам полинома секонд(меньшего)
        uint64_t value = second->data[word];//тек слово
        while (value != 0)
        {
            unsigned int bit = 0;
            uint64_t temp = value;
            while ((temp & 1ULL) == 0)//находим поз младшего установленного бита в value
            {
                temp >>= 1;//сдвигаем вправо пока мл бит не станет 1
                ++bit;//кол-во сдвигов
            }
            result.xor_shifted(*first, word * 64 + bit);//соотв умн на x^(word*64 + bit)
            value &= value - 1;//сбрасывает мл установл бит в 
        }
    }
    return result;
}
polynomial polynomial::mod(const polynomial& modulus) const
{
    if (modulus.is_zero())
    {
        throw std::invalid_argument("zero modulus");
    }
    polynomial result = *this;
    size_t modulus_degree = modulus.degree();
    while (!result.is_zero() && result.degree() >= modulus_degree)
    {
        result.xor_shifted(modulus, result.degree() - modulus_degree);
    }//хор резалт с модулем сдвинутым на разность степений, после хор степень резалт уменьшится, цикл повт пока не получим остаток
    return result;
}
polynomial polynomial::multiply_mod(
    const polynomial& other,
    const polynomial& modulus
) const
{
    return multiply(other).mod(modulus);
}
polynomial polynomial::power_mod(
    uint64_t exponent,
    const polynomial& modulus
) const
{
    if (modulus.is_zero())
    {
        throw std::invalid_argument("zero modulus");
    }
    polynomial result = polynomial::one();
    polynomial base = mod(modulus);
    while (exponent != 0)
    {
        if (exponent & 1ULL)//мл бит 1
        {
            result = result.multiply_mod(base, modulus);
        }
        exponent >>= 1;
        if (exponent != 0)
        {
            base = base.multiply_mod(base, modulus);
        }
    }
    return result;
}
polynomial polynomial::gcd(const polynomial& other) const
{
    polynomial first = *this;
    polynomial second = other;
    while (!second.is_zero())
    {
        polynomial remainder = first.mod(second);
        first = second;
        second = remainder;
    }
    return first;
}
polynomial polynomial::divide_remainder(
    const polynomial& divisor,
    polynomial& remainder
) const
{
    if (divisor.is_zero())
    {
        throw std::invalid_argument("division by zero polynomial");
    }
    polynomial quotient;
    remainder = *this;
    size_t divisor_degree = divisor.degree();//степень делителя
    while (!remainder.is_zero() && remainder.degree() >= divisor_degree)
    {
        size_t shift = remainder.degree() - divisor_degree;
        quotient.set_bit(shift, true);//устанавливаем бит в чатсном на поз шифт
        remainder.xor_shifted(divisor, shift);
    }
    return quotient;
}
polynomial polynomial::inverse(const polynomial& modulus) const
{
    if (is_zero())
    {
        throw std::invalid_argument("zero polynomial has no inverse");
    }
    if (modulus.is_zero())
    {
        throw std::invalid_argument("zero modulus");
    }
    polynomial r0 = modulus;
    polynomial r1 = mod(modulus);
    polynomial t0;
    polynomial t1 = polynomial::one();
    while (!r1.is_zero())
    {
        polynomial remainder;
        polynomial quotient = r0.divide_remainder(r1, remainder);
        polynomial next_t = t0.add(quotient.multiply(t1));
        r0 = r1;
        r1 = remainder;
        t0 = t1;
        t1 = next_t;
    }
    if (!r0.is_one())
    {
        throw std::runtime_error("inverse does not exist");
    }
    return t0.mod(modulus);
}
bool polynomial::is_irreducible() const//полином степени n неприводим, если он взаимно прост с (x^2^i)-x для всех i от 1 до n/2 или для всех простых делителей n
{
    if (is_zero() || is_one())
    {
        return false;
    }
    size_t n = degree();
    if (n == 0)
    {
        return false;
    }
    std::vector<size_t> factors;
    size_t value = n;
    for (size_t p = 2; p * p <= value; ++p)
    {
        if (value % p == 0)
        {
            factors.push_back(p);
            while (value % p == 0)
            {
                value /= p;
            }
        }
    }
    if (value > 1)
    {
        factors.push_back(value);
    }
    polynomial base = polynomial::x();
    for (size_t q : factors)
    {
        size_t count = n / q;
        polynomial current = base;
        for (size_t i = 0; i < count; ++i)
        {//х послед возводится в квадрат каунт раз
            current = current.multiply_mod(current, *this);//х возводится в квадрат умн на себя по модулю f, x^(2^count) mod f 
        }
        if (!current.add(base).gcd(*this).is_one())//x^(2^(n/q)) + x затем нод если не равен 1 то полином приводоим
        {
            return false;
        }
    }
    polynomial current = base;//вторая ч проверки вычисление x^(2^n) mod f
    for (size_t i = 0; i < n; ++i)
    {
        current = current.multiply_mod(current, *this);
    }
    return current == base;//x^(2^n) ≡ x (mod f)
}
finite_field::finite_field(const polynomial& modulus_value)
    : modulus(modulus_value)
{
    if (modulus.is_zero() || !modulus.is_irreducible())
    {
        throw std::invalid_argument("modulus must be irreducible");
    }
}
polynomial finite_field::add(
    const polynomial& first,
    const polynomial& second
) const
{
    return first.add(second);
}
polynomial finite_field::multiply(
    const polynomial& first,
    const polynomial& second
) const
{
    return first.multiply_mod(second, modulus);
}
polynomial finite_field::inverse(const polynomial& value) const
{
    return value.inverse(modulus);
}
polynomial finite_field::power(
    const polynomial& value,
    uint64_t exponent
) const
{
    return value.power_mod(exponent, modulus);
}
const polynomial& finite_field::get_modulus() const
{
    return modulus;
}
