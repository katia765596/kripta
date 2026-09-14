#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <cassert>
using boost::multiprecision::cpp_int;
using boost::multiprecision::pow;
const double pi = 3.14159265358979323846;
typedef cpp_int big_int;
big_int mod_pow(big_int base, big_int exp, big_int mod) {
    if (mod == 0) throw std::invalid_argument("mod cannot be zero");
    big_int result = 1;//нейтр эл умн
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = (result * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;//сдвиг показателя вправо на 1 бит (делим на 2)
    }
    return result;
}
big_int gcd_euclid(big_int a, big_int b) {
    if (a < 0) a = -a;//абсл знач
    if (b < 0) b = -b;
    while (b != 0) {
        big_int r = a % b;
        a = b;
        b = r;
    }
    return a;
}
struct extended_gcd_result {
    big_int gcd;
    big_int x;
    big_int y;
};
extended_gcd_result extended_gcd(big_int a, big_int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    big_int x0 = 1, x1 = 0;
    big_int y0 = 0, y1 = 1;
    big_int a0 = a, b0 = b;
    while (b0 != 0) {
        big_int q = a0 / b0;
        big_int r = a0 - q * b0;
        big_int x = x0 - q * x1;
        big_int y = y0 - q * y1;
        a0 = b0;
        b0 = r;
        x0 = x1;
        x1 = x;
        y0 = y1;
        y1 = y;
    }
    extended_gcd_result res;
    res.gcd = a0;
    res.x = x0;
    res.y = y0;
    return res;
}
big_int mod_inverse(big_int a, big_int n) {
    if (n <= 0) throw std::invalid_argument("modulus must be positive");
    a %= n; //присваивание с остатком
    if (a < 0) a += n;//тк остаток от деления отриц прибавляем n чтобы не был отриц
    extended_gcd_result eg = extended_gcd(a, n);
    if (eg.gcd != 1) throw std::invalid_argument("inverse does not exist");//обр сущ при нод()=1
    big_int inv = eg.x % n;//берем остаток от x
    if (inv < 0) inv += n;
    return inv;
}
int legendre_symbol(big_int a, big_int p) {//p-нечетное число простое число
    if (p <= 2 || p % 2 == 0) throw std::invalid_argument("p must be odd prime");
    a %= p;
    if (a < 0) a += p;
    if (a == 0) return 0;//если а делится на р
    big_int exp = (p - 1) / 2;
    big_int res = mod_pow(a, exp, p);//критерий эйлера (a/p) ≡ a^((p-1)/2) (mod p)
    if (res == 1) return 1;//а -квадратичный вычет
    if (res == p - 1) return -1;
    return 0;
}
int jacobi_symbol(big_int a, big_int n) {//n-нечетное положительное
    if (n <= 0 || n % 2 == 0) throw std::invalid_argument("n must be positive odd");
    a %= n;
    if (a < 0) a += n;
    int result = 1;//переменная для накопления знака
    while (a != 0) {
        while (a % 2 == 0) {//когда а четное
            a /= 2;
            big_int n_mod8 = n % 8;
            if (n_mod8 == 3 || n_mod8 == 5) result = -result;//(2/n) = 1 если n ≡ ±1 mod 8, иначе -1 если +-3mod8
        }
        big_int temp = a;
        a = n;//з-н взаимности меняем a и n местами (факт из теории гаусса о квадр полях)
        n = temp;
        if (a % 4 == 3 && n % 4 == 3) result = -result;//если оа числа дают ост 3 по мод 4,знак меняется
        a %= n;//приводим новое а к остатку от деления на новое n
    }
    if (n == 1) return result;//если n!=1 значит нод не равен 1 и символ 0
    return 0;
}
big_int phi_definition(big_int n) {//ф.эйлера по опр
    if (n <= 0) throw std::invalid_argument("n must be positive");
    big_int count = 0;
    for (big_int k = 1; k <= n; ++k) {
        if (gcd_euclid(k, n) == 1) ++count;//нод()=1
    }
    return count;
}
std::vector<big_int> factorize(big_int n) {//возвр вектор уник простых делителей (без учета степеней)
    std::vector<big_int> factors;
    if (n <= 1) return factors;
    big_int d = 2;//начинаем перебор делителей с 2
    while (d * d <= n) {//достаточно проверить делители до корня из n
        if (n % d == 0) {
            factors.push_back(d);
            while (n % d == 0) n /= d;//удал все вхождения д из n
        }
        if (d == 2) d = 3;
        else d += 2;//если д=2 то след 3 иначе увел на 2 пропускаем четные кроме 2
    }
    if (n > 1) factors.push_back(n);//если осталось то простой множ больший чем корень из исх числа
    return factors;
}
big_int phi_factorization(big_int n) {
    if (n <= 0) throw std::invalid_argument("n must be positive");
    if (n == 1) return 1;//φ(1) = 1
    auto factors = factorize(n);//уник простые делители n
    big_int result = n;
    for (auto p : factors) {//для каждого простого p из вектора
        result = result / p * (p - 1);//n * (1 - 1/p) для каждого простого дел
    }
    return result;
}
big_int phi_dft(big_int n) {
    if (n <= 0) throw std::invalid_argument("n must be positive");
    double sum = 0.0;
    for (big_int k = 1; k <= n; ++k) {
        big_int g = gcd_euclid(k, n);
        double angle = 2.0 * pi * (double)k / (double)n;//угол
        sum += (double)g * cos(angle);
    }
    return (big_int)std::llround(sum);//округляем сам до ближ целого long long затем к big_int
}
class i_prime_test {//интерфейс для вероятностного теста простоты
public:
    virtual bool is_prime(const big_int& n, double min_probability) = 0;
    virtual ~i_prime_test() {}//обеспечивает корректное удал объектов наследников через указ на баз класс
};
class prime_test_base : public i_prime_test {//баз класс с шабл методом
public://наследник переопределяет метод, реализует паттерн шабл метод - общий алгоритм проверки пустоты, который вызывает абстрактынй метод one_iter 
    bool is_prime(const big_int& n, double min_probability) override {
        if (min_probability < 0.5 || min_probability >= 1.0)
            throw std::invalid_argument("probability must be in [0.5, 1)");
        if (n < 2) return false;
        if (n == 2 || n == 3) return true;
        if (n % 2 == 0) return false;
        for (big_int p : {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {//проверка делимости на простые числа для отсева сост чисел
            if (n % p == 0) return n == p;//простое только если n==p
        }
        int iterations = (int)std::ceil(-std::log2(1.0 - min_probability));//итер теста фера умен вер-ть ошибки в 2р отриц тк лог от числа мен 1 отриц округ вверх итер целое число приводим к инт
        if (iterations < 1) iterations = 1;
        for (int i = 0; i < iterations; ++i) {
            if (!one_iteration(n)) return false;//вызов вирт метода (реал 1 проверку малой теор) если метод вернул фолс сразу выходим с рез фолс
        }
        return true;
    }
protected://ток наследникам обязаны реал метод
    virtual bool one_iteration(const big_int& n) = 0;//нельзя создать объект этого класса только наследие
};
class fermat_prime_test : public prime_test_base {
protected:
    bool one_iteration(const big_int& n) override {
        if (n < 4) return true;//простое 0 и 1 уже отсеяны в из прайм
        std::random_device rd;
        std::mt19937_64 gen(rd());
        big_int a;
        do {//цикл идет как мин 1р и повт пока усл истинно
            uint64_t r = std::uniform_int_distribution<uint64_t>(2, 1000000)(gen);
            a = r % (n - 2);
            if (a < 2) a = 2;
        } while (a < 2 || a > n - 2);//генер пока a не будет от 2 до n-2 
        big_int mod = mod_pow(a, n - 1, n);
        return mod == 1;//прошло итер теста ферма (вер простое)
    }
};
void test_legendre() {
    assert(legendre_symbol(2, 7) == 1);
    assert(legendre_symbol(3, 7) == -1);
    assert(legendre_symbol(0, 7) == 0);
    std::cout << "Legendre tests passed\n";
}
void test_jacobi() {
    assert(jacobi_symbol(2, 9) == 1);
    assert(jacobi_symbol(3, 5) == -1);
    assert(jacobi_symbol(0, 5) == 0);
    std::cout << "Jacobi tests passed\n";
}
void test_gcd() {
    assert(gcd_euclid(48, 18) == 6);
    assert(gcd_euclid(-48, 18) == 6);
    std::cout << "GCD tests passed\n";
}
void test_extended_gcd() {
    auto eg = extended_gcd(48, 18);
    assert(eg.gcd == 6);
    assert(48 * eg.x + 18 * eg.y == 6);
    std::cout << "Extended GCD tests passed\n";
}
void test_mod_pow() {
    assert(mod_pow(2, 10, 1000) == 24);
    assert(mod_pow(5, 3, 13) == 8);
    std::cout << "Mod pow tests passed\n";
}
void test_mod_inverse() {
    assert(mod_inverse(3, 11) == 4); // 3*4=12≡1
    assert(mod_inverse(17, 3120) == 2753); // known
    std::cout << "Mod inverse tests passed\n";
}
void test_phi() {
    big_int n = 10;
    assert(phi_definition(n) == 4);
    assert(phi_factorization(n) == 4);
    assert(phi_dft(n) == 4);
    n = 1;
    assert(phi_definition(n) == 1);
    assert(phi_factorization(n) == 1);
    assert(phi_dft(n) == 1);
    n = 7;
    assert(phi_definition(n) == 6);
    assert(phi_factorization(n) == 6);
    assert(phi_dft(n) == 6);
    std::cout << "Phi tests passed\n";
}
void test_prime_test() {
    fermat_prime_test test;
    assert(test.is_prime(2, 0.9) == true);
    assert(test.is_prime(3, 0.9) == true);
    assert(test.is_prime(4, 0.9) == false);
    assert(test.is_prime(17, 0.9) == true);
    assert(test.is_prime(561, 0.9) == false); 
    std::cout << "Prime test tests passed\n";
}
void demo() {
    big_int a = 7, b = 11;
    std::cout << "Legendre(2,7) = " << legendre_symbol(2, 7) << "\n";
    std::cout << "Jacobi(2,9) = " << jacobi_symbol(2, 9) << "\n";
    std::cout << "GCD(48,18) = " << gcd_euclid(48, 18) << "\n";
    auto eg = extended_gcd(48, 18);
    std::cout << "Extended GCD: gcd=" << eg.gcd << ", x=" << eg.x << ", y=" << eg.y << "\n";
    std::cout << "mod_pow(2,10,1000) = " << mod_pow(2, 10, 1000) << "\n";
    std::cout << "mod_inverse(3,11) = " << mod_inverse(3, 11) << "\n";
    big_int n = 10;
    std::cout << "phi(10) by definition = " << phi_definition(n) << "\n";
    std::cout << "phi(10) by factorization = " << phi_factorization(n) << "\n";
    std::cout << "phi(10) by DFT = " << phi_dft(n) << "\n";
    fermat_prime_test test;
    std::cout << "Is 17 prime? " << (test.is_prime(17, 0.99) ? "yes" : "no") << "\n";
    std::cout << "Is 561 prime? " << (test.is_prime(561, 0.99) ? "yes" : "no") << "\n";
}
int main() {
    test_legendre();
    test_jacobi();
    test_gcd();
    test_extended_gcd();
    test_mod_pow();
    test_mod_inverse();
    test_phi();
    test_prime_test();
    demo();
    std::cout << "All tests passed\n";
    return 0;
}