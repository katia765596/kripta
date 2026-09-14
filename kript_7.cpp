#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>//для возврата пар
#include <vector>
using boost::multiprecision::cpp_int;
using big_int = cpp_int;
big_int abs_big_int(big_int x) { return x < 0 ? -x : x; }
big_int abs_diff(const big_int& a, const big_int& b) { return a >= b ? a - b : b - a; }
big_int mod_pow(big_int base, big_int exp, const big_int& mod) {
    if (mod <= 0) throw std::invalid_argument("modulus must be positive");
    if (exp < 0) throw std::invalid_argument("exponent must be non-negative");
    big_int result = 1;
    base %= mod;
    while (exp > 0) {
        if ((exp & 1) != 0) result = (result * base) % mod;//мл бит равен 1 то
        exp >>= 1;//вправо на 1 бит(делим на 2)
        if (exp != 0) base = (base * base) % mod;
    }
    return result;
}
big_int gcd_euclid(big_int a, big_int b) {
    a = abs_big_int(a);
    b = abs_big_int(b);
    while (b != 0) {
        big_int r = a % b;
        a = b;
        b = r;
    }
    return a;
}
struct extended_gcd_result { big_int gcd; big_int x; big_int y; };
extended_gcd_result extended_gcd(big_int a, big_int b) {
    big_int old_r = a, r = b;
    big_int old_s = 1, s = 0;
    big_int old_t = 0, t = 1;
    while (r != 0) {
        big_int q = old_r / r;
        big_int tmp = old_r - q * r; old_r = r; r = tmp;//обновляем остатки
        tmp = old_s - q * s; old_s = s; s = tmp;//коэф для а
        tmp = old_t - q * t; old_t = t; t = tmp;//коэф для b
    }
    if (old_r < 0) { old_r = -old_r; old_s = -old_s; old_t = -old_t; }
    return { old_r, old_s, old_t };
}
big_int mod_inverse(big_int a, const big_int& n) {
    if (n <= 0) throw std::invalid_argument("modulus must be positive");
    a %= n;
    if (a < 0) a += n;
    extended_gcd_result eg = extended_gcd(a, n);
    if (eg.gcd != 1) throw std::invalid_argument("inverse does not exist");
    big_int result = eg.x % n;
    if (result < 0) result += n;
    return result;
}
big_int sqrt_big_int(const big_int& n) {
    if (n < 0) throw std::invalid_argument("sqrt of negative number");
    if (n == 0) return 0;
    unsigned bit_count = static_cast<unsigned>(boost::multiprecision::msb(n)) + 1;//индекс самого стар бита +1 -колво бит в числе (старишй -конечный бит)
    big_int x = big_int(1) << ((bit_count + 1) / 2);//приближ для корня берем 1, сдвигаем на сток бит - дает число, которое явл степенью двойки примерно равной половине бит длины n, испол для ускорения сходимости алг ньютона нач приближ больше или равно корню
    while (true) {
        big_int y = (x + n / x) >> 1;//сдвиг вправо на 1бит, деление на 2, у след приближ
        if (y >= x) break;//достигли сходимости иначе расход (приближ перестало умен) тек х явл целой частью корня
        x = y;
    }
    while ((x + 1) * (x + 1) <= n) ++x;//х целая часть корня, может быть на 1 меньге из-за делений, корректировка сверху
    while (x * x > n) --x;//x^2<=n х-наиб целое квадрат которого не преврросходит n
    return x;
}
bool is_perfect_square(const big_int& n, big_int& root) {
    if (n < 0) return false;
    root = sqrt_big_int(n);
    return root * root == n;
}
struct rsa_public_key { big_int e; big_int n; };
struct rsa_private_key { big_int d; big_int n; };
class rsa_weak_key_generator {
public:
    rsa_weak_key_generator() : gen_(std::random_device{}()) {}
    std::pair<rsa_public_key, rsa_private_key> generate_fermat_vulnerable(int bits, double min_probability) {
        validate_bits(bits);
        const int p_bits = bits / 2;//бит длина для каждого простого числа чтобы получить модуль N длиной битс, если битс нечетное, то одно из чисел будет на бит длинее
        const big_int max_difference = big_int(1) << 20;//макс разница до 2^20
        big_int p = generate_prime(p_bits, min_probability);
        big_int q;
        while (true) {
            uint64_t delta = random_uint64(1, (uint64_t(1) << 20) - 1);
            q = p + delta;//q больше p на небольшую велечину разность дельта меньше 2^20
            if (boost::multiprecision::msb(q) + 1 != static_cast<unsigned>(p_bits)) continue;//если бит длина q!=p_bits ищем новую дельут (т.е старший бит на позиции p_bits-1)
            if ((q & 1) == 0) ++q;//если четное
            if (q == p) continue;
            if (abs_diff(p, q) >= max_difference) continue;
            if (is_prime(q, min_probability)) break;
        }
        return make_rsa_key_pair(p, q);
    }
    std::pair<rsa_public_key, rsa_private_key> generate_wiener_vulnerable(int bits, double /*min_probability*/) {//вер-ть комментарий для совместимости с интерфейсом тк испол детерминированный тест миллера-рабина с фикс основаниями
        validate_bits(bits);
        const int p_bits = bits / 2;
        big_int p, q;
        do {
            p = generate_prime(p_bits, 0.8);
            q = generate_prime(p_bits, 0.8);
        } while (p == q);
        const big_int phi = (p - 1) * (q - 1);
        const uint64_t max_d = (uint64_t(1) << 20) - 1;//опред макс знач для закрытой экс d 2^20-1 (всего 20 бит) усл d<N^(1/4) для 256-битного N выполняется
        big_int d;
        do { d = random_uint64(2, max_d); } while (gcd_euclid(d, phi) != 1);//d должно быть взаимно простым с phi иначе обр эл не сущ
        big_int e = mod_inverse(d, phi);//вычисляем отк эксп e как обр к d по модулю phi e=d^(-1)mod phi,тк d маленькое e будет близко к phi
        rsa_public_key pub{ e, p * q };
        rsa_private_key priv{ d, p * q };
        return { pub, priv };
    }
private:
    std::mt19937_64 gen_;
    static void validate_bits(int bits) { if (bits < 16) throw std::invalid_argument("RSA key size must be at least 16 bits"); }
    uint64_t random_uint64(uint64_t min_value, uint64_t max_value) {
        std::uniform_int_distribution<uint64_t> dist(min_value, max_value);
        return dist(gen_);
    }
    big_int random_big_int(const big_int& min_value, const big_int& max_value) {
        if (min_value > max_value) throw std::invalid_argument("invalid random range");
        const big_int range = max_value - min_value + 1;
        const unsigned bits = static_cast<unsigned>(boost::multiprecision::msb(range)) + 1;
        const unsigned blocks = (bits + 63) / 64;//вычисляем сколько 64битных блоков нужно чтобы покрыть битс бит округляем вверх
        big_int result = 0;
        for (unsigned i = 0; i < blocks; ++i) {
            result <<= 64;
            result |= random_uint64(0, std::numeric_limits<uint64_t>::max());
        }
        return min_value + result % range;
    }
    big_int generate_prime(int bits, double min_probability) {
        if (bits < 2) throw std::invalid_argument("prime size is too small");
        big_int min_value = big_int(1) << (bits - 1);//старший бит =1
        big_int max_value = (big_int(1) << bits) - 1;//все биты 1
        while (true) {
            big_int candidate = random_big_int(min_value, max_value);
            if ((candidate & 1) == 0) ++candidate;//если четное
            if (candidate > max_value) continue;
            if (is_prime(candidate, min_probability)) return candidate;
        }
    }
    bool is_prime(const big_int& n, double /*min_probability*/) {
        if (n < 2) return false;
        static const uint32_t small_primes[] = { 2,3,5,7,11,13,17,19,23,29,31,37 };
        for (uint32_t p : small_primes) {
            if (n == p) return true;
            if (n % p == 0) return false;
        }
        return miller_rabin(n);
    }
    bool miller_rabin(const big_int& n) {
        big_int d = n - 1;
        unsigned s = 0;
        while ((d & 1) == 0) { d >>= 1; ++s; }//n-1=d*2^s где d нечетное, пока оно четное сдвигаем вправо на 1 бит (деление на 2) и увелич счетчик степени двойки
        static const uint64_t bases[] = { 2,3,5,7,11,13,17 };//фикс основания для теста
        for (uint64_t base : bases) {
            if (big_int(base) >= n) continue;
            big_int x = mod_pow(big_int(base), d, n);//вычисляем base^d mod n
            if (x == 1 || x == n - 1) continue;
            bool passed = false;
            for (unsigned r = 1; r < s; ++r) {
                x = (x * x) % n;
                if (x == n - 1) { passed = true; break; }
            }
            if (!passed) return false;
        }
        return true;//если все основания прошли-число простое
    }
    std::pair<rsa_public_key, rsa_private_key> make_rsa_key_pair(const big_int& p, const big_int& q) {
        if (p <= 1 || q <= 1) throw std::invalid_argument("invalid RSA primes");
        if (p == q) throw std::invalid_argument("RSA primes must be different");
        const big_int n = p * q;
        const big_int phi = (p - 1) * (q - 1);
        big_int e = 65537;
        if (gcd_euclid(e, phi) != 1) {
            e = 3;
            while (e < phi && gcd_euclid(e, phi) != 1) e += 2;
            if (e >= phi) throw std::runtime_error("unable to find RSA public exponent");
        }
        const big_int d = mod_inverse(e, phi);
        return { {e, n}, {d, n} };
    }
};
std::pair<big_int, big_int> fermat_attack(const rsa_public_key& pub) {
    const big_int& n = pub.n;
    if (n <= 0) throw std::invalid_argument("invalid RSA modulus");
    big_int a = sqrt_big_int(n);
    if (a * a < n) ++a;
    big_int b2 = a * a - n;//если p,q близки, то b2 должно быть полным квадратом
    while (true) {
        big_int b;
        if (is_perfect_square(b2, b)) {//функция проверяет явл ли b2 полнам квадратом если да в b записывается полный квадрат
            big_int p = a - b, q = a + b;
            if (p > 1 && q > 1 && p * q == n) return { p, q };
        }
        ++a;
        b2 += 2 * a - 1;//при переходе от a к а+1  разность (а+1)-а^2=2а+1 тк a увел то добавляем 2а+1 новое b2 соотв новому а
    }
}
std::vector<big_int> continued_fraction(big_int a, big_int b) {//возвр вектор неполных частных цепной дроби
    if (a < 0 || b <= 0) throw std::invalid_argument("a must be non-negative and b positive");
    std::vector<big_int> cf;
    while (b != 0) {
        cf.push_back(a / b);
        big_int r = a % b;
        a = b;
        b = r;
    }
    return cf;
}
std::pair<big_int, big_int> continued_fraction_to_fraction(const std::vector<big_int>& cf) {//ищем дробь по вектору не к исх а к нужной
    if (cf.empty()) return { 0, 1 };
    big_int numerator = 1, denominator = 0;

    for (auto it = cf.rbegin(); it != cf.rend(); ++it) {//идем по неполным частным в обр порядке
        big_int new_numerator = (*it) * numerator + denominator;//p'=a*p+q,q'=p, *it-тек неполное частное
        denominator = numerator;//старый числ новый знам
        numerator = new_numerator;//у числ новое знач
    }
    return { numerator, denominator };//возвр полученную дробь
}
big_int continued_fraction_to_value(const std::vector<big_int>& cf) {
    auto fraction = continued_fraction_to_fraction(cf);
    return fraction.first / fraction.second;//целая часть от деления числ на знам (первый член цепной дроби)
}
std::vector<std::pair<big_int, big_int>> convergents(const std::vector<big_int>& cf) {//все подходящие дроби возвр вектор пар(числ,знам) для каждой подходящей дроби
    std::vector<std::pair<big_int, big_int>> result;
    big_int p_minus_2 = 0, p_minus_1 = 1;//иниц нач значения для реккурентых формул
    big_int q_minus_2 = 1, q_minus_1 = 0;
    for (const big_int& a : cf) {
        big_int p = a * p_minus_1 + p_minus_2;//Р_k = a_k * P_{k-1} + P_{k-2}
        big_int q = a * q_minus_1 + q_minus_2;//Q_k = a_k * Q_{k-1} + Q_{k-2}
        result.push_back({ p, q });
        p_minus_2 = p_minus_1; p_minus_1 = p;//сдвиг
        q_minus_2 = q_minus_1; q_minus_1 = q;
    }
    return result;
}
std::string calkin_wilf_path(big_int a, big_int b) {//путь от корня к дроби 
    if (a <= 0 || b <= 0) throw std::invalid_argument("positive required");
    std::string path;
    while (a != b) {
        if (a < b) {
            big_int count = (b - 1) / a;//переход к родителю лев потомка a/(b-a) делаем несколько гагов пока b не станет меньше или равно а
            path.append(count.convert_to<std::size_t>(), 'L');//добавляем каунт раз символ L преобразуем бигинт в сайзт
            b -= count * a;//обновляем знам значением: сколько вычитали числитель из знаменателя для получения родителя
        }
        else {
            big_int count = (a - 1) / b;//от правого потомка к родителю (a-b)/b
            path.append(count.convert_to<std::size_t>(), 'R');
            a -= count * b;//обн числ знач: сколько раз вычлили знам из числ для получения родителя прав потомка
        }
    }
    return path;//путь от дроби к корню порядок от послед шага к первому
}
std::pair<big_int, big_int> calkin_wilf_from_path(const std::string& path) {
    big_int a = 1, b = 1;
    for (char c : path) {
        if (c == 'L' || c == '0') b = a + b;
        else if (c == 'R' || c == '1') a = a + b;
        else throw std::invalid_argument("invalid path");
    }
    return { a, b };
}
std::string stern_brocot_path(big_int a, big_int b) {//путь к дроби
    if (a <= 0 || b <= 0) throw std::invalid_argument("positive required");
    big_int l_num = 0, l_den = 1;//границы левая 0/1 правая 1/0
    big_int r_num = 1, r_den = 0;
    std::string path;//строка для пути
    while (true) {
        big_int m_num = l_num + r_num;//медианта 
        big_int m_den = l_den + r_den;
        if (m_num * b == a * m_den) { //равна ли медианта иск дроби
            break;
        }
        if (a * m_den < m_num * b) { //если иск дробь меньше медианты - левое поддерво
            path.push_back('L');
            r_num = m_num;//правая граница медианта
            r_den = m_den;
        }
        else {
            path.push_back('R');
            l_num = m_num;//левая граница медианта
            l_den = m_den;
        }
    }
    return path;//от корня к дроби
}
std::pair<big_int, big_int> stern_brocot_from_path(const std::string& path) {
    return calkin_wilf_from_path(path); 
}
std::vector<std::pair<big_int, big_int>> stern_brocot_convergents_by_path(const std::string& path) {//пос-ть всех дробей на пути от корня до дроби
    std::vector<std::pair<big_int, big_int>> result;
    big_int a = 1, b = 1;
    result.push_back({ a, b });//корень первая подход дробь
    for (char c : path) {
        if (c == 'L' || c == '0') b = a + b;
        else if (c == 'R' || c == '1') a = a + b;
        else throw std::invalid_argument("invalid path");
        result.push_back({ a, b });
    }
    return result;
}
std::string stern_brocot_bit_path(big_int a, big_int b) {
    std::string p = stern_brocot_path(a, b);
    for (char& c : p) c = (c == 'L') ? '0' : '1';
    return p;
}
std::pair<big_int, big_int> wiener_attack(const rsa_public_key& pub) {//атакаааа вин принимает отк ключ и восстанавл простые
    const big_int& e = pub.e;
    const big_int& n = pub.n;
    if (e <= 0 || n <= 0) throw std::invalid_argument("invalid RSA public key");
    std::vector<big_int> cf = continued_fraction(e, n);//функция разлагает e/n в цепную дробь и возвр вектор неполных частных
    std::vector<std::pair<big_int, big_int>> convs = convergents(cf);//
    for (const auto& fraction : convs) {
        const big_int& k = fraction.first;//числ
        const big_int& d = fraction.second;//знам
        if (k == 0 || d == 0) continue;
        big_int ed_minus_1 = e * d - 1;//e*d-1 по теории должно быть кратно ф(n) если (k,d) верная пара
        if (ed_minus_1 <= 0) continue;
        if (ed_minus_1 % k != 0) continue;//e*d-1/k нацело
        big_int phi = ed_minus_1 / k;
        big_int s = n - phi + 1;//φ(N) = p*q - p - q + 1 = N - (p + q) + 1 =>p + q = N - φ(N) + 1
        if (s <= 0) continue;
        big_int disc = s * s - 4 * n;//x^2 - S * x + N = 0
        if (disc < 0) continue;
        big_int sqrt_disc;
        if (!is_perfect_square(disc, sqrt_disc)) continue;
        if ((s - sqrt_disc) % 2 != 0) continue;
        big_int p = (s - sqrt_disc) / 2;//корни уравнения
        big_int q = (s + sqrt_disc) / 2;
        if (p <= 1 || q <= 1) continue;
        if (p * q != n) continue;
        if ((e * d) % ((p - 1) * (q - 1)) != 1) continue;//e*d ≡ 1 (mod φ(N))
        return { p, q };
    }
    throw std::runtime_error("Wiener attack failed");
}
bool same_factors(big_int p1, big_int q1, big_int p2, big_int q2) {//сравнивает как неупор пары
    if (p1 > q1) std::swap(p1, q1);
    if (p2 > q2) std::swap(p2, q2);
    return p1 == p2 && q1 == q2;//если обе пары после сортировки равны
}
void require(bool condition, const std::string& message) {//принимает булево усл и строку соо
    if (!condition) throw std::runtime_error("TEST FAILED: " + message);//передаем упр обработчки искл
}
void test_fermat_attack() {
    rsa_weak_key_generator gen;
    auto keys = gen.generate_fermat_vulnerable(256, 0.8);
    const rsa_public_key& pub = keys.first;
    auto factors = fermat_attack(pub);
    big_int p = factors.first, q = factors.second;
    require(p * q == pub.n, "Fermat: p*q != n");
    require(p > 1 && q > 1, "Fermat: invalid factors");
    std::cout << "Fermat attack test passed\n";
}
void test_wiener_attack() {
    rsa_weak_key_generator gen;
    auto keys = gen.generate_wiener_vulnerable(256, 0.8);
    const rsa_public_key& pub = keys.first;
    auto factors = wiener_attack(pub);
    big_int p = factors.first, q = factors.second;
    require(p * q == pub.n, "Wiener: p*q != n");
    require(p > 1 && q > 1, "Wiener: invalid factors");
    std::cout << "Wiener attack test passed\n";
}
void test_extended_gcd() {
    big_int a = 240, b = 46;
    auto result = extended_gcd(a, b);
    require(result.gcd == 2, "extended_gcd: wrong gcd");
    require(a * result.x + b * result.y == result.gcd, "extended_gcd: Bezout identity");//беееезу
    std::cout << "Extended GCD tests passed\n";
}
void test_mod_inverse() {
    big_int a = 3, n = 11;
    big_int inv = mod_inverse(a, n);
    require(inv == 4, "mod_inverse: 3^-1 mod 11 != 4");
    require((a * inv) % n == 1, "mod_inverse: verification failed");//реал обр эл
    std::cout << "Mod inverse tests passed\n";
}
void test_continued_fraction() {
    big_int a = 355, b = 113;
    auto cf = continued_fraction(a, b);//вектор неполных частных
    require(cf.size() == 3, "CF: unexpected size");
    require(cf[0] == 3 && cf[1] == 7 && cf[2] == 16, "CF: incorrect coefficients");
    auto fraction = continued_fraction_to_fraction(cf);//воссстан дробь из цепной
    require(fraction.first == a && fraction.second == b, "CF: fraction reconstruction failed");
    auto convs = convergents(cf);//подходящие дроби
    require(!convs.empty(), "CF: convergents are empty");
    require(convs.back().first == a && convs.back().second == b, "CF: final convergent incorrect");//посл подход дробь равна исх дроби
    std::cout << "Continued fraction tests passed\n";
}
void test_trees() {
    big_int a = 3, b = 5;
    std::string cw_path = calkin_wilf_path(a, b);//строит путь от корня к  дроби
    auto cw_fraction = calkin_wilf_from_path(cw_path);//восст дробь по пути принимает страку и возв пару
    require(cw_fraction.first == a && cw_fraction.second == b, "Calkin-Wilf reconstruction failed");//вост дробь равна исх из строки
    std::string cw_bits = calkin_wilf_path(a, b);
    for (char& c : cw_bits) c = (c == 'L') ? '0' : '1';
    auto cw_fraction_bits = calkin_wilf_from_path(cw_bits);
    require(cw_fraction_bits.first == a && cw_fraction_bits.second == b, "Calkin-Wilf bit path failed");//из бит строки в исх дробь
    std::string sb_path = stern_brocot_path(a, b);//путь от корня к дроби
    auto sb_fraction = stern_brocot_from_path(sb_path);
    require(sb_fraction.first == a && sb_fraction.second == b, "Stern-Brocot reconstruction failed");//строка проверка
    std::string sb_bits = stern_brocot_bit_path(a, b);
    auto sb_fraction_bits = stern_brocot_from_path(sb_bits);
    require(sb_fraction_bits.first == a && sb_fraction_bits.second == b, "Stern-Brocot bit path failed");//бит проверка
    auto values = stern_brocot_convergents_by_path(sb_path);//возвр все дроби(узлы) на этом пути с корня. посл подх дробей
    require(!values.empty(), "Stern-Brocot values are empty");
    require(values.front().first == 1 && values.front().second == 1, "Stern-Brocot root incorrect");
    require(values.back().first == a && values.back().second == b, "Stern-Brocot final value incorrect");
    std::cout << "Tree tests passed\n";
}
void test_sqrt() {//проверяет функцию извлечения целой части квадратного корня
    big_int n = 123456789;
    big_int r = sqrt_big_int(n);
    require(r * r <= n, "sqrt: r^2 > n");//r-наиб целое квадрат которого не превосходит n
    require((r + 1) * (r + 1) > n, "sqrt: r is not floor(sqrt(n))");
    big_int perfect = big_int(12345) * 12345;
    require(sqrt_big_int(perfect) == 12345, "sqrt: perfect square");
    std::cout << "Sqrt tests passed\n";
}
void demonstrate_fermat() {
    std::cout << "\nFermat attack demonstration\n";
    rsa_weak_key_generator gen;
    auto keys = gen.generate_fermat_vulnerable(256, 0.8);
    const rsa_public_key& pub = keys.first;
    std::cout << "Public key:\n";
    std::cout << "e = " << pub.e << "\n";
    std::cout << "n = " << pub.n << "\n";
    auto start = std::chrono::high_resolution_clock::now();
    auto factors = fermat_attack(pub);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    std::cout << "p = " << factors.first << "\n";
    std::cout << "q = " << factors.second << "\n";
    std::cout << "p*q == n: " << (factors.first * factors.second == pub.n ? "true" : "false") << "\n";
    std::cout << "Attack time: " << elapsed.count() << " ms\n";
}
void demonstrate_wiener() {
    std::cout << "\nWiener attack demonstration\n";
    rsa_weak_key_generator gen;
    auto keys = gen.generate_wiener_vulnerable(256, 0.8);
    const rsa_public_key& pub = keys.first;
    std::cout << "Public key:\n";
    std::cout << "e = " << pub.e << "\n";
    std::cout << "n = " << pub.n << "\n";
    auto start = std::chrono::high_resolution_clock::now();
    auto factors = wiener_attack(pub);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    std::cout << "p = " << factors.first << "\n";
    std::cout << "q = " << factors.second << "\n";
    std::cout << "p*q == n: " << (factors.first * factors.second == pub.n ? "true" : "false") << "\n";
    std::cout << "Attack time: " << elapsed.count() << " ms\n";
}
void demonstrate_continued_fractions() {
    std::cout << "\nContinued fraction demo\n";
    big_int a = 355, b = 113;
    auto cf = continued_fraction(a, b);//возвр вектор неполных частных
    std::cout << "CF of " << a << "/" << b << ": [";
    for (size_t i = 0; i < cf.size(); ++i) {
        if (i != 0) std::cout << ", ";
        std::cout << cf[i];
    }
    std::cout << "]\n";
    auto fraction = continued_fraction_to_fraction(cf);//восст по вектору
    std::cout << "Recovered fraction: " << fraction.first << "/" << fraction.second << "\n";
    auto convs = convergents(cf);
    std::cout << "Convergents:\n";
    for (const auto& f : convs) std::cout << f.first << "/" << f.second << " ";
    std::cout << "\n";
}
void demonstrate_trees() {
    std::cout << "\nTree demo\n";
    big_int a = 3, b = 5;
    std::string cw_path = calkin_wilf_path(a, b);
    std::string cw_bits = calkin_wilf_path(a, b);
    for (char& c : cw_bits) c = (c == 'L') ? '0' : '1';
    std::cout << "Calkin-Wilf:\n";
    std::cout << "fraction = " << a << "/" << b << "\n";
    std::cout << "L/R path = " << cw_path << "\n";
    std::cout << "bit path = " << cw_bits << "\n";
    auto cw_value = calkin_wilf_from_path(cw_path);//восстанавл дробь по пути
    std::cout << "restored = " << cw_value.first << "/" << cw_value.second << "\n";
    std::string sb_path = stern_brocot_path(a, b);
    std::string sb_bits = stern_brocot_bit_path(a, b);
    std::cout << "\nStern-Brocot:\n";
    std::cout << "L/R path = " << sb_path << "\n";
    std::cout << "bit path = " << sb_bits << "\n";
    auto sb_value = stern_brocot_from_path(sb_path);
    std::cout << "restored = " << sb_value.first << "/" << sb_value.second << "\n";
    auto values = stern_brocot_convergents_by_path(sb_path);
    std::cout << "Path values:\n";
    for (const auto& f : values) std::cout << f.first << "/" << f.second << " ";
    std::cout << "\n";
}
int main() {
    try {
        std::cout << " RSA / Continued Fractions / Trees\n";
        std::cout << "Unit tests\n";
        test_extended_gcd();
        test_mod_inverse();
        test_sqrt();
        test_continued_fraction();
        test_trees();
        test_fermat_attack();
        test_wiener_attack();
        std::cout << "\nAll unit tests passed\n";
        demonstrate_fermat();
        demonstrate_wiener();
        demonstrate_continued_fractions();
        demonstrate_trees();
        std::cout << "All tests and demonstrations passed\n";
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "\nERROR: " << ex.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "\nERROR: unknown exception\n";
        return 1;
    }
}