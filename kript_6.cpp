#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <chrono>
using boost::multiprecision::cpp_int;
using boost::multiprecision::pow;
typedef cpp_int big_int;
big_int mod_pow(big_int base, big_int exp, big_int mod) {
    if (mod == 0) throw std::invalid_argument("mod cannot be zero");
    big_int result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = (result * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return result;
}
big_int gcd_euclid(big_int a, big_int b) {
    if (a < 0) a = -a;
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
        big_int x = x0 - q * x1;//коэф для нового остатка
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
    a %= n;
    if (a < 0) a += n;
    extended_gcd_result eg = extended_gcd(a, n);
    if (eg.gcd != 1) throw std::invalid_argument("inverse does not exist");
    big_int inv = eg.x % n;
    if (inv < 0) inv += n;
    return inv;
}
int jacobi_symbol(big_int a, big_int n) {
    if (n <= 0 || n % 2 == 0) throw std::invalid_argument("n must be positive odd");
    a %= n;
    if (a < 0) a += n;
    int result = 1;
    while (a != 0) {
        while (a % 2 == 0) {
            a /= 2;
            big_int n_mod8 = n % 8;
            if (n_mod8 == 3 || n_mod8 == 5) result = -result;
        }
        big_int temp = a;
        a = n;
        n = temp;
        if (a % 4 == 3 && n % 4 == 3) result = -result;
        a %= n;//а к мен знач
    }
    if (n == 1) return result;
    return 0;
}
class i_prime_test {
public:
    virtual bool is_prime(const big_int& n, double min_probability) = 0;
    virtual ~i_prime_test() {}
};
class prime_test_base : public i_prime_test {
public:
    bool is_prime(const big_int& n, double min_probability) override {
        if (min_probability < 0.5 || min_probability >= 1.0)
            throw std::invalid_argument("probability must be in [0.5, 1)");
        if (n < 2) return false;
        if (n == 2 || n == 3) return true;
        if (n % 2 == 0) return false;
        for (big_int p : {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
            if (n % p == 0) return n == p;
        }
        int iterations = (int)std::ceil(-std::log2(1.0 - min_probability));
        if (iterations < 1) iterations = 1;
        for (int i = 0; i < iterations; ++i) {
            if (!one_iteration(n)) return false;
        }
        return true;
    }
protected:
    virtual bool one_iteration(const big_int& n) = 0;
};
class solovay_strassen_prime_test : public prime_test_base {
protected:
    bool one_iteration(const big_int& n) override {
        if (n < 4) return true;
        std::random_device rd;
        std::mt19937_64 gen(rd());
        big_int a;
        do {//генерирует случ а [2,n-2]
            uint64_t r = std::uniform_int_distribution<uint64_t>(2, 1000000)(gen);
            a = r % (n - 2);//берем осаток от деления r на n-2 чтобы число в [0,n-3]
            if (a < 2) a = 2;
        } while (a < 2 || a > n - 2);
        if (gcd_euclid(a, n) != 1) return false;//нод=1
        big_int left = mod_pow(a, (n - 1) / 2, n);
        int right = jacobi_symbol(a, n);
        if (right < 0) right = (int)(n - 1);//если якоби=-1 то в кольце вычетов n-1
        if (left != right) return false;
        return true;
    }
};
class miller_rabin_prime_test : public prime_test_base {
protected:
    bool one_iteration(const big_int& n) override {
        if (n < 4) return true;
        std::random_device rd;
        std::mt19937_64 gen(rd());
        big_int a;
        do {
            uint64_t r = std::uniform_int_distribution<uint64_t>(2, 1000000)(gen);
            a = r % (n - 2);
            if (a < 2) a = 2;
        } while (a < 2 || a > n - 2);
        if (gcd_euclid(a, n) != 1) return false;
        big_int d = n - 1;
        int s = 0;
        while (d % 2 == 0) {
            d /= 2;
            ++s;
        }//для нчетного n n-1=d*2^s d-нечетное
        big_int x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) return true;
        for (int r = 1; r < s; ++r) {//от r до s-1 это проверяет усл сущ ли r такое, что a^(d*2^r)=-1(modn)
            x = (x * x) % n;
            if (x == n - 1) return true;
        }
        return false;
    }
};
struct rsa_public_key {//для открытого ключа
    big_int e;//откр экспонента
    big_int n;
};
struct rsa_private_key {//для закрытого
    big_int d;//закр экспонента
    big_int n;
};
class rsa_key_generator {
public:
    rsa_key_generator(i_prime_test* test) : prime_test(test) {}//констр принимает указ на объект, реализующий интерфейс шзпраймтест, прив поле иниц переданным указ - внедрение зав
    std::pair<rsa_public_key, rsa_private_key> generate_keys(int bits, double min_probability) {//битс-общ битовая длина модуля N
        if (bits < 16) throw std::invalid_argument("bits too small");//рса работает с числами которые должны быть больши шифр данных, макс число будет меньше 65536 1 байт не поместиться
        big_int p, q;
        do {
            p = generate_prime(bits / 2 + (bits % 2), min_probability);//если битс нечетное то p получает на 1байт больше q чтобы общ длина была точно bits
            q = generate_prime(bits / 2, min_probability);
        } while (abs_diff(p, q) < (big_int(1) << (bits / 2 - 100)) || p == q);//если p=q,разность меньше 2^(bits/2-100) -защита от атаки ферма если p,q близки то N модуль можно разложить
        big_int n = p * q;//сдвиг 1 влево на битс\2-100 бит, создает мин допуст знач разности
        big_int phi = (p - 1) * (q - 1);//для вычисл закрытой экспоненты d
        big_int e = 65537;//отк эксп число ферма взаимно просто с ф(н)
        if (gcd_euclid(e, phi) != 1) {
            for (big_int candidate = 65537; candidate < phi; ++candidate) {
                if (gcd_euclid(candidate, phi) == 1) { e = candidate; break; }
            }//при нод=1 гарантия сущ обр элем d
        }
        big_int d = mod_inverse(e, phi);//закр эксп
        big_int n_4th_root = big_int(1) << (bits / 4);
        if (d < n_4th_root) {//если d меньше 4корня из N(2^(bits/4)) то ключ уязвим для атаки Винера
            return generate_keys(bits, min_probability);
        }
        rsa_public_key pub{ e, n };
        rsa_private_key priv{ d, n };
        return { pub, priv };
    }
private:
    i_prime_test* prime_test;//указ на объект теста простоты
    big_int generate_prime(int bits, double min_probability) {
        if (bits < 2) throw std::invalid_argument("bits too small for prime");
        std::random_device rd;
        std::mt19937_64 gen(rd());
        big_int min_val = big_int(1) << (bits - 1);//старший бит 1
        big_int max_val = (big_int(1) << bits) - 1;//все биты 1
        while (true) {
            big_int candidate = random_big_int(min_val, max_val, gen);//случ простое
            if (candidate % 2 == 0) ++candidate;
            if (prime_test->is_prime(candidate, min_probability)) return candidate;
        }
    }
    big_int random_big_int(const big_int& min_val, const big_int& max_val, std::mt19937_64& gen) {
        big_int range = max_val - min_val + 1;
        big_int result = 0;
        int blocks = (int)std::ceil((double)boost::multiprecision::msb(range) / 64) + 1;//скок 64-битных блоков нужно сгенерировать чтобы покртыь все
        for (int i = 0; i < blocks; ++i) {
            uint64_t r = std::uniform_int_distribution<uint64_t>(0, UINT64_MAX)(gen);
            result = (result << 64) | r;//добавляем в рез (сдвиг влево на 64 бита и побитовое ИЛИ)
        }
        result = min_val + (result % range);//приводим к диапозону
        return result;
    }
    big_int abs_diff(const big_int& a, const big_int& b) {
        return (a > b) ? a - b : b - a;
    }
};
class rsa_cipher {
public:
    rsa_cipher(const rsa_public_key& pub_key) : e(pub_key.e), n(pub_key.n), d(0) {
        init_block_size();//вычисляет размер блока на основе n
    }
    rsa_cipher(const rsa_private_key& priv_key) : d(priv_key.d), n(priv_key.n), e(0) {
        init_block_size();
    }
    std::vector<uint8_t> encrypt(const std::vector<uint8_t>& data) {
        if (e == 0) throw std::runtime_error("no public key");
        return process(data, true);
    }
    std::vector<uint8_t> decrypt(const std::vector<uint8_t>& data) {
        if (d == 0) throw std::runtime_error("no private key");
        return process(data, false);
    }
private:
    big_int e, d, n;
    int block_size_bits;//блок меньше n
    int block_size_bytes;
    void init_block_size() {
        block_size_bits = boost::multiprecision::msb(n) - 1;//возвр индекс самого старшего бита в числе n,-1 чтобы получить макс число бит для блока
        block_size_bytes = (block_size_bits + 7) / 8;//с окр вверх число байт, которое занимает блок
        if (block_size_bytes < 1) block_size_bytes = 1;
    }
    std::vector<uint8_t> process(const std::vector<uint8_t>& data, bool encrypt) {
        std::vector<uint8_t> result;
        size_t len = data.size();
        size_t block_len = block_size_bytes;
        for (size_t i = 0; i < len; i += block_len) {
            size_t cur_len = std::min(block_len, len - i);//фактическая длина тек блока (последний мб крч)
            std::vector<uint8_t> block(data.begin() + i, data.begin() + i + cur_len);
            std::vector<uint8_t> padded_block(block_len, 0);
            std::copy(block.begin(), block.end(), padded_block.end() - cur_len);//содержимое блок в конец паддедблок помещает блок в младшие байты,старшие нулевые, так если блок короткий он доп нулями слева (чтобы не получить число большее n)
            big_int m = bytes_to_big_int(padded_block);//преобр дополн байтовый массив в число
            big_int c;//рез
            if (encrypt) c = mod_pow(m, e, n);
            else c = mod_pow(m, d, n);
            std::vector<uint8_t> out_block = big_int_to_bytes(c, block_len);
            result.insert(result.end(), out_block.begin(), out_block.end());
        }
        return result;
    }
    big_int bytes_to_big_int(const std::vector<uint8_t>& bytes) {
        big_int val = 0;
        for (unsigned char b : bytes) {
            val = (val << 8) | b;
        }
        return val;
    }
    std::vector<uint8_t> big_int_to_bytes(const big_int& val, size_t len) {//преоб большое целое обр в вектор байтов фикс длинв len
        std::vector<uint8_t> bytes(len, 0);
        big_int temp = val;
        for (int i = (int)len - 1; i >= 0; --i) {
            bytes[i] = (uint8_t)(temp & 0xFF);
            temp >>= 8;//вправо на 8 бит
        }
        return bytes;
    }
};
void test_solovay_strassen() {
    solovay_strassen_prime_test test;
    assert(test.is_prime(2, 0.9) == true);
    assert(test.is_prime(3, 0.9) == true);
    assert(test.is_prime(4, 0.9) == false);
    assert(test.is_prime(17, 0.9) == true);
    assert(test.is_prime(561, 0.9) == false);
    std::cout << "Solovay-Strassen tests passed\n";
}
void test_miller_rabin() {
    miller_rabin_prime_test test;
    assert(test.is_prime(2, 0.9) == true);
    assert(test.is_prime(3, 0.9) == true);
    assert(test.is_prime(4, 0.9) == false);
    assert(test.is_prime(17, 0.9) == true);
    assert(test.is_prime(561, 0.9) == false);
    std::cout << "Miller-Rabin tests passed\n";
}
void test_rsa() {
    miller_rabin_prime_test prime_test;
    rsa_key_generator key_gen(&prime_test);
    auto start = std::chrono::high_resolution_clock::now();
    auto keys = key_gen.generate_keys(256, 0.99);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Key generation time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms\n";
    rsa_public_key pub = keys.first;
    rsa_private_key priv = keys.second;
    rsa_cipher enc(pub);
    rsa_cipher dec(priv);
    std::string message = "Hello RSA!";
    std::vector<uint8_t> plain(message.begin(), message.end());
    std::vector<uint8_t> ciphertext = enc.encrypt(plain);
    std::vector<uint8_t> decrypted = dec.decrypt(ciphertext);
    size_t start_pos = 0;
    while (start_pos < decrypted.size() && decrypted[start_pos] == 0) ++start_pos;
    std::vector<uint8_t> trimmed(decrypted.begin() + start_pos, decrypted.end());
    assert(plain == trimmed);
    std::cout << "RSA encryption/decryption test passed\n";
}
int main() {
    test_solovay_strassen();
    test_miller_rabin();
    test_rsa();
    miller_rabin_prime_test prime_test;
    rsa_key_generator key_gen(&prime_test);
    auto keys = key_gen.generate_keys(256, 0.99);
    rsa_public_key pub = keys.first;
    rsa_private_key priv = keys.second;
    std::cout << "Public key: e=" << pub.e << ", n=" << pub.n << "\n";
    std::cout << "Private key: d=" << priv.d << ", n=" << priv.n << "\n";
    rsa_cipher enc(pub);
    rsa_cipher dec(priv);
    std::string msg = "RSA is working!";
    std::vector<uint8_t> plain(msg.begin(), msg.end());
    std::vector<uint8_t> cipher = enc.encrypt(plain);
    std::vector<uint8_t> recovered = dec.decrypt(cipher);
    size_t start_pos = 0;
    while (start_pos < recovered.size() && recovered[start_pos] == 0) ++start_pos;
    std::vector<uint8_t> trimmed(recovered.begin() + start_pos, recovered.end());
    std::string recovered_str(trimmed.begin(), trimmed.end());
    std::cout << "Original: " << msg << "\n";
    std::cout << "Decrypted: " << recovered_str << "\n";
    assert(msg == recovered_str);
    std::cout << "All tests passed\n";
    return 0;
}