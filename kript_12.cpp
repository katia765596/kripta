#include <vector>
#include <cstdint>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <string>
class BigInt {
public:
    using Digit = uint32_t;
    static constexpr int BITS = 32;
    static constexpr uint64_t BASE = 1ULL << BITS;
    std::vector<Digit> digits;
    BigInt() : digits(1, 0) {}
    BigInt(uint64_t val) {
        if (val == 0) digits = { 0 };
        else {
            while (val > 0) {
                digits.push_back(val & (BASE - 1));
                val >>= BITS;
            }
        }
    }
    BigInt(const std::string& s, int base = 10) { *this = fromString(s, base); }
    bool isZero() const { return digits.size() == 1 && digits[0] == 0; }
    size_t bitLength() const {
        if (isZero()) return 0;
        size_t bits = (digits.size() - 1) * BITS;
        Digit top = digits.back();
        while (top) { top >>= 1; ++bits; }
        return bits;
    }
    bool operator==(const BigInt& rhs) const { return digits == rhs.digits; }
    bool operator!=(const BigInt& rhs) const { return !(*this == rhs); }
    bool operator<(const BigInt& rhs) const {
        if (digits.size() != rhs.digits.size())
            return digits.size() < rhs.digits.size();
        for (int i = (int)digits.size() - 1; i >= 0; --i)
            if (digits[i] != rhs.digits[i]) return digits[i] < rhs.digits[i];
        return false;
    }
    bool operator>(const BigInt& rhs) const { return rhs < *this; }
    bool operator<=(const BigInt& rhs) const { return !(rhs < *this); }
    bool operator>=(const BigInt& rhs) const { return !(*this < rhs); }
    BigInt operator+(const BigInt& rhs) const {
        BigInt res;
        size_t n = std::max(digits.size(), rhs.digits.size());
        res.digits.resize(n + 1, 0);
        uint64_t carry = 0;
        for (size_t i = 0; i < n; ++i) {
            uint64_t a = i < digits.size() ? digits[i] : 0;
            uint64_t b = i < rhs.digits.size() ? rhs.digits[i] : 0;
            uint64_t sum = a + b + carry;
            res.digits[i] = (Digit)(sum % BASE);
            carry = sum / BASE;
        }
        if (carry) res.digits[n] = (Digit)carry;
        else res.digits.pop_back();
        return res;
    }
    BigInt operator-(const BigInt& rhs) const {
        if (*this < rhs) throw std::runtime_error("Negative result unsupported");
        BigInt res(*this);
        int64_t borrow = 0;
        for (size_t i = 0; i < rhs.digits.size() || borrow; ++i) {
            int64_t sub = (i < rhs.digits.size() ? rhs.digits[i] : 0) + borrow;
            if (res.digits[i] < sub) {
                res.digits[i] = (Digit)(res.digits[i] + BASE - sub);
                borrow = 1;
            }
            else {
                res.digits[i] = (Digit)(res.digits[i] - sub);
                borrow = 0;
            }
        }
        while (res.digits.size() > 1 && res.digits.back() == 0)
            res.digits.pop_back();
        return res;
    }
    BigInt operator*(const BigInt& rhs) const {
        BigInt res;
        res.digits.resize(digits.size() + rhs.digits.size(), 0);
        for (size_t i = 0; i < digits.size(); ++i) {
            uint64_t carry = 0;
            for (size_t j = 0; j < rhs.digits.size(); ++j) {
                uint64_t cur = res.digits[i + j] + (uint64_t)digits[i] * rhs.digits[j] + carry;
                res.digits[i + j] = (Digit)(cur % BASE);
                carry = cur / BASE;
            }
            res.digits[i + rhs.digits.size()] = (Digit)carry;
        }
        while (res.digits.size() > 1 && res.digits.back() == 0)
            res.digits.pop_back();
        return res;
    }
    BigInt operator/(const BigInt& divisor) const {
        BigInt q, r;
        divide(*this, divisor, q, r);
        return q;
    }
    BigInt operator%(const BigInt& divisor) const {
        BigInt q, r;
        divide(*this, divisor, q, r);
        return r;
    }
    BigInt operator<<(int shift) const {
        if (shift <= 0) return *this;
        BigInt res;
        int wordShift = shift / BITS;
        int bitShift = shift % BITS;
        res.digits.resize(digits.size() + wordShift + (bitShift ? 1 : 0), 0);
        Digit carry = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            uint64_t cur = ((uint64_t)digits[i] << bitShift) | carry;
            res.digits[i + wordShift] = (Digit)(cur % BASE);
            carry = (Digit)(cur / BASE);
        }
        if (bitShift > 0 && carry) {
            res.digits[digits.size() + wordShift] = carry;
        }
        while (res.digits.size() > 1 && res.digits.back() == 0)
            res.digits.pop_back();
        return res;
    }
    BigInt operator>>(int shift) const {
        if (shift <= 0) return *this;
        BigInt res = *this;
        int wordShift = shift / BITS;
        int bitShift = shift % BITS;
        if (wordShift > 0) {
            if (wordShift >= (int)res.digits.size()) return BigInt(0);
            res.digits.erase(res.digits.begin(), res.digits.begin() + wordShift);
        }
        if (bitShift > 0) {
            Digit carry = 0;
            for (int i = (int)res.digits.size() - 1; i >= 0; --i) {
                Digit new_carry = res.digits[i] << (BITS - bitShift);
                res.digits[i] = (res.digits[i] >> bitShift) | carry;
                carry = new_carry;
            }
            while (res.digits.size() > 1 && res.digits.back() == 0)
                res.digits.pop_back();
        }
        return res;
    }
    BigInt powMod(BigInt exp, const BigInt& mod) const {
        if (mod.isZero()) throw std::runtime_error("powMod with modulo zero");
        BigInt base = *this % mod;
        BigInt result(1);
        while (!exp.isZero()) {
            if ((exp.digits[0] & 1) == 1)
                result = (result * base) % mod;
            exp = exp >> 1;
            base = (base * base) % mod;
        }
        return result;
    }
    BigInt modInverse(const BigInt& mod) const {
        if (mod.isZero()) throw std::runtime_error("modInverse with zero modulus");
        BigInt a = *this % mod;
        if (a.isZero()) throw std::runtime_error("No inverse");
        BigInt b = mod;
        BigInt x0(0), x1(1);
        while (!a.isZero()) {
            if (b < a) {
                std::swap(a, b);
                std::swap(x0, x1);
            }
            BigInt q = b / a;
            BigInt r = b % a;
            b = a;
            a = r;

            BigInt qx1 = q * x1;
            if (x0 >= qx1)
                x0 = x0 - qx1;
            else
                x0 = mod - ((qx1 - x0) % mod);

            std::swap(x0, x1);
        }
        if (b != BigInt(1)) throw std::runtime_error("No inverse");
        return x0 % mod;
    }
    static BigInt gcd(BigInt a, BigInt b) {
        while (!b.isZero()) {
            BigInt t = b;
            b = a % b;
            a = t;
        }
        return a;
    }
    static bool isPrime(const BigInt& n, int rounds = 3) {
        if (n < BigInt(2)) return false;
        if (n == BigInt(2) || n == BigInt(3)) return true;
        if ((n.digits[0] & 1) == 0) return false;
        static const int smallPrimes[] = { 3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97 };
        for (int p : smallPrimes) {
            BigInt bp(p);
            if (n == bp) return true;
            if (n % bp == BigInt(0)) return false;
        }
        BigInt d = n - BigInt(1);
        int s = 0;
        while ((d.digits[0] & 1) == 0) {
            d = d >> 1;
            ++s;
        }
        static std::random_device rd;
        static std::mt19937_64 gen(rd());
        BigInt n_minus_3 = n - BigInt(3);

        for (int i = 0; i < rounds; ++i) {
            BigInt a = BigInt::random((int)n.bitLength(), gen) % n_minus_3 + BigInt(2);
            BigInt x = a.powMod(d, n);
            if (x == BigInt(1) || x == n - BigInt(1)) continue;
            bool composite = true;
            for (int r = 1; r < s; ++r) {
                x = (x * x) % n;
                if (x == n - BigInt(1)) { composite = false; break; }
            }
            if (composite) return false;
        }
        return true;
    }
    static BigInt random(int bits, std::mt19937_64& rng) {
        if (bits <= 0) return BigInt(0);
        int fullWords = bits / BITS;
        int rem = bits % BITS;
        std::uniform_int_distribution<uint64_t> dist;
        BigInt res;
        res.digits.resize(fullWords + (rem > 0 ? 1 : 0), 0);
        for (int i = 0; i < fullWords; ++i)
            res.digits[i] = (Digit)(dist(rng) & (BASE - 1));
        if (rem > 0) {
            Digit mask = (1ULL << rem) - 1;
            res.digits.back() = (Digit)(dist(rng) & mask);
            res.digits.back() |= (1ULL << (rem - 1));
        }
        else {
            res.digits.back() |= (1ULL << (BITS - 1));
        }
        while (res.digits.size() > 1 && res.digits.back() == 0)
            res.digits.pop_back();
        return res;
    }
    std::vector<uint8_t> toBytes() const {
        if (isZero()) return { 0 };
        size_t bytes = (bitLength() + 7) / 8;
        std::vector<uint8_t> out(bytes, 0);
        for (size_t i = 0; i < digits.size(); ++i) {
            for (int b = 0; b < 4; ++b) {
                size_t idx = i * 4 + b;
                if (idx < bytes)
                    out[bytes - 1 - idx] = (digits[i] >> (8 * b)) & 0xFF;
            }
        }
        return out;
    }
    static BigInt fromBytes(const std::vector<uint8_t>& bytes) {
        BigInt res(0);
        for (uint8_t b : bytes)
            res = res * BigInt(256) + BigInt(b);
        return res;
    }
    std::string toString(int base = 10) const {
        if (isZero()) return "0";
        static const char chars[] = "0123456789ABCDEF";
        BigInt n = *this;
        std::string s;
        BigInt bbase(base);
        while (!n.isZero()) {
            BigInt rem = n % bbase;
            s += chars[rem.digits[0]];
            n = n / bbase;
        }
        std::reverse(s.begin(), s.end());
        return s;
    }
    friend std::ostream& operator<<(std::ostream& os, const BigInt& bi) {
        os << bi.toString();
        return os;
    }
private:
    static BigInt fromString(const std::string& s, int base) {
        BigInt res(0);
        BigInt bbase(base);
        for (char c : s) {
            int val;
            if (c >= '0' && c <= '9') val = c - '0';
            else if (c >= 'A' && c <= 'F') val = c - 'A' + 10;
            else if (c >= 'a' && c <= 'f') val = c - 'a' + 10;
            else throw std::runtime_error("Invalid character");
            res = res * bbase + BigInt(val);
        }
        return res;
    }
    static void divide(const BigInt& dividend, const BigInt& divisor,
        BigInt& quotient, BigInt& remainder) {
        if (divisor.isZero()) throw std::runtime_error("Division by zero");
        if (dividend < divisor) {
            quotient = BigInt(0);
            remainder = dividend;
            return;
        }
        if (dividend.digits.size() <= 2 && divisor.digits.size() <= 2) {
            uint64_t a = dividend.digits[0];
            if (dividend.digits.size() > 1) a |= (uint64_t)dividend.digits[1] << BITS;
            uint64_t b = divisor.digits[0];
            if (divisor.digits.size() > 1) b |= (uint64_t)divisor.digits[1] << BITS;
            quotient = BigInt(a / b);
            remainder = BigInt(a % b);
            return;
        }
        int s = 0;
        BigInt v = divisor;
        while (v.digits.back() < (BASE / 2)) {
            v = v << 1;
            ++s;
        }
        BigInt u = dividend << s;
        u.digits.push_back(0);
        size_t m = u.digits.size() - v.digits.size();
        size_t n = v.digits.size();
        quotient.digits.assign(m, 0);
        Digit v_n_1 = v.digits.back();
        Digit v_n_2 = n >= 2 ? v.digits[v.digits.size() - 2] : 0;
        for (int j = (int)m - 1; j >= 0; --j) {
            uint64_t u_jn = ((uint64_t)u.digits[j + n] << BITS) | u.digits[j + n - 1];
            uint64_t q_hat = u_jn / v_n_1;
            uint64_t r_hat = u_jn % v_n_1;
            while (q_hat >= BASE ||
                (q_hat * v_n_2 > ((r_hat << BITS) |
                    (j + n - 2 < u.digits.size() ? u.digits[j + n - 2] : 0)))) {
                --q_hat;
                r_hat += v_n_1;
                if (r_hat >= BASE) break;
            }
            Digit carry = 0;
            for (size_t i = 0; i < n; ++i) {
                uint64_t prod = (uint64_t)q_hat * v.digits[i] + carry;
                Digit sub = (Digit)(prod % BASE);
                carry = (Digit)(prod / BASE);
                if (u.digits[j + i] < sub) {
                    u.digits[j + i] = (Digit)(u.digits[j + i] + BASE - sub);
                    ++carry;
                }
                else {
                    u.digits[j + i] -= sub;
                }
            }
            if (j + n < u.digits.size()) {
                if (u.digits[j + n] < carry) {
                    u.digits[j + n] = (Digit)(u.digits[j + n] + BASE - carry);
                    --q_hat;
                    carry = 0;
                    for (size_t i = 0; i < n; ++i) {
                        uint64_t sum = (uint64_t)u.digits[j + i] + v.digits[i] + carry;
                        u.digits[j + i] = (Digit)(sum % BASE);
                        carry = (Digit)(sum / BASE);
                    }
                    u.digits[j + n] += carry;
                }
                else {
                    u.digits[j + n] -= carry;
                }
            }
            else {
                if (carry) {
                    --q_hat;
                    carry = 0;
                    for (size_t i = 0; i < n; ++i) {
                        uint64_t sum = (uint64_t)u.digits[j + i] + v.digits[i] + carry;
                        u.digits[j + i] = (Digit)(sum % BASE);
                        carry = (Digit)(sum / BASE);
                    }
                }
            }

            quotient.digits[j] = (Digit)q_hat;
        }
        remainder.digits.assign(u.digits.begin(), u.digits.begin() + n);
        remainder = remainder >> s;

        while (quotient.digits.size() > 1 && quotient.digits.back() == 0)
            quotient.digits.pop_back();
        while (remainder.digits.size() > 1 && remainder.digits.back() == 0)
            remainder.digits.pop_back();
    }
};
class ElGamalSignature {
public:
    static void generateKeyPair(int bits, BigInt& p, BigInt& g, BigInt& y, BigInt& x) {
        static std::random_device rd;
        static std::mt19937_64 gen(rd());
        do {
            p = BigInt::random(bits, gen);
            p.digits[0] |= 1;              
        } while (!BigInt::isPrime(p, 3));  
        BigInt pMinus1 = p - BigInt(1);
        g = 2;
        while (true) {
            if (g.powMod(pMinus1, p) != BigInt(1)) {
                g = g + BigInt(1);
                continue;
            }
            if (g.powMod(pMinus1 >> 1, p) == BigInt(1)) {
                g = g + BigInt(1);
                continue;
            }
            break;
        }
        x = BigInt::random(bits - 1, gen) % (p - BigInt(2)) + BigInt(1);
        y = g.powMod(x, p);
    }
    static void sign(const std::vector<uint8_t>& message,
        const BigInt& p, const BigInt& g, const BigInt& x,
        BigInt& r, BigInt& s) {
        static std::random_device rd;
        static std::mt19937_64 gen(rd());
        BigInt h = BigInt::fromBytes(message);
        h = h % (p - BigInt(1));
        BigInt pMinus1 = p - BigInt(1);
        BigInt k;
        do {
            k = BigInt::random((int)pMinus1.bitLength(), gen) % (pMinus1 - BigInt(1)) + BigInt(1);
        } while (BigInt::gcd(k, pMinus1) != BigInt(1));
        r = g.powMod(k, p);
        BigInt xr = (x * r) % pMinus1;
        BigInt diff;
        if (h >= xr) diff = h - xr;
        else diff = pMinus1 - ((xr - h) % pMinus1);
        BigInt kInv = k.modInverse(pMinus1);
        s = (diff * kInv) % pMinus1;
    }
    static bool verify(const std::vector<uint8_t>& message,
        const BigInt& p, const BigInt& g, const BigInt& y,
        const BigInt& r, const BigInt& s) {
        if (r <= BigInt(0) || r >= p) return false;
        if (s <= BigInt(0) || s >= p - BigInt(1)) return false;
        BigInt h = BigInt::fromBytes(message);
        h = h % (p - BigInt(1));
        BigInt v1 = (y.powMod(r, p) * r.powMod(s, p)) % p;
        BigInt v2 = g.powMod(h, p);
        return v1 == v2;
    }
};
void runTests() {
    BigInt p, g, y, x;
    ElGamalSignature::generateKeyPair(64, p, g, y, x);
    assert(BigInt::isPrime(p, 3));
    assert(y == g.powMod(x, p));
    std::vector<uint8_t> msg = { 'H','e','l','l','o' };
    BigInt r, s;
    ElGamalSignature::sign(msg, p, g, x, r, s);
    assert(ElGamalSignature::verify(msg, p, g, y, r, s));
    std::vector<uint8_t> msg2 = { 'W','o','r','l','d' };
    assert(!ElGamalSignature::verify(msg2, p, g, y, r, s));
    BigInt p2, g2, y2, x2;
    ElGamalSignature::generateKeyPair(64, p2, g2, y2, x2);
    ElGamalSignature::sign(msg, p, g, x, r, s);
    assert(!ElGamalSignature::verify(msg, p2, g2, y2, r, s));
    std::cout << "PASS: All tests passed.\n";
}
int main() {
    runTests();
    BigInt p, g, y, x;
    ElGamalSignature::generateKeyPair(64, p, g, y, x);
    std::cout << "p = " << p << "\ng = " << g << "\ny = " << y << "\nx = " << x << "\n";
    std::string text = "ElGamal signature";
    std::vector<uint8_t> msg(text.begin(), text.end());
    BigInt r, s;
    ElGamalSignature::sign(msg, p, g, x, r, s);
    std::cout << "Signature: r = " << r << "\n          s = " << s << "\n";
    bool ok = ElGamalSignature::verify(msg, p, g, y, r, s);
    std::cout << "Verification: " << (ok ? "OK" : "FAIL") << "\n";
    msg.push_back('!');
    ok = ElGamalSignature::verify(msg, p, g, y, r, s);
    std::cout << "Verification of altered message: " << (ok ? "OK" : "FAIL") << "\n";
    return 0;
}