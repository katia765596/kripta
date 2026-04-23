#include "bigint.hpp"
#include <cassert>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <limits>
static const uint64_t BASE = (uint64_t(1) << 32);
void BigInt::trim() {
    while (!digits_.empty() && digits_.back() == 0)
        digits_.pop_back();
    if (digits_.empty()) negative_ = false;
}
int BigInt::compare_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    if (a.size() != b.size())
        return (a.size() < b.size()) ? -1 : 1;
    for (size_t i = a.size(); i-- > 0; ) {
        if (a[i] != b[i])
            return (a[i] < b[i]) ? -1 : 1;
    }
    return 0;
}
std::vector<uint32_t> BigInt::add_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    std::vector<uint32_t> res;
    res.reserve(std::max(a.size(), b.size()) + 1);
    uint64_t carry = 0;
    size_t i = 0;
    for (; i < a.size() || i < b.size() || carry; ++i) {
        uint64_t sum = carry;
        if (i < a.size()) sum += a[i];
        if (i < b.size()) sum += b[i];
        res.push_back(static_cast<uint32_t>(sum & 0xFFFFFFFF));
        carry = sum >> 32;
    }
    return res;
}
std::vector<uint32_t> BigInt::sub_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    std::vector<uint32_t> res;
    res.reserve(a.size());
    uint64_t borrow = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        uint64_t sub = a[i];
        sub -= borrow;
        if (i < b.size()) sub -= b[i];
        if (static_cast<int64_t>(sub) < 0) {
            sub += BASE;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res.push_back(static_cast<uint32_t>(sub));
    }
    while (!res.empty() && res.back() == 0)
        res.pop_back();
    return res;
}
std::vector<uint32_t> BigInt::mul_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    std::vector<uint32_t> res(a.size() + b.size(), 0);
    for (size_t i = 0; i < a.size(); ++i) {
        uint64_t carry = 0;
        for (size_t j = 0; j < b.size(); ++j) {
            uint64_t prod = uint64_t(a[i]) * b[j] + res[i + j] + carry;
            res[i + j] = static_cast<uint32_t>(prod & 0xFFFFFFFF);
            carry = prod >> 32;
        }
        if (carry)
            res[i + b.size()] += static_cast<uint32_t>(carry);
    }
    while (!res.empty() && res.back() == 0)
        res.pop_back();
    return res;
}
std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
BigInt::divmod_abs(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    if (b.empty() || (b.size() == 1 && b[0] == 0))
        throw std::domain_error("Division by zero");
    if (b.size() == 1) {
        uint64_t divisor = b[0];
        std::vector<uint32_t> q;
        uint64_t rem = 0;
        for (size_t i = a.size(); i-- > 0; ) {
            uint64_t cur = (rem << 32) | a[i];
            q.push_back(static_cast<uint32_t>(cur / divisor));
            rem = cur % divisor;
        }
        std::reverse(q.begin(), q.end());
        while (!q.empty() && q.back() == 0) q.pop_back();
        std::vector<uint32_t> r;
        if (rem != 0) r.push_back(static_cast<uint32_t>(rem));
        return { q, r };
    }
    if (compare_abs(a, b) < 0)
        return { std::vector<uint32_t>(), a };
    uint64_t norm = BASE / (uint64_t(b.back()) + 1);
    std::vector<uint32_t> norm_a = mul_abs(a, { static_cast<uint32_t>(norm) });
    std::vector<uint32_t> norm_b = mul_abs(b, { static_cast<uint32_t>(norm) });
    norm_a.push_back(0);
    size_t n = norm_b.size();
    if (n < 2) {
        auto res = divmod_abs(norm_a, { static_cast<uint32_t>(norm) });
        return { res.first, res.second };
    }
    size_t m = norm_a.size() - n - 1;
    std::vector<uint32_t> q(m + 1, 0);
    for (size_t i = m; i != (size_t)-1; --i) {
        uint64_t dividend = (uint64_t(norm_a[i + n]) << 32) | norm_a[i + n - 1];
        uint64_t q_hat = dividend / norm_b.back();
        uint64_t r_hat = dividend % norm_b.back();
        while (q_hat >= BASE ||
               (q_hat * uint64_t(norm_b[n - 2]) > (r_hat << 32) + norm_a[i + n - 2])) {
            --q_hat;
            r_hat += norm_b.back();
            if (r_hat >= BASE) break;
        }
        uint64_t carry = 0, borrow = 0;
        for (size_t j = 0; j < n; ++j) {
            uint64_t prod = q_hat * norm_b[j] + carry;
            carry = prod >> 32;
            prod &= 0xFFFFFFFF;
            uint64_t sub = norm_a[i + j];
            sub -= prod + borrow;
            if (static_cast<int64_t>(sub) < 0) {
                sub += BASE;
                borrow = 1;
            } else {
                borrow = 0;
            }
            norm_a[i + j] = static_cast<uint32_t>(sub);
        }
        if (borrow) {
            --q_hat;
            carry = 0;
            for (size_t j = 0; j < n; ++j) {
                uint64_t sum = uint64_t(norm_a[i + j]) + norm_b[j] + carry;
                norm_a[i + j] = static_cast<uint32_t>(sum & 0xFFFFFFFF);
                carry = sum >> 32;
            }
        }
        q[i] = static_cast<uint32_t>(q_hat);
    }
    norm_a.pop_back();
    auto div_res = divmod_abs(norm_a, { static_cast<uint32_t>(norm) });
    std::vector<uint32_t> r = div_res.second;
    while (!q.empty() && q.back() == 0) q.pop_back();
    while (!r.empty() && r.back() == 0) r.pop_back();
    return { q, r };
}
BigInt::BigInt() : negative_(false) {}
BigInt::BigInt(int64_t value) {
    negative_ = value < 0;
    if (negative_) value = -value;
    do {
        digits_.push_back(static_cast<uint32_t>(value & 0xFFFFFFFF));
        value >>= 32;
    } while (value > 0);
    trim();
}
BigInt::BigInt(const std::string& str) { from_string(str); }
void BigInt::from_string(const std::string& str) {
    negative_ = false;
    digits_.clear();
    if (str.empty()) return;
    size_t pos = 0;
    if (str[0] == '-') { negative_ = true; pos = 1; }
    BigInt res(0);
    BigInt ten(10);
    for (; pos < str.size(); ++pos) {
        char c = str[pos];
        if (c < '0' || c > '9') throw std::invalid_argument("Invalid char");
        res = res * ten + BigInt(c - '0');
    }
    *this = negative_ ? -res : res;
}
BigInt BigInt::operator+(const BigInt& other) const {
    if (negative_ == other.negative_) {
        BigInt res; res.digits_ = add_abs(digits_, other.digits_); res.negative_ = negative_; res.trim(); return res;
    }
    int cmp = compare_abs(digits_, other.digits_);
    if (cmp == 0) return BigInt(0);
    BigInt res;
    if (cmp > 0) { res.digits_ = sub_abs(digits_, other.digits_); res.negative_ = negative_; }
    else { res.digits_ = sub_abs(other.digits_, digits_); res.negative_ = other.negative_; }
    res.trim(); return res;
}
BigInt BigInt::operator-(const BigInt& other) const { BigInt neg = other; neg.negative_ = !neg.negative_; return *this + neg; }
BigInt BigInt::operator*(const BigInt& other) const {
    if (is_zero() || other.is_zero()) return BigInt(0);
    BigInt res; res.digits_ = mul_abs(digits_, other.digits_); res.negative_ = negative_ ^ other.negative_; res.trim(); return res;
}
BigInt BigInt::operator/(const BigInt& other) const { return divmod(*this, other).first; }
BigInt BigInt::operator%(const BigInt& other) const {
    auto p = divmod(*this, other);
    BigInt rem = p.second;
    if (rem.negative_ && !other.negative_) rem = rem + other;
    return rem;
}
BigInt BigInt::operator-() const { if (is_zero()) return *this; BigInt res = *this; res.negative_ = !negative_; return res; }
BigInt& BigInt::operator+=(const BigInt& other) { *this = *this + other; return *this; }
BigInt& BigInt::operator-=(const BigInt& other) { *this = *this - other; return *this; }
BigInt& BigInt::operator*=(const BigInt& other) { *this = *this * other; return *this; }
BigInt& BigInt::operator/=(const BigInt& other) { *this = *this / other; return *this; }
BigInt& BigInt::operator%=(const BigInt& other) { *this = *this % other; return *this; }

bool BigInt::operator==(const BigInt& other) const { return negative_ == other.negative_ && digits_ == other.digits_; }
bool BigInt::operator!=(const BigInt& other) const { return !(*this == other); }
bool BigInt::operator<(const BigInt& other) const {
    if (negative_ != other.negative_) return negative_;
    if (negative_) return compare_abs(digits_, other.digits_) > 0;
    return compare_abs(digits_, other.digits_) < 0;
}
bool BigInt::operator<=(const BigInt& other) const { return !(other < *this); }
bool BigInt::operator>(const BigInt& other) const { return other < *this; }
bool BigInt::operator>=(const BigInt& other) const { return !(*this < other); }
BigInt& BigInt::operator++() { *this += BigInt(1); return *this; }
BigInt BigInt::operator++(int) { BigInt tmp = *this; ++*this; return tmp; }
BigInt& BigInt::operator--() { *this -= BigInt(1); return *this; }
BigInt BigInt::operator--(int) { BigInt tmp = *this; --*this; return tmp; }
BigInt BigInt::operator<<(size_t bits) const {
    if (is_zero()) return *this;
    size_t word_shift = bits / 32;
    size_t bit_shift = bits % 32;
    std::vector<uint32_t> res(digits_.size() + word_shift + (bit_shift ? 1 : 0), 0);
    uint64_t carry = 0;
    for (size_t i = 0; i < digits_.size(); ++i) {
        uint64_t val = (uint64_t(digits_[i]) << bit_shift) | carry;
        carry = val >> 32;
        res[i + word_shift] = static_cast<uint32_t>(val & 0xFFFFFFFF);
    }
    if (carry) res.back() = static_cast<uint32_t>(carry);
    BigInt result; result.digits_ = std::move(res); result.negative_ = negative_; result.trim(); return result;
}
BigInt BigInt::operator>>(size_t bits) const {
    if (is_zero()) return BigInt(0);
    size_t word_shift = bits / 32;
    size_t bit_shift = bits % 32;
    if (word_shift >= digits_.size())
        return BigInt(0);
    std::vector<uint32_t> res(digits_.size() - word_shift, 0);
    uint64_t carry = 0;
    for (size_t i = digits_.size(); i-- > word_shift;) {
        uint64_t cur = digits_[i];
        uint64_t val = (cur >> bit_shift) | carry;
        carry = (cur << (32 - bit_shift)) & 0xFFFFFFFF;
        res[i - word_shift] = static_cast<uint32_t>(val);
    }
    BigInt result;
    result.digits_ = res;
    result.negative_ = negative_;
    result.trim();
    return result;
}
BigInt& BigInt::operator<<=(size_t bits) { *this = *this << bits; return *this; }
BigInt& BigInt::operator>>=(size_t bits) { *this = *this >> bits; return *this; }
std::string BigInt::to_string() const {
    if (digits_.empty()) return "0";
    BigInt ten(10);
    BigInt num = this->abs();
    std::string digits;
    while (!num.is_zero()) {
        auto div_res = divmod(num, ten);
        int64_t digit = div_res.second.to_int64();
        digits.push_back('0' + static_cast<char>(digit));
        num = div_res.first;
    }
    if (negative_) digits.push_back('-');
    std::reverse(digits.begin(), digits.end());
    return digits;
}
int64_t BigInt::to_int64() const {
    if (digits_.empty()) return 0;
    uint64_t val = 0;
    for (size_t i = 0; i < digits_.size() && i < 2; ++i) {
        if (i == 1 && (digits_[1] & 0x80000000)) throw std::overflow_error("overflow");
        val |= uint64_t(digits_[i]) << (32 * i);
    }
    if (negative_) {
        if (val > uint64_t(std::numeric_limits<int64_t>::max()) + 1) throw std::overflow_error("overflow");
        return -static_cast<int64_t>(val);
    } else {
        if (val > uint64_t(std::numeric_limits<int64_t>::max())) throw std::overflow_error("overflow");
        return static_cast<int64_t>(val);
    }
}
uint64_t BigInt::to_uint64() const {
    if (negative_ || digits_.size() > 2) throw std::overflow_error("not convertible to uint64_t");
    uint64_t val = 0;
    for (size_t i = 0; i < digits_.size(); ++i) val |= uint64_t(digits_[i]) << (32 * i);
    return val;
}
bool BigInt::is_zero() const { return digits_.empty(); }
bool BigInt::is_odd() const { return !digits_.empty() && (digits_[0] & 1); }
bool BigInt::is_negative() const { return negative_; }
int BigInt::sign() const { return is_zero() ? 0 : (negative_ ? -1 : 1); }
BigInt BigInt::abs() const { BigInt res = *this; res.negative_ = false; return res; }
std::pair<BigInt, BigInt> BigInt::divmod(const BigInt& a, const BigInt& b) {
    if (b.is_zero()) throw std::domain_error("Division by zero");
    if (a.is_zero()) return { BigInt(0), BigInt(0) };
    bool neg_q = a.negative_ ^ b.negative_;
    bool neg_r = a.negative_;
    auto res_abs = divmod_abs(a.digits_, b.digits_);
    BigInt q, r;
    q.digits_ = std::move(res_abs.first);
    q.negative_ = neg_q && !q.digits_.empty();
    r.digits_ = std::move(res_abs.second);
    r.negative_ = neg_r && !r.digits_.empty();
    q.trim(); r.trim();
    if (!r.is_zero() && r.negative_) {
        r = r + b.abs();
        q = q - 1;
    }
    if (r.abs() >= b.abs()) {
        r = r - b.abs();
        q = q + 1;
    }
    return { q, r };
}
BigInt BigInt::mod_pow(BigInt base, BigInt exp, const BigInt& mod) {
    if (mod == 1) return BigInt(0);
    BigInt result(1);
    base = base % mod;
    while (exp > 0) {
        if (exp.is_odd()) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp = exp >> 1;
    }
    return result;
}
BigInt BigInt::mod_inverse(const BigInt& a, const BigInt& mod) {
    BigInt t(0), newt(1), r = mod, newr = a % mod;
    if (newr < 0) newr += mod;
    while (!newr.is_zero()) {
        BigInt q = r / newr;
        BigInt tmp = t; t = newt; newt = tmp - q * newt;
        tmp = r; r = newr; newr = tmp - q * newr;
    }
    if (r > 1) throw std::runtime_error("Not invertible");
    if (t < 0) t += mod;
    return t;
}
static std::mt19937_64& rng() {
    static std::mt19937_64 gen(
        static_cast<uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count()));
    return gen;
}
BigInt BigInt::random(const BigInt& min, const BigInt& max) {
    if (min > max) throw std::invalid_argument("min > max");
    BigInt range = max - min + 1;
    size_t bits = 0;
    BigInt tmp = range;
    while (tmp > 0) { ++bits; tmp >>= 1; }
    BigInt result;
    auto& gen = rng();
    do {
        result = 0;
        for (size_t i = 0; i < bits; i += 64) {
            uint64_t rnd = gen();
            result = (result << 64) + BigInt(static_cast<int64_t>(rnd));
        }
        result = result % range;
    } while (result + min > max);
    return result + min;
}
bool BigInt::is_prime(const BigInt& n, int certainty) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;
    const std::vector<int> small_primes = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};
    for (int p : small_primes) {
        if (n == BigInt(p)) return true;
        if (n % p == 0) return false;
    }
    BigInt d = n - 1;
    int s = 0;
    while (d % 2 == 0) { d = d / 2; ++s; }
    for (int i = 0; i < certainty; ++i) {
        BigInt a = random(BigInt(2), n - 2);
        BigInt x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s; ++r) {
            x = (x * x) % n;
            if (x == n - 1) { composite = false; break; }
            if (x == 1) break;
        }
        if (composite) return false;
    }
    return true;
}
BigInt BigInt::generate_prime(int bits) {
    if (bits < 2) throw std::invalid_argument("bits >= 2");
    BigInt min = BigInt(1) << (bits - 1);
    BigInt max = (BigInt(1) << bits) - 1;
    while (true) {
        BigInt candidate = random(min, max);
        if (candidate % 2 == 0) candidate += 1;
        if (is_prime(candidate, 20)) return candidate;
    }
}
BigInt BigInt::center_lift(const BigInt& x, const BigInt& mod) {
    BigInt r = x % mod;
    if (r < 0) r += mod;
    BigInt half = mod / 2;
    if (r > half) r -= mod;
    return r;
}