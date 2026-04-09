#include <iostream>
#include <bitset>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <cassert>
namespace poly_utils {
    template <size_t N>
    size_t degree(const std::bitset<N>& poly) {
        for (size_t i = N; i-- > 0; ) {
            if (poly[i]) return i;
        }
        return 0;
    }
    template <size_t N, size_t M>
    std::bitset<M> shift_left(const std::bitset<N>& poly, size_t shift) {
        std::bitset<M> result;
        for (size_t i = 0; i < N && i + shift < M; ++i) {
            result[i + shift] = poly[i];
        }
        return result;
    }
    template <size_t N, size_t M>
    std::pair<std::bitset<N>, std::bitset<N>> div_mod(
        const std::bitset<N>& a,
        const std::bitset<M>& b)
    {
        if (b.none()) throw std::domain_error("Division by zero polynomial");
        size_t deg_b = degree(b);
        std::bitset<N> q, r = a;
        while (r.any() && degree(r) >= deg_b) {
            size_t shift = degree(r) - deg_b;
            q.set(shift);
            r ^= shift_left<M, N>(b, shift);
        }
        return { q, r };
    }
    template <size_t N, size_t M, size_t Res>
    std::bitset<Res> multiply_raw(const std::bitset<N>& a, const std::bitset<M>& b) {
        std::bitset<Res> result;
        for (size_t i = 0; i < N; ++i) {
            if (a[i]) {
                result ^= shift_left<M, Res>(b, i);
            }
        }
        return result;
    }
    template <size_t N>
    std::bitset<N> gcd(std::bitset<N> a, std::bitset<N> b) {
        while (b.any()) {
            auto [_, r] = div_mod(a, b);
            a = b;
            b = r;
        }
        return a;
    }
}
template <size_t N>
class GF2mPoly {
public:
    static std::bitset<N + 1> modulus;

private:
    std::bitset<N> poly_;
public:
    GF2mPoly() : poly_(0) {}
    explicit GF2mPoly(const std::bitset<N>& bits) : poly_(bits) {}
    explicit GF2mPoly(unsigned long long val) {
        for (size_t i = 0; i < N && val; ++i, val >>= 1) {
            poly_[i] = val & 1;
        }
    }
    const std::bitset<N>& bits() const { return poly_; }
    GF2mPoly operator+(const GF2mPoly& other) const {
        return GF2mPoly(poly_ ^ other.poly_);
    }
    GF2mPoly& operator+=(const GF2mPoly& other) {
        poly_ ^= other.poly_;
        return *this;
    }
    GF2mPoly operator*(const GF2mPoly& other) const {
        constexpr size_t ProdSize = 2 * N - 1;
        auto prod = poly_utils::multiply_raw<N, N, ProdSize>(poly_, other.poly_);
        std::bitset<ProdSize> r = prod;
        size_t deg_mod = poly_utils::degree(modulus);
        while (poly_utils::degree(r) >= deg_mod) {
            size_t shift = poly_utils::degree(r) - deg_mod;
            r ^= poly_utils::shift_left<N + 1, ProdSize>(modulus, shift);
        }
        std::bitset<N> reduced;
        for (size_t i = 0; i < N; ++i) reduced[i] = r[i];
        return GF2mPoly(reduced);
    }
    GF2mPoly& operator*=(const GF2mPoly& other) {
        *this = *this * other;
        return *this;
    }
    bool operator==(const GF2mPoly& other) const {
        return poly_ == other.poly_;
    }
    bool operator!=(const GF2mPoly& other) const {
        return poly_ != other.poly_;
    }
    GF2mPoly pow(unsigned long long exponent) const {
        GF2mPoly result(1ULL);
        GF2mPoly base = *this;
        while (exponent > 0) {
            if (exponent & 1) result = result * base;
            base = base * base;
            exponent >>= 1;
        }
        return result;
    }
    GF2mPoly inverse() const {
        if (poly_.none()) throw std::domain_error("Zero has no inverse");
        constexpr size_t BigSize = 2 * N + 1;
        using BigPoly = std::bitset<BigSize>;
        BigPoly a;
        for (size_t i = 0; i < N + 1; ++i) a[i] = modulus[i];
        BigPoly b;
        for (size_t i = 0; i < N; ++i) b[i] = poly_[i];
        BigPoly old_r = a, r = b;
        BigPoly old_s{ 1 }, s{ 0 };
        BigPoly old_t{ 0 }, t{ 1 };
        while (r.any()) {
            auto [q, rem] = poly_utils::div_mod(old_r, r);
            old_r = r;
            r = rem;
            auto qs = poly_utils::multiply_raw<BigSize, BigSize, 2 * BigSize - 1>(q, s);
            BigPoly new_old_s = old_s;
            for (size_t i = 0; i < BigSize; ++i)
                new_old_s[i] = old_s[i] ^ qs[i];
            old_s = s;
            s = new_old_s;
            auto qt = poly_utils::multiply_raw<BigSize, BigSize, 2 * BigSize - 1>(q, t);
            BigPoly new_old_t = old_t;
            for (size_t i = 0; i < BigSize; ++i)
                new_old_t[i] = old_t[i] ^ qt[i];
            old_t = t;
            t = new_old_t;
        }
        auto [_, rem_final] = poly_utils::div_mod(old_t, a);
        std::bitset<N> result;
        for (size_t i = 0; i < N; ++i) result[i] = rem_final[i];
        return GF2mPoly(result);
    }
    static bool is_irreducible(const std::bitset<N + 1>& p, size_t degree) {
        if (degree == 0) return false;
        if (degree == 1) return true;
        auto old_mod = modulus;
        modulus = p;
        GF2mPoly x(2ULL);
        auto x_pow = x.pow(1ULL << degree);
        if ((x_pow + x).bits().any()) {
            modulus = old_mod;
            return false;
        }
        std::vector<size_t> prime_divs;
        size_t n = degree;
        for (size_t d = 2; d * d <= n; ++d) {
            if (n % d == 0) {
                prime_divs.push_back(d);
                while (n % d == 0) n /= d;
            }
        }
        if (n > 1) prime_divs.push_back(n);
        for (size_t d : prime_divs) {
            size_t q = degree / d;
            auto x_pow_q = x.pow(1ULL << q);
            auto candidate = (x_pow_q + x).bits();
            std::bitset<N + 1> cand_poly;
            cand_poly.reset();
            for (size_t i = 0; i < N; ++i) cand_poly[i] = candidate[i];
            auto g = poly_utils::gcd(p, cand_poly);
            if (g.count() != 1 || !g[0]) {
                modulus = old_mod;
                return false;
            }
        }
        modulus = old_mod;
        return true;
    }
};
template <size_t N>
std::bitset<N + 1> GF2mPoly<N>::modulus;
template <size_t N>
void print_poly(const GF2mPoly<N>& p) {
    auto bits = p.bits();
    std::cout << "0b";
    for (size_t i = N; i-- > 0; ) std::cout << bits[i];
}
void run_tests() {
    using Poly = GF2mPoly<8>;
    Poly::modulus = std::bitset<9>(0x11B);
    std::cout << "Running unit tests...\n";
    Poly a1(0x89), b1(0x26);
    Poly sum = a1 + b1;
    assert(sum.bits().to_ulong() == (0x89 ^ 0x26));
    Poly a2(0x57), b2(0x83);
    Poly prod = a2 * b2;
    assert(prod.bits().to_ulong() == 0xC1);
    Poly a3(0x0B);
    Poly inv = a3.inverse();
    Poly check = a3 * inv;
    assert(check.bits().to_ulong() == 1);
    Poly one(1);
    assert(one.inverse().bits().to_ulong() == 1);
    bool caught = false;
    try {
        Poly zero(0);
        zero.inverse();
    }
    catch (const std::domain_error&) {
        caught = true;
    }
    assert(caught);
    assert(Poly::is_irreducible(Poly::modulus, 8) == true);
    std::bitset<9> reducible_poly(0x1FF);
    assert(Poly::is_irreducible(reducible_poly, 8) == false);
    std::cout << "All tests passed!\n\n";
}
int main() {
    run_tests();
    const size_t N = 8;
    using Poly = GF2mPoly<N>;
    Poly::modulus = std::bitset<N + 1>(0x11B);
    std::cout << "Demonstration of operations in GF(2^8) \n";
    std::cout << "Modulus m(x) = x^8 + x^4 + x^3 + x + 1 (0x11B)\n\n";
    Poly a(0x89), b(0x26);
    Poly c = a + b;
    std::cout << "a = "; print_poly(a); std::cout << " (0x" << std::hex << a.bits().to_ulong() << ")\n";
    std::cout << "b = "; print_poly(b); std::cout << " (0x" << std::hex << b.bits().to_ulong() << ")\n";
    std::cout << "a + b = "; print_poly(c); std::cout << " (0x" << std::hex << c.bits().to_ulong() << ")\n\n";
    Poly p1(0x53), p2(0x83);
    Poly prod = p1 * p2;
    std::cout << "p1 = "; print_poly(p1); std::cout << " (0x" << std::hex << p1.bits().to_ulong() << ")\n";
    std::cout << "p2 = "; print_poly(p2); std::cout << " (0x" << std::hex << p2.bits().to_ulong() << ")\n";
    std::cout << "p1 * p2 mod m(x) = "; print_poly(prod);
    std::cout << " (0x" << std::hex << prod.bits().to_ulong() << ")\n\n";
    Poly inv_test(0x0B);
    Poly inv = inv_test.inverse();
    Poly check = inv * inv_test;
    std::cout << "Element a = "; print_poly(inv_test); std::cout << " (0x" << std::hex << inv_test.bits().to_ulong() << ")\n";
    std::cout << "Inverse of a = "; print_poly(inv); std::cout << " (0x" << std::hex << inv.bits().to_ulong() << ")\n";
    std::cout << "Check a * a^{-1} = "; print_poly(check);
    std::cout << " (0x" << std::hex << check.bits().to_ulong() << ")\n\n";
    bool irr_mod = Poly::is_irreducible(Poly::modulus, N);
    std::cout << "Is modulus 0x11B irreducible? " << std::boolalpha << irr_mod << "\n";
    std::bitset<N + 1> test_poly(0x1FF);
    bool irr_test = Poly::is_irreducible(test_poly, 8);
    std::cout << "Is polynomial 0x1FF irreducible? " << irr_test << " (should be false)\n";
    return 0;
}