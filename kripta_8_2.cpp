#include <iostream>
#include <bitset>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <array>
#include <iomanip>
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
uint8_t affine_transform(uint8_t b) {
    uint8_t result = 0;
    for (int i = 0; i < 8; ++i) {
        uint8_t bit = ((b >> i) & 1) ^
            ((b >> ((i + 4) % 8)) & 1) ^
            ((b >> ((i + 5) % 8)) & 1) ^
            ((b >> ((i + 6) % 8)) & 1) ^
            ((b >> ((i + 7) % 8)) & 1);
        result |= (bit << i);
    }
    return result ^ 0x63;
}
uint8_t inv_affine_transform(uint8_t b) {
    uint8_t result = 0;
    for (int i = 0; i < 8; ++i) {
        uint8_t bit = ((b >> ((i + 2) % 8)) & 1) ^
            ((b >> ((i + 5) % 8)) & 1) ^
            ((b >> ((i + 7) % 8)) & 1);
        result |= (bit << i);
    }
    return result ^ 0x05;
}
using SBox = std::array<uint8_t, 256>;
SBox build_sbox(const std::bitset<9>& modulus) {
    using Poly = GF2mPoly<8>;
    Poly::modulus = modulus;
    SBox sbox;
    for (unsigned int x = 0; x < 256; ++x) {
        if (x == 0) {
            sbox[x] = affine_transform(0);
        }
        else {
            Poly p(x);
            Poly inv = p.inverse();
            uint8_t inv_byte = static_cast<uint8_t>(inv.bits().to_ulong());
            sbox[x] = affine_transform(inv_byte);
        }
    }
    return sbox;
}
SBox build_inv_sbox(const std::bitset<9>& modulus) {
    using Poly = GF2mPoly<8>;
    Poly::modulus = modulus;
    SBox inv_sbox;
    for (unsigned int y = 0; y < 256; ++y) {
        uint8_t y_aff = inv_affine_transform(static_cast<uint8_t>(y));
        if (y_aff == 0) {
            inv_sbox[y] = 0;
        }
        else {
            Poly p(y_aff);
            Poly inv = p.inverse();
            inv_sbox[y] = static_cast<uint8_t>(inv.bits().to_ulong());
        }
    }
    return inv_sbox;
}
void run_tests() {
    using Poly = GF2mPoly<8>;
    Poly::modulus = std::bitset<9>(0x11B); 
    std::cout << "Running unit tests...\n";
    {
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
    }
    {
        auto sbox = build_sbox(std::bitset<9>(0x11B));
        auto inv_sbox = build_inv_sbox(std::bitset<9>(0x11B));
        for (int i = 0; i < 256; ++i) {
            assert(inv_sbox[sbox[i]] == static_cast<uint8_t>(i));
            assert(sbox[inv_sbox[i]] == static_cast<uint8_t>(i));
        }
        assert(sbox[0x00] == 0x63);
        assert(sbox[0x01] == 0x7C);
        assert(sbox[0x02] == 0x77);
        assert(sbox[0x03] == 0x7B);
        assert(sbox[0x04] == 0xF2);
        assert(sbox[0x05] == 0x6B);
        assert(sbox[0x89] == 0xA7);
    }
    std::cout << "All tests passed!\n\n";
}
void print_sbox_partial(const SBox& sbox, const std::string& name, int rows = 8, int cols = 16) {
    std::cout << name << " (first " << rows << " rows):\n";
    for (int r = 0; r < rows; ++r) {
        std::cout << "  ";
        for (int c = 0; c < cols; ++c) {
            std::cout << std::hex << std::setw(2) << std::setfill('0')
                << static_cast<int>(sbox[r * 16 + c]) << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::dec;
}
int main() {
    run_tests();
    const size_t N = 8;
    using Poly = GF2mPoly<N>;
    Poly::modulus = std::bitset<N + 1>(0x11B);
    std::cout << "Demonstration of operations in GF(2^8)\n";
    std::cout << "Modulus m(x) = x^8 + x^4 + x^3 + x + 1 (0x11B)\n\n";
    auto sbox = build_sbox(Poly::modulus);
    auto inv_sbox = build_inv_sbox(Poly::modulus);
    print_sbox_partial(sbox, "AES S-Box", 8, 16);
    std::cout << "\n";
    print_sbox_partial(inv_sbox, "AES Inverse S-Box", 8, 16);
    std::cout << "\nExample: SBox[0x53] = 0x" << std::hex << static_cast<int>(sbox[0x53])
        << ", InvSBox[0x" << std::hex << static_cast<int>(sbox[0x53]) << "] = 0x"
        << std::hex << static_cast<int>(inv_sbox[sbox[0x53]]) << "\n";
    std::cout << "\nIs modulus 0x11B irreducible? " << std::boolalpha << Poly::is_irreducible(Poly::modulus, N) << "\n";
    std::bitset<N + 1> test_poly(0x1FF);
    std::cout << "Is polynomial 0x1FF irreducible? " << std::boolalpha << Poly::is_irreducible(test_poly, 8) << "\n";

    return 0;
}
