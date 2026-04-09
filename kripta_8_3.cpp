#include <iostream>
#include <bitset>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <array>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <thread>
#include <future>
#include <mutex>
#include <random>
#include <cstring>
#include <memory>
#include <functional>
#include <queue>
#include <condition_variable>
namespace poly_utils {
    template <size_t N>
    size_t degree(const std::bitset<N>& poly) {
        for (size_t i = N; i-- > 0; ) if (poly[i]) return i;
        return 0;
    }
    template <size_t N, size_t M>
    std::bitset<M> shift_left(const std::bitset<N>& poly, size_t shift) {
        std::bitset<M> result;
        for (size_t i = 0; i < N && i + shift < M; ++i) result[i + shift] = poly[i];
        return result;
    }
    template <size_t N, size_t M>
    std::pair<std::bitset<N>, std::bitset<N>> div_mod(const std::bitset<N>& a, const std::bitset<M>& b) {
        if (b.none()) throw std::domain_error("Division by zero polynomial");
        size_t deg_b = degree(b);
        std::bitset<N> q, r = a;
        while (r.any() && degree(r) >= deg_b) {
            size_t shift = degree(r) - deg_b;
            q.set(shift);
            r ^= shift_left<M, N>(b, shift);
        }
        return std::make_pair(q, r);
    }
    template <size_t N, size_t M, size_t Res>
    std::bitset<Res> multiply_raw(const std::bitset<N>& a, const std::bitset<M>& b) {
        std::bitset<Res> result;
        for (size_t i = 0; i < N; ++i) if (a[i]) result ^= shift_left<M, Res>(b, i);
        return result;
    }
    template <size_t N>
    std::bitset<N> gcd(std::bitset<N> a, std::bitset<N> b) {
        while (b.any()) {
            std::pair<std::bitset<N>, std::bitset<N>> div_res = div_mod(a, b);
            std::bitset<N> r = div_res.second;
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
        for (size_t i = 0; i < N && val; ++i, val >>= 1) poly_[i] = val & 1;
    }
    const std::bitset<N>& bits() const { return poly_; }
    GF2mPoly operator+(const GF2mPoly& o) const { return GF2mPoly(poly_ ^ o.poly_); }
    GF2mPoly& operator+=(const GF2mPoly& o) { poly_ ^= o.poly_; return *this; }
    GF2mPoly operator*(const GF2mPoly& o) const {
        constexpr size_t ProdSize = 2 * N - 1;
        auto prod = poly_utils::multiply_raw<N, N, ProdSize>(poly_, o.poly_);
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
    GF2mPoly pow(unsigned long long exp) const {
        GF2mPoly res(1), base = *this;
        while (exp) { if (exp & 1) res = res * base; base = base * base; exp >>= 1; }
        return res;
    }
    GF2mPoly inverse() const {
        if (poly_.none()) throw std::domain_error("Zero has no inverse");
        constexpr size_t BigSize = 2 * N + 1;
        using BigPoly = std::bitset<BigSize>;
        BigPoly a, b;
        for (size_t i = 0; i < N + 1; ++i) a[i] = modulus[i];
        for (size_t i = 0; i < N; ++i) b[i] = poly_[i];
        BigPoly old_r = a, r = b, old_s{ 1 }, s{ 0 }, old_t{ 0 }, t{ 1 };
        while (r.any()) {
            std::pair<BigPoly, BigPoly> div_res = poly_utils::div_mod(old_r, r);
            BigPoly q = div_res.first;
            BigPoly rem = div_res.second;
            old_r = r;
            r = rem;
            auto qs = poly_utils::multiply_raw<BigSize, BigSize, 2 * BigSize - 1>(q, s);
            BigPoly new_old_s = old_s;
            for (size_t i = 0; i < BigSize; ++i) new_old_s[i] = old_s[i] ^ qs[i];
            old_s = s;
            s = new_old_s;
            auto qt = poly_utils::multiply_raw<BigSize, BigSize, 2 * BigSize - 1>(q, t);
            BigPoly new_old_t = old_t;
            for (size_t i = 0; i < BigSize; ++i) new_old_t[i] = old_t[i] ^ qt[i];
            old_t = t;
            t = new_old_t;
        }
        std::pair<BigPoly, BigPoly> div_res2 = poly_utils::div_mod(old_t, a);
        BigPoly rem_final = div_res2.second;
        std::bitset<N> res;
        for (size_t i = 0; i < N; ++i) res[i] = rem_final[i];
        return GF2mPoly(res);
    }
    static bool is_irreducible(const std::bitset<N + 1>& p, size_t degree) {
        if (degree == 0) return false;
        if (degree == 1) return true;
        auto old = modulus; modulus = p;
        GF2mPoly x(2);
        auto xp = x.pow(1ULL << degree);
        if ((xp + x).bits().any()) { modulus = old; return false; }
        std::vector<size_t> divs; size_t n = degree;
        for (size_t d = 2; d * d <= n; ++d) if (n % d == 0) { divs.push_back(d); while (n % d == 0) n /= d; }
        if (n > 1) divs.push_back(n);
        for (size_t d : divs) {
            size_t q = degree / d;
            auto xq = x.pow(1ULL << q);
            auto cand = (xq + x).bits();
            std::bitset<N + 1> cp; cp.reset();
            for (size_t i = 0; i < N; ++i) cp[i] = cand[i];
            auto g = poly_utils::gcd(p, cp);
            if (g.count() != 1 || !g[0]) { modulus = old; return false; }
        }
        modulus = old; return true;
    }
};
template <size_t N> std::bitset<N + 1> GF2mPoly<N>::modulus;

class Rijndael {
public:
    enum class BlockSize { B128 = 16, B192 = 24, B256 = 32 };
    enum class KeySize { K128 = 16, K192 = 24, K256 = 32 };
private:
    BlockSize Nb_;
    KeySize Nk_;
    size_t Nr_;
    std::vector<uint8_t> key_;
    std::vector<uint32_t> w_;
    static const std::array<uint8_t, 256> sbox;
    static const std::array<uint8_t, 256> inv_sbox;
    static const std::array<uint8_t, 30> rcon;
    static std::array<uint8_t, 256> build_sbox() {
        GF2mPoly<8>::modulus = std::bitset<9>(0x11B);
        std::array<uint8_t, 256> box;
        for (unsigned x = 0; x < 256; ++x) {
            if (x == 0) box[x] = 0x63;
            else {
                GF2mPoly<8> p(x);
                uint8_t inv = static_cast<uint8_t>(p.inverse().bits().to_ulong());
                uint8_t res = 0;
                for (int i = 0; i < 8; ++i) {
                    uint8_t bit = ((inv >> i) & 1) ^ ((inv >> ((i + 4) % 8)) & 1) ^
                        ((inv >> ((i + 5) % 8)) & 1) ^ ((inv >> ((i + 6) % 8)) & 1) ^
                        ((inv >> ((i + 7) % 8)) & 1);
                    res |= (bit << i);
                }
                box[x] = res ^ 0x63;
            }
        }
        return box;
    }
    static std::array<uint8_t, 256> build_inv_sbox() {
        GF2mPoly<8>::modulus = std::bitset<9>(0x11B);
        std::array<uint8_t, 256> box;
        for (unsigned y = 0; y < 256; ++y) {
            uint8_t yaff = 0;
            for (int i = 0; i < 8; ++i) {
                uint8_t bit = ((y >> ((i + 2) % 8)) & 1) ^ ((y >> ((i + 5) % 8)) & 1) ^ ((y >> ((i + 7) % 8)) & 1);
                yaff |= (bit << i);
            }
            yaff ^= 0x05;
            if (yaff == 0) box[y] = 0;
            else {
                GF2mPoly<8> p(yaff);
                box[y] = static_cast<uint8_t>(p.inverse().bits().to_ulong());
            }
        }
        return box;
    }
    uint8_t xtime(uint8_t x) const { return (x << 1) ^ ((x & 0x80) ? 0x1B : 0x00); }
    void key_expansion() {
        size_t Nk_words = static_cast<size_t>(Nk_) / 4;
        size_t Nb_words = static_cast<size_t>(Nb_) / 4;
        w_.resize(Nb_words * (Nr_ + 1));
        for (size_t i = 0; i < Nk_words; ++i) {
            w_[i] = (key_[4 * i] << 24) | (key_[4 * i + 1] << 16) | (key_[4 * i + 2] << 8) | key_[4 * i + 3];
        }
        for (size_t i = Nk_words; i < Nb_words * (Nr_ + 1); ++i) {
            uint32_t temp = w_[i - 1];
            if (i % Nk_words == 0) {
                temp = (temp << 8) | (temp >> 24);
                temp = (sbox[(temp >> 24) & 0xFF] << 24) |
                    (sbox[(temp >> 16) & 0xFF] << 16) |
                    (sbox[(temp >> 8) & 0xFF] << 8) |
                    sbox[temp & 0xFF];
                temp ^= (rcon[i / Nk_words] << 24);
            }
            else if (Nk_words > 6 && i % Nk_words == 4) {
                temp = (sbox[(temp >> 24) & 0xFF] << 24) |
                    (sbox[(temp >> 16) & 0xFF] << 16) |
                    (sbox[(temp >> 8) & 0xFF] << 8) |
                    sbox[temp & 0xFF];
            }
            w_[i] = w_[i - Nk_words] ^ temp;
        }
    }
public:
    Rijndael(BlockSize bs, KeySize ks, const std::vector<uint8_t>& key)
        : Nb_(bs), Nk_(ks), key_(key)
    {
        if (key.size() != static_cast<size_t>(ks))
            throw std::invalid_argument("Key size mismatch");
        size_t Nb_words = static_cast<size_t>(Nb_) / 4;
        size_t Nk_words = static_cast<size_t>(Nk_) / 4;
        if (Nk_words == 4) Nr_ = 10;
        else if (Nk_words == 6) Nr_ = 12;
        else Nr_ = 14;
        key_expansion();
    }
    std::vector<uint8_t> encrypt_block(const std::vector<uint8_t>& block) const {
        if (block.size() != static_cast<size_t>(Nb_))
            throw std::invalid_argument("Block size mismatch");
        std::vector<uint8_t> state = block;
        size_t Nb = static_cast<size_t>(Nb_);
        size_t Nb_words = Nb / 4;
        for (size_t c = 0; c < Nb_words; ++c) {
            uint32_t w = w_[c];
            state[4 * c] ^= (w >> 24) & 0xFF;
            state[4 * c + 1] ^= (w >> 16) & 0xFF;
            state[4 * c + 2] ^= (w >> 8) & 0xFF;
            state[4 * c + 3] ^= w & 0xFF;
        }
        for (size_t round = 1; round < Nr_; ++round) {
            for (size_t i = 0; i < Nb; ++i) state[i] = sbox[state[i]];
            std::vector<uint8_t> temp(Nb);
            for (size_t r = 0; r < 4; ++r) {
                size_t shift = r;
                if (Nb == 32) {
                    if (r == 0) shift = 0;
                    else if (r == 1) shift = 1;
                    else if (r == 2) shift = 3;
                    else shift = 4;
                }
                else if (Nb == 24) {
                    if (r == 0) shift = 0;
                    else if (r == 1) shift = 1;
                    else if (r == 2) shift = 2;
                    else shift = 3;
                }
                for (size_t c = 0; c < Nb_words; ++c)
                    temp[4 * c + r] = state[4 * ((c + shift) % Nb_words) + r];
            }
            state = temp;
            for (size_t c = 0; c < Nb_words; ++c) {
                uint8_t a0 = state[4 * c], a1 = state[4 * c + 1], a2 = state[4 * c + 2], a3 = state[4 * c + 3];
                state[4 * c] = xtime(a0) ^ (xtime(a1) ^ a1) ^ a2 ^ a3;
                state[4 * c + 1] = a0 ^ xtime(a1) ^ (xtime(a2) ^ a2) ^ a3;
                state[4 * c + 2] = a0 ^ a1 ^ xtime(a2) ^ (xtime(a3) ^ a3);
                state[4 * c + 3] = (xtime(a0) ^ a0) ^ a1 ^ a2 ^ xtime(a3);
            }
            for (size_t c = 0; c < Nb_words; ++c) {
                uint32_t w = w_[round * Nb_words + c];
                state[4 * c] ^= (w >> 24) & 0xFF;
                state[4 * c + 1] ^= (w >> 16) & 0xFF;
                state[4 * c + 2] ^= (w >> 8) & 0xFF;
                state[4 * c + 3] ^= w & 0xFF;
            }
        }
        for (size_t i = 0; i < Nb; ++i) state[i] = sbox[state[i]];
        std::vector<uint8_t> temp(Nb);
        for (size_t r = 0; r < 4; ++r) {
            size_t shift = r;
            if (Nb == 32) {
                if (r == 0) shift = 0;
                else if (r == 1) shift = 1;
                else if (r == 2) shift = 3;
                else shift = 4;
            }
            else if (Nb == 24) {
                if (r == 0) shift = 0;
                else if (r == 1) shift = 1;
                else if (r == 2) shift = 2;
                else shift = 3;
            }
            for (size_t c = 0; c < Nb_words; ++c)
                temp[4 * c + r] = state[4 * ((c + shift) % Nb_words) + r];
        }
        state = temp;
        for (size_t c = 0; c < Nb_words; ++c) {
            uint32_t w = w_[Nr_ * Nb_words + c];
            state[4 * c] ^= (w >> 24) & 0xFF;
            state[4 * c + 1] ^= (w >> 16) & 0xFF;
            state[4 * c + 2] ^= (w >> 8) & 0xFF;
            state[4 * c + 3] ^= w & 0xFF;
        }
        return state;
    }
    std::vector<uint8_t> decrypt_block(const std::vector<uint8_t>& block) const {
        if (block.size() != static_cast<size_t>(Nb_))
            throw std::invalid_argument("Block size mismatch");
        std::vector<uint8_t> state = block;
        size_t Nb = static_cast<size_t>(Nb_);
        size_t Nb_words = Nb / 4;
        for (size_t c = 0; c < Nb_words; ++c) {
            uint32_t w = w_[Nr_ * Nb_words + c];
            state[4 * c] ^= (w >> 24) & 0xFF;
            state[4 * c + 1] ^= (w >> 16) & 0xFF;
            state[4 * c + 2] ^= (w >> 8) & 0xFF;
            state[4 * c + 3] ^= w & 0xFF;
        }
        for (size_t round = Nr_ - 1; round >= 1; --round) {
            std::vector<uint8_t> temp(Nb);
            for (size_t r = 0; r < 4; ++r) {
                size_t shift = r;
                if (Nb == 32) {
                    if (r == 0) shift = 0;
                    else if (r == 1) shift = 1;
                    else if (r == 2) shift = 3;
                    else shift = 4;
                }
                else if (Nb == 24) {
                    if (r == 0) shift = 0;
                    else if (r == 1) shift = 1;
                    else if (r == 2) shift = 2;
                    else shift = 3;
                }
                shift = (Nb_words - (shift % Nb_words)) % Nb_words;
                for (size_t c = 0; c < Nb_words; ++c)
                    temp[4 * c + r] = state[4 * ((c + shift) % Nb_words) + r];
            }
            state = temp;
            for (size_t i = 0; i < Nb; ++i) state[i] = inv_sbox[state[i]];
            for (size_t c = 0; c < Nb_words; ++c) {
                uint32_t w = w_[round * Nb_words + c];
                state[4 * c] ^= (w >> 24) & 0xFF;
                state[4 * c + 1] ^= (w >> 16) & 0xFF;
                state[4 * c + 2] ^= (w >> 8) & 0xFF;
                state[4 * c + 3] ^= w & 0xFF;
            }
            for (size_t c = 0; c < Nb_words; ++c) {
                uint8_t a0 = state[4 * c], a1 = state[4 * c + 1], a2 = state[4 * c + 2], a3 = state[4 * c + 3];
                auto mul = [this](uint8_t a, uint8_t b) -> uint8_t {
                    if (b == 0x0E) { uint8_t t = xtime(a); return xtime(t) ^ xtime(xtime(t)) ^ t; }
                    if (b == 0x0B) { uint8_t t = xtime(a); return a ^ t ^ xtime(xtime(t)); }
                    if (b == 0x0D) { uint8_t t = xtime(a); return a ^ xtime(t) ^ xtime(xtime(t)); }
                    if (b == 0x09) { uint8_t t = xtime(a); return a ^ xtime(xtime(t)); }
                    return 0;
                    };
                state[4 * c] = mul(a0, 0x0E) ^ mul(a1, 0x0B) ^ mul(a2, 0x0D) ^ mul(a3, 0x09);
                state[4 * c + 1] = mul(a0, 0x09) ^ mul(a1, 0x0E) ^ mul(a2, 0x0B) ^ mul(a3, 0x0D);
                state[4 * c + 2] = mul(a0, 0x0D) ^ mul(a1, 0x09) ^ mul(a2, 0x0E) ^ mul(a3, 0x0B);
                state[4 * c + 3] = mul(a0, 0x0B) ^ mul(a1, 0x0D) ^ mul(a2, 0x09) ^ mul(a3, 0x0E);
            }
        }
        std::vector<uint8_t> temp(Nb);
        for (size_t r = 0; r < 4; ++r) {
            size_t shift = r;
            if (Nb == 32) {
                if (r == 0) shift = 0;
                else if (r == 1) shift = 1;
                else if (r == 2) shift = 3;
                else shift = 4;
            }
            else if (Nb == 24) {
                if (r == 0) shift = 0;
                else if (r == 1) shift = 1;
                else if (r == 2) shift = 2;
                else shift = 3;
            }
            shift = (Nb_words - (shift % Nb_words)) % Nb_words;
            for (size_t c = 0; c < Nb_words; ++c)
                temp[4 * c + r] = state[4 * ((c + shift) % Nb_words) + r];
        }
        state = temp;
        for (size_t i = 0; i < Nb; ++i) state[i] = inv_sbox[state[i]];
        for (size_t c = 0; c < Nb_words; ++c) {
            uint32_t w = w_[c];
            state[4 * c] ^= (w >> 24) & 0xFF;
            state[4 * c + 1] ^= (w >> 16) & 0xFF;
            state[4 * c + 2] ^= (w >> 8) & 0xFF;
            state[4 * c + 3] ^= w & 0xFF;
        }
        return state;
    }
    size_t block_size() const { return static_cast<size_t>(Nb_); }
};
const std::array<uint8_t, 256> Rijndael::sbox = Rijndael::build_sbox();
const std::array<uint8_t, 256> Rijndael::inv_sbox = Rijndael::build_inv_sbox();
const std::array<uint8_t, 30> Rijndael::rcon = {
    0x00, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36,
    0x6C, 0xD8, 0xAB, 0x4D, 0x9A, 0x2F, 0x5E, 0xBC, 0x63, 0xC6, 0x97,
    0x35, 0x6A, 0xD4, 0xB3, 0x7D, 0xFA, 0xEF, 0xC5
};
class Padding {
public:
    enum class Scheme { Zeros, PKCS7, ISO10126, ANSIX923 };
    static std::vector<uint8_t> pad(const std::vector<uint8_t>& data, size_t block_size, Scheme scheme) {
        size_t pad_len = block_size - (data.size() % block_size);
        if (pad_len == 0) pad_len = block_size;
        std::vector<uint8_t> padded = data;
        padded.resize(data.size() + pad_len);
        if (scheme == Scheme::Zeros) {
        }
        else if (scheme == Scheme::PKCS7) {
            std::fill(padded.end() - pad_len, padded.end(), static_cast<uint8_t>(pad_len));
        }
        else if (scheme == Scheme::ISO10126) {
            static thread_local std::mt19937 rng(std::random_device{}());
            for (size_t i = 0; i < pad_len - 1; ++i)
                padded[data.size() + i] = rng() & 0xFF;
            padded.back() = static_cast<uint8_t>(pad_len);
        }
        else if (scheme == Scheme::ANSIX923) {
            std::fill(padded.end() - pad_len, padded.end() - 1, 0);
            padded.back() = static_cast<uint8_t>(pad_len);
        }
        return padded;
    }
    static std::vector<uint8_t> unpad(const std::vector<uint8_t>& data, Scheme scheme) {
        if (data.empty()) return {};
        uint8_t pad_len = data.back();
        if (pad_len == 0 || pad_len > data.size()) throw std::runtime_error("Invalid padding");
        if (scheme == Scheme::PKCS7) {
            for (size_t i = data.size() - pad_len; i < data.size(); ++i)
                if (data[i] != pad_len) throw std::runtime_error("PKCS7 padding invalid");
        }
        else if (scheme == Scheme::ANSIX923) {
            for (size_t i = data.size() - pad_len; i < data.size() - 1; ++i)
                if (data[i] != 0) throw std::runtime_error("ANSI X9.23 padding invalid");
        }
        else if (scheme == Scheme::ISO10126) {
        }
        else if (scheme == Scheme::Zeros) {
            return data; 
        }
        std::vector<uint8_t> unpadded(data.begin(), data.end() - pad_len);
        return unpadded;
    }
};
class CipherMode {
public:
    virtual ~CipherMode() = default;
    virtual std::vector<uint8_t> encrypt(const Rijndael& cipher, const std::vector<uint8_t>& plain, const std::vector<uint8_t>& iv) = 0;
    virtual std::vector<uint8_t> decrypt(const Rijndael& cipher, const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>& iv) = 0;
    virtual size_t iv_size(size_t block_size) const = 0;
};

class ECBMode : public CipherMode {
public:
    std::vector<uint8_t> encrypt(const Rijndael& cipher, const std::vector<uint8_t>& plain, const std::vector<uint8_t>&) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result;
        for (size_t i = 0; i < plain.size(); i += bs) {
            std::vector<uint8_t> block(plain.begin() + i, plain.begin() + i + bs);
            auto enc = cipher.encrypt_block(block);
            result.insert(result.end(), enc.begin(), enc.end());
        }
        return result;
    }
    std::vector<uint8_t> decrypt(const Rijndael& cipher, const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>&) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result;
        for (size_t i = 0; i < ciphertext.size(); i += bs) {
            std::vector<uint8_t> block(ciphertext.begin() + i, ciphertext.begin() + i + bs);
            auto dec = cipher.decrypt_block(block);
            result.insert(result.end(), dec.begin(), dec.end());
        }
        return result;
    }
    size_t iv_size(size_t) const override { return 0; }
};
class CBCMode : public CipherMode {
public:
    std::vector<uint8_t> encrypt(const Rijndael& cipher, const std::vector<uint8_t>& plain, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result, prev = iv;
        for (size_t i = 0; i < plain.size(); i += bs) {
            std::vector<uint8_t> block(plain.begin() + i, plain.begin() + i + bs);
            for (size_t j = 0; j < bs; ++j) block[j] ^= prev[j];
            auto enc = cipher.encrypt_block(block);
            result.insert(result.end(), enc.begin(), enc.end());
            prev = enc;
        }
        return result;
    }
    std::vector<uint8_t> decrypt(const Rijndael& cipher, const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result, prev = iv;
        for (size_t i = 0; i < ciphertext.size(); i += bs) {
            std::vector<uint8_t> block(ciphertext.begin() + i, ciphertext.begin() + i + bs);
            auto dec = cipher.decrypt_block(block);
            for (size_t j = 0; j < bs; ++j) dec[j] ^= prev[j];
            result.insert(result.end(), dec.begin(), dec.end());
            prev = block;
        }
        return result;
    }
    size_t iv_size(size_t bs) const override { return bs; }
};
class PCBCMode : public CipherMode {
public:
    std::vector<uint8_t> encrypt(const Rijndael& cipher, const std::vector<uint8_t>& plain, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result, prev_plain = iv, prev_cipher = iv;
        for (size_t i = 0; i < plain.size(); i += bs) {
            std::vector<uint8_t> block(plain.begin() + i, plain.begin() + i + bs);
            std::vector<uint8_t> feed(bs);
            for (size_t j = 0; j < bs; ++j) feed[j] = prev_plain[j] ^ prev_cipher[j];
            for (size_t j = 0; j < bs; ++j) block[j] ^= feed[j];
            auto enc = cipher.encrypt_block(block);
            result.insert(result.end(), enc.begin(), enc.end());
            prev_plain = block;
            prev_cipher = enc;
        }
        return result;
    }
    std::vector<uint8_t> decrypt(const Rijndael& cipher, const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result, prev_plain = iv, prev_cipher = iv;
        for (size_t i = 0; i < ciphertext.size(); i += bs) {
            std::vector<uint8_t> block(ciphertext.begin() + i, ciphertext.begin() + i + bs);
            auto dec = cipher.decrypt_block(block);
            std::vector<uint8_t> feed(bs);
            for (size_t j = 0; j < bs; ++j) feed[j] = prev_plain[j] ^ prev_cipher[j];
            for (size_t j = 0; j < bs; ++j) dec[j] ^= feed[j];
            result.insert(result.end(), dec.begin(), dec.end());
            prev_plain = dec;
            prev_cipher = block;
        }
        return result;
    }
    size_t iv_size(size_t bs) const override { return bs; }
};
class CFBMode : public CipherMode {
public:
    std::vector<uint8_t> encrypt(const Rijndael& cipher, const std::vector<uint8_t>& plain, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result(plain.size());
        std::vector<uint8_t> feedback = iv;
        for (size_t i = 0; i < plain.size(); i += bs) {
            auto keystream = cipher.encrypt_block(feedback);
            size_t chunk = std::min(bs, plain.size() - i);
            for (size_t j = 0; j < chunk; ++j) {
                result[i + j] = plain[i + j] ^ keystream[j];
                feedback[j] = result[i + j];
            }
        }
        return result;
    }
    std::vector<uint8_t> decrypt(const Rijndael& cipher, const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result(ciphertext.size());
        std::vector<uint8_t> feedback = iv;
        for (size_t i = 0; i < ciphertext.size(); i += bs) {
            auto keystream = cipher.encrypt_block(feedback);
            size_t chunk = std::min(bs, ciphertext.size() - i);
            for (size_t j = 0; j < chunk; ++j) {
                result[i + j] = ciphertext[i + j] ^ keystream[j];
                feedback[j] = ciphertext[i + j];
            }
        }
        return result;
    }
    size_t iv_size(size_t bs) const override { return bs; }
};
class OFBMode : public CipherMode {
public:
    std::vector<uint8_t> encrypt(const Rijndael& cipher, const std::vector<uint8_t>& plain, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result(plain.size());
        std::vector<uint8_t> feedback = iv;
        for (size_t i = 0; i < plain.size(); i += bs) {
            feedback = cipher.encrypt_block(feedback);
            size_t chunk = std::min(bs, plain.size() - i);
            for (size_t j = 0; j < chunk; ++j)
                result[i + j] = plain[i + j] ^ feedback[j];
        }
        return result;
    }
    std::vector<uint8_t> decrypt(const Rijndael& cipher, const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>& iv) override {
        return encrypt(cipher, ciphertext, iv);
    }
    size_t iv_size(size_t bs) const override { return bs; }
};
class CTRMode : public CipherMode {
public:
    std::vector<uint8_t> encrypt(const Rijndael& cipher, const std::vector<uint8_t>& plain, const std::vector<uint8_t>& nonce) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result(plain.size());
        std::vector<uint8_t> counter = nonce;
        size_t processed = 0;
        while (processed < plain.size()) {
            auto keystream = cipher.encrypt_block(counter);
            for (size_t j = 0; j < bs && processed + j < plain.size(); ++j)
                result[processed + j] = plain[processed + j] ^ keystream[j];
            for (int j = bs - 1; j >= 0; --j)
                if (++counter[j] != 0) break;
            processed += bs;
        }
        return result;
    }
    std::vector<uint8_t> decrypt(const Rijndael& cipher, const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>& nonce) override {
        return encrypt(cipher, ciphertext, nonce);
    }
    size_t iv_size(size_t bs) const override { return bs; }
};
class RandomDeltaMode : public CipherMode {
public:
    std::vector<uint8_t> encrypt(const Rijndael& cipher, const std::vector<uint8_t>& plain, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result;
        std::vector<uint8_t> prev = iv;
        static thread_local std::mt19937 rng(std::random_device{}());
        for (size_t i = 0; i < plain.size(); i += bs) {
            std::vector<uint8_t> block(plain.begin() + i, plain.begin() + i + bs);
            std::vector<uint8_t> delta(bs);
            for (auto& b : delta) b = rng() & 0xFF;
            auto enc_delta = cipher.encrypt_block(delta);
            for (size_t j = 0; j < bs; ++j) block[j] ^= enc_delta[j] ^ prev[j];
            auto enc = cipher.encrypt_block(block);
            result.insert(result.end(), delta.begin(), delta.end());
            result.insert(result.end(), enc.begin(), enc.end());
            prev = enc;
        }
        return result;
    }
    std::vector<uint8_t> decrypt(const Rijndael& cipher, const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>& iv) override {
        size_t bs = cipher.block_size();
        std::vector<uint8_t> result;
        std::vector<uint8_t> prev = iv;
        for (size_t i = 0; i < ciphertext.size(); i += 2 * bs) {
            std::vector<uint8_t> delta(ciphertext.begin() + i, ciphertext.begin() + i + bs);
            std::vector<uint8_t> block(ciphertext.begin() + i + bs, ciphertext.begin() + i + 2 * bs);
            auto dec = cipher.decrypt_block(block);
            auto enc_delta = cipher.encrypt_block(delta);
            for (size_t j = 0; j < bs; ++j) dec[j] ^= enc_delta[j] ^ prev[j];
            result.insert(result.end(), dec.begin(), dec.end());
            prev = block;
        }
        return result;
    }
    size_t iv_size(size_t bs) const override { return bs; }
};
class RijndaelEngine {
    Rijndael cipher_;
    std::unique_ptr<CipherMode> mode_;
    Padding::Scheme padding_;
    size_t num_threads_;
    class ThreadPool {
    public:
        explicit ThreadPool(size_t threads) : stop(false) {
            for (size_t i = 0; i < threads; ++i)
                workers.emplace_back([this] { worker_loop(); });
        }
        ~ThreadPool() {
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                stop = true;
            }
            condition.notify_all();
            for (std::thread& w : workers) w.join();
        }
        template<class F>
        auto enqueue(F&& f) -> std::future<typename std::result_of<F()>::type> {
            using return_type = typename std::result_of<F()>::type;
            auto task = std::make_shared<std::packaged_task<return_type()>>(std::forward<F>(f));
            std::future<return_type> res = task->get_future();
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                if (stop) throw std::runtime_error("enqueue on stopped ThreadPool");
                tasks.emplace([task]() { (*task)(); });
            }
            condition.notify_one();
            return res;
        }
    private:
        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;
        std::mutex queue_mutex;
        std::condition_variable condition;
        bool stop;
        void worker_loop() {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    condition.wait(lock, [this] { return stop || !tasks.empty(); });
                    if (stop && tasks.empty()) return;
                    task = std::move(tasks.front());
                    tasks.pop();
                }
                task();
            }
        }
    };
public:
    RijndaelEngine(Rijndael::BlockSize bs, Rijndael::KeySize ks, const std::vector<uint8_t>& key,
        std::unique_ptr<CipherMode> mode, Padding::Scheme pad, size_t threads = 1)
        : cipher_(bs, ks, key), mode_(std::move(mode)), padding_(pad), num_threads_(threads) {
        if (!mode_) throw std::invalid_argument("Cipher mode cannot be null");
    }
    std::vector<uint8_t> encrypt_array(const std::vector<uint8_t>& plain, const std::vector<uint8_t>& iv) {
        auto padded = Padding::pad(plain, cipher_.block_size(), padding_);
        if (num_threads_ <= 1 || padded.size() < cipher_.block_size() * 2) {
            return mode_->encrypt(cipher_, padded, iv);
        }
        CipherMode* mode_ptr = mode_.get();
        if (dynamic_cast<ECBMode*>(mode_ptr) || dynamic_cast<CTRMode*>(mode_ptr) || dynamic_cast<OFBMode*>(mode_ptr)) {
            size_t bs = cipher_.block_size();
            size_t num_blocks = padded.size() / bs;
            ThreadPool pool(num_threads_);
            std::vector<std::future<std::vector<uint8_t>>> futures;
            for (size_t i = 0; i < num_blocks; ++i) {
                futures.push_back(pool.enqueue([this, &padded, &iv, bs, i]() {
                    std::vector<uint8_t> block(padded.begin() + i * bs, padded.begin() + (i + 1) * bs);
                    if (auto* ecb = dynamic_cast<ECBMode*>(mode_.get())) {
                        return ecb->encrypt(cipher_, block, {});
                    }
                    else if (auto* ctr = dynamic_cast<CTRMode*>(mode_.get())) {
                        std::vector<uint8_t> nonce = iv;
                        uint64_t counter = 0;
                        for (size_t j = 0; j < bs; ++j) counter = (counter << 8) | nonce[j];
                        counter += i;
                        for (int j = bs - 1; j >= 0; --j) {
                            nonce[j] = counter & 0xFF;
                            counter >>= 8;
                        }
                        return ctr->encrypt(cipher_, block, nonce);
                    }
                    else if (auto* ofb = dynamic_cast<OFBMode*>(mode_.get())) {
                        return ofb->encrypt(cipher_, block, iv);
                    }
                    return std::vector<uint8_t>();
                    }));
            }
            std::vector<uint8_t> result;
            for (auto& f : futures) {
                auto part = f.get();
                result.insert(result.end(), part.begin(), part.end());
            }
            return result;
        }
        else {
            return mode_->encrypt(cipher_, padded, iv);
        }
    }
    std::vector<uint8_t> decrypt_array(const std::vector<uint8_t>& ciphertext, const std::vector<uint8_t>& iv) {
        if (num_threads_ <= 1 || ciphertext.size() < cipher_.block_size() * 2) {
            auto dec = mode_->decrypt(cipher_, ciphertext, iv);
            return Padding::unpad(dec, padding_);
        }
        CipherMode* mode_ptr = mode_.get();
        if (dynamic_cast<ECBMode*>(mode_ptr) || dynamic_cast<CTRMode*>(mode_ptr) || dynamic_cast<OFBMode*>(mode_ptr)) {
            size_t bs = cipher_.block_size();
            size_t num_blocks = ciphertext.size() / bs;
            ThreadPool pool(num_threads_);
            std::vector<std::future<std::vector<uint8_t>>> futures;
            for (size_t i = 0; i < num_blocks; ++i) {
                futures.push_back(pool.enqueue([this, &ciphertext, &iv, bs, i]() {
                    std::vector<uint8_t> block(ciphertext.begin() + i * bs, ciphertext.begin() + (i + 1) * bs);
                    if (auto* ecb = dynamic_cast<ECBMode*>(mode_.get())) {
                        return ecb->decrypt(cipher_, block, {});
                    }
                    else if (auto* ctr = dynamic_cast<CTRMode*>(mode_.get())) {
                        std::vector<uint8_t> nonce = iv;
                        uint64_t counter = 0;
                        for (size_t j = 0; j < bs; ++j) counter = (counter << 8) | nonce[j];
                        counter += i;
                        for (int j = bs - 1; j >= 0; --j) {
                            nonce[j] = counter & 0xFF;
                            counter >>= 8;
                        }
                        return ctr->decrypt(cipher_, block, nonce);
                    }
                    return std::vector<uint8_t>();
                    }));
            }
            std::vector<uint8_t> result;
            for (auto& f : futures) {
                auto part = f.get();
                result.insert(result.end(), part.begin(), part.end());
            }
            return Padding::unpad(result, padding_);
        }
        else {
            auto dec = mode_->decrypt(cipher_, ciphertext, iv);
            return Padding::unpad(dec, padding_);
        }
    }
    std::future<void> encrypt_file_async(const std::string& in_path, const std::string& out_path, const std::vector<uint8_t>& iv) {
        return std::async(std::launch::async, [this, in_path, out_path, iv]() {
            std::ifstream in(in_path, std::ios::binary);
            if (!in) throw std::runtime_error("Cannot open input file");
            std::ofstream out(out_path, std::ios::binary);
            std::vector<uint8_t> buffer((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            auto encrypted = encrypt_array(buffer, iv);
            out.write(reinterpret_cast<const char*>(encrypted.data()), encrypted.size());
            });
    }
    std::future<void> decrypt_file_async(const std::string& in_path, const std::string& out_path, const std::vector<uint8_t>& iv) {
        return std::async(std::launch::async, [this, in_path, out_path, iv]() {
            std::ifstream in(in_path, std::ios::binary);
            if (!in) throw std::runtime_error("Cannot open input file");
            std::ofstream out(out_path, std::ios::binary);
            std::vector<uint8_t> buffer((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            auto decrypted = decrypt_array(buffer, iv);
            out.write(reinterpret_cast<const char*>(decrypted.data()), decrypted.size());
            });
    }
};
void run_tests() {
    using namespace std;
    GF2mPoly<8>::modulus = bitset<9>(0x11B);
    GF2mPoly<8> a(0x57), b(0x83);
    assert((a * b).bits().to_ulong() == 0xC1);
    assert(a.inverse().bits().to_ulong() != 0);
    vector<uint8_t> key = { 0x2b,0x7e,0x15,0x16,0x28,0xae,0xd2,0xa6,0xab,0xf7,0x15,0x88,0x09,0xcf,0x4f,0x3c };
    vector<uint8_t> plain = { 0x32,0x43,0xf6,0xa8,0x88,0x5a,0x30,0x8d,0x31,0x31,0x98,0xa2,0xe0,0x37,0x07,0x34 };
    Rijndael aes(Rijndael::BlockSize::B128, Rijndael::KeySize::K128, key);
    auto cipher = aes.encrypt_block(plain);
    vector<uint8_t> expected = { 0x39,0x25,0x84,0x1d,0x02,0xdc,0x09,0xfb,0xdc,0x11,0x85,0x97,0x19,0x6a,0x0b,0x32 };
    assert(cipher == expected);
    auto dec = aes.decrypt_block(cipher);
    assert(dec == plain);
    vector<uint8_t> data = { 1,2,3,4,5 };
    vector<Padding::Scheme> schemes = { Padding::Scheme::PKCS7, Padding::Scheme::ISO10126, Padding::Scheme::ANSIX923 };
    for (auto scheme : schemes) {
        auto padded = Padding::pad(data, 8, scheme);
        assert(padded.size() == 8);
        auto unpadded = Padding::unpad(padded, scheme);
        assert(unpadded == data);
    }
    vector<uint8_t> iv(16, 0x42);
    vector<unique_ptr<CipherMode>> modes;
    modes.push_back(make_unique<ECBMode>());
    modes.push_back(make_unique<CBCMode>());
    modes.push_back(make_unique<CTRMode>());
    for (auto& mode : modes) {
        RijndaelEngine engine(Rijndael::BlockSize::B128, Rijndael::KeySize::K128, key,
            move(mode), Padding::Scheme::PKCS7, 1);
        vector<uint8_t> test_data = { 'T','e','s','t',' ','d','a','t','a' };
        auto encrypted = engine.encrypt_array(test_data, iv);
        auto decrypted = engine.decrypt_array(encrypted, iv);
        assert(decrypted == test_data);
    }
    vector<uint8_t> key192(24, 0x55);
    vector<uint8_t> key256(32, 0xAA);
    Rijndael aes192(Rijndael::BlockSize::B192, Rijndael::KeySize::K192, key192);
    Rijndael aes256(Rijndael::BlockSize::B256, Rijndael::KeySize::K256, key256);
    vector<uint8_t> block192(24, 0x11);
    vector<uint8_t> block256(32, 0x22);
    auto enc192 = aes192.encrypt_block(block192);
    auto dec192 = aes192.decrypt_block(enc192);
    assert(dec192 == block192);
    auto enc256 = aes256.encrypt_block(block256);
    auto dec256 = aes256.decrypt_block(enc256);
    assert(dec256 == block256);
    {
        auto ecb = make_unique<ECBMode>();
        RijndaelEngine engine(Rijndael::BlockSize::B128, Rijndael::KeySize::K128, key,
            move(ecb), Padding::Scheme::PKCS7, 4);
        vector<uint8_t> big_data(10000, 'X');
        auto encrypted = engine.encrypt_array(big_data, iv);
        auto decrypted = engine.decrypt_array(encrypted, iv);
        assert(decrypted == big_data);
    }
    cout << "All tests passed!\n";
}
int main() {
    run_tests();
    std::cout << "\nAll tests completed successfully.\n";
    return 0;
}