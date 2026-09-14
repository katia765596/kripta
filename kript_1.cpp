#include <vector>
#include <cstdint>
#include <stdexcept>
#include <cassert>
#include <iostream>
enum class bit_order { lsb_first, msb_first };
enum class index_base { zero, one };
class bit_permutation {
public:
    static std::vector<uint8_t> permute(
        const std::vector<uint8_t>& input,
        const std::vector<int>& permutation,
        bit_order order,
        index_base base
    ) {
        size_t num_bits = input.size() * 8;//кол-во бит
        if (permutation.size() != num_bits) {
            throw std::invalid_argument("permutation size mismatch");
        }
        std::vector<uint8_t> output(input.size(), 0);
        for (size_t pos = 0; pos < num_bits; ++pos) {//pos-номер поз вых бита
            int raw = permutation[pos];//вых бит с номером pos берется из вход бита с таким-то номером
            if (base == index_base::one) {
                raw -= 1;//индексы в правиле на 1 больше чем рил номера
            }
            if (raw < 0 || static_cast<size_t>(raw) >= num_bits) {
                throw std::out_of_range("index out of range");
            }
            size_t src = static_cast<size_t>(raw);//номер бита во вход массиве но в нумер которая в правиле + база
            if (order == bit_order::msb_first) {
                size_t byte = src / 8;
                size_t bit = src % 8;
                src = byte * 8 + (7 - bit);//пересчитывает физ номер
            }
            uint8_t bit_val = (input[src / 8] >> (src % 8)) & 1;//сдвиг байта вправо на нужное число поз, чтобы нужный бит оказался на месте младшего и обнуляем все кроме младешго
            size_t dst = pos;//номер вых бита + преобр в физ номер
            if (order == bit_order::msb_first) {
                size_t byte = dst / 8;
                size_t bit = dst % 8;
                dst = byte * 8 + (7 - bit);
            }//запись бита в вых массив
            if (bit_val) {
                output[dst / 8] |= (1 << (dst % 8));//маска с 1 на нужной позиции + побит или устанавливает этот бит в 1
            }
            else {
                output[dst / 8] &= ~(1 << (dst % 8));//инверт маску (везде 1 кроме нужной там 0) и побит и обнуляет только указаггый бит
            }
        }
        return output;
    }
};
static void test_identity() {
    std::vector<uint8_t> in = { 0b10110010, 0b01010101 };
    size_t n = in.size() * 8;
    std::vector<int> p(n);
    for (size_t i = 0; i < n; ++i) p[i] = static_cast<int>(i);
    auto out = bit_permutation::permute(in, p, bit_order::lsb_first, index_base::zero);
    assert(out == in);
}
static void test_reverse_lsb() {
    std::vector<uint8_t> in = { 0b10110010 };
    std::vector<int> p(8);
    for (int i = 0; i < 8; ++i) p[i] = 7 - i;
    auto out = bit_permutation::permute(in, p, bit_order::lsb_first, index_base::zero);
    assert(out.size() == 1);
    assert(out[0] == 0b01001101);
}
static void test_reverse_msb() {
    std::vector<uint8_t> in = { 0b10110010 };
    std::vector<int> p(8);
    for (int i = 0; i < 8; ++i) p[i] = 7 - i;
    auto out = bit_permutation::permute(in, p, bit_order::msb_first, index_base::zero);
    assert(out.size() == 1);
    assert(out[0] == 0b01001101);
}
static void test_base_one() {
    std::vector<uint8_t> in = { 0b00000001 };
    std::vector<int> p(8);
    for (int i = 0; i < 8; ++i) p[i] = i + 1;
    auto out = bit_permutation::permute(in, p, bit_order::lsb_first, index_base::one);
    assert(out.size() == 1);
    assert(out[0] == 0b00000001);
}
static void test_wrong_size() {
    std::vector<uint8_t> in = { 0xaa };
    std::vector<int> p = { 0, 1 };
    bool ok = false;
    try {
        bit_permutation::permute(in, p, bit_order::lsb_first, index_base::zero);
    }
    catch (const std::invalid_argument&) {
        ok = true;
    }
    assert(ok);
}
static void test_out_of_range() {
    std::vector<uint8_t> in = { 0x00 };
    std::vector<int> p(8);
    for (int i = 0; i < 8; ++i) p[i] = 0;
    p[0] = 8;
    bool ok = false;
    try {
        bit_permutation::permute(in, p, bit_order::lsb_first, index_base::zero);
    }
    catch (const std::out_of_range&) {
        ok = true;
    }
    assert(ok);
}
static void test_two_bytes_msb() {
    std::vector<uint8_t> in = { 0x01, 0x80 };
    std::vector<int> p(16);
    for (int i = 0; i < 16; ++i) p[i] = 15 - i;
    auto out = bit_permutation::permute(in, p, bit_order::msb_first, index_base::zero);
    auto out2 = bit_permutation::permute(out, p, bit_order::msb_first, index_base::zero);
    assert(out2 == in);
}
static void demo() {
    std::vector<uint8_t> data = { 0b10110010 };
    std::cout << "original byte: " << std::hex << (int)data[0] << std::dec << std::endl;
    std::vector<int> rev(8);
    for (int i = 0; i < 8; ++i) rev[i] = 7 - i;
    auto res = bit_permutation::permute(data, rev, bit_order::lsb_first, index_base::zero);
    std::cout << "reversed (lsb): " << std::hex << (int)res[0] << std::dec << std::endl;
    std::vector<int> identity(8);
    for (int i = 0; i < 8; ++i) identity[i] = i;
    auto id_res = bit_permutation::permute(data, identity, bit_order::lsb_first, index_base::zero);
    std::cout << "identity: " << std::hex << (int)id_res[0] << std::dec << std::endl;
}
int main() {
    test_identity();
    test_reverse_lsb();
    test_reverse_msb();
    test_base_one();
    test_wrong_size();
    test_out_of_range();
    test_two_bytes_msb();
    demo();
    std::cout << "all tests passed" << std::endl;
    return 0;
}