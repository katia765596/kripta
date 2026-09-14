//#include <vector>
//#include <cstdint>
//#include <stdexcept>
//#include <cassert>
//#include <iostream>
//class bit_operations {
//private:
//    static bool get_bit(const std::vector<uint8_t>& data, size_t n, size_t pos) {
//        if (pos >= n) {
//            throw std::out_of_range("bit position out of range");
//        }
//        size_t byte_idx = pos / 8;
//        size_t bit_idx = pos % 8;
//        return (data[byte_idx] >> bit_idx) & 1;//сдвиг байта вправо чтобы нужный бит был на поз младшего разряда + побитовое И с 1 - получаем 0 или 1
//    }
//    static void set_bit(std::vector<uint8_t>& data, size_t n, size_t pos, int val) {
//        if (pos >= n) {//бит pos в знач val (0 или 1) в массиве data 
//            throw std::out_of_range("bit position out of range");
//        }
//        size_t byte_idx = pos / 8;
//        size_t bit_idx = pos % 8;
//        if (val) {//vall=1 устанавить  бит в 1. маска с 1 на поз bit_idx + побитовое или с текущим байтом устанавл этот бит в 1, остальные не меняются
//            data[byte_idx] |= (1 << bit_idx);
//        }
//        else {
//            data[byte_idx] &= ~(1 << bit_idx);//если vall=0, обнуляем бит, маска с нулем на бит_идх и единицами на других разрдах, побит И с маской сбрасывает нужный бит в 9
//        }
//    }
//    static size_t byte_count(size_t bits) {
//        return (bits + 7) / 8;//округ вверх получаем колво байт для всех битов
//    }
//    static void trim_extra_bits(std::vector<uint8_t>& data, size_t n) {//обнул лишнние биты в посл байте
//        size_t extra = n % 8;
//        if (extra != 0) {
//            data.back() &= (1 << extra) - 1;//в посл байте побит И с маской (маска из extra единиц в мл разрядах) обнуляем все биты, с ехстра и выше
//        }
//    }
//public:
//    static std::vector<uint8_t> rotate_left(const std::vector<uint8_t>& data, size_t n, size_t k) {
//        if (n == 0) return std::vector<uint8_t>();
//        k = k % n;
//        if (k == 0) return data;
//        size_t bc = byte_count(n);
//        if (data.size() < bc) throw std::invalid_argument("data too short");
//        std::vector<uint8_t> out(bc, 0);
//        for (size_t pos = 0; pos < n; ++pos) {
//            size_t new_pos = (pos + k) % n;//куда перейдет бит после сдвига влево на к
//            int val = get_bit(data, n, pos);//получаем val - знач бита на pos из вход массива
//            set_bit(out, n, new_pos, val);//записывает знач в вых массив
//        }
//        trim_extra_bits(out, n);
//        return out;
//    }
//    static std::vector<uint8_t> rotate_right(const std::vector<uint8_t>& data, size_t n, size_t k) {
//        if (n == 0) return std::vector<uint8_t>();
//        k = k % n;
//        if (k == 0) return data;
//        size_t bc = byte_count(n);
//        if (data.size() < bc) throw std::invalid_argument("data too short");
//        std::vector<uint8_t> out(bc, 0);
//        for (size_t pos = 0; pos < n; ++pos) {
//            size_t new_pos = (pos + n - k) % n;//цикл сдвиг вправо на к
//            int val = get_bit(data, n, pos);
//            set_bit(out, n, new_pos, val);
//        }
//        trim_extra_bits(out, n);
//        return out;
//    }
//    static std::vector<uint8_t> apply_mask(const std::vector<uint8_t>& data, size_t n,
//        const std::vector<uint8_t>& mask, size_t k) {//обнуляет все биты, кроме тех где в маске 1
//        if (n == 0) return std::vector<uint8_t>();
//        if (k > n) throw std::invalid_argument("mask length exceeds value length");
//        size_t bc = byte_count(n);
//        if (data.size() < bc) throw std::invalid_argument("data too short");
//        size_t mask_bc = byte_count(k);
//        if (mask.size() < mask_bc) throw std::invalid_argument("mask too short");
//        std::vector<uint8_t> out(data.begin(), data.begin() + bc);
//        for (size_t pos = 0; pos < k; ++pos) {//проходим по всем битам маски
//            int mask_bit = get_bit(mask, k, pos);//получаем 0 или 1
//            if (mask_bit == 0) {//если бит маски =0, то соот бит вых знач обнуляется
//                set_bit(out, n, pos, 0);
//            }
//        }
//        for (size_t pos = k; pos < n; ++pos) {
//            set_bit(out, n, pos, 0);//все биты от к до n-1 обнуляем (маска их не покрыла)
//        }
//        trim_extra_bits(out, n);
//        return out;
//    }
//    static std::vector<uint8_t> extract_bits(const std::vector<uint8_t>& data, size_t n,
//        size_t i, size_t j) {//извлечение подстроки битов от i до j
//        if (i > j) throw std::invalid_argument("i must be <= j");
//        if (j >= n) throw std::out_of_range("j out of range");
//        size_t len = j - i + 1;
//        size_t bc = byte_count(n);
//        if (data.size() < bc) throw std::invalid_argument("data too short");
//        std::vector<uint8_t> out(byte_count(len), 0);
//        for (size_t pos = 0; pos < len; ++pos) {
//            int val = get_bit(data, n, i + pos);//для каждой pos в вых знач читаем бит из исх по индексу i+pos и записываем в вых на pos
//            set_bit(out, len, pos, val);
//        }
//        trim_extra_bits(out, len);
//        return out;
//    }
//    static std::vector<uint8_t> swap_bits(const std::vector<uint8_t>& data, size_t n,
//        size_t i, size_t j) {//обмен местами двух бит
//        if (i >= n || j >= n) throw std::out_of_range("bit index out of range");
//        size_t bc = byte_count(n);
//        if (data.size() < bc) throw std::invalid_argument("data too short");
//        std::vector<uint8_t> out(data.begin(), data.begin() + bc);
//        int bi = get_bit(out, n, i);
//        int bj = get_bit(out, n, j);
//        if (bi != bj) {
//            set_bit(out, n, i, bj);//в поз i знач bj
//            set_bit(out, n, j, bi);//в поз j знач bi
//        }
//        trim_extra_bits(out, n);
//        return out;
//    }
//    static std::vector<uint8_t> set_bit_value(const std::vector<uint8_t>& data, size_t n,
//        size_t i, int value) {//установка бита с номером i в заданное значение 0 или 1
//        if (i >= n) throw std::out_of_range("bit index out of range");
//        if (value != 0 && value != 1) throw std::invalid_argument("value must be 0 or 1");
//        size_t bc = byte_count(n);
//        if (data.size() < bc) throw std::invalid_argument("data too short");
//        std::vector<uint8_t> out(data.begin(), data.begin() + bc);
//        set_bit(out, n, i, value);//устанав нужный бит
//        trim_extra_bits(out, n);
//        return out;
//    }
//};
//static void test_rotate_left() {
//    std::vector<uint8_t> in = { 0b11010011 };//211 0xD3
//    auto out = bit_operations::rotate_left(in, 8, 3);
//    assert(out.size() == 1 && out[0] == 0b10011110);//0x9E
//}
//static void test_rotate_right() {
//    std::vector<uint8_t> in = { 0b11010011 };//OxD3
//    auto out = bit_operations::rotate_right(in, 8, 3);
//    assert(out.size() == 1 && out[0] == 0b01111010);//0x7A
//}
//static void test_apply_mask() {
//    std::vector<uint8_t> in = { 0b10101100 };
//    std::vector<uint8_t> mask = { 0b1010 };
//    auto out = bit_operations::apply_mask(in, 8, mask, 4);
//    assert(out.size() == 1 && out[0] == 0b00001000);//младшие энд с маской+ старшие обнул
//}
//static void test_extract_bits() {
//    std::vector<uint8_t> in = { 0b10110101 };//1101 с 2 по 5 бит извлечь
//    auto out = bit_operations::extract_bits(in, 8, 2, 5);
//    assert(out.size() == 1 && out[0] == 0b00001101);
//}
//static void test_swap_bits() {
//    std::vector<uint8_t> in = { 0b10101010 };
//    auto out = bit_operations::swap_bits(in, 8, 0, 7);
//    assert(out.size() == 1 && out[0] == 0b00101011);
//}
//static void test_set_bit() {
//    std::vector<uint8_t> in = { 0b11110000 };
//    auto out = bit_operations::set_bit_value(in, 8, 3, 0);
//    assert(out[0] == 0b11110000);
//    auto out2 = bit_operations::set_bit_value(in, 8, 3, 1);
//    assert(out2[0] == 0b11111000);
//}
//static void demo() {
//    std::vector<uint8_t> data = { 0b10101100, 0b11110000 };
//    auto rl = bit_operations::rotate_left(data, 16, 5);
//    std::cout << "rotate left 5: ";
//    for (auto b : rl) std::cout << std::hex << (int)b << " ";
//    std::cout << std::dec << std::endl;
//    auto rr = bit_operations::rotate_right(data, 16, 5);
//    std::cout << "rotate right 5: ";
//    for (auto b : rr) std::cout << std::hex << (int)b << " ";
//    std::cout << std::dec << std::endl;
//    std::vector<uint8_t> mask = { 0b1100, 0b1010 };
//    auto masked = bit_operations::apply_mask(data, 16, mask, 8);
//    std::cout << "apply mask (8 bits): ";
//    for (auto b : masked) std::cout << std::hex << (int)b << " ";
//    std::cout << std::dec << std::endl;
//    auto ext = bit_operations::extract_bits(data, 16, 3, 10);
//    std::cout << "extract bits 3..10: ";
//    for (auto b : ext) std::cout << std::hex << (int)b << " ";
//    std::cout << std::dec << std::endl;
//    auto swapped = bit_operations::swap_bits(data, 16, 0, 15);
//    std::cout << "swap bit 0 and 15: ";
//    for (auto b : swapped) std::cout << std::hex << (int)b << " ";
//    std::cout << std::dec << std::endl;
//    auto set0 = bit_operations::set_bit_value(data, 16, 7, 0);
//    std::cout << "set bit 7 to 0: ";
//    for (auto b : set0) std::cout << std::hex << (int)b << " ";
//    std::cout << std::dec << std::endl;
//}
//int main() {
//    test_rotate_left();
//    test_rotate_right();
//    test_apply_mask();
//    test_extract_bits();
//    test_swap_bits();
//    test_set_bit();
//    demo();
//    std::cout << "All bit operations tests passed" << std::endl;
//    return 0;
//}