//#include <vector>
//#include <cstdint>
//#include <map>
//#include <functional>
//#include <stdexcept>
//#include <cassert>
//#include <iostream>
//class substitution_box {
//public:
//    static std::vector<uint8_t> substitute(
//        const std::vector<uint8_t>& input,
//        const std::map<uint8_t, uint8_t>& sbox
//    ) {
//        if (input.empty()) {
//            return std::vector<uint8_t>();
//        }
//        std::vector<uint8_t> output;
//        output.reserve(input.size());
//        for (size_t i = 0; i < input.size(); ++i) {
//            uint8_t byte = input[i];
//            uint8_t high = (byte >> 4) & 0x0f;//сдвиг на 4 бита получаем старший полубайт в мл разрядах, + маска (оставляем 4 млад бита)
//            uint8_t low = byte & 0x0f;//млад полубайт
//            auto it_high = sbox.find(high);//ищем в мап элемент с ключом хайт возвр итератор
//            if (it_high == sbox.end()) {
//                throw std::out_of_range("sbox missing entry for high nibble");
//            }
//            auto it_low = sbox.find(low);
//            if (it_low == sbox.end()) {
//                throw std::out_of_range("sbox missing entry for low nibble");
//            }
//            uint8_t new_byte = (it_high->second << 4) | (it_low->second & 0x0f);//значение для замены старш полкбайта и сдвиг влево - старшие 4 бита вых байта + побит или + млад полубайт с маской если число >15
//            output.push_back(new_byte);
//        }
//        return output;
//    }
//    static std::vector<uint8_t> substitute(
//        const std::vector<uint8_t>& input,
//        std::function<uint8_t(uint8_t)> func//перед по значению
//    ) {
//        if (input.empty()) {
//            return std::vector<uint8_t>();
//        }
//        std::vector<uint8_t> output;
//        output.reserve(input.size());
//        for (size_t i = 0; i < input.size(); ++i) {
//            uint8_t byte = input[i];
//            uint8_t high = (byte >> 4) & 0x0f;
//            uint8_t low = byte & 0x0f;
//            uint8_t new_high = func(high) & 0x0f;//чтоб не было лишних битов (знач в пределах полубайта)
//            uint8_t new_low = func(low) & 0x0f;
//            output.push_back((new_high << 4) | new_low);
//        }
//        return output;
//    }
//};
//static void test_substitution_map() {
//    std::vector<uint8_t> in = { 0xab, 0x12 };
//    std::map<uint8_t, uint8_t> sbox;
//    sbox[0xa] = 0x5;//заполн таблицу замены
//    sbox[0xb] = 0xe;
//    sbox[0x1] = 0x9;
//    sbox[0x2] = 0x3;
//    auto out = substitution_box::substitute(in, sbox);//комплитяр выбир 1 вар тк 2 арг - мап
//    assert(out.size() == 2);
//    assert(out[0] == 0x5e); 
//    assert(out[1] == 0x93); 
//}
//static void test_substitution_func() {
//    std::vector<uint8_t> in = { 0xab, 0x12 };
//    auto func = [](uint8_t x) -> uint8_t {//лям ф принимает инт8 и возвр инт8
//        return (x + 1) & 0x0f;//тело увел знач на 1 и маскирует мл 4 бита (каждый полубайт заменяется на след по модулю 16)
//        };
//    auto out = substitution_box::substitute(in, func);//лямбда неявно преобразуется в стдфанкшн
//    assert(out.size() == 2);
//    assert(out[0] == 0xbc);
//    assert(out[1] == 0x23);
//}
//static void test_substitution_map_missing_key() {//в мап нет для б проверка на искл
//    std::vector<uint8_t> in = { 0xab };
//    std::map<uint8_t, uint8_t> sbox;
//    sbox[0xa] = 0x5; 
//    bool caught = false;
//    try {
//        substitution_box::substitute(in, sbox);
//    }
//    catch (const std::out_of_range&) {
//        caught = true;
//    }
//    assert(caught);
//}
//static void test_substitution_empty_input() {//проверка на пустой массив
//    std::vector<uint8_t> in;
//    std::map<uint8_t, uint8_t> sbox;
//    auto out = substitution_box::substitute(in, sbox);
//    assert(out.empty());
//    auto func = [](uint8_t x) { return x; };
//    auto out2 = substitution_box::substitute(in, func);
//    assert(out2.empty());
//}
//static void demo() {
//    std::vector<uint8_t> data = { 0x3f, 0x8a };
//    std::map<uint8_t, uint8_t> sbox = {
//        {0x3, 0x7}, {0xf, 0x0}, {0x8, 0x2}, {0xa, 0x9}
//    };
//    auto res1 = substitution_box::substitute(data, sbox);
//    std::cout << "substitution with map: ";
//    for (auto b : res1) {
//        std::cout << std::hex << (int)b << " ";
//    }
//    std::cout << std::dec << std::endl;
//    auto func = [](uint8_t x) -> uint8_t {
//        return (x * 2) & 0x0f;//эквивал умнож по модулю 16 но с потерей старших бит
//        };
//    auto res2 = substitution_box::substitute(data, func);
//    std::cout << "substitution with function: ";
//    for (auto b : res2) {
//        std::cout << std::hex << (int)b << " ";
//    }
//    std::cout << std::dec << std::endl;
//}
//int main() {
//    test_substitution_map();
//    test_substitution_func();
//    test_substitution_map_missing_key();
//    test_substitution_empty_input();
//    demo();
//    std::cout << "all substitution tests passed" << std::endl;
//    return 0;
//}