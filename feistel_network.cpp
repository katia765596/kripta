#include "feistel_network.h"
#include <cstring>
#include <stdexcept>
feistel_network::feistel_network(i_key_schedule* ks, i_feistel_round* fr, int rounds)//1д
    : key_sched(ks), round_func(fr), num_rounds(rounds) {//сохр переданные объекты в приватные поля класса
}
std::vector<uint8_t> feistel_network::process_block(const std::vector<uint8_t>& block, bool encrypt, const std::vector<uint8_t>& key) {
    if (block.size() != 8) throw std::invalid_argument("block size must be 8 bytes");
    std::vector<std::vector<uint8_t>> round_keys = key_sched->expand_key(key);//вызываем метод экспенд у переданного объект кей шед (ген вектор раунд ключей)
    if (round_keys.size() != num_rounds) throw std::runtime_error("wrong number of round keys");
    std::vector<uint8_t> left(4), right(4);//л и п половина блока по 4 б
    std::memcpy(left.data(), block.data(), 4);
    std::memcpy(right.data(), block.data() + 4, 4);//с конца, п половина
    for (int r = 0; r < num_rounds; ++r) {
        int idx = encrypt ? r : (num_rounds - 1 - r);//опр номер раунд ключа если шифр по порядку если нет то в обр
        std::vector<uint8_t> round_key = round_keys[idx];
        std::vector<uint8_t> new_right = round_func->round_function(right, round_key);//раунд ф у объекта раунд фанк
        for (int i = 0; i < 4; ++i) new_right[i] ^= left[i];
        left = right;
        right = new_right;//меняем + это новая левая половина для след раунда еще
    }
    std::vector<uint8_t> result(8);
    std::memcpy(result.data(), right.data(), 4);//в конце раунда мы меняли местами, поэтому итоговый (райт,левт)
    std::memcpy(result.data() + 4, left.data(), 4);
    return result;
}
feistel_wrapper::feistel_wrapper(i_key_schedule* ks, i_feistel_round* fr, int rounds)
    : network(ks, fr, rounds) {//констр обертки передает в констр влож объекта параметры
}
void feistel_wrapper::set_key(const std::vector<uint8_t>& key) {
    current_key = key;//сохр переданный ключ в поле
}
std::vector<uint8_t> feistel_wrapper::encrypt_block(const std::vector<uint8_t>& block) {
    std::vector<uint8_t> b = block;
    b = pre_processing(b);//переопр в наследниках
    b = network.process_block(b, true, current_key);
    b = post_processing(b);
    return b;
}
std::vector<uint8_t> feistel_wrapper::decrypt_block(const std::vector<uint8_t>& block) {
    std::vector<uint8_t> b = block;
    b = pre_processing(b);
    b = network.process_block(b, false, current_key);
    b = post_processing(b);
    return b;
}
std::vector<uint8_t> feistel_wrapper::pre_processing(const std::vector<uint8_t>& block) {
    return block;
}
std::vector<uint8_t> feistel_wrapper::post_processing(const std::vector<uint8_t>& block) {
    return block;
}