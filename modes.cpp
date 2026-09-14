#include "modes.h"
#include <cstring>
#include <stdexcept>
void ecb_mode::init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) {
    this->iv = iv;//вис чтобы не было конфликта имен параметр и поле одинаково наз-ся
}
std::vector<uint8_t> ecb_mode::process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) {
    if (encrypt) return cipher->encrypt_block(block);
    else return cipher->decrypt_block(block);
}
bool ecb_mode::can_parallel_encrypt() const { return true; }
bool ecb_mode::can_parallel_decrypt() const { return true; }
bool ecb_mode::needs_padding() const { return true; }
void ecb_mode::reset() { iv.clear(); }
void cbc_mode::init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) {
    prev_block = iv;
    initialized = true;
}
std::vector<uint8_t> cbc_mode::process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) {
    if (!initialized) throw std::runtime_error("cbc not initialized");
    size_t bs = cipher->block_size();
    std::vector<uint8_t> result(bs);
    if (encrypt) {
        for (size_t i = 0; i < bs; ++i) result[i] = block[i] ^ prev_block[i];//хорим тек блок с пред шифротекстом
        result = cipher->encrypt_block(result);
        prev_block = result;//зашифр блок в прев_блок для след шага
    }
    else {
        std::vector<uint8_t> dec = cipher->decrypt_block(block);
        for (size_t i = 0; i < bs; ++i) result[i] = dec[i] ^ prev_block[i];//хорим дек с пред шифр
        prev_block = block;//сохраняем тек блок шифротекста как прев_блок для след шага
    }
    return result;
}
bool cbc_mode::can_parallel_encrypt() const { return false; }
bool cbc_mode::can_parallel_decrypt() const { return true; }//можно тк каждый блок расшифр нез испол предыдущий шифр уже известный
bool cbc_mode::needs_padding() const { return true; }
void cbc_mode::reset() { prev_block.clear(); initialized = false; }
void pcbc_mode::init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) {
    prev_cipher = iv;
    prev_plain = iv;
    initialized = true;
}
std::vector<uint8_t> pcbc_mode::process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) {
    if (!initialized) throw std::runtime_error("pcbc not initialized");
    size_t bs = cipher->block_size();
    std::vector<uint8_t> result(bs);
    if (encrypt) {
        for (size_t i = 0; i < bs; ++i) result[i] = block[i] ^ prev_cipher[i] ^ prev_plain[i];
        result = cipher->encrypt_block(result);//испол тек блок прнд шифрот пред открыт текст
        prev_cipher = result;
        prev_plain = block;//пред откр текст исх и т д
    }
    else {
        std::vector<uint8_t> dec = cipher->decrypt_block(block);
        for (size_t i = 0; i < bs; ++i) result[i] = dec[i] ^ prev_cipher[i] ^ prev_plain[i];
        prev_cipher = block;
        prev_plain = result;
    }
    return result;
}
bool pcbc_mode::can_parallel_encrypt() const { return false; }
bool pcbc_mode::can_parallel_decrypt() const { return false; }
bool pcbc_mode::needs_padding() const { return true; }
void pcbc_mode::reset() { prev_cipher.clear(); prev_plain.clear(); initialized = false; }
void cfb_mode::init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) {
    shift_reg = iv;//потоковый вместо шифрования блоков откр текста, шифр генерирует гамму (ключ поток) на основе пред блока шифро затем гамма хорится с откр текстом
}//гамма зав от шифротекста и ошибка в одном блоке шифрокеста влияет на все
std::vector<uint8_t> cfb_mode::process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) {
    size_t bs = cipher->block_size();
    std::vector<uint8_t> enc = cipher->encrypt_block(shift_reg);
    std::vector<uint8_t> result(bs);
    for (size_t i = 0; i < bs; ++i) result[i] = block[i] ^ enc[i];
    shift_reg = result;
    return result;
}
bool cfb_mode::can_parallel_encrypt() const { return false; }
bool cfb_mode::can_parallel_decrypt() const { return false; }
bool cfb_mode::needs_padding() const { return false; }
void cfb_mode::reset() { shift_reg.clear(); }//шифр симметр одна и та же опер
void ofb_mode::init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) {
    state = iv;//хранит тек знач, которое шифруется на каждом шаге для генерации гаммы, уник для каждого соо
}//на каждом шаге подается предыд вых блок шифра(гамма) для первого iv,шифр дает новую гамму которая хорится с откр текстом,ошибка в шифротексте портит лишь один блок.,гамма не зав от шифротекста
std::vector<uint8_t> ofb_mode::process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) {
    size_t bs = cipher->block_size();
    state = cipher->encrypt_block(state);
    std::vector<uint8_t> result(bs);
    for (size_t i = 0; i < bs; ++i) result[i] = block[i] ^ state[i];
    return result;
}
bool ofb_mode::can_parallel_encrypt() const { return false; }
bool ofb_mode::can_parallel_decrypt() const { return false; }
bool ofb_mode::needs_padding() const { return false; }
void ofb_mode::reset() { state.clear(); }//опер идентичная
void ctr_mode::init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) {
    if (iv.size() != 8) throw std::invalid_argument("ctr iv must be 8 bytes");
    counter = iv;
    counter_value = 0;//счетчик
    for (size_t i = 0; i < 8; ++i) counter_value = (counter_value << 8) | iv[i];
}//цикл по 8 байтам преобраз байты iv в 64 битное целое это будет начальное знач счётчика для первого броска
std::vector<uint8_t> ctr_mode::process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) {
    size_t bs = cipher->block_size();
    std::vector<uint8_t> counter_bytes(bs);//вектор для байтового представления тек счетчика
    uint64_t val = counter_value;
    for (int i = bs - 1; i >= 0; --i) {//заполняем от млад байта до старшего порядок об большего
        counter_bytes[i] = val & 0xFF;
        val >>= 8;
    }
    std::vector<uint8_t> keystream = cipher->encrypt_block(counter_bytes);//шифруем байты счетчика получаем гамму
    std::vector<uint8_t> result(bs);
    for (size_t i = 0; i < bs; ++i) result[i] = block[i] ^ keystream[i];
    counter_value++;
    return result;//счетчика известны зараннее и незав, можно параллелить, шифр и дешифр симметричные 
}
bool ctr_mode::can_parallel_encrypt() const { return true; }
bool ctr_mode::can_parallel_decrypt() const { return true; }
bool ctr_mode::needs_padding() const { return false; }
void ctr_mode::reset() { counter_value = 0; }
void random_delta_mode::init(const std::vector<uint8_t>& iv, const std::vector<uint8_t>& params) {
    if (iv.size() != 8) throw std::invalid_argument("rd iv must be 8 bytes (initial)");
    if (params.size() != 8) throw std::invalid_argument("rd params must be 8 bytes (delta)");
    initial = iv;
    delta = params;
    current_counter = 0;
    for (size_t i = 0; i < 8; ++i) current_counter = (current_counter << 8) | iv[i];
    initialized = true;
}
std::vector<uint8_t> random_delta_mode::process_block(const std::vector<uint8_t>& block, bool encrypt, i_symmetric_cipher* cipher) {
    if (!initialized) throw std::runtime_error("random delta not initialized");
    size_t bs = cipher->block_size();
    std::vector<uint8_t> counter_bytes(bs);
    uint64_t val = current_counter;
    for (int i = bs - 1; i >= 0; --i) {
        counter_bytes[i] = val & 0xFF;
        val >>= 8;
    }
    std::vector<uint8_t> keystream = cipher->encrypt_block(counter_bytes);
    std::vector<uint8_t> result(bs);
    for (size_t i = 0; i < bs; ++i) result[i] = block[i] ^ keystream[i];
    uint64_t delta_val = 0;
    for (size_t i = 0; i < 8; ++i) delta_val = (delta_val << 8) | delta[i];
    current_counter += delta_val;//счетчик на дельту
    return result;
}
bool random_delta_mode::can_parallel_encrypt() const { return true; }
bool random_delta_mode::can_parallel_decrypt() const { return true; }
bool random_delta_mode::needs_padding() const { return false; }
void random_delta_mode::reset() { current_counter = 0; initialized = false; }