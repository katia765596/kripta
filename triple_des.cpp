#include "triple_des.h"//пункт 4 реал 3дес с поддержкой интерфейса 1с (класс обязан следовать интерфейсу (напрямую реализует симетри шифр но внутри использует объекты которые тоже реализуют этот интерфейс_
#include <stdexcept>
triple_des::triple_des(des_mode mode) : mode(mode) {}
void triple_des::set_key(const std::vector<uint8_t>& key) {//реализует метод сеткей из интерф 1с (разобрать переданный ключ на составляющие)
    if (mode == des_mode::eee3 || mode == des_mode::ede3) {
        if (key.size() != 24) throw std::invalid_argument("3des needs 24 bytes key");
        std::vector<uint8_t> k1(key.begin(), key.begin() + 8);
        std::vector<uint8_t> k2(key.begin() + 8, key.begin() + 16);
        std::vector<uint8_t> k3(key.begin() + 16, key.begin() + 24);//конструктур вектора извлекает 8 байтовых ключа для каждого внутр дес
        des1.set_key(k1);//каждый объект дес инициал своим ключом, внутри дес вызывается генерация раундовых ключей
        des2.set_key(k2);
        des3.set_key(k3);
    }
    else {
        if (key.size() != 16) throw std::invalid_argument("3des 2-key needs 16 bytes");
        std::vector<uint8_t> k1(key.begin(), key.begin() + 8);
        std::vector<uint8_t> k2(key.begin() + 8, key.begin() + 16);
        des1.set_key(k1);
        des2.set_key(k2);
        des3.set_key(k1);
    }
}
std::vector<uint8_t> triple_des::encrypt_block(const std::vector<uint8_t>& block) {//переор вирт метод интерф
    if (mode == des_mode::eee3 || mode == des_mode::eee2) {
        return des3.encrypt_block(des2.encrypt_block(des1.encrypt_block(block)));
    }//еее3 еее2
    else {
        return des3.encrypt_block(des2.decrypt_block(des1.encrypt_block(block)));
    }//ede3 ede2
}
std::vector<uint8_t> triple_des::decrypt_block(const std::vector<uint8_t>& block) {
    if (mode == des_mode::eee3 || mode == des_mode::eee2) {
        return des1.decrypt_block(des2.decrypt_block(des3.decrypt_block(block)));
    }
    else {
        return des1.decrypt_block(des2.encrypt_block(des3.decrypt_block(block)));
    }
}