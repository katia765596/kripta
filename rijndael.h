#ifndef RIJNDAEL_H
#define RIJNDAEL_H
#include "byte_array.h"
#include "block_cipher.h"//содержит объявление абстрактного баз класса блок_шифер, класс рейнджел наследует интерфейс и реализует его вирт методы
#include <cstddef>
#include <cstdint>
class rijndael : public block_cipher//рейнджел кокретная реализация блочного шифра и обязан предоставить реализации всех чисто вирт методов базового класса 
{
public:
    rijndael(size_t block_size, const byte_array& key);
    byte_array encrypt_block(const byte_array& block) const override;
    byte_array decrypt_block(const byte_array& block) const override;
    const byte_array& get_s_box() const;//геттер для таблицы замен s-box (256 байт)
    const byte_array& get_inv_s_box() const;
    size_t get_block_size() const override;
    size_t get_key_size() const;
    size_t get_round_count() const;//получение кол-ва раундов шифрования
private:
    size_t block_size;//стандарт для AES 16(128БИТ),24(192 бит),32(256 бит)
    size_t key_size;
    size_t nb;//кол-во столбцов в матрице состояния (размер блока на 4)
    size_t nk;//кол-во столбцов ключа на 4 делим
    size_t nr;//кол-во раундов (комбинация nb и nk)
    byte_array key;//используется при генерации раунд ключей
    byte_array round_keys;
    byte_array s_box;
    byte_array inv_s_box;
    static uint8_t multiply(uint8_t first, uint8_t second);//умн двуз байтов в поле gf(2^8) с исполь неприводимого полинома aes
    static uint8_t gf_inverse(uint8_t value);//вычисляет мультипл обратный элемент байта value в поле, используется при построении s-box
    static uint8_t rotl8(uint8_t value, unsigned int shift);//испол в аффинном преобраз для построения s-box
    void generate_s_boxes();
    void generate_round_keys();
    void sub_bytes(byte_array& state) const;//замена каждого байта состояния с помощью s-box принимает ссылку на состояние
    void inv_sub_bytes(byte_array& state) const;
    void shift_rows(byte_array& state) const;//циклический сдвиг строк состояния (представлкно как иассив байтов -матрица 4хnb)
    void inv_shift_rows(byte_array& state) const;
    void mix_columns(byte_array& state) const;//умн каждого столбца сост на фикс матрица в поле
    void inv_mix_columns(byte_array& state) const;
    void add_round_key(byte_array& state, size_t round) const;//XOR состояния с раунд ключом для заданного номера раунд
    uint8_t get_round_constant(size_t index) const;//возвр раунд константу для расширения ключа по индексу
};
#endif
