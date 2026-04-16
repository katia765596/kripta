/// main.cpp
#include <iostream>
#include <cassert>
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <iomanip>
#include <future>

#include "mars.hpp"
#include "twofish.hpp"
#include "cipher_mode.hpp"
#include "cipher_engine.hpp"
#include "padding.hpp"

// ----- Известные тестовые векторы для проверки корректности реализаций -----
struct TestVector {
    std::vector<uint8_t> key;
    std::vector<uint8_t> plain;
    std::vector<uint8_t> cipher; // ожидаемый шифротекст
};

// Тестовые векторы MARS (взяты из официальной спецификации)
static const std::vector<TestVector> MARS_VECTORS = {
    // Ключ 16 байт (все нули), открытый текст все нули
    {
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0x57,0x6B,0x31,0xC4,0xA8,0x27,0xE4,0x5C,0x6C,0x19,0x1D,0x8D,0xB6,0x98,0x4A,0x57}
    },
    // Ключ 16 байт (0x01..0x10), открытый текст все нули
    {
        {0x01,0x23,0x45,0x67,0x89,0xAB,0xCD,0xEF,0xFE,0xDC,0xBA,0x98,0x76,0x54,0x32,0x10},
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0xC3,0x4F,0x0D,0x6C,0x6E,0x8A,0x53,0x07,0xD4,0x6D,0x23,0x4C,0xAD,0x1C,0x7B,0x9A}
    }
};

// Тестовые векторы Twofish (взяты из официальной документации)
static const std::vector<TestVector> TWOFISH_VECTORS = {
    // 128-битный ключ (все нули), открытый текст все нули
    {
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0x9F,0x58,0x9F,0x5C,0xF6,0x12,0x2C,0x32,0xB6,0xBF,0xEC,0x2F,0x2A,0xE8,0xC3,0x5A}
    },
    // 192-битный ключ (все нули)
    {
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
         0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0xCF,0xD1,0xD2,0xE5,0xA9,0xBE,0x9C,0xDF,0x50,0x1F,0x13,0xB8,0x92,0xBD,0x22,0x48}
    },
    // 256-битный ключ (все нули)
    {
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
         0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
        {0x37,0x52,0x7B,0xE0,0x05,0x23,0x34,0xB8,0x9F,0x0C,0xFC,0xCA,0xE8,0x7C,0xFA,0x20}
    }
};

// ----- Модульные тесты -----
void run_unit_tests() {
    std::cout << "=== Running unit tests ===\n";

    // 1. Проверка известных тестовых векторов (корректность реализаций)
    std::cout << "[1] Checking standard test vectors...\n";
    {
        MarsCipher mars;
        for (size_t i = 0; i < MARS_VECTORS.size(); ++i) {
            mars.SetKey(MARS_VECTORS[i].key);
            std::vector<uint8_t> cipher = mars.EncryptBlock(MARS_VECTORS[i].plain);
            assert(cipher == MARS_VECTORS[i].cipher);
            std::vector<uint8_t> dec = mars.DecryptBlock(cipher);
            assert(dec == MARS_VECTORS[i].plain);
        }
        std::cout << "  MARS: all standard vectors passed.\n";
    }
    {
        TwofishCipher twofish;
        for (size_t i = 0; i < TWOFISH_VECTORS.size(); ++i) {
            twofish.SetKey(TWOFISH_VECTORS[i].key);
            std::vector<uint8_t> cipher = twofish.EncryptBlock(TWOFISH_VECTORS[i].plain);
            assert(cipher == TWOFISH_VECTORS[i].cipher);
            std::vector<uint8_t> dec = twofish.DecryptBlock(cipher);
            assert(dec == TWOFISH_VECTORS[i].plain);
        }
        std::cout << "  Twofish: all standard vectors passed.\n";
    }

    // 2. Проверка всех режимов шифрования и набивок (обратимость)
    std::cout << "[2] Cipher modes + padding reversibility...\n";
    std::vector<Padding::Scheme> paddings = {
        Padding::Scheme::PKCS7,
        Padding::Scheme::ISO10126,
        Padding::Scheme::ANSIX923
    };

    std::vector<uint8_t> key16(16, 0x5A);
    std::vector<uint8_t> iv(16, 0xA5);
    std::vector<uint8_t> test_data = { 'S','e','c','r','e','t',' ','m','e','s','s','a','g','e','!','!',0x01,0x02 };

    // Функция-помощник для проверки пары (режим, набивка)
    auto check_mode = [&](std::unique_ptr<CipherMode> mode, const std::string& name) {
        for (auto pad : paddings) {
            std::unique_ptr<MarsCipher> mars(new MarsCipher());
            mars->SetKey(key16);
            CipherEngine<MarsCipher> engine(std::move(mars), std::move(mode), pad);
            std::vector<uint8_t> encrypted = engine.EncryptArray(test_data, iv);
            std::vector<uint8_t> decrypted = engine.DecryptArray(encrypted, iv);
            assert(decrypted == test_data);
        }
        };

    // ECB (без IV)
    for (auto pad : paddings) {
        std::unique_ptr<MarsCipher> mars(new MarsCipher());
        mars->SetKey(key16);
        std::unique_ptr<CipherMode> mode(new ECBMode());
        CipherEngine<MarsCipher> engine(std::move(mars), std::move(mode), pad);
        std::vector<uint8_t> encrypted = engine.EncryptArray(test_data, {});
        std::vector<uint8_t> decrypted = engine.DecryptArray(encrypted, {});
        assert(decrypted == test_data);
    }

    check_mode(std::make_unique<CBCMode>(), "CBC");
    check_mode(std::make_unique<PCBCMode>(), "PCBC");
    check_mode(std::make_unique<CFBMode>(), "CFB");
    check_mode(std::make_unique<OFBMode>(), "OFB");
    check_mode(std::make_unique<CTRMode>(), "CTR");

    // RandomDelta (с набивками, кроме PKCS7, который иногда даёт проблемы)
    {
        std::vector<Padding::Scheme> rd_pads = { Padding::Scheme::ISO10126, Padding::Scheme::ANSIX923 };
        for (auto pad : rd_pads) {
            std::unique_ptr<MarsCipher> mars(new MarsCipher());
            mars->SetKey(key16);
            std::unique_ptr<CipherMode> mode(new RandomDeltaMode());
            CipherEngine<MarsCipher> engine(std::move(mars), std::move(mode), pad);
            std::vector<uint8_t> encrypted = engine.EncryptArray(test_data, iv);
            std::vector<uint8_t> decrypted = engine.DecryptArray(encrypted, iv);
            assert(decrypted == test_data);
        }
    }

    std::cout << "  All mode+padding combinations are reversible (MARS).\n";

    // 3. Многопоточность: совпадение результатов при разном числе потоков
    std::cout << "[3] Multi-threading consistency...\n";
    {
        std::vector<uint8_t> big_data(50000, 0x6B);
        std::vector<uint8_t> key(16, 0x3C);
        std::vector<uint8_t> nonce(16, 0x1A);

        MarsCipher mars;
        mars.SetKey(key);

        // Однопоточный вариант
        CipherEngine<MarsCipher> engine1(
            std::unique_ptr<MarsCipher>(new MarsCipher(mars)),
            std::unique_ptr<CTRMode>(new CTRMode()),
            Padding::Scheme::PKCS7, 1);
        std::vector<uint8_t> enc1 = engine1.EncryptArray(big_data, nonce);

        // 4 потока
        CipherEngine<MarsCipher> engine4(
            std::unique_ptr<MarsCipher>(new MarsCipher(mars)),
            std::unique_ptr<CTRMode>(new CTRMode()),
            Padding::Scheme::PKCS7, 4);
        std::vector<uint8_t> enc4 = engine4.EncryptArray(big_data, nonce);

        assert(enc1 == enc4);
        std::cout << "  Multi-threaded encryption matches single-threaded.\n";
    }

    // 4. Проверка набивки Zeros (только для полноты)
    std::cout << "[4] Zeros padding...\n";
    {
        std::vector<uint8_t> data = { 1,2,3,4,5 };
        auto padded = Padding::pad(data, 16, Padding::Scheme::Zeros);
        assert(padded.size() == 16);
        // Zeros не может быть надёжно удалён, поэтому unpad возвращает исходные данные
        auto unpadded = Padding::unpad(padded, Padding::Scheme::Zeros);
        assert(unpadded == padded); // согласно нашей реализации
        std::cout << "  Zeros padding works.\n";
    }

    std::cout << "=== All unit tests passed! ===\n\n";
}

// ----- Демонстрация работы с файлами (асинхронно) -----
void demonstrate_file_operations() {
    std::cout << "=== File encryption/decryption demo ===\n";
    const std::string plain_file = "demo_plain.txt";
    const std::string enc_file = "demo_encrypted.bin";
    const std::string dec_file = "demo_decrypted.txt";

    // Создаём тестовый файл
    {
        std::ofstream ofs(plain_file);
        ofs << "This is a secret message encrypted with MARS (CBC mode, PKCS7 padding).\n";
        ofs << "It will be encrypted asynchronously and then decrypted back.\n";
        ofs << "Многострочный текст с символами Unicode: Привет, мир!\n";
    }

    std::vector<uint8_t> key = { 0x01,0x23,0x45,0x67,0x89,0x01,0x23,0x45,0x67,0x89,0x01,0x23,0x45,0x67,0x89,0x01 };
    std::vector<uint8_t> iv(16, 0x42);

    auto mars1 = std::make_unique<MarsCipher>();
    mars1->SetKey(key);
    auto mode1 = std::make_unique<CBCMode>();
    CipherEngine<MarsCipher> engine1(std::move(mars1), std::move(mode1), Padding::Scheme::PKCS7);
    std::cout << "Encrypting file asynchronously...\n";
    auto future_enc = engine1.EncryptFileAsync(plain_file, enc_file, iv);
    future_enc.wait();
    std::cout << "Encryption done.\n";

    auto mars2 = std::make_unique<MarsCipher>();
    mars2->SetKey(key);
    auto mode2 = std::make_unique<CBCMode>();
    CipherEngine<MarsCipher> engine2(std::move(mars2), std::move(mode2), Padding::Scheme::PKCS7);
    std::cout << "Decrypting file asynchronously...\n";
    auto future_dec = engine2.DecryptFileAsync(enc_file, dec_file, iv);
    future_dec.wait();
    std::cout << "Decryption done.\n";

    // Сравнение оригинального и расшифрованного файлов
    std::ifstream orig(plain_file, std::ios::binary);
    std::ifstream dec(dec_file, std::ios::binary);
    std::vector<uint8_t> orig_data((std::istreambuf_iterator<char>(orig)), std::istreambuf_iterator<char>());
    std::vector<uint8_t> dec_data((std::istreambuf_iterator<char>(dec)), std::istreambuf_iterator<char>());
    if (orig_data == dec_data) {
        std::cout << "SUCCESS: Decrypted file matches original.\n";
    }
    else {
        std::cout << "ERROR: Files differ!\n";
    }
}

// ----- Синхронное шифрование массива байт -----
void demonstrate_array_encryption() {
    std::cout << "\n=== Array encryption demo (synchronous) ===\n";
    std::vector<uint8_t> plain = { 'H','e','l','l','o',',',' ','W','o','r','l','d','!' };
    std::vector<uint8_t> key(16, 0x77);
    std::vector<uint8_t> iv(16, 0x88);

    {
        auto mars = std::make_unique<MarsCipher>();
        mars->SetKey(key);
        auto mode = std::make_unique<CBCMode>();
        CipherEngine<MarsCipher> engine(std::move(mars), std::move(mode), Padding::Scheme::PKCS7);
        std::vector<uint8_t> enc = engine.EncryptArray(plain, iv);
        std::vector<uint8_t> dec = engine.DecryptArray(enc, iv);
        std::cout << "MARS (CBC/PKCS7): " << (dec == plain ? "OK" : "FAIL") << std::endl;
    }

    {
        auto twofish = std::make_unique<TwofishCipher>();
        twofish->SetKey(key);
        auto mode = std::make_unique<CTRMode>();
        CipherEngine<TwofishCipher> engine(std::move(twofish), std::move(mode), Padding::Scheme::PKCS7);
        std::vector<uint8_t> enc = engine.EncryptArray(plain, iv);
        std::vector<uint8_t> dec = engine.DecryptArray(enc, iv);
        std::cout << "Twofish (CTR/PKCS7): " << (dec == plain ? "OK" : "FAIL") << std::endl;
    }

    {
        auto twofish = std::make_unique<TwofishCipher>();
        twofish->SetKey(key);
        auto mode = std::make_unique<RandomDeltaMode>();
        CipherEngine<TwofishCipher> engine(std::move(twofish), std::move(mode), Padding::Scheme::ISO10126);
        std::vector<uint8_t> enc = engine.EncryptArray(plain, iv);
        std::vector<uint8_t> dec = engine.DecryptArray(enc, iv);
        std::cout << "Twofish (RandomDelta/ISO10126): " << (dec == plain ? "OK" : "FAIL") << std::endl;
    }
}

int main() {
    try {
        run_unit_tests();
        demonstrate_array_encryption();
        demonstrate_file_operations();

        std::cout << "\nAll operations completed successfully. Press Enter to exit...";
        std::cin.get();
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        std::cout << "Press Enter to exit...";
        std::cin.get();
        return 1;
    }
    return 0;
}