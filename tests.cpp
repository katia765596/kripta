#include "polynomial.h"
#include "padding.h"
#include "rijndael.h"
#include "twofish.h"
#include "mars.h"
#include "crypto_modes.h"
#include "file_processor.h"
#include "twofish.h"
#include "mars.h"
#include "des.h"
#include "primitive_roots.h"
#include "diffie_hellman.h"
#include "lan_key_exchange.h"
#include <cstdint>
#include <cstdio>
#include <future>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <thread>
#include <chrono>
void check(bool value, const std::string& message)
{
    if (!value)
    {
        throw std::runtime_error(message);
    }
}
byte_array from_hex(const std::string& value)
{
    byte_array result;
    for (size_t i = 0; i < value.size(); i += 2)
    {
        result.push_back(static_cast<uint8_t>(std::stoul(value.substr(i, 2), 0, 16)));//извлекает подстроку из двух символов
    }
    return result;
}
byte_array make_data(size_t size)
{
    byte_array result(size);
    for (size_t i = 0; i < size; ++i)
    {
        result[i] = static_cast<uint8_t>((i * 37 + 11) & 0xff);
    }
    return result;
}
void test_polynomial()
{
    polynomial modulus(0x11b);
    polynomial first(0x57);
    polynomial second(0x83);
    check(first.add(second).to_uint64() == 0xd4, "polynomial addition failed");
    check(first.multiply_mod(second, modulus).to_uint64() == 0xc1, "polynomial multiplication failed");
    check(first.inverse(modulus).to_uint64() == 0xbf, "polynomial inverse failed");
    check(first.multiply_mod(first.inverse(modulus), modulus).to_uint64() == 1, "inverse verification failed");
    check(modulus.is_irreducible(), "irreducibility test failed");
    polynomial large = polynomial::from_binary("10000000000000001");
    check(large.degree() == 16, "large polynomial degree failed");
    check(large.to_binary() == "10000000000000001", "binary polynomial failed");
}
void test_padding_type(padding_type type)
{
    byte_array source = { 1, 2, 3, 4, 5, 6, 7 };
    byte_array padded = padding::add(source, 16, type);
    check(!padded.empty(), "padding is empty");
    check(padded.size() % 16 == 0, "padding size failed");
    check(padding::remove(padded, 16, type) == source, "padding restore failed");
}
void test_s_box()
{
    rijndael algorithm(16, byte_array(16, 0));
    const byte_array& s_box = algorithm.get_s_box();
    const byte_array& inv_s_box = algorithm.get_inv_s_box();
    check(s_box.size() == 256, "s box size failed");
    check(inv_s_box.size() == 256, "inverse s box size failed");
    check(s_box[0x00] == 0x63, "s box failed");
    check(s_box[0x53] == 0xed, "s box value failed");
    check(inv_s_box[0x63] == 0x00, "inverse s box failed");
    for (size_t i = 0; i < 256; ++i)
    {
        check(inv_s_box[s_box[i]] == i, "s box inverse relation failed");
    }
}
void test_aes_vector()
{
    byte_array key = from_hex("000102030405060708090a0b0c0d0e0f");
    byte_array source = from_hex("00112233445566778899aabbccddeeff");
    byte_array expected = from_hex("69c4e0d86a7b0430d8cdb78070b4c55a");
    rijndael algorithm(16, key);
    byte_array encrypted = algorithm.encrypt_block(source);
    check(encrypted == expected, "AES-128 test vector failed");
    check(algorithm.decrypt_block(encrypted) == source, "AES-128 decryption failed");
    check(algorithm.get_block_size() == 16, "block size getter failed");
    check(algorithm.get_key_size() == 16, "key size getter failed");
    check(algorithm.get_round_count() == 10, "round count failed");
}
void test_rijndael(size_t block_size, size_t key_size)
{
    byte_array key(key_size);
    byte_array source(block_size);
    for (size_t i = 0; i < key.size(); ++i)
    {
        key[i] = static_cast<uint8_t>(i);
    }
    for (size_t i = 0; i < source.size(); ++i)
    {
        source[i] = static_cast<uint8_t>(i + 1);
    }
    rijndael algorithm(block_size, key);
    check(algorithm.decrypt_block(algorithm.encrypt_block(source)) == source, "Rijndael encryption failed");
}
void test_mode(block_cipher& algorithm, crypto_mode mode, padding_type padding_mode)
{
    byte_array iv(16);
    for (size_t i = 0; i < iv.size(); ++i)
    {
        iv[i] = static_cast<uint8_t>(16 - i);
    }
    byte_array source = make_data(100);
    crypto_modes cryptor(algorithm, mode, padding_mode, iv);
    byte_array encrypted = cryptor.encrypt(source);
    check(cryptor.decrypt(encrypted) == source, "crypto mode failed");
    check(cryptor.get_mode() == mode, "mode getter failed");
    check(cryptor.get_padding_type() == padding_mode, "padding getter failed");
    check(cryptor.get_initialization_vector() == iv, "iv getter failed");
    if (mode != crypto_mode::ecb && mode != crypto_mode::ctr)
    {
        byte_array fallback_encrypted = cryptor.encrypt_parallel(source, 4);
        check(cryptor.decrypt_parallel(fallback_encrypted, 4) == source, "parallel fallback decryption failed");
    }
    else
    {
        byte_array encrypted_parallel = cryptor.encrypt_parallel(source, 4);
        check(cryptor.decrypt_parallel(encrypted_parallel, 4) == source, "parallel decryption failed");
        byte_array single_thread_encrypted = cryptor.encrypt_parallel(source, 1);
        check(cryptor.decrypt_parallel(single_thread_encrypted, 1) == source, "single thread decryption failed");
        if (padding_mode != padding_type::iso_10126)
        {//для паддинга этого рез шифр может различ из-за случ байтов поэтому сравнение с синхр шифр опускается
            check(single_thread_encrypted == encrypted, "single thread encryption failed");
        }
        if (padding_mode != padding_type::iso_10126)
        {
            check(encrypted_parallel == encrypted, "parallel encryption failed");
        }
    }
}
void test_algorithm_modes(block_cipher& algorithm)
{
    crypto_mode modes[] = {
        crypto_mode::ecb,
        crypto_mode::cbc,
        crypto_mode::pcbc,
        crypto_mode::cfb,
        crypto_mode::ofb,
        crypto_mode::ctr,
        crypto_mode::random_delta
    };
    padding_type paddings[] = {
        padding_type::zeros,
        padding_type::pkcs7,
        padding_type::iso_10126,
        padding_type::ansi_x923
    };
    for (crypto_mode mode : modes)
    {
        for (padding_type padding_mode : paddings)
        {
            test_mode(algorithm, mode, padding_mode);
        }
    }
}
void test_twofish_vectors()
{
    struct test_case
    {
        std::string key;
        std::string plain;
        std::string expected;
    };
    test_case cases[] = {//ключи разной длины 128,196,256 бит
        {"00000000000000000000000000000000", "00000000000000000000000000000000", "9f589f5cf6122c32b6bfec2f2ae8c35a"},
        {"0123456789abcdeffedcba98765432100011223344556677", "00000000000000000000000000000000", "cfd1d2e5a9be9cdf501f13b892bd2248"},
        {"0123456789abcdeffedcba987654321000112233445566778899aabbccddeeff", "00000000000000000000000000000000", "37527be0052334b89f0cfccae87cfa20"}
    };
    for (const test_case& item : cases)
    {
        twofish algorithm(from_hex(item.key));
        byte_array plain = from_hex(item.plain);
        byte_array encrypted = algorithm.encrypt_block(plain);
        check(encrypted == from_hex(item.expected), "Twofish test vector failed");
        check(algorithm.decrypt_block(encrypted) == plain, "Twofish decryption failed");
        check(algorithm.get_block_size() == 16, "Twofish block size failed");
        check(algorithm.get_key_size() == item.key.size() / 2, "Twofish key size failed");//каждый байт кодир двумя символами
        check(algorithm.get_round_count() == 16, "Twofish round count failed");
    }
}
void test_mars_vectors()
{
    struct test_case
    {
        std::string key;
        std::string plain;
        std::string expected;
    };
    test_case cases[] = {
        {"00000000000000000000000000000000", "00000000000000000000000000000000", "dcc07b8dfb0738d6e30a22dfcf27e886"},
        {"80000000000000000000000000000000", "00000000000000000000000000000000", "b3e2ad5608ac1b6733a7cb4fdf8f9952"},
        {"000000000000000000000000000000000000000000000000", "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "97778747d60e425c2b4202599db856fb"},
        {"0000000000000000000000000000000000000000000000000000000000000000", "62e45b4cf3477f1dd65063729d9aba8f", "0f4b897ea014d21fbc20f1054a42f719"}
    };
    for (const test_case& item : cases)
    {
        mars algorithm(from_hex(item.key));
        byte_array plain = from_hex(item.plain);
        byte_array encrypted = algorithm.encrypt_block(plain);
        check(encrypted == from_hex(item.expected), "MARS test vector failed");
        check(algorithm.decrypt_block(encrypted) == plain, "MARS decryption failed");
        check(algorithm.get_block_size() == 16, "MARS block size failed");
        check(algorithm.get_key_size() == item.key.size() / 2, "MARS key size failed");
        check(algorithm.get_round_count() == 32, "MARS round count failed");
    }
}
void test_new_algorithms()
{
    byte_array key128(16, 0);
    byte_array key192(24, 0);
    byte_array key256(32, 0);
    twofish twofish_128(key128);
    twofish twofish_192(key192);
    twofish twofish_256(key256);
    mars mars_128(key128);
    mars mars_192(key192);
    mars mars_256(key256);
    block_cipher* algorithms[] = {
        &twofish_128,
        &twofish_192,
        &twofish_256,
        &mars_128,
        &mars_192,
        &mars_256
    };
    for (block_cipher* algorithm : algorithms)
    {
        test_algorithm_modes(*algorithm);
    }
}
void test_file_processor_for_algorithm(block_cipher& algorithm, const std::string& name)
{
    const std::string input = name + "_input.bin";
    const std::string encrypted = name + "_encrypted.bin";
    const std::string decrypted = name + "_decrypted.bin";
    byte_array source = make_data(257);
    {
        std::ofstream file(input, std::ios::binary);
        check(static_cast<bool>(file), "file open failed");
        file.write(reinterpret_cast<const char*>(source.data()), static_cast<std::streamsize>(source.size()));
    }
    crypto_modes cryptor(algorithm, crypto_mode::ctr, padding_type::pkcs7, byte_array(16, 3));
    file_processor processor(cryptor, 4);
    check(processor.get_thread_count() == 4, "thread count getter failed");
    check(processor.encrypt_file(input, encrypted), "sync file encryption failed");
    check(processor.decrypt_file(encrypted, decrypted), "sync file decryption failed");
    std::ifstream file(decrypted, std::ios::binary);
    byte_array restored((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    check(restored == source, "sync file content failed");
    check(processor.encrypt_file_async(input, encrypted).get(), "async file encryption failed");
    check(processor.decrypt_file_async(encrypted, decrypted).get(), "async file decryption failed");
    std::ifstream async_file(decrypted, std::ios::binary);
    byte_array async_restored((std::istreambuf_iterator<char>(async_file)), std::istreambuf_iterator<char>());
    check(async_restored == source, "async file content failed");
    std::remove(input.c_str());
    std::remove(encrypted.c_str());
    std::remove(decrypted.c_str());
}
void test_invalid_thread_count(block_cipher& algorithm)
{
    crypto_modes cryptor(algorithm, crypto_mode::ecb, padding_type::pkcs7);
    bool encryption_failed = false;
    bool decryption_failed = false;
    try
    {
        cryptor.encrypt_parallel(make_data(32), 0);
    }
    catch (const std::invalid_argument&)
    {
        encryption_failed = true;
    }
    try
    {
        cryptor.decrypt_parallel(make_data(32), 0);
    }
    catch (const std::invalid_argument&)
    {
        decryption_failed = true;
    }
    check(encryption_failed, "invalid encryption thread count was accepted");
    check(decryption_failed, "invalid decryption thread count was accepted");
}
void test_file_processor()
{
    byte_array key(16, 7);
    rijndael rijndael_algorithm(16, key);
    twofish twofish_algorithm(key);
    mars mars_algorithm(key);
    test_file_processor_for_algorithm(rijndael_algorithm, "rijndael");
    test_file_processor_for_algorithm(twofish_algorithm, "twofish");
    test_file_processor_for_algorithm(mars_algorithm, "mars");
}

void test_twofish_and_mars()
{
    byte_array key = from_hex("000102030405060708090a0b0c0d0e0f");
    byte_array source = from_hex("00112233445566778899aabbccddeeff");
    twofish twofish_algorithm(key);
    mars mars_algorithm(key);
    check(twofish_algorithm.get_block_size() == 16, "Twofish block size failed");
    check(twofish_algorithm.get_key_size() == 16, "Twofish key size failed");
    check(twofish_algorithm.get_round_count() == 16, "Twofish round count failed");
    check(mars_algorithm.get_block_size() == 16, "MARS block size failed");
    check(mars_algorithm.get_key_size() == 16, "MARS key size failed");
    check(mars_algorithm.get_round_count() == 32, "MARS round count failed");
    check(twofish_algorithm.decrypt_block(twofish_algorithm.encrypt_block(source)) == source, "Twofish failed");
    check(mars_algorithm.decrypt_block(mars_algorithm.encrypt_block(source)) == source, "MARS failed");
}
void test_des()
{
    byte_array key = from_hex("133457799bbcdff1");
    byte_array source = from_hex("0123456789abcdef");
    byte_array expected = from_hex("85e813540f0ab405");
    des algorithm(key);
    byte_array encrypted = algorithm.encrypt_block(source);
    check(encrypted == expected, "DES test vector failed");
    check(algorithm.decrypt_block(encrypted) == source, "DES decryption failed");
    check(algorithm.get_block_size() == 8, "DES block size failed");
    check(algorithm.get_key_size() == 8, "DES key size failed");
}
void test_primitive_roots()
{
    check(primitive_roots::exists(2), "primitive roots n=2 failed");
    check(primitive_roots::exists(4), "primitive roots n=4 failed");
    check(primitive_roots::exists(7), "primitive roots n=7 failed");
    check(primitive_roots::exists(14), "primitive roots n=14 failed");
    check(!primitive_roots::exists(8), "primitive roots n=8 failed");
    check(!primitive_roots::exists(12), "primitive roots n=12 failed");
    std::vector<uint64_t> roots = primitive_roots::get_all(7);
    std::vector<uint64_t> expected = { 3, 5 };
    check(roots == expected, "primitive roots values failed");
    roots = primitive_roots::get_all(9);
    check(roots.size() == 2, "primitive roots count failed");
    check(roots[0] != roots[1], "primitive roots uniqueness failed");
}
void test_diffie_hellman()
{
    diffie_hellman first(23, 5, 6);
    diffie_hellman second(23, 5, 15);
    check(first.get_public_key() == 8, "DH first public key failed");
    check(second.get_public_key() == 19, "DH second public key failed");
    check(first.make_shared_key(second.get_public_key()) == 2, "DH first shared key failed");
    check(second.make_shared_key(first.get_public_key()) == 2, "DH second shared key failed");
    check(first.make_key(second.get_public_key(), 16) == second.make_key(first.get_public_key(), 16), "DH derived key failed");
}
void test_distributed_keys()
{
    diffie_hellman first(467, 2, 127);
    diffie_hellman second(467, 2, 211);
    byte_array des_key_first = first.make_key(second.get_public_key(), 8);
    byte_array des_key_second = second.make_key(first.get_public_key(), 8);
    byte_array aes_key_first = first.make_key(second.get_public_key(), 16);
    byte_array aes_key_second = second.make_key(first.get_public_key(), 16);
    byte_array mars_key_first = first.make_key(second.get_public_key(), 16);
    byte_array mars_key_second = second.make_key(first.get_public_key(), 16);
    check(des_key_first == des_key_second, "DES distributed key failed");
    check(aes_key_first == aes_key_second, "AES distributed key failed");
    check(mars_key_first == mars_key_second, "MARS distributed key failed");
    des des_algorithm(des_key_first);
    rijndael aes_algorithm(16, aes_key_first);
    mars mars_algorithm(mars_key_first);
    byte_array des_source(8, 0x11);
    byte_array aes_source(16, 0x22);
    byte_array mars_source(16, 0x33);
    check(des_algorithm.decrypt_block(des_algorithm.encrypt_block(des_source)) == des_source, "DES distributed encryption failed");
    check(aes_algorithm.decrypt_block(aes_algorithm.encrypt_block(aes_source)) == aes_source, "AES distributed encryption failed");
    check(mars_algorithm.decrypt_block(mars_algorithm.encrypt_block(mars_source)) == mars_source, "MARS distributed encryption failed");
}

void test_lan_key_exchange()
{
    diffie_hellman first(467, 2, 127);
    diffie_hellman second(467, 2, 211);
    const uint16_t port = 5057;
    std::future<byte_array> server = std::async(std::launch::async, [&first, port]()
        {
            return lan_key_exchange::run_server(first, port, 16);
        });
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    byte_array client = lan_key_exchange::run_client(second, "127.0.0.1", port, 16);
    byte_array server_key = server.get();
    check(server_key == client, "LAN key exchange failed");
}
void run_tests()
{
    test_polynomial();
    test_padding_type(padding_type::zeros);
    test_padding_type(padding_type::pkcs7);
    test_padding_type(padding_type::iso_10126);
    test_padding_type(padding_type::ansi_x923);
    test_s_box();
    test_aes_vector();
    test_rijndael(16, 16);
    test_rijndael(16, 24);
    test_rijndael(16, 32);
    test_rijndael(24, 16);
    test_rijndael(24, 24);
    test_rijndael(24, 32);
    test_rijndael(32, 16);
    test_rijndael(32, 24);
    test_rijndael(32, 32);
    rijndael rijndael_algorithm(16, byte_array(16, 1));
    test_algorithm_modes(rijndael_algorithm);
    test_twofish_vectors();
    test_mars_vectors();
    test_new_algorithms();
    twofish thread_test_algorithm(byte_array(16, 0));
    test_invalid_thread_count(thread_test_algorithm);
    test_file_processor();
    test_twofish_and_mars();
    test_des();
    test_primitive_roots();
    test_diffie_hellman();
    test_distributed_keys();
    test_lan_key_exchange();
    std::cout << "all tests passed" << std::endl;
}
