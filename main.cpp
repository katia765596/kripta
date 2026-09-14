#include "polynomial.h"
#include "rijndael.h"
#include "twofish.h"
#include "mars.h"
#include "crypto_modes.h"
#include "file_processor.h"
#include "primitive_roots.h"
#include "diffie_hellman.h"
#include "des.h"
#include "lan_key_exchange.h"
#include <future>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
void run_tests();
byte_array string_to_bytes(const std::string& value)
{
    return byte_array(value.begin(), value.end());
}
void print_bytes(const byte_array& data)
{
    for (size_t i = 0; i < data.size(); ++i)
    {
        std::cout
            << std::hex
            << std::setw(2)
            << std::setfill('0')
            << static_cast<int>(data[i])
            << " ";
    }
    std::cout << std::dec << std::endl;
}

void demonstrate_lab_11()
{
    std::cout << "primitive roots" << std::endl;
    std::vector<uint64_t> roots = primitive_roots::get_all(7);
    for (uint64_t value : roots)
        std::cout << value << " ";
    std::cout << std::endl;
    diffie_hellman first(467, 2, 127);
    diffie_hellman second(467, 2, 211);
    uint64_t first_shared = first.make_shared_key(second.get_public_key());//вычисляет общий секрет для каждого участника испол открытые кл.чи друг друга
    uint64_t second_shared = second.make_shared_key(first.get_public_key());
    std::cout << "diffie hellman" << std::endl;
    std::cout << "public a: " << first.get_public_key() << std::endl;
    std::cout << "public b: " << second.get_public_key() << std::endl;
    std::cout << "shared a: " << first_shared << std::endl;
    std::cout << "shared b: " << second_shared << std::endl;
    byte_array des_key = first.make_key(second.get_public_key(), 8);
    byte_array aes_key = first.make_key(second.get_public_key(), 16);
    byte_array mars_key = first.make_key(second.get_public_key(), 16);
    des des_algorithm(des_key);//создает объекты шифров с полученными ключами
    rijndael aes_algorithm(16, aes_key);
    mars mars_algorithm(mars_key);
    byte_array des_source(8, 0x11);//тест блоки для шифрования
    byte_array aes_source(16, 0x22);
    byte_array mars_source(16, 0x33);
    if (des_algorithm.decrypt_block(des_algorithm.encrypt_block(des_source)) != des_source)
        throw std::runtime_error("DES DH demonstration failed");
    if (aes_algorithm.decrypt_block(aes_algorithm.encrypt_block(aes_source)) != aes_source)
        throw std::runtime_error("AES DH demonstration failed");
    if (mars_algorithm.decrypt_block(mars_algorithm.encrypt_block(mars_source)) != mars_source)
        throw std::runtime_error("MARS DH demonstration failed");
    std::cout << "DES key distributed" << std::endl;
    std::cout << "AES key distributed" << std::endl;
    std::cout << "MARS key distributed" << std::endl;
}
void demonstrate_polynomial()
{
    polynomial modulus(0x11b);//x^8 + x^4 + x^3 + x + 1
    polynomial first(0x57);
    polynomial second(0x83);
    finite_field field(modulus);
    std::cout << "polynomial demonstration" << std::endl;
    std::cout << "addition: "
        << std::hex
        << field.add(first, second).to_uint64()
        << std::endl;
    std::cout << "multiplication: "
        << std::hex
        << field.multiply(first, second).to_uint64()
        << std::endl;
    std::cout << "inverse: "
        << std::hex
        << field.inverse(first).to_uint64()
        << std::endl;
    std::cout << "irreducible: "
        << std::boolalpha
        << modulus.is_irreducible()
        << std::endl;
}
void demonstrate_s_box()
{
    rijndael algorithm(16, byte_array(16, 0));
    const byte_array& s_box = algorithm.get_s_box();
    const byte_array& inv_s_box = algorithm.get_inv_s_box();
    std::cout << "s box" << std::endl;
    for (size_t i = 0; i < 16; ++i)
    {
        std::cout
            << std::hex
            << std::setw(2)
            << std::setfill('0')
            << static_cast<int>(s_box[i])
            << " ";
    }
    std::cout << std::endl;
    std::cout << "inverse s box" << std::endl;
    for (size_t i = 0; i < 16; ++i)
    {
        std::cout
            << std::hex
            << std::setw(2)
            << std::setfill('0')
            << static_cast<int>(inv_s_box[i])
            << " ";
    }
    std::cout << std::dec << std::endl;
}
void demonstrate_memory_encryption()
{
    byte_array key = string_to_bytes("1234567890123456");
    byte_array initialization_vector = string_to_bytes("abcdefghijklmnop");
    byte_array source = string_to_bytes("hello rijndael encryption demonstration");
    rijndael algorithm(16, key);
    crypto_modes cryptor(
        algorithm,
        crypto_mode::cbc,
        padding_type::pkcs7,
        initialization_vector
    );
    byte_array encrypted = cryptor.encrypt(source);
    byte_array decrypted = cryptor.decrypt(encrypted);
    std::cout << "source:" << std::endl;
    print_bytes(source);
    std::cout << "encrypted:" << std::endl;
    print_bytes(encrypted);
    std::cout << "decrypted:" << std::endl;
    print_bytes(decrypted);
}
void demonstrate_new_algorithms()
{
    byte_array key = string_to_bytes("1234567890123456");
    byte_array initialization_vector = string_to_bytes("abcdefghijklmnop");
    byte_array source = string_to_bytes("hello twofish mars encryption demonstration");
    twofish twofish_algorithm(key);
    mars mars_algorithm(key);
    crypto_modes twofish_cryptor(
        twofish_algorithm,
        crypto_mode::cbc,
        padding_type::pkcs7,
        initialization_vector
    );
    crypto_modes mars_cryptor(
        mars_algorithm,
        crypto_mode::cbc,
        padding_type::pkcs7,
        initialization_vector
    );
    byte_array twofish_encrypted = twofish_cryptor.encrypt(source);
    byte_array mars_encrypted = mars_cryptor.encrypt(source);
    std::cout << "twofish decrypted:" << std::endl;
    print_bytes(twofish_cryptor.decrypt(twofish_encrypted));
    std::cout << "mars decrypted:" << std::endl;
    print_bytes(mars_cryptor.decrypt(mars_encrypted));
}
void demonstrate_new_files()
{
    byte_array key = string_to_bytes("1234567890123456");
    byte_array initialization_vector = string_to_bytes("abcdefghijklmnop");
    twofish twofish_algorithm(key);
    mars mars_algorithm(key);
    crypto_modes twofish_cryptor(
        twofish_algorithm,
        crypto_mode::ctr,
        padding_type::pkcs7,
        initialization_vector
    );
    crypto_modes mars_cryptor(
        mars_algorithm,
        crypto_mode::ctr,
        padding_type::pkcs7,
        initialization_vector
    );
    file_processor twofish_processor(twofish_cryptor, 4);
    file_processor mars_processor(mars_cryptor, 4);
    if (!twofish_processor.encrypt_file_async("input.txt", "twofish_encrypted.bin").get() ||
        !twofish_processor.decrypt_file_async("twofish_encrypted.bin", "twofish_decrypted.txt").get() ||
        !mars_processor.encrypt_file_async("input.txt", "mars_encrypted.bin").get() ||
        !mars_processor.decrypt_file_async("mars_encrypted.bin", "mars_decrypted.txt").get())
    {
        throw std::runtime_error("new algorithm file processing failed");
    }
    std::cout << "twofish file encrypted and decrypted" << std::endl;
    std::cout << "mars file encrypted and decrypted" << std::endl;
}
void demonstrate_file_encryption()
{
    byte_array key = string_to_bytes("1234567890123456");
    byte_array initialization_vector = string_to_bytes("abcdefghijklmnop");
    rijndael algorithm(16, key);
    crypto_modes cryptor(
        algorithm,
        crypto_mode::ctr,
        padding_type::pkcs7,
        initialization_vector
    );
    file_processor processor(cryptor, 4);
    std::future<bool> encrypt_result =
        processor.encrypt_file_async(
            "input.txt",
            "encrypted.bin"
        );
    bool encrypted = encrypt_result.get();
    if (!encrypted)
    {
        std::cout << "file encryption error" << std::endl;
        return;
    }
    std::cout << "file encrypted" << std::endl;
    std::future<bool> decrypt_result =
        processor.decrypt_file_async(
            "encrypted.bin",
            "decrypted.txt"
        );
    bool decrypted = decrypt_result.get();
    if (decrypted)
    {
        std::cout << "file decrypted" << std::endl;
    }
    else
    {
        std::cout << "file decryption error" << std::endl;
    }
}
int main(int argc, char* argv[])
{
    try
    {
        if (argc >= 2 && std::string(argv[1]) == "dh_server")
        {
            if (argc != 3)
                throw std::invalid_argument("usage: kript_8.exe dh_server port");
            diffie_hellman participant(467, 2, 127);
            byte_array key = lan_key_exchange::run_server(participant, static_cast<uint16_t>(std::stoul(argv[2])), 16);
            std::cout << "received shared DES key" << std::endl;
            print_bytes(byte_array(key.begin(), key.begin() + 8));
            std::cout << "received shared AES key" << std::endl;
            print_bytes(key);
            std::cout << "received shared MARS key" << std::endl;
            print_bytes(key);
            return 0;
        }
        if (argc >= 2 && std::string(argv[1]) == "dh_client")
        {
            if (argc != 4)
                throw std::invalid_argument("usage: kript_8.exe dh_client address port");
            diffie_hellman participant(467, 2, 211);
            byte_array key = lan_key_exchange::run_client(participant, argv[2], static_cast<uint16_t>(std::stoul(argv[3])), 16);
            std::cout << "received shared AES key" << std::endl;
            print_bytes(key);
            return 0;
        }
        run_tests();
        demonstrate_polynomial();
        demonstrate_s_box();
        demonstrate_memory_encryption();
        demonstrate_new_algorithms();
        demonstrate_file_encryption();
        demonstrate_new_files();
        demonstrate_lab_11();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "error: " << error.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}
