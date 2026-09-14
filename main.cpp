#include "des.h"
#include "triple_des.h"
#include "processor.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
static void test_des() {
    des d;
    std::vector<uint8_t> key = { 0x13,0x34,0x57,0x79,0x9B,0xBC,0xDF,0xF1 };
    d.set_key(key);
    std::vector<uint8_t> plain = { 0x01,0x23,0x45,0x67,0x89,0xAB,0xCD,0xEF };
    std::vector<uint8_t> enc = d.encrypt_block(plain);
    std::vector<uint8_t> dec = d.decrypt_block(enc);
    assert(plain == dec);
}
static void test_triple_des() {
    triple_des t(des_mode::ede3);
    std::vector<uint8_t> key(24, 0x01);//24 байта 
    t.set_key(key);
    std::vector<uint8_t> plain(8, 0xAA);
    std::vector<uint8_t> enc = t.encrypt_block(plain);
    std::vector<uint8_t> dec = t.decrypt_block(enc);
    assert(plain == dec);
}
static void test_processor() {
    des d;
    std::vector<uint8_t> key(8, 0x11);
    d.set_key(key);
    cipher_processor proc(&d);// проц принимает любой алг реализ симм шифр
    std::vector<uint8_t> data = { 1,2,3,4,5,6,7,8,9 };
    std::vector<uint8_t> enc = proc.process_data(data, true);
    std::vector<uint8_t> dec = proc.process_data(enc, false);
    assert(data == dec);
}
static void test_files() {
    des d;
    std::vector<uint8_t> key(8, 0x22);
    d.set_key(key);
    cipher_processor proc(&d);
    std::string in = "test_in.txt";
    std::string enc = "test_enc.bin";
    std::string dec = "test_dec.txt";
    std::remove(in.c_str());//удал старые файлы с этими именами
    std::remove(enc.c_str());
    std::remove(dec.c_str());
    std::ofstream f(in, std::ios::binary);
    f << "hello world";
    f.close();
    proc.process_file(in, enc, true);//шифр дешифр с файлами
    proc.process_file(enc, dec, false);
    assert(compare_files(in, dec));
    std::remove(in.c_str());
    std::remove(enc.c_str());//удал временные файлы
    std::remove(dec.c_str());
}
int main() {
    test_des();
    test_triple_des();
    test_processor();
    test_files();
    std::cout << "All tests passed\n";
    des d1;
    std::vector<uint8_t> key1(8, 0x33);
    d1.set_key(key1);
    triple_des t1(des_mode::ede3);
    std::vector<uint8_t> key3(24, 0x44);
    t1.set_key(key3);
    cipher_processor proc_des(&d1);
    cipher_processor proc_3des(&t1);
    std::vector<std::string> files = { "text.txt", "image.png", "music.mp3", "video.mp4", "empty.bin" };
    for (auto& f : files) {
        try {
            proc_des.process_file(f, f + ".des.enc", true);
            proc_des.process_file(f + ".des.enc", f + ".des.dec", false);
            if (compare_files(f, f + ".des.dec"))
                std::cout << "DES: " << f << " OK\n";
            else
                std::cout << "DES: " << f << " FAIL\n";
        }
        catch (...) {
            std::cout << "DES: " << f << " error\n";
        }
        try {
            proc_3des.process_file(f, f + ".3des.enc", true);
            proc_3des.process_file(f + ".3des.enc", f + ".3des.dec", false);
            if (compare_files(f, f + ".3des.dec"))
                std::cout << "3DES: " << f << " OK\n";
            else
                std::cout << "3DES: " << f << " FAIL\n";
        }
        catch (...) {
            std::cout << "3DES: " << f << " error\n";
        }
    }
    return 0;
}