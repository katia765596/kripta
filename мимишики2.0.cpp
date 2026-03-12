#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <thread>
#include <future>
#include <mutex>
#include <random>
#include <stdexcept>
#include <filesystem>
#include <algorithm>
#include <locale>
#include <iterator>
namespace fs = std::filesystem;
using bytearr = std::vector<uint8_t>;
std::mutex g_logmutex;
void loginfo(const std::string& msg) { std::lock_guard<std::mutex> lock(g_logmutex); std::cout << "[INFO] " << msg << std::endl; }
void logerror(const std::string& msg) { std::lock_guard<std::mutex> lock(g_logmutex); std::cerr << "[ERROR] " << msg << std::endl; }
class Isymmshifr {
public:
    virtual size_t bl_size() const = 0;
    virtual void encryptb(const uint8_t* in, uint8_t* out, const bytearr& key) = 0;
    virtual void decryptb(const uint8_t* in, uint8_t* out, const bytearr& key) = 0;
    virtual ~Isymmshifr() = default;
};
class Imode {
public:
    virtual void encrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t thr) = 0;
    virtual void decrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t thr) = 0;
    virtual ~Imode() = default;
};
class Ipad {
public:
    virtual void addp(bytearr& data, size_t bsize) = 0;
    virtual void rempad(bytearr& data, size_t bsize) = 0;
    virtual ~Ipad() = default;
};
class Zerosp : public Ipad {
public:
    void addp(bytearr& data, size_t bsize) override {
        size_t pad = bsize - (data.size() % bsize);
        if (pad != bsize) data.insert(data.end(), pad, 0);
    }
    void rempad(bytearr& data, size_t) override {
        while (!data.empty() && data.back() == 0) data.pop_back();
    }
};
class ANSIX923p : public Ipad {
public:
    void addp(bytearr& data, size_t bsize) override {
        size_t pad = bsize - (data.size() % bsize);
        if (pad == 0) pad = bsize;
        data.insert(data.end(), pad - 1, 0);
        data.push_back(static_cast<uint8_t>(pad));
    }
    void rempad(bytearr& data, size_t) override {
        if (data.empty()) return;
        uint8_t pad = data.back();
        if (pad > data.size()) throw std::runtime_error("Неверный паддинг ANSI X.923");
        data.erase(data.end() - pad, data.end());
    }
};
class PKCS7p : public Ipad {
public:
    void addp(bytearr& data, size_t bsize) override {
        size_t pad = bsize - (data.size() % bsize);
        data.insert(data.end(), pad, static_cast<uint8_t>(pad));
    }
    void rempad(bytearr& data, size_t) override {
        if (data.empty()) return;
        uint8_t pad = data.back();
        if (pad > data.size()) throw std::runtime_error("Неверный паддинг PKCS7");
        data.erase(data.end() - pad, data.end());
    }
};
class ISO10126p : public Ipad {
public:
    void addp(bytearr& data, size_t bsize) override {
        size_t pad = bsize - (data.size() % bsize);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(0, 255);
        for (size_t i = 0; i < pad - 1; i++)
            data.push_back(static_cast<uint8_t>(dis(gen)));
        data.push_back(static_cast<uint8_t>(pad));
    }
    void rempad(bytearr& data, size_t) override {
        if (data.empty()) return;
        uint8_t pad = data.back();
        if (pad > data.size()) throw std::runtime_error("Неверный паддинг ISO10126");
        data.erase(data.end() - pad, data.end());
    }
};
class ECBm : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr&, size_t ths) override {
        out = data;
        size_t blocks = data.size() / alg->bl_size();
        std::vector<std::thread> thre;
        for (size_t t = 0; t < ths; t++) {
            thre.emplace_back([=, &out]() {
                for (size_t i = t; i < blocks; i += ths)
                    alg->encryptb(&data[i * alg->bl_size()], &out[i * alg->bl_size()], key);
                });
        }
        for (auto& th : thre) th.join();
    }
    void decrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr&, size_t thr) override {
        out = data;
        size_t blocks = data.size() / alg->bl_size();
        std::vector<std::thread> ths;
        for (size_t t = 0; t < thr; t++) {
            ths.emplace_back([=, &out]() {
                for (size_t i = t; i < blocks; i += thr)
                    alg->decryptb(&data[i * alg->bl_size()], &out[i * alg->bl_size()], key);
                });
        }
        for (auto& th : ths) th.join();
    }
};
class CBCm : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        if (iv.size() != block) throw std::runtime_error("Неверный IV для CBC");
        out.resize(data.size());
        bytearr prev = iv;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr buf(data.begin() + i, data.begin() + i + block);
            for (size_t j = 0; j < block; j++) buf[j] ^= prev[j];
            alg->encryptb(buf.data(), &out[i], key);
            prev.assign(out.begin() + i, out.begin() + i + block);
        }
    }
    void decrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        if (iv.size() != block) throw std::runtime_error("Неверный IV для CBC");
        out.resize(data.size());
        bytearr prev = iv;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr temp(data.begin() + i, data.begin() + i + block);
            bytearr buf(block);
            alg->decryptb(&data[i], buf.data(), key);
            for (size_t j = 0; j < block; j++) out[i + j] = buf[j] ^ prev[j];
            prev = temp;
        }
    }
};
class PCBCm : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        if (iv.size() != block) throw std::runtime_error("Неверный IV для PCBC");
        out.resize(data.size());
        bytearr prev = iv;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr buf(data.begin() + i, data.begin() + i + block);
            for (size_t j = 0; j < block; j++) buf[j] ^= prev[j];
            alg->encryptb(buf.data(), &out[i], key);
            for (size_t j = 0; j < block; j++) prev[j] = buf[j] ^ out[i + j];
        }
    }
    void decrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        if (iv.size() != block) throw std::runtime_error("Неверный IV для PCBC");
        out.resize(data.size());
        bytearr prev = iv;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr cipher(data.begin() + i, data.begin() + i + block);
            bytearr buf(block);
            alg->decryptb(&data[i], buf.data(), key);
            for (size_t j = 0; j < block; j++) out[i + j] = buf[j] ^ prev[j];
            for (size_t j = 0; j < block; j++) prev[j] = out[i + j] ^ cipher[j];
        }
    }
};
class CFBm : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        out.resize(data.size());
        bytearr prev = iv;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr buf(block);
            alg->encryptb(prev.data(), buf.data(), key);
            size_t sz = std::min(block, data.size() - i);
            for (size_t j = 0; j < sz; j++) { out[i + j] = data[i + j] ^ buf[j]; prev[j] = out[i + j]; }
        }
    }
    void decrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        out.resize(data.size());
        bytearr prev = iv;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr buf(block);
            alg->encryptb(prev.data(), buf.data(), key);
            size_t sz = std::min(block, data.size() - i);
            for (size_t j = 0; j < sz; j++) { out[i + j] = data[i + j] ^ buf[j]; prev[j] = data[i + j]; }
        }
    }
};
class OFBm : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        out.resize(data.size());
        bytearr feedback = iv;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr buf(block);
            alg->encryptb(feedback.data(), buf.data(), key);
            feedback = buf;
            size_t sz = std::min(block, data.size() - i);
            for (size_t j = 0; j < sz; j++) out[i + j] = data[i + j] ^ buf[j];
        }
    }
    void decrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t thr) override {
        encrypt(data, out, alg, key, iv, thr);
    }
};
class CTRm : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        out.resize(data.size());
        uint64_t counter = 0;
        bytearr counterblock(block);
        for (size_t i = 0; i < data.size(); i += block) {
            for (size_t j = 0; j < block; j++) counterblock[j] = iv[j] ^ ((counter >> (8 * j)) & 0xFF);
            bytearr buf(block);
            alg->encryptb(counterblock.data(), buf.data(), key);
            size_t sz = std::min(block, data.size() - i);
            for (size_t j = 0; j < sz; j++) out[i + j] = data[i + j] ^ buf[j];
            counter++;
        }
    }
    void decrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t thr) override {
        encrypt(data, out, alg, key, iv, thr);
    }
};
class RandomDeltam : public Imode {
public:
    void encrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        out.resize(data.size());
        bytearr delta = iv;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(0, 255);;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr buf(block);
            for (size_t j = 0; j < block; j++) buf[j] = delta[j] ^ data[i + j];
            alg->encryptb(buf.data(), &out[i], key);
            for (size_t j = 0; j < block; j++) delta[j] = buf[j] ^ static_cast<uint8_t>(dis(gen));
        }
    }
    void decrypt(const bytearr& data, bytearr& out, Isymmshifr* alg,
        const bytearr& key, const bytearr& iv, size_t) override {
        size_t block = alg->bl_size();
        out.resize(data.size());
        bytearr delta = iv;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(0, 255);;
        for (size_t i = 0; i < data.size(); i += block) {
            bytearr buf(block);
            alg->decryptb(&data[i], buf.data(), key);
            for (size_t j = 0; j < block; j++) out[i + j] = buf[j] ^ delta[j];
            for (size_t j = 0; j < block; j++) delta[j] = buf[j] ^ static_cast<uint8_t>(dis(gen));;
        }
    }
};
enum class shifrmode { ECB, CBC, PCBC, CFB, OFB, CTR, RandomDelta };
enum class paddingmode { Zeros, ANSI_X923, PKCS7, ISO10126 };
class Symm {
public:
    Symm(std::shared_ptr<Isymmshifr> algo,
        shifrmode encmode,
        paddingmode padmode,
        const bytearr& iv = {})
        : algo_(algo), iv_(iv) {
        padding_ = createpadding(padmode);
        mode_ = createmode(encmode, iv_);
    }
    void encrypt(const bytearr& input, bytearr& out, const bytearr& key, size_t thr = 1) {
        bytearr data = input;
        padding_->addp(data, algo_->bl_size());
        mode_->encrypt(data, out, algo_.get(), key, iv_, thr);
    }
    void decrypt(const bytearr& input, bytearr& out, const bytearr& key, size_t thr = 1) {
        mode_->decrypt(input, out, algo_.get(), key, iv_, thr);
        padding_->rempad(out, algo_->bl_size());
    }
    std::future<void> encryptfasy(const std::string& inf, const std::string& outf,
        const bytearr& key, size_t thr = 1) {
        return std::async(std::launch::async, [=]() {
            std::ifstream in(inf, std::ios::binary);
            if (!in) throw std::runtime_error("Ошибка открытия файла: " + inf);
            bytearr buf((std::istreambuf_iterator<char>(in)), {});
            bytearr enc;
            encrypt(buf, enc, key, thr);
            std::ofstream out(outf, std::ios::binary);
            out.write((char*)enc.data(), enc.size());
            });
    }
    void encryptf(const std::string& f, const std::string& outf,
        const bytearr& key, size_t thr = 1)
    {
        fs::create_directories(outf);
        for (auto& entry : fs::recursive_directory_iterator(f))
        {
            if (!entry.is_regular_file())
                continue;
            try
            {
                std::string relp = fs::relative(entry.path(), f).string();
                std::string outf= (fs::path(outf) / relp).string();
                fs::create_directories(fs::path(outf).parent_path());
                std::ifstream in(entry.path(), std::ios::binary);
                if (!in)
                {
                    std::cout << "Не удалось открыть файл: " << entry.path().string() << std::endl;
                    continue;
                }
                bytearr buf((std::istreambuf_iterator<char>(in)), {});
                bytearr enc;
                encrypt(buf, enc, key, thr);
                std::ofstream out(outf, std::ios::binary);
                out.write((char*)enc.data(), enc.size());
                std::cout << "Файл зашифрован: " << entry.path().filename().string() << std::endl;
            }
            catch (...)
            {
                std::cout << "Не удалось зашифровать файл: "
                    << entry.path().filename().string() << std::endl;
            }
        }
        std::cout << std::endl;
        std::cout << "   Ваши файлы зашифрованы." << std::endl;
    }
private:
    std::shared_ptr<Isymmshifr> algo_;
    std::shared_ptr<Imode> mode_;
    std::shared_ptr<Ipad> padding_;
    bytearr iv_;
    std::shared_ptr<Imode> createmode(shifrmode m, const bytearr& iv) {
        switch (m) {
        case shifrmode::ECB: return std::make_shared<ECBm>();
        case shifrmode::CBC: return std::make_shared<CBCm>();
        case shifrmode::PCBC: return std::make_shared<PCBCm>();
        case shifrmode::CFB: return std::make_shared<CFBm>();
        case shifrmode::OFB: return std::make_shared<OFBm>();
        case shifrmode::CTR: return std::make_shared<CTRm>();
        case shifrmode::RandomDelta: return std::make_shared<RandomDeltam>();
        default: throw std::runtime_error("Unsupported mode");
        }
    }
    std::shared_ptr<Ipad> createpadding(paddingmode p) {
        switch (p) {
        case paddingmode::Zeros: return std::make_shared<Zerosp>();
        case paddingmode::ANSI_X923: return std::make_shared<ANSIX923p>();
        case paddingmode::PKCS7: return std::make_shared<PKCS7p>();
        case paddingmode::ISO10126: return std::make_shared<ISO10126p>();
        default: throw std::runtime_error("Unsupported padding");
        }
    }
};
class XORalg : public Isymmshifr {
public:
    size_t bl_size() const override { return 8; }
    void encryptb(const uint8_t* in, uint8_t* out, const bytearr& key) override {
        for (size_t i = 0; i < bl_size(); i++) out[i] = in[i] ^ key[i % key.size()];
    }
    void decryptb(const uint8_t* in, uint8_t* out, const bytearr& key) override { encryptb(in, out, key); }
};
bytearr randkey(size_t length) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> dis(0, 255); 
    bytearr key(length);
    for (auto& b : key)
        b = static_cast<uint8_t>(dis(gen)); 
    return key;
}
void unittests() {
    loginfo("Запуск unit-тестов...");
    auto algo = std::make_shared<XORalg>();
    bytearr key = randkey(8);
    bytearr iv(8, 0);
    std::vector<paddingmode> paddings = { paddingmode::Zeros, paddingmode::ANSI_X923, paddingmode::PKCS7, paddingmode::ISO10126 };
    std::vector<shifrmode> modes = { shifrmode::ECB,shifrmode::CBC,shifrmode::PCBC,
                                         shifrmode::CFB,shifrmode::OFB,shifrmode::CTR,shifrmode::RandomDelta };
    for (auto pad : paddings) {
        for (auto mode : modes) {
            Symm ctx(algo, mode, pad, iv);
            bytearr data = { 1,2,3,4,5,6,7 };
            bytearr encrypted, decrypted;
            ctx.encrypt(data, encrypted, key);
            ctx.decrypt(encrypted, decrypted, key);
            if (data != decrypted) throw std::runtime_error("Unit test failed for mode " + std::to_string((int)mode) + " pad " + std::to_string((int)pad));
        }
    }
    loginfo("Все unit-тесты пройдены!");
}
int main() {
    setlocale(LC_ALL, "Russian");
    try {
        auto algo = std::make_shared<XORalg>();
        bytearr key = randkey(8);
        std::cout << "Сгенерирован ключ: ";
        for (auto b : key) std::cout << (int)b << " ";
        std::cout << "\n";
        int modeinput;
        std::cout << "Выберите режим (0-ECB,1-CBC,2-PCBC,3-CFB,4-OFB,5-CTR,6-RandomDelta): ";
        std::cin >> modeinput;
        int padinput;
        std::cout << "Выберите паддинг (0-Zeros,1-ANSI X.923,2-PKCS7,3-ISO10126): ";
        std::cin >> padinput;
        bytearr iv(8, 0);
        Symm ctx(algo, (shifrmode)modeinput, (paddingmode)padinput, iv);
        std::string fpath, outf;
        std::cout << "Введите путь к папке для шифрования: ";
        std::cin >> fpath;
        std::cout << "Папка для сохранения зашифрованных файлов: ";
        std::cin >> outf;
        ctx.encryptf(fpath, outf, key, 2);
        unittests();
    }
    catch (const std::exception& e) { logerror(e.what()); return 1; }
    return 0;
}

