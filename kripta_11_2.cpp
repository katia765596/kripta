#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "ws2_32.lib")
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#define SOCKET int
#define INVALID_SOCKET (SOCKET)(~0)
#define SOCKET_ERROR (-1)
#define closesocket close
#endif
class DiffieHellman {
public:
    using Integer = uint64_t;

    DiffieHellman(Integer p, Integer g) : p(p), g(g) {
        if (p < 2 || g < 2 || g >= p)
            throw std::invalid_argument("Invalid DH parameters");
    }
    Integer generatePrivateKey() {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<Integer> dist(2, p - 2);
        privateKey = dist(gen);
        return privateKey;
    }
    Integer computePublicKey() {
        if (privateKey == 0)
            throw std::logic_error("Private key not generated");
        publicKey = powm(g, privateKey, p);
        return publicKey;
    }
    Integer computeSharedSecret(Integer peerPublic) {
        if (privateKey == 0)
            throw std::logic_error("Private key not generated");
        sharedSecret = powm(peerPublic, privateKey, p);
        return sharedSecret;
    }
    Integer getSharedSecret() const {
        if (sharedSecret == 0)
            throw std::logic_error("Shared secret not computed");
        return sharedSecret;
    }
    std::vector<uint8_t> deriveKey(size_t bits) const {
        if (sharedSecret == 0)
            throw std::logic_error("Shared secret not computed");
        std::vector<uint8_t> secretBytes;
        Integer tmp = sharedSecret;
        while (tmp > 0) {
            secretBytes.push_back(static_cast<uint8_t>(tmp & 0xFF));
            tmp >>= 8;
        }
        if (secretBytes.empty()) secretBytes.push_back(0);
        size_t bytes = (bits + 7) / 8;
        std::vector<uint8_t> key(bytes, 0);
        for (size_t i = 0; i < secretBytes.size() && i < bytes; ++i)
            key[i] = secretBytes[i];
        size_t copied = secretBytes.size();
        while (copied < bytes) {
            for (size_t i = 0; i < secretBytes.size() && copied < bytes; ++i)
                key[copied++] = secretBytes[i];
        }
        return key;
    }
private:
    Integer p, g;
    Integer privateKey = 0;
    Integer publicKey = 0;
    Integer sharedSecret = 0;
    static Integer powm(Integer base, Integer exp, Integer mod) {
        Integer result = 1;
        base %= mod;
        while (exp > 0) {
            if (exp & 1) result = (result * base) % mod;
            base = (base * base) % mod;
            exp >>= 1;
        }
        return result;
    }
};
const uint64_t DH_P = 23;
const uint64_t DH_G = 5;
void initWinsock() {
#ifdef _WIN32
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0)
        throw std::runtime_error("WSAStartup failed");
#endif
}
void cleanupWinsock() {
#ifdef _WIN32
    WSACleanup();
#endif
}
SOCKET createServerSocket(int port) {
    SOCKET listenSock = socket(AF_INET, SOCK_STREAM, 0);
    if (listenSock == INVALID_SOCKET)
        throw std::runtime_error("socket failed");
    sockaddr_in addr;
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(port);
    if (bind(listenSock, (sockaddr*)&addr, sizeof(addr)) == SOCKET_ERROR) {
        closesocket(listenSock);
        throw std::runtime_error("bind failed");
    }
    if (listen(listenSock, 1) == SOCKET_ERROR) {
        closesocket(listenSock);
        throw std::runtime_error("listen failed");
    }
    return listenSock;
}
SOCKET connectToServer(const std::string& serverIp, int port) {
    SOCKET sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == INVALID_SOCKET)
        throw std::runtime_error("socket failed");
    sockaddr_in addr;
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    addr.sin_addr.s_addr = inet_addr(serverIp.c_str());
    if (addr.sin_addr.s_addr == INADDR_NONE)
        throw std::runtime_error("invalid IP address");
    if (connect(sock, (sockaddr*)&addr, sizeof(addr)) == SOCKET_ERROR) {
        closesocket(sock);
        throw std::runtime_error("connect failed");
    }
    return sock;
}
void sendUint64(SOCKET sock, uint64_t value) {
    unsigned char buf[8];
    for (int i = 0; i < 8; ++i) {
        buf[i] = (value >> (56 - 8 * i)) & 0xFF;
    }
    send(sock, (char*)buf, 8, 0);
}
uint64_t recvUint64(SOCKET sock) {
    unsigned char buf[8];
    int recvd = 0;
    while (recvd < 8) {
        int r = recv(sock, (char*)(buf + recvd), 8 - recvd, 0);
        if (r <= 0) throw std::runtime_error("recv failed");
        recvd += r;
    }
    uint64_t value = 0;
    for (int i = 0; i < 8; ++i) {
        value = (value << 8) | buf[i];
    }
    return value;
}
void runServer(int port) {
    initWinsock();
    SOCKET listenSock = createServerSocket(port);
    std::cout << "Server: waiting for client on port " << port << " ...\n";
    SOCKET clientSock = accept(listenSock, nullptr, nullptr);
    if (clientSock == INVALID_SOCKET) {
        closesocket(listenSock);
        cleanupWinsock();
        throw std::runtime_error("accept failed");
    }
    closesocket(listenSock);
    std::cout << "Server: client connected.\n";
    uint64_t clientPub = recvUint64(clientSock);
    DiffieHellman dh(DH_P, DH_G);
    dh.generatePrivateKey();
    uint64_t serverPub = dh.computePublicKey();
    sendUint64(clientSock, serverPub);
    uint64_t shared = dh.computeSharedSecret(clientPub);
    std::cout << "Server: shared secret = " << shared << "\n\n";
    std::vector<size_t> keySizes = { 56, 128, 128 };
    std::vector<std::string> algNames = { "DES", "AES", "MARS" };
    for (size_t i = 0; i < keySizes.size(); ++i) {
        auto key = dh.deriveKey(keySizes[i]);
        std::cout << algNames[i] << " key (" << keySizes[i] << " bits): ";
        for (auto b : key) printf("%02X", b);
        std::cout << std::endl;
    }

    closesocket(clientSock);
    cleanupWinsock();
}
void runClient(const std::string& serverIp, int port) {
    initWinsock();
    SOCKET sock = connectToServer(serverIp, port);
    std::cout << "Client: connected to server.\n";
    DiffieHellman dh(DH_P, DH_G);
    dh.generatePrivateKey();
    uint64_t clientPub = dh.computePublicKey();
    sendUint64(sock, clientPub);
    uint64_t serverPub = recvUint64(sock);
    uint64_t shared = dh.computeSharedSecret(serverPub);
    std::cout << "Client: shared secret = " << shared << "\n\n";

    std::vector<size_t> keySizes = { 56, 128, 128 };
    std::vector<std::string> algNames = { "DES", "AES", "MARS" };
    for (size_t i = 0; i < keySizes.size(); ++i) {
        auto key = dh.deriveKey(keySizes[i]);
        std::cout << algNames[i] << " key (" << keySizes[i] << " bits): ";
        for (auto b : key) printf("%02X", b);
        std::cout << std::endl;
    }

    closesocket(sock);
    cleanupWinsock();
}
void testDHMath() {
    DiffieHellman alice(DH_P, DH_G);
    DiffieHellman bob(DH_P, DH_G);
    alice.generatePrivateKey();
    bob.generatePrivateKey();
    auto pubA = alice.computePublicKey();
    auto pubB = bob.computePublicKey();
    auto secretA = alice.computeSharedSecret(pubB);
    auto secretB = bob.computeSharedSecret(pubA);
    assert(secretA == secretB);
    std::cout << "TEST: DH math test: PASSED\n";
    auto key56a = alice.deriveKey(56);
    auto key56b = bob.deriveKey(56);
    assert(key56a == key56b);
    std::cout << "TEST: Key derivation test: PASSED\n";
}
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage:\n"
            << "  " << argv[0] << " server port\n"
            << "  " << argv[0] << " client port\n"
            << "  " << argv[0] << " test\n";
        return 1;
    }
    std::string mode = argv[1];
    int port = 8080;
    if (argc > 2) port = std::stoi(argv[2]);
    if (mode == "test") {
        testDHMath();
        return 0;
    }
    try {
        if (mode == "server")
            runServer(port);
        else if (mode == "client")
            runClient("127.0.0.1", port);
        else
            std::cerr << "Unknown mode\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}