#include "lan_key_exchange.h"//обеспечивает обмен ключами по локалке с протоколом,испол усл компилчяцию
#include <stdexcept>
#include <string>
#ifdef _WIN32
#include <winsock2.h>//для работы с соекатами
#include <ws2tcpip.h>//для работы с протоколами TCP/IP
#pragma comment(lib, "Ws2_32.lib")//указ компоновщику включить библ
using socket_type = SOCKET;//тип определенный в винсок для дескриптора сокета
const socket_type invalid_socket_value = INVALID_SOCKET;//для проверки ошибок при создании сокета
#else
#include <arpa/inet.h>//для работы с сетевыми адресами
#include <netinet/in.h>//определение структур для интернет-адресов
#include <sys/socket.h>
#include <unistd.h>//для системных вызово unix
using socket_type = int;//дескриптор файла
const socket_type invalid_socket_value = -1;
#endif
namespace
{
    class winsock_guard//для иниц и очистки винсок
    {
    public:
        winsock_guard()
        {
#ifdef _WIN32
            WSADATA data;//для хранения инфы и реал винсок
            if (WSAStartup(MAKEWORD(2, 2), &data) != 0)
                throw std::runtime_error("WSAStartup failed");
#endif
        }
        ~winsock_guard()
        {
#ifdef _WIN32
            WSACleanup();//деструктор для освобож ресурсов винсок
#endif
        }
    };
    void close_socket(socket_type socket)
    {
#ifdef _WIN32
        closesocket(socket);
#else
        close(socket);
#endif
    }
    void send_all(socket_type socket, const char* data, size_t size)//декср сокета, указ на данные и кол-во ба  т
    {//отпрака всех данных через сокет
        size_t sent = 0;//скок отправлено
        while (sent < size)
        {
#ifdef _WIN32
            int result = send(socket, data + sent, static_cast<int>(size - sent), 0);
#else
            ssize_t result = send(socket, data + sent, size - sent, 0);
#endif
            if (result <= 0)
                throw std::runtime_error("send failed");
            sent += static_cast<size_t>(result);
        }
    }
    void receive_all(socket_type socket, char* data, size_t size)//получает все сайз байт из сокета
    {
        size_t received = 0;
        while (received < size)
        {
#ifdef _WIN32
            int result = recv(socket, data + received, static_cast<int>(size - received), 0);
#else//recv для чтения данных
            ssize_t result = recv(socket, data + received, size - received, 0);
#endif
            if (result <= 0)
                throw std::runtime_error("receive failed");
            received += static_cast<size_t>(result);
        }
    }
    void send_u64(socket_type socket, uint64_t value)//функция отправ 64битное через сокет в сетевом порядке байтов бигэнджен
    {
        char data[8];
        for (int i = 7; i >= 0; --i)
        {
            data[i] = static_cast<char>(value & 0xff);
            value >>= 8;
        }
        send_all(socket, data, 8);
    }
    uint64_t receive_u64(socket_type socket)//читает 64битное из сокета
    {
        char data[8];
        receive_all(socket, data, 8);
        uint64_t value = 0;
        for (int i = 0; i < 8; ++i)
            value = (value << 8) | static_cast<unsigned char>(data[i]);//число в лок порядке байтов
        return value;
    }
    socket_type make_server(uint16_t port)//созд серверный сокет принимает порт, возвр дескриптор сокета
    {
        socket_type server = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);//созд соект с помощью системного вызова сокет, af_inet-семейство адресов IPv4, сок_стрим-тип соке для потоковой передачи (tcp),ипрото_тсп-протокол TCP
        if (server == invalid_socket_value)
            throw std::runtime_error("socket creation failed");
        sockaddr_in address{};//объявляет и иниц нулями структуру (адрес сокета для IPv4)
        address.sin_family = AF_INET;
        address.sin_addr.s_addr = htonl(INADDR_ANY);//конст, означает принимать соед на любом сетевом интерфейсе,htonl преобразует это значение в сетевой порядок байтов (big-endian)
        address.sin_port = htons(port);//преоб порт в сетевой порядок байтов 
        if (bind(server, reinterpret_cast<sockaddr*>(&address), sizeof(address)) != 0 || listen(server, 1) != 0)
        {//привязывает сокет к адресу и порту,принимает сокет,указ на структуру sockaddrЮ и размер структуру), переводит с помощью листен в режим прослушивая (ожид вход соедин), второй параметр 1 - макс длина очереди ожид соед
            close_socket(server);
            throw std::runtime_error("server socket setup failed");
        }
        return server;//возвр дескриптор серверного сокета
    }
    socket_type make_client(const std::string& address, uint16_t port)
    {//для созд клиентского сокета и подклю к серверу
        socket_type socket_value = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
        if (socket_value == invalid_socket_value)
            throw std::runtime_error("socket creation failed");
        sockaddr_in server{};
        server.sin_family = AF_INET;
        server.sin_port = htons(port);
        if (inet_pton(AF_INET, address.c_str(), &server.sin_addr) != 1)
        {//преобразует строковый ip-адрес в бин формат и сохр в син_адр, возвр 1 при успехе
            close_socket(socket_value);
            throw std::runtime_error("invalid address");
        }
        if (connect(socket_value, reinterpret_cast<sockaddr*>(&server), sizeof(server)) != 0)
        {//подкл к серверу
            close_socket(socket_value);
            throw std::runtime_error("connect failed");
        }
        return socket_value;
    }//возвр подкл сокет
}
byte_array lan_key_exchange::run_server(const diffie_hellman& participant, uint16_t port, size_t key_size)
{//возвр общий секретный ключ в виде массива
    winsock_guard guard;
    socket_type server = make_server(port);
    sockaddr_in client_address{};//структура для хранения адреса клиента
#ifdef _WIN32
    int client_size = sizeof(client_address);//размер адреса клиента
#else
    socklen_t client_size = sizeof(client_address);
#endif
    socket_type client = accept(server, reinterpret_cast<sockaddr*>(&client_address), &client_size);//асепт принимает вход соед на серверном сокете,заполняет клиент адрес, адресом клиента, возвр дескриптоп сокета для общения с клиентом, блок выполн до появления соединения
    close_socket(server);
    if (client == invalid_socket_value)
        throw std::runtime_error("accept failed");
    try
    {
        send_u64(client, participant.get_public_key());//отпр отк ключ участника через сокет 
        uint64_t other = receive_u64(client);//получает откр ключ другой стороны
        byte_array result = participant.make_key(other, key_size);//вычисляет общий секрет
        close_socket(client);
        return result;
    }
    catch (...)
    {
        close_socket(client);
        throw;
    }
}
byte_array lan_key_exchange::run_client(const diffie_hellman& participant, const std::string& address, uint16_t port, size_t key_size)
{//клиентская часть обмена, принимает объект участника, адрес сервера, порт и размер ключа, возвр общий ключ
    winsock_guard guard;
    socket_type socket_value = make_client(address, port);//подкл клиентски сокет к серверу
    try
    {
        uint64_t other = receive_u64(socket_value);//клиент первый получает открытый ключ сервера
        send_u64(socket_value, participant.get_public_key());//отправляет откр ключ свой
        byte_array result = participant.make_key(other, key_size);//вычисляет общиц секрет
        close_socket(socket_value);
        return result;
    }
    catch (...)
    {
        close_socket(socket_value);
        throw;
    }
}
