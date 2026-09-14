#ifndef ENUMS_H
#define ENUMS_H
enum class cipher_mode { ecb, cbc, pcbc, cfb, ofb, ctr, random_delta };
enum class padding_mode { zeros, ansi_x923, pkcs7, iso10126 };
#endif