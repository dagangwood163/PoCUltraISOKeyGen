#include "UltraISOValidator.h"
#include <iostream>
#include <cstring>
#include "hash.h"
#include "lip.h"
#include "constdef.h"

static char modulus_n[] = "A70F8F62AA6E97D1";
static char d_str[] = "A7CAD9177AE95A9";
static char op_str5[] = "10001";
static char check_str6[] = "55776B";
static char UTRISO_str[] = "UTRISO";
static char REGCODE_str[] = "4BA90D54214AC938";
void RSA_decrypt(char* encrypted_data)
{

    std::string regcode;
    regcode.resize(64);

    long* m = nullptr;
    long* n = nullptr;
    long* result = nullptr;
    long* d = nullptr;
    zhsread(encrypted_data, &m);
    zhsread(modulus_n, &n);
    zhsread(d_str, &d);
    zexpmod(m, d, n, &result);
    zhswrite(const_cast<char*>(regcode.c_str()), result);
    std::cout << regcode << " len:" << strlen(regcode.c_str()) << std::endl;
}

void RSA_encrypt(char* plaintext)
{
    long* a = nullptr;
    long* b = nullptr;
    long* c = nullptr;
    long* d = nullptr;
    zhsread(plaintext, &a);
    zhsread(modulus_n, &c);
    zhsread(op_str5, &b);
    zexpmod(a, b, c, &d);//d=(a^b)%c
    std::string e;
    e.resize(64);
    zhswrite(const_cast<char*>(e.c_str()), d);

}
void validate()
{
    long* a = nullptr;
    long* b = nullptr;
    long* c = nullptr;
    long* d = nullptr;
    int flag = 0;
    std::string operated_regcode;
    operated_regcode.resize(56);
    int seat_num = 0;
    int seat_c = 0;

    std::string regcode;
    std::string name;
    name.resize(32);
    std::cout << "Registration Name: " << std::endl;
    std::cin.getline(const_cast<char *>(name.c_str()),name.size());
    name.resize(strlen(name.c_str()));
    std::cout << "Registration Code: " << std::endl;
    std::cin >> regcode;
    size_t pos;
    while ( pos = regcode.find("-"), pos!= std::string::npos)
    {
        regcode.replace(pos,1,"");
    }

    HASH hasher;

    zhsread(const_cast<char*>(regcode.c_str()), &a);
    zhsread(modulus_n, &c);
    zhsread(op_str5, &b);
    zexpmod(a, b, c, &d);//d=(a^b)%c
    zhswrite(const_cast<char*>(operated_regcode.c_str()), d);
    std::cout << "RSA-64 encrypted registration code: " << operated_regcode << std::endl;
    hasher.lowerStr(const_cast<char*>(operated_regcode.c_str()));

    if (operated_regcode[0] != check_str6[0] &&
        operated_regcode[0] != '4')
    {
        if (operated_regcode[0] != check_str6[2] &&
            operated_regcode[1] != check_str6[3])
        {
            if (operated_regcode[0] == check_str6[4] ||
                operated_regcode[1] == check_str6[5])
            {
                flag += 4;
                if (operated_regcode[14] < 'a')
                    seat_c = operated_regcode[14] - '0';
                else
                    seat_c = operated_regcode[14] - 'W';
                seat_c *= 16;
                if (operated_regcode[15] < 'a')
                    seat_c += operated_regcode[15] - '0';
                else
                    seat_c += operated_regcode[15] - 'W';
                seat_c -= 32;
                if (operated_regcode[8] < 'a')
                    seat_num = operated_regcode[8] - '0';
                else
                    seat_num = operated_regcode[8] - 'W';
                seat_num *= 16;
                if (operated_regcode[9] < 'a')
                    seat_num += operated_regcode[9] - '0';
                else
                    seat_num += operated_regcode[9] - 'W';
                seat_num -= 32;
                seat_num += seat_c << 6;
            }
        }
        else
        {
            flag += 4;
            if (operated_regcode[8] < 97)
                seat_num = operated_regcode[8] - '0';
            else
                seat_num = operated_regcode[8] - 'W';
            seat_num *= 16;
            if (operated_regcode[9] < 97)
                seat_num += operated_regcode[9] - '0';
            else
                seat_num += operated_regcode[9] - 'W';
            seat_num -= 32;
        }
    }
    else
    {

        ++flag;
        if (operated_regcode[1] == check_str6[1] || operated_regcode[1] == '6')
            ++flag;
        if (operated_regcode[8] == '5')
            ++flag;
        if (operated_regcode[9] == '3')
            ++flag;
    }

    std::string composite_str = "UTRISO" + name + regcode;

    std::string composite_md5 = hasher(composite_str.c_str(), "MD5");
    hasher.upperStr(const_cast<char*>(composite_md5.c_str()));
    std::cout << "composite UTRISO,name,regcode md5: " << composite_md5 << std::endl;

    unsigned int comp_md5_len;
    const unsigned char* hex_composite_md5 = hasher.hash(composite_str.c_str(), "MD5", comp_md5_len);
    if ((operated_regcode[6] < '2' || operated_regcode[6] > '2') && (operated_regcode[6] < 'd' || operated_regcode[6] > 'd'))
        --flag;
    if (operated_regcode[7] != 'a' && operated_regcode[7] != 'c' && operated_regcode[7] != 'b')
        --flag;
    flag -= 2;
    int len_data = 52413;
    int i = 0;
    int index = 0;
    do
    {

        index = (len_data + i) / 2;
        int result = memcmp(&data[6 * index], hex_composite_md5, 6);
        if (result <= 0)
        {
            if (result >= 0)
            {
                flag -= 2;
                if (flag == 0)
                {
                    //success
                    std::cout << "matched key check data at offset " << 6 * index << std::endl;
                }
                break;
            }
            i = index + 1;
        }
        else
        {
            len_data = index - 1;
        }
    } while (i <= len_data);
    std::cout << "last compared data: ";
    for (int i = 0; i < 6; i++)
    {
        printf("%02X", data[6 * index + i]);
    }
    std::cout << std::endl;




    std::string check_data_md5 = hasher(reinterpret_cast<const char*>(data), "MD5", 314484);
    hasher.upperStr(const_cast<char*>(check_data_md5.c_str()));
    std::cout << "key check data md5: " << check_data_md5 << std::endl;

    if (flag == 0)
    {
        std::cout << "correct registration!" << std::endl;
    }
    else
    {
        std::cout << "wrong registration!" << std::endl;
    }
    FREESPACE(a);
    FREESPACE(b);
    FREESPACE(c);
    FREESPACE(d);
}