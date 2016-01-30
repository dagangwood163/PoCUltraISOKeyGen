#include "hash.h"
#include <cstring>

const char* HASH::operator()(const char* data, const char* type, unsigned int data_len)
{
    unsigned int md_len;
    const unsigned char * md_value = hash(data, type, md_len, data_len);
    char* hashstr = const_cast<char*>(m_mtxb_hash.block_acquire());
    for (unsigned int i = 0; i < md_len; i++)
        hex2str(&hashstr[i * 2], md_value[i]);
    hashstr[md_len * 2] = '\0';
    m_mtxb_mdvalue.release(md_value);
    return hashstr;
}

const unsigned char* HASH::hash(const char* data, const char* type, unsigned int &len, unsigned int data_len)
{
    if (data_len == 0)
    {
        data_len = strlen(data);
    }
    unsigned char * md_value = const_cast<unsigned char*>(m_mtxb_mdvalue.block_acquire());
    EVP_MD_CTX *mdctx = EVP_MD_CTX_create();
    const EVP_MD* md = EVP_get_digestbyname(type);
    EVP_DigestInit_ex(mdctx, md, NULL);
    EVP_DigestUpdate(mdctx, data, data_len);
    EVP_DigestFinal_ex(mdctx, md_value, &len);
    EVP_MD_CTX_destroy(mdctx);
    return md_value;
}

void HASH::release_hash_str(const char* p_hash)
{
    m_mtxb_hash.release(p_hash);
}

void HASH::release_hash_mdvalue(unsigned const char* p_mdvalue)
{
    m_mtxb_mdvalue.release(p_mdvalue);
}

void HASH::hex2str(char* dest, unsigned int hex)
{
    if (hex < 0x0f)
    {
        dest[0] = '0';
        if (hex < 0x0a)
        {
            dest[1] = hex + 0x30;
        }
        else
        {
            dest[1] = hex + 0x57;
        }
    }
    else
    {
        unsigned int sndBit = 0;
        unsigned int fstBit = 0;

        sndBit = hex / 0x10;
        fstBit = hex - sndBit * 0x10;
        dest[0] = sndBit < 0x0a ? sndBit + 0x30 : sndBit + 0x57;
        dest[1] = fstBit < 0x0a ? fstBit + 0x30 : fstBit + 0x57;
    }
}

void HASH::upperStr(char* str)
{
    for (int i = strlen(str) - 1; i >= 0; i--)
    {
        str[i] = toupper(str[i]);
    }
}

void HASH::lowerStr(char* str)
{
    for (int i = strlen(str) - 1; i >= 0; i--)
    {
        str[i] = tolower(str[i]);
    }
}

