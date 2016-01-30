
#ifndef HASH_H
#define HASH_H
#include <openssl/evp.h>
#include <string>
#include "MutexBlock.h"
class HASH
{
public:
    explicit HASH()
    {
       
        OpenSSL_add_all_digests();
    }

    ~HASH()
    {
        
        /* Call this once before exit. */
        EVP_cleanup();
    }
    void upperStr(char* str);
    void lowerStr(char* str);
    void hex2str(char* dest, unsigned int hex);
    const char* operator()(const char* data, const char* type, unsigned int data_len = 0);
    const unsigned char* hash(const char* data, const char* type, unsigned int &len, unsigned int data_len = 0);
    void release_hash_str(const char*);
    void release_hash_mdvalue(unsigned  const char*);
private:
    MutexBlock<unsigned char, EVP_MAX_MD_SIZE,16> m_mtxb_mdvalue;
    MutexBlock<char, 17, 64> m_mtxb_hash;
};
#endif