#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <string.h>
#include <stdio.h>
#define m_res  0x58F0709D5591682F
typedef unsigned long long ull;

__device__ __inline__ void add_mod(ull &a, ull b, ull m)
{
    b = a + b;
    if (b < a)
    {
        a = (b%m + m_res) % m;
    }
    else
    {
        a = b%m;
    }
}

__device__ __inline__ void  mul_mod(ull &a, ull b, ull m)
{
    ull tmp = a%m;
    a = 0;
    while (b)
    {
        if (b & 1)
            add_mod(a, tmp, m);

        add_mod(tmp, tmp, m);
        b >>= 1;
    }
}

__device__ __inline__ ull exp_mod(ull c, ull e, ull n)
{
    c %= n;
    ull result = 1;
    while (e > 0)
    {
        if (e & 1)
        {
            mul_mod(result, c, n);//result = (result*c) % n;
        }
        e >>= 1;
        mul_mod(c, c, n);//c = (c*c) % n;
    }

    return result;
}



typedef unsigned long long ull;
typedef unsigned int uint;
typedef unsigned char uchar;
#define shift(x, n) (((x) << (n)) | ((x) >> (32-(n))))
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))    
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) ((y) ^ ((x) | (~z)))

#define FF(a ,b ,c ,d ,Mj ,s ,ti)  a = b + (shift((a + F(b, c, d) + Mj + ti) , s))
#define GG(a, b, c, d, Mj, s, ti)  a = b + (shift((a + G(b, c, d) + Mj + ti) , s))
#define HH(a, b, c, d, Mj, s, ti)  a = b + (shift((a + H(b, c, d) + Mj + ti) , s))
#define II(a, b, c, d, Mj, s, ti)  a = b + (shift((a + I(b, c, d) + Mj + ti) , s))

#define A 0x67452301
#define B 0xefcdab89
#define C 0x98badcfe
#define D 0x10325476

//msg must be padded 512 bit.
__device__ __inline__ ull md5_512(uint* msg)
{
    uint a = A, b = B, c = C, d = D;

    FF(a, b, c, d, msg[0], 7, 0xd76aa478);
    FF(d, a, b, c, msg[1], 12, 0xe8c7b756);
    FF(c, d, a, b, msg[2], 17, 0x242070db);
    FF(b, c, d, a, msg[3], 22, 0xc1bdceee);
    FF(a, b, c, d, msg[4], 7, 0xf57c0faf);
    FF(d, a, b, c, msg[5], 12, 0x4787c62a);
    FF(c, d, a, b, msg[6], 17, 0xa8304613);
    FF(b, c, d, a, msg[7], 22, 0xfd469501);
    FF(a, b, c, d, msg[8], 7, 0x698098d8);
    FF(d, a, b, c, msg[9], 12, 0x8b44f7af);
    FF(c, d, a, b, msg[10], 17, 0xffff5bb1);
    FF(b, c, d, a, msg[11], 22, 0x895cd7be);
    FF(a, b, c, d, msg[12], 7, 0x6b901122);
    FF(d, a, b, c, msg[13], 12, 0xfd987193);
    FF(c, d, a, b, msg[14], 17, 0xa679438e);
    FF(b, c, d, a, msg[15], 22, 0x49b40821);

    GG(a, b, c, d, msg[1], 5, 0xf61e2562);
    GG(d, a, b, c, msg[6], 9, 0xc040b340);
    GG(c, d, a, b, msg[11], 14, 0x265e5a51);
    GG(b, c, d, a, msg[0], 20, 0xe9b6c7aa);
    GG(a, b, c, d, msg[5], 5, 0xd62f105d);
    GG(d, a, b, c, msg[10], 9, 0x02441453);
    GG(c, d, a, b, msg[15], 14, 0xd8a1e681);
    GG(b, c, d, a, msg[4], 20, 0xe7d3fbc8);
    GG(a, b, c, d, msg[9], 5, 0x21e1cde6);
    GG(d, a, b, c, msg[14], 9, 0xc33707d6);
    GG(c, d, a, b, msg[3], 14, 0xf4d50d87);
    GG(b, c, d, a, msg[8], 20, 0x455a14ed);
    GG(a, b, c, d, msg[13], 5, 0xa9e3e905);
    GG(d, a, b, c, msg[2], 9, 0xfcefa3f8);
    GG(c, d, a, b, msg[7], 14, 0x676f02d9);
    GG(b, c, d, a, msg[12], 20, 0x8d2a4c8a);

    HH(a, b, c, d, msg[5], 4, 0xfffa3942);
    HH(d, a, b, c, msg[8], 11, 0x8771f681);
    HH(c, d, a, b, msg[11], 16, 0x6d9d6122);
    HH(b, c, d, a, msg[14], 23, 0xfde5380c);
    HH(a, b, c, d, msg[1], 4, 0xa4beea44);
    HH(d, a, b, c, msg[4], 11, 0x4bdecfa9);
    HH(c, d, a, b, msg[7], 16, 0xf6bb4b60);
    HH(b, c, d, a, msg[10], 23, 0xbebfbc70);
    HH(a, b, c, d, msg[13], 4, 0x289b7ec6);
    HH(d, a, b, c, msg[0], 11, 0xeaa127fa);
    HH(c, d, a, b, msg[3], 16, 0xd4ef3085);
    HH(b, c, d, a, msg[6], 23, 0x04881d05);
    HH(a, b, c, d, msg[9], 4, 0xd9d4d039);
    HH(d, a, b, c, msg[12], 11, 0xe6db99e5);
    HH(c, d, a, b, msg[15], 16, 0x1fa27cf8);
    HH(b, c, d, a, msg[2], 23, 0xc4ac5665);

    II(a, b, c, d, msg[0], 6, 0xf4292244);
    II(d, a, b, c, msg[7], 10, 0x432aff97);
    II(c, d, a, b, msg[14], 15, 0xab9423a7);
    II(b, c, d, a, msg[5], 21, 0xfc93a039);
    II(a, b, c, d, msg[12], 6, 0x655b59c3);
    II(d, a, b, c, msg[3], 10, 0x8f0ccc92);
    II(c, d, a, b, msg[10], 15, 0xffeff47d);
    II(b, c, d, a, msg[1], 21, 0x85845dd1);
    II(a, b, c, d, msg[8], 6, 0x6fa87e4f);
    II(d, a, b, c, msg[15], 10, 0xfe2ce6e0);
    II(c, d, a, b, msg[6], 15, 0xa3014314);
    II(b, c, d, a, msg[13], 21, 0x4e0811a1);
    II(a, b, c, d, msg[4], 6, 0xf7537e82);
    II(d, a, b, c, msg[11], 10, 0xbd3af235);
    II(c, d, a, b, msg[2], 15, 0x2ad7d2bb);
    II(b, c, d, a, msg[9], 21, 0xeb86d391);

    a = a + A;
    b = b + B;
    //c = c + C;
    //d = d + D;
    return (static_cast<ull>(b) << 32) | a;
}


__constant__ char hex2str_map[] = "0123456789ABCDEF";
char* device_msg;
int device_msg_len;
unsigned long long * device_data;
void check(cudaError ret)
{
    if (ret != cudaError::cudaSuccess)
        printf("cuda failed!\n");
}

__global__ void KGkernel(unsigned long long* data, char* msg_, int len)
{
    int i =blockIdx.x*blockDim.x + threadIdx.x;
    ull pat = data[i];
    ull num = exp_mod(pat, 0x58C5D3F7, 0xF513E783);
    pat = exp_mod(pat, 0xAC3A102B, 0xAE818F1B);
    mul_mod(num, 0x5C4AF104DA37C96D, 0xA70F8F62AA6E97D1);
    mul_mod(pat, 0x4AC49E5DD036CE65, 0xA70F8F62AA6E97D1);
    add_mod(num,pat,0xA70F8F62AA6E97D1);
    uchar msg[64]{ 0 };
    memset(msg, 0, 64);
    *reinterpret_cast<uint*>(&msg[56]) = (len + 16) * 8;
    memcpy(msg, msg_, len);
    msg[len + 16] = 0x80;

    //hex to str start.
    (&msg[len])[0] = hex2str_map[(num >> 15 * 4) & 0xF];
    (&msg[len])[1] = hex2str_map[(num >> 14 * 4) & 0xF];
    (&msg[len])[2] = hex2str_map[(num >> 13 * 4) & 0xF];
    (&msg[len])[3] = hex2str_map[(num >> 12 * 4) & 0xF];
    (&msg[len])[4] = hex2str_map[(num >> 11 * 4) & 0xF];
    (&msg[len])[5] = hex2str_map[(num >> 10 * 4) & 0xF];
    (&msg[len])[6] = hex2str_map[(num >> 9 * 4) & 0xF];
    (&msg[len])[7] = hex2str_map[(num >> 8 * 4) & 0xF];
    (&msg[len])[8] = hex2str_map[(num >> 7 * 4) & 0xF];
    (&msg[len])[9] = hex2str_map[(num >> 6 * 4) & 0xF];
    (&msg[len])[10] = hex2str_map[(num >> 5 * 4) & 0xF];
    (&msg[len])[11] = hex2str_map[(num >> 4 * 4) & 0xF];
    (&msg[len])[12] = hex2str_map[(num >> 3 * 4) & 0xF];
    (&msg[len])[13] = hex2str_map[(num >> 2 * 4) & 0xF];
    (&msg[len])[14] = hex2str_map[(num >> 1 * 4) & 0xF];
    (&msg[len])[15] = hex2str_map[num & 0xF];
    //hex to str end.

    num = md5_512(reinterpret_cast<uint*>(msg));
    data[i] = num;
}

void initKG(char* msg)
{
    device_msg_len = strlen(msg);
    check(cudaMalloc(&device_msg, device_msg_len*sizeof(char)));
    check(cudaMemcpy(device_msg, msg, device_msg_len*sizeof(char), cudaMemcpyHostToDevice));
}

void exitKG()
{
    check(cudaFree(device_msg));
}


//the data also hold the calculated data.
//data len is blocks*threads.
void runKGkernel(unsigned long long* data, int blocks, int threads)
{
    check(cudaMalloc((void**)&device_data, blocks*threads*sizeof(unsigned long long)));
    check(cudaMemcpy(device_data, data, blocks*threads*sizeof(unsigned long long), cudaMemcpyHostToDevice));

    KGkernel << <blocks, threads >> >(device_data, device_msg, device_msg_len);
    check(cudaMemcpy(data, device_data, blocks*threads*sizeof(unsigned long long), cudaMemcpyDeviceToHost));
    check(cudaDeviceSynchronize());
    

    check(cudaFree(device_data));
}