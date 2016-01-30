//
// Created by bhzyf on 2016/1/30.
//

#include "FASTRSA64.h"


void add_mod(ull &a, ull b, ull m)
{
    b = a + b;
    if (b < a)
    {
        a = (b % m + m_res) % m;
    }
    else
    {
        a = b % m;
    }
}

void mul_mod(ull &a, ull b, ull m)
{
    ull tmp = a % m;
    a = 0;
    while (b)
    {
        if (b & 1)
            add_mod(a, tmp, m);

        add_mod(tmp, tmp, m);
        b >>= 1;
    }
}

ull exp_mod(ull c, ull e, ull n)
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