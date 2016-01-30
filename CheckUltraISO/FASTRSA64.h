//
// Created by bhzyf on 2016/1/30.
//

#ifndef CHECKULTRAISO_FASTRSA64_H
#define CHECKULTRAISO_FASTRSA64_H
#define m_res  0x58F0709D5591682F
typedef unsigned long long ull;
void add_mod(ull &a, ull b, ull m);
void  mul_mod(ull &a, ull b, ull m);
ull exp_mod(ull c, ull e, ull n);
#endif //CHECKULTRAISO_FASTRSA64_H
