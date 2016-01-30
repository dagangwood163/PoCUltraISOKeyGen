#include "Helper.h"


Helper::factor Helper::factoring(int num, int max_factor_b)
{
    int factor_a = 1, factor_b = 1;
    while (factor_a * factor_b != num)
    {
        if (num%factor_a != 0)
        {
            factor_a += 1;
        }
        else
        {
            factor_b += 1;
        }
        if (factor_b > max_factor_b)
        {
            factor_b = 1;
            factor_a += 1;
        }

    }
    return{ factor_a, factor_b };
}

Helper::Helper()
{
}


Helper::~Helper()
{
}
