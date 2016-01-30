#pragma once
class Helper
{
public:
    struct factor
    {
        int a;
        int b;
    };

    /*
    * always return the minimum factor_a and the maximum under
    * max_factor_b.
    * */
    static factor factoring(int num, int max_factor_b);
    Helper();
    ~Helper();
};

