
#include "third_party\gmp-6.0.0\gmp.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "KGCPU.h"
#include "KGCUDA.h"

int main()
{
    KGCPU cpu;
    cpu.drive();
    char input;
    while (scanf("%c",&input) != EOF);
    return 0;
}

