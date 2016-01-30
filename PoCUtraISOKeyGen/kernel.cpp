
#include "third_party\gmp-6.0.0\gmp.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "KGCPU.h"
#include "KGCUDA.h"

int main()
{
    printf(
        R"(THIS IS A PROOF OF CONCEPT KEYGEN OF ULTRAISO.
IT ONLY HAS THE POSSIBILITIES TO CALCULATE THE CORRECT KEY.
IT USES BOTH CPU AND GPU TO ACCELERATE CALCULATION.

SRC ON GITHUB: https://github.com/yufengzjj/PoCUtraISOKeyGen.git


GOOD LUCK!


)");
    KGCPU cpu;
    cpu.drive();
    char input;
    while (scanf("%c",&input) != EOF);
    return 0;
}

