#include "KGCUDA.h"
#include "device_launch_parameters.h"
#include "Task.h"
#include "third_party\gmp-6.0.0\gmp.h"
#include "Helper.h"
#include <cstdio>
static char modulus_n[] = "A70F8F62AA6E97D1";
static char d_str[] = "A7CAD9177AE95A9";
static mpz_t d, n;
KGCUDA::KGCUDA()
{
    cudaGetDeviceProperties(&m_gpu_prop, 0);
    compatible();
    mpz_init_set_str(n, modulus_n, 16);
    mpz_init_set_str(d, d_str, 16);
}

void exitKG();
KGCUDA::~KGCUDA()
{
    mpz_clears(n, d, NULL);
    exitKG();
}

void runKGkernel(unsigned long long* data, int blocks, int threads);
void KGCUDA::dispatch(Task* _task)
{
    Helper::factor factor;
    if (_task->count == capability())
    {
        factor = { num_blocks, num_threads_preblock };
    }
    else
    {
        factor = Helper::factoring(_task->count, num_threads_preblock);
    }
    unsigned long long* data = generate(_task);
    runKGkernel(data, factor.a, factor.b);
    find(data, factor.a*factor.b);
    delete data;
}

void initKG(char* msg);
void KGCUDA::set_name(std::string name)
{
    m_name = name;
    initKG(const_cast<char*>(m_name.c_str()));
}

void KGCUDA::print_gpu_properties()
{
    printf("GPU:%s\n", m_gpu_prop.name);
    printf("multiProcessor Count:%d\n", m_gpu_prop.multiProcessorCount);
    printf("launched threads:%d\n", num_threads_preblock);
    printf("launched blocks:%d\n", num_blocks);
}

void KGCUDA::test()
{
    test_capability();
}

void KGCUDA::generate_test()
{

}

void KGCUDA::find(unsigned long long* md5, int len)
{
    extern unsigned char data[314484];
    for (int step = len - 1; step >= 0; step--)
    {
        int len_data = 52413;
        int i = 0;
        int index;
        do
        {
            index = (len_data + i) / 2;
            int result = memcmp(&data[6 * index], reinterpret_cast<unsigned char*>(&md5[step]), 6);
            if (result <= 0)
            {
                if (result >= 0)
                {
                    //success
                    locate(step);
                    break;

                }
                i = index + 1;
            }
            else
            {
                len_data = index - 1;
            }
        } while (i <= len_data);
    }


}

unsigned long long* KGCUDA::generate(Task* _task)
{
    using ull = unsigned long long;
    ull pat;

    short* selected_pos = nullptr;
    short curr_line = 15; // 0-15
    // there shouldn't minus 1 as we do not want to roll selected position with 1.
    // thus we start after (globe) selected position.
    selected_pos = _task->start_range;
    curr_task = *_task;
    ull* data_p = new ull[_task->count];
    unsigned int count_i = 0;
    // generate strings start
    while (curr_line >= 0 && _task->count > 0)
    {
        if (++selected_pos[curr_line] > Task::m_globe_range[curr_line * 2 + 1])
        {
            selected_pos[curr_line] = Task::m_globe_range[curr_line * 2];
            --curr_line;
        }
        else
        {
            _task->count -= 1;
            curr_line = 15;
            pat = 0;
            pat |= static_cast<ull>(m_dict[0][selected_pos[0]]) << 15 * 4;
            pat |= static_cast<ull>(m_dict[1][selected_pos[1]]) << 14 * 4;
            pat |= static_cast<ull>(m_dict[2][selected_pos[2]]) << 13 * 4;
            pat |= static_cast<ull>(m_dict[3][selected_pos[3]]) << 12 * 4;
            pat |= static_cast<ull>(m_dict[4][selected_pos[4]]) << 11 * 4;
            pat |= static_cast<ull>(m_dict[5][selected_pos[5]]) << 10 * 4;
            pat |= static_cast<ull>(m_dict[6][selected_pos[6]]) << 9 * 4;
            pat |= static_cast<ull>(m_dict[7][selected_pos[7]]) << 8 * 4;
            pat |= static_cast<ull>(m_dict[8][selected_pos[8]]) << 7 * 4;
            pat |= static_cast<ull>(m_dict[9][selected_pos[9]]) << 6 * 4;
            pat |= static_cast<ull>(m_dict[10][selected_pos[10]]) << 5 * 4;
            pat |= static_cast<ull>(m_dict[11][selected_pos[11]]) << 4 * 4;
            pat |= static_cast<ull>(m_dict[12][selected_pos[12]]) << 3 * 4;
            pat |= static_cast<ull>(m_dict[13][selected_pos[13]]) << 2 * 4;
            pat |= static_cast<ull>(m_dict[14][selected_pos[14]]) << 1 * 4;
            pat |= static_cast<ull>(m_dict[15][selected_pos[15]]);
            data_p[count_i++] = pat;
        }
    }
    return data_p;
}

const char map[] = { "0123456789ABCDEF" };
void KGCUDA::locate(int step)
{
    short selected_pos[16]{0};
    memset(selected_pos, 0, 16 * sizeof(short));
    short curr_line = 15;//0-15
    memcpy(selected_pos, curr_task.start_range, sizeof(selected_pos));
    int local_step = -1;
    while (curr_line >= 0)
    {
        if (++selected_pos[curr_line] > Task::m_globe_range[curr_line * 2 + 1])
        {
            selected_pos[curr_line] = Task::m_globe_range[curr_line * 2];
            --curr_line;
        }
        else
        {
            curr_line = 15;
            if (++local_step >= step)
            {
                printf("registration code:");
                std::string pat;
                for (int i = 0; i < 16; i++)
                {
                    pat += map[m_dict[i][selected_pos[i]]];

                }
                mpz_t m;
                mpz_init_set_str(m, pat.c_str(), 16);
                mpz_powm(m, m, d, n);
                char* regcode = mpz_get_str(nullptr, 16, m);
                mpz_clear(m);
                printf("%s\n", regcode);
                delete regcode;
                break;
            }
        }
    }
}

void KGCUDA::compatible()
{
    int dev_n = 0;
    cudaGetDeviceCount(&dev_n);
    if (dev_n <= 0)
        num_blocks = num_threads_preblock = 0;
    else
    {
        num_threads_preblock = m_gpu_prop.maxThreadsPerBlock;
        num_blocks = m_gpu_prop.maxThreadsPerMultiProcessor / num_threads_preblock * m_gpu_prop.multiProcessorCount;
    }

}

void KGCUDA::test_capability()
{
    if (capability() > 0)
        printf("capability test passed\n");
    print_gpu_properties();

}
