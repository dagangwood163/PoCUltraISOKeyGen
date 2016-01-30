#include "KGCPU.h"
#include "constdef.h"
#include "third_party\gmp-6.0.0\gmp.h"
#include "Task.h"
#include "Helper.h"
#include <thread>
#include <iostream>
#include <vector>
#include <cstdio>
extern int get_CPU_core_num();
void runKGkernel(unsigned long long* data, int blocks, int threads);
static char modulus_n[] = "A70F8F62AA6E97D1";
static char d_str[] = "A7CAD9177AE95A9";
static mpz_t d, n;

KGCPU::KGCPU()
    :m_finished(false),
    capability(1000000)
{
    produce_dict();
    m_gpu.set_dict(m_gpu_dict);

    m_name.resize(64);
    printf("enter user name:");
    std::cin.getline(const_cast<char*>(m_name.c_str()), m_name.size());
    m_name = "UTRISO" + m_name;
    m_name.resize(strlen(m_name.c_str()));
    m_gpu.set_name(m_name);

    memset(m_globe_selected_pos, 0, 16 * sizeof(short));
    m_globe_selected_pos[15] = -1; // it won't miss the first pat.

    mpz_init_set_str(n, modulus_n, 16);
    mpz_init_set_str(d, d_str, 16);
}

KGCPU::~KGCPU()
{
    mpz_clears(n, d, NULL);
}

void KGCPU::drive()
{
    m_time_probe.tick_start("keygen");
    std::thread td(&KGCPU::drive_gpu, this);

    /*int core_num = get_CPU_core_num();
    core_num /= 2;
    if (core_num == 0)
    {
    core_num = 1;
    }
    std::thread** cpu_threads = new std::thread *[core_num];
    while (!m_finished)
    {
    for (int i = core_num - 1; i >= 0; i--)
    {
    cpu_threads[i] =
    new std::thread(&KGCPU::drive_cpu, this, request_task(capability));
    }

    for (int i = core_num - 1; i >= 0; i--)
    {
    if (cpu_threads[i]->joinable())
    cpu_threads[i]->join();
    delete cpu_threads[i];
    }
    }*/
    if (td.joinable())
    {
        td.join();
    }
    m_time_probe.tick_stop("keygen");
    m_time_probe.print();
}

void KGCPU::test()
{
    //product_dict_test();
    // request_task_test();
    // drive_cpu_test();
    //drive_gpu_test();
    drive_test();

}

void KGCPU::product_dict_test()
{
    char cpu_char_temp[2];
    short range_i[] = { 2, 2, 15, 15, 15, 15, 2, 3, 1, 1, 15, 15, 15, 15, 15, 15 };

    printf("generated dictionary:\n");
    for (auto i = 0; i < 16; i++)
    {
        for (int k = 0; k < range_i[i]; k++)
        {
            printf("%c", m_cpu_dict[i][k]);
        }
        printf("\n");
    }

    std::string unique_str;

    // check whether m_cpu_dict ==  m_gpu_dict
    for (auto i = 0; i < 16; i++)
    {
        unique_str.clear();
        for (int j = range_i[i] - 1; j >= 0; j--)
        {
            sprintf_s(cpu_char_temp, "%1X", m_gpu_dict[i][j]);
            if (m_cpu_dict[i][j] != cpu_char_temp[0])
            {
                printf("dictionary generation failed!\n");
                for (int k = 0; k < range_i[i] - 1; k++)
                {
                    printf("%c", m_cpu_dict[i][k]);
                }
                printf("\n");

                for (int k = 0; k < range_i[i] - 1; k++)
                {
                    printf("%01X", m_gpu_dict[i][k]);
                }
                printf("\n");
                return;
            }
            unique_str += m_cpu_dict[i][j];
        }

        // check duplicate keys.
        for (int j = unique_str.length() - 1; j >= 0; j--)
        {
            if (unique_str.find(unique_str[j]) != unique_str.rfind(unique_str[j]))
            {
                printf("duplicate char found!\n");
                printf("dictionary generation failed!\n");
                return;
            }
        }
    }
    printf("dictionary generation test passed.\n");
}

void KGCPU::request_task_test()
{
    for (int i = 1; i < 999999; i++)
    {
        auto range = request_task(i);
        int r = 0;
        auto comsume = [&r, i, range, this]() -> bool
        {
            short selected_pos[16]{0};
            memset(selected_pos, 0, 16 * sizeof(short));
            short curr_line = 15;
            selected_pos[curr_line] = range->start_range[curr_line];
            while (curr_line >= 0)
            {
                if (++selected_pos[curr_line] >
                    Task::m_globe_range[curr_line * 2 + 1])
                {
                    selected_pos[curr_line] = Task::m_globe_range[curr_line * 2];
                    --curr_line;
                }
                else
                {
                    r += 1;
                    if (r == i)
                        break;
                }
            }
            return memcmp(m_globe_selected_pos, selected_pos,
                sizeof(m_globe_selected_pos)) == 0;
        };

        if (!comsume())
        {
            printf("request task test failed with request:%d with count:%d\n", i, r);
            delete range;
            return;
        }
        delete range;
    }
    printf("task generation test passed.\n");
}

void KGCPU::drive_cpu_test()
{
    char* _dict[]{
        "45",
            "56",
            "0123456789ABCDEF",
            "0123456789ABCDEF",
            "0123456789ABCDEF",
            "0123456789ABCDEF",
            "2D",
            "ABC",
            "5",
            "3",
            "2",
            "0",
            "2",
            "0",
            "4",
            "DEF"};
    char* dict[]{
        "5",
            "5",
            "4",
            "8",
            "6",
            "D",
            "2",
            "C",
            "5",
            "3",
            "2",
            "0",
            "2",
            "0",
            "4",
            "F"};
    short range[32] = {
        0, 1,
        0, 1,
        0, 15,
        0, 15,
        0, 15,
        0, 15,
        0, 1,
        0, 2,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 2 };
    auto exhange = [&](int n, int m)
    {
        for (int i = n; i < m + 1; i++)
        {
            m_cpu_dict[i] = _dict[i];
        }
    };
    auto fixrange = [&]()
    {
        for (int i = 0; i < 32; i += 1)
        {
            Task::m_globe_range[i] = range[i];
        }
    };

    exhange(0, 15);
    fixrange();

    Task* task = request_task(1000000);
    while (task->count > 0)
    {
        drive_cpu(task);
        delete task;
        task = request_task(1000000);
    }
}

void KGCPU::drive_gpu_test()
{
    char* _dict[]{
        "\x4\x5",
            "\x5\x6",
            "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x2\xD",
            "\xA\xB\xC",
            "\x5",
            "\x3",
            "\x2",
            "\x0",
            "\x2",
            "\x0",
            "\x4",
            "\xD\xE\xF"};
    char* dict[]{
        "\x4\x6\x5",
            "\x5",
            "\x4",
            "\x8",
            "\x6",
            "\xD",
            "\x2",
            "\xC",
            "\x5",
            "\x3",
            "\x2",
            "\x0",
            "\x2",
            "\x0",
            "\x4",
            "\xF"};
    short _range[32] = {
        0, 1,
        0, 1,
        0, 15,
        0, 15,
        0, 15,
        0, 15,
        0, 1,
        0, 2,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 2 };
    short range[32] = {
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0 };
    auto exhange = [&](int n, int m)
    {
        for (int i = n; i < m + 1; i++)
        {
            m_gpu_dict[i] = _dict[i];
        }
    };
    auto fixrange = [&]()
    {
        for (int i = 0; i < 32; i += 1)
        {
            Task::m_globe_range[i] = _range[i];
        }
    };

    exhange(0, 15);
    fixrange();
    Task* task = request_task(m_gpu.capability());
    while (task->count > 0)
    {
        Helper::factor factor;
        if (task->count == m_gpu.capability())
        {
            factor = { m_gpu.num_blocks, m_gpu.num_threads_preblock };
        }
        else
        {
            factor = Helper::factoring(task->count, m_gpu.num_threads_preblock);
        }
        auto data = m_gpu.generate(task);
        runKGkernel(data, factor.a, factor.b);
        m_gpu.find(data, factor.a*factor.b);
        delete data;
        delete task;
        task = request_task(m_gpu.capability());
    }

}

void KGCPU::drive_test()
{
    char* _gpu_dict[]{
        "\x4\x5",
            "\x5\x6",
            "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x0\x1\x2\x3\x4\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x2\xD",
            "\xA\xB\xC",
            "\x5",
            "\x3",
            "\x6",
            "\x5",
            "\x4",
            "\xF",
            "\x4",
            "\xF"};
    char* gpu_dict[]{
        "\x5",
            "\x5",
            "\x5",
            "\x3",
            "\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x5\x6\x7\x8\x9\xA\xB\xC\xD\xE\xF",
            "\x2\xD",
            "\xC",
            "\x5",
            "\x3",
            "\x6",
            "\x5",
            "\x4",
            "\xF",
            "\x4",
            "\xF"};
    char* cpu_dict[]{
        "45",
            "56",
            "0123456789ABCDEF",
            "0123456789ABCDEF",
            "0123456789ABCDEF",
            "0123456789ABCDEF",
            "2D",
            "ABC",
            "5",
            "3",
            "6",
            "5",
            "4",
            "F",
            "4",
            "F"};
    short _range[32] = {
        0, 1,
        0, 1,
        0, 15,
        0, 15,
        0, 15,
        0, 15,
        0, 1,
        0, 2,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0 };
    short range[32] = {
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0,
        0, 0 };
    auto exhange = [&](int n, int m)
    {
        for (int i = n; i < m + 1; i++)
        {
            m_gpu_dict[i] = gpu_dict[i];                                                      
            m_cpu_dict[i] = cpu_dict[i];
        }
    };
    auto fixrange = [&]()
    {
        for (int i = 0; i < 32; i += 1)
        {
            Task::m_globe_range[i] = range[i];
        }
    };
    exhange(0, 15);
    fixrange();
    drive();
}

Task* KGCPU::request_task(size_t num)
{
    Task* range_result_p = new Task();
    m_mtx.lock();
    short selected_pos[16];
    // we start from the position after the begin of globe position as
    // it used in previous production.
    memcpy(selected_pos, m_globe_selected_pos, sizeof(selected_pos));
    short curr_line = 15; // 0-15

    unsigned count = 0;
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
            if (++count >= num)
            {
                break;
            }
        }
    }
    if (count >= num || curr_line < 0)
    {
        for (int i = 0; i < 16; i += 1)
        {
            range_result_p->start_range[i] = m_globe_selected_pos[i];
        }
        range_result_p->count = count;
        memcpy(m_globe_selected_pos, selected_pos, 16 * sizeof(short));
    }
    if (m_finished)
    {
        range_result_p->count = 0;
    }
    if (curr_line < 0)
    {
        m_finished = true;
    }
    m_mtx.unlock();
    return range_result_p;
}

void KGCPU::produce_sqnce()
{
    std::vector<short> s;
    m_rand.change_seed(0);
    s.assign({ 0, 1 });
    m_sqnce[0] = m_rand.unique_randomize(s, s.size());

    m_rand.change_seed(1);
    m_sqnce[1] = m_rand.unique_randomize(s, s.size());

    m_rand.change_seed(6);
    m_sqnce[6] = m_rand.unique_randomize(s, s.size());

    m_rand.change_seed(7);
    s.assign({ 0, 1, 2 });
    m_sqnce[7] = m_rand.unique_randomize(s, s.size());

    s.assign({ 0 });
    m_sqnce[8] = m_rand.unique_randomize(s, s.size());
    m_sqnce[9] = m_rand.unique_randomize(s, s.size());

    s.assign({ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 });
    for (int i = 2; i < 6; i++)
    {
        m_rand.change_seed(i);
        m_sqnce[i] = m_rand.unique_randomize(s, s.size());
    }

    for (int i = 10; i < 16; i++)
    {
        m_rand.change_seed(i);
        m_sqnce[i] = m_rand.unique_randomize(s, s.size());
    }
}

void KGCPU::produce_cpu_dict()
{
    std::string s;
    s = "45";
    m_cpu_dict[0] = new char[2];
    for (int i = 0; i < 2; i++)
    {
        m_cpu_dict[0][i] = s[m_sqnce[0][i]];
    }

    s = "56";
    m_cpu_dict[1] = new char[2];
    for (int i = 0; i < 2; i++)
    {
        m_cpu_dict[1][i] = s[m_sqnce[1][i]];
    }

    s = "2D";
    m_cpu_dict[6] = new char[2];
    for (int i = 0; i < 2; i++)
    {
        m_cpu_dict[6][i] = s[m_sqnce[6][i]];
    }

    s = "ABC";
    m_cpu_dict[7] = new char[2];
    for (int i = 0; i < 3; i++)
    {
        m_cpu_dict[7][i] = s[m_sqnce[7][i]];
    }

    m_cpu_dict[8] = new char[1];
    m_cpu_dict[9] = new char[1];
    m_cpu_dict[8][0] = '5';
    m_cpu_dict[9][0] = '3';

    s = "0123456789ABCDEF";
    for (int i = 2; i < 6; i++)
    {
        m_cpu_dict[i] = new char[16];
        for (int j = 0; j < 16; j++)
        {
            m_cpu_dict[i][j] = s[m_sqnce[i][j]];
        }
    }

    for (int i = 10; i < 16; i++)
    {
        m_cpu_dict[i] = new char[16];
        for (int j = 0; j < 16; j++)
        {
            m_cpu_dict[i][j] = s[m_sqnce[i][j]];
        }
    }
}

void KGCPU::produce_gpu_dict()
{
    std::vector<char> s;
    s = { '\x4', '\x5' };
    m_gpu_dict[0] = new char[2];
    for (int i = 0; i < 2; i++)
    {
        m_gpu_dict[0][i] = s[m_sqnce[0][i]];
    }

    s = { '\x5', '\x6' };
    m_gpu_dict[1] = new char[2];
    for (int i = 0; i < 2; i++)
    {
        m_gpu_dict[1][i] = s[m_sqnce[1][i]];
    }

    s = { '\x2', '\xD' };
    m_gpu_dict[6] = new char[2];
    for (int i = 0; i < 2; i++)
    {
        m_gpu_dict[6][i] = s[m_sqnce[6][i]];
    }

    s = { '\xA', '\xB', '\xC' };
    m_gpu_dict[7] = new char[2];
    for (int i = 0; i < 3; i++)
    {
        m_gpu_dict[7][i] = s[m_sqnce[7][i]];
    }

    m_gpu_dict[8] = new char[1];
    m_gpu_dict[9] = new char[1];
    m_gpu_dict[8][0] = 5;
    m_gpu_dict[9][0] = 3;

    s = { '\x0', '\x1', '\x2', '\x3', '\x4', '\x5', '\x6', '\x7',
        '\x8', '\x9', '\xA', '\xB', '\xC', '\xD', '\xE', '\xF' };
    for (int i = 2; i < 6; i++)
    {
        m_gpu_dict[i] = new char[16];
        for (int j = 0; j < 16; j++)
        {
            m_gpu_dict[i][j] = s[m_sqnce[i][j]];
        }
    }

    for (int i = 10; i < 16; i++)
    {
        m_gpu_dict[i] = new char[16];
        for (int j = 0; j < 16; j++)
        {
            m_gpu_dict[i][j] = s[m_sqnce[i][j]];
        }
    }
}

void KGCPU::produce_dict()
{
    produce_sqnce();
    produce_cpu_dict();
    produce_gpu_dict();
}

void KGCPU::drive_gpu()
{
    m_gpu.print_gpu_properties();
    unsigned int gpu_capability = m_gpu.capability();
    Task* task = nullptr;
    while (task = request_task(gpu_capability), task->count > 0)
    {
        m_gpu.dispatch(task);
        delete task;
    }
    delete task;
}

typedef unsigned long long ull;
typedef unsigned int uint;
typedef unsigned char uchar;
ull md5_512(uint* msg);

void KGCPU::drive_cpu(Task* range)
{
    mpz_t m;
    char pat[17]{0};

    short* selected_pos = nullptr;
    short curr_line = 15; // 0-15
    // there shouldn't minus 1 as we do not want to roll selected position with 1.
    // thus we start after (globe) selected position.
    selected_pos = range->start_range;

    uchar msg[64]{0};
    memset(msg, 0, 64);

    auto len = m_name.length();
    // I believe len won't beyond unsigned int type representation range.
    *reinterpret_cast<uint *>(&msg[56]) = (len + 16) * 8;
    msg[len + 16] = 0x80;
    memcpy(msg, m_name.c_str(), len);

    ull md5;

    // generate strings start
    while (curr_line >= 0 && range->count > 0)
    {
        if (++selected_pos[curr_line] > Task::m_globe_range[curr_line * 2 + 1])
        {
            selected_pos[curr_line] = Task::m_globe_range[curr_line * 2];
            --curr_line;
        }
        else
        {
            range->count -= 1;
            curr_line = 15;
            pat[0] = m_cpu_dict[0][selected_pos[0]];
            pat[1] = m_cpu_dict[1][selected_pos[1]];
            pat[2] = m_cpu_dict[2][selected_pos[2]];
            pat[3] = m_cpu_dict[3][selected_pos[3]];
            pat[4] = m_cpu_dict[4][selected_pos[4]];
            pat[5] = m_cpu_dict[5][selected_pos[5]];
            pat[6] = m_cpu_dict[6][selected_pos[6]];
            pat[7] = m_cpu_dict[7][selected_pos[7]];
            pat[8] = m_cpu_dict[8][selected_pos[8]];
            pat[9] = m_cpu_dict[9][selected_pos[9]];
            pat[10] = m_cpu_dict[10][selected_pos[10]];
            pat[11] = m_cpu_dict[11][selected_pos[11]];
            pat[12] = m_cpu_dict[12][selected_pos[12]];
            pat[13] = m_cpu_dict[13][selected_pos[13]];
            pat[14] = m_cpu_dict[14][selected_pos[14]];
            pat[15] = m_cpu_dict[15][selected_pos[15]];
            mpz_init_set_str(m, pat, 16);
            mpz_powm(m, m, d, n);
            mpz_get_str(reinterpret_cast<char *>(&msg[len]), 16, m);
            mpz_clear(m);

            msg[len + 16] = 0x80; // md5 padding

            md5 = md5_512(reinterpret_cast<uint *>(msg));

            // finder start
            int len_data = 52413;
            int i = 0;
            int index;
            do
            {
                index = (len_data + i) / 2;
                int result =
                    memcmp(&data[6 * index], reinterpret_cast<uchar *>(&md5), 6);
                if (result <= 0)
                {
                    if (result >= 0)
                    {
                        // success
                        printf("registration code:");
                        printf("%s\n", reinterpret_cast<char *>(&msg[len]));
                    }
                    i = index + 1;
                }
                else
                {
                    len_data = index - 1;
                }
            } while (i <= len_data);
            // finder end
        }
    }
    // generate strings end
}
