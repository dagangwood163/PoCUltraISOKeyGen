#pragma once
#ifndef KGCUDA_h__
#define KGCUDA_h__
#include "cuda_runtime.h"
#include <string>
#include "Task.h"
class KGCUDA
{
public:
    KGCUDA();
    ~KGCUDA();
    void dispatch(Task* _task);
    void set_dict(char** dict){ m_dict = dict; }
    void set_name(std::string name);
    unsigned int capability() const
    {
        return num_blocks*num_threads_preblock;
    }
    void print_gpu_properties();

    void test();
    void generate_test();
    //note: the md5 array length is num_blocks*num_threads_preblock.
    void find(unsigned long long* md5, int len);

    unsigned long long* generate(Task* _task);
    void locate(int step);

    //calculate kernel blocks and threads.
    void compatible();

    void init_globe_range();


    //test functions
    void test_capability();

    //char(of each line of dict) position from start to end
    short m_globe_range[32];
    cudaDeviceProp m_gpu_prop;
    unsigned int num_threads_preblock;
    unsigned int num_blocks;
    char** m_dict;
    std::string m_name;
    Task curr_task;
};
#endif // KGCUDA_h__

