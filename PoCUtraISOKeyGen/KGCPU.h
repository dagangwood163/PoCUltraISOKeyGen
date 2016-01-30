#pragma once
#ifndef KGCPU_h__
#define KGCPU_h__
#include "KGCUDA.h"
#include "randomizer.h"
#include "TimeProbe.h"
#include <string>
#include <mutex>
class Task;
class KGCPU
{
public:
    KGCPU();
    ~KGCPU();
    void drive();

    void test();
private:
    void product_dict_test();
    void request_task_test();
    void drive_cpu_test();
    void drive_gpu_test();
    void drive_test();

private:

    Task* request_task(size_t num);
    void produce_sqnce();
    void produce_cpu_dict();
    void produce_gpu_dict();
    //generate both cpu and gpu dictionary.they must look like the same one.
    void produce_dict();

    void drive_gpu();
    //high cpu occupancy rate function.
    void drive_cpu(Task* range);

private:
    char* m_cpu_dict[16];
    char* m_gpu_dict[16];
    short* m_sqnce[16];
    short m_globe_selected_pos[16];

    std::string m_name;
    KGCUDA m_gpu;
    Randomizer m_rand;
    bool m_finished;

    std::mutex m_mtx;

    int capability;

    TimeProbe m_time_probe;
};
#endif // KGCPU_h__

