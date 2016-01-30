#ifndef TIMEPROBE_H__
#define TIMEPROBE_H__
#include <string>
#include <unordered_map>
#include <chrono>
#include <tuple>
#include <vector>
#include <mutex>
#pragma once
class TimeProbe
{
public:
    enum Status
    {
        started = 0,
        stoped
    };
    TimeProbe();
    ~TimeProbe();
    void tick_start(std::string probe_name);
    void tick_stop(std::string probe_name);
    void print();
private:
    typedef std::tuple <
        std::chrono::duration<long double>,/*total time*/
        std::chrono::steady_clock::time_point/* start time*/,
        long/* call times*/,
        Status
    > probe_data;
    std::unordered_map < std::string,/*name of probe*/
        probe_data>
        m_ticker;

    std::mutex m_mtx;
};

#endif // TimeProbe_h__
