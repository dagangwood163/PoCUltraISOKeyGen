#include "TimeProbe.h"
#include <iostream>

TimeProbe::TimeProbe()
{
}


TimeProbe::~TimeProbe()
{
}

void TimeProbe::tick_start(std::string probe_name)
{
    m_mtx.lock();
    if (m_ticker.find(probe_name) == m_ticker.end())
    {
        std::chrono::steady_clock::time_point t = std::chrono::steady_clock::now();
        m_ticker[probe_name] = probe_data(std::chrono::duration<long double>(0), t, 0, Status::started);
    }
    m_mtx.unlock();

}

void TimeProbe::tick_stop(std::string probe_name)
{
    m_mtx.lock();
    std::chrono::system_clock::time_point t = std::chrono::system_clock::now();
    if (m_ticker.find(probe_name) != m_ticker.end())
    {
        if (std::get<3>(m_ticker[probe_name]) == Status::started)
        {
            std::get<3>(m_ticker[probe_name]) = Status::stoped;

        }
        std::get<0>(m_ticker[probe_name]) += t - std::get<1>(m_ticker[probe_name]);
        std::get<2>(m_ticker[probe_name]) += 1;
    }
    m_mtx.unlock();

}

void TimeProbe::print()
{
    m_mtx.lock();
    for (auto probe : m_ticker)
    {
        std::cout << probe.first << " COSTS:";
        std::cout << std::get<0>(probe.second).count() << " SECONDS ";
        std::cout << std::endl;
    }
    m_mtx.unlock();
}
