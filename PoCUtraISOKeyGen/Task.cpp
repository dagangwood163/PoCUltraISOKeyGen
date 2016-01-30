#include "Task.h"
#include <cstring>

short Task::m_globe_range[32];

Task::Task()
    :
    count(0)
{
    memset(start_range, 0, sizeof(start_range));
}


Task::~Task()
{
}

void Task::init_globe_range()
{
    for (int i = 0; i < 32; i += 2)
    {
        m_globe_range[i] = 0;
    }
    m_globe_range[1] = 2;
    m_globe_range[3] = 2;
    m_globe_range[13] = 2;
    m_globe_range[15] = 3;
    m_globe_range[17] = 1;
    m_globe_range[19] = 1;

    for (int i = 2; i < 6; i++)
    {
        m_globe_range[i * 2 + 1] = 15;
    }

    for (int i = 10; i < 16; i++)
    {
        m_globe_range[i * 2 + 1] = 15;
    }
}
