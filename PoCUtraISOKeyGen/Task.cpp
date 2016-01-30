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
        0, 15,
        0, 15,
        0, 15,
        0, 15,
        0, 15,
        0, 15 };
    for (int i = 0; i < 32; i += 1)
    {
        m_globe_range[i] = range[i];
    }
    
}
