#pragma once
class Task
{
public:
    Task();
    ~Task();
private:
    static void init_globe_range();
public:
    size_t count;
    short start_range[16];
    //char(of each line of dict) position from start to end
    static short m_globe_range[32];
};

