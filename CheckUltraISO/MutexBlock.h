#ifndef MutexBlock_h__
#define MutexBlock_h__

#pragma once
#include <array>
#include <type_traits>
#include <mutex>
#include <thread>
#include <chrono>    
#include <atomic>
template<typename _Type, size_t piec_len, size_t count>
class MutexBlock
{
public:
    typedef size_t pos_t;
    typedef typename std::add_const<_Type>::type* block_pointer;
    enum class Status
    {
        free = 0,
        locked = 1,
        inuse
    };
    struct Block
    {
        _Type data[piec_len];
        std::atomic<Status> m_stat;
    };
    MutexBlock();
    ~MutexBlock();
    __forceinline typename std::add_const<_Type>::type* acquire(Status stat = Status::free);

    //will be blocked when there is no free data block.
    __forceinline typename std::add_const<_Type>::type* block_acquire();

    __forceinline void release(typename std::add_const<_Type>::type* p_data);

    void change_stat(typename std::add_const<_Type>::type* p_data,Status stat);
private:
    std::array<Block, count> m_data;
    pos_t m_curr_pos;

    //thread safety
    std::mutex m_mtx;
public:
    const size_t m_block_count = count;
    const size_t m_block_size = count*piec_len*sizeof(_Type);
};

template<typename _Type, size_t piec_len, size_t count>
void MutexBlock<_Type, piec_len, count>::change_stat(typename std::add_const<_Type>::type* p_data,Status stat)
{
    m_mtx.lock();
    if (p_data)
    {
        pos_t pos = (p_data - reinterpret_cast<_Type*>(m_data.data())) / sizeof(Block);
        if (pos < count)
        {
            m_data[pos].m_stat.store(stat);
        }
    }
    m_mtx.unlock();
}

template<typename _Type, size_t piec_len, size_t count>
__forceinline typename std::add_const<_Type>::type* MutexBlock<_Type, piec_len, count>::block_acquire()
{
    typename std::add_const<_Type>::type* temp = nullptr;
    while (!(temp = acquire()))
    {

        std::this_thread::sleep_for(std::chrono::milliseconds(300));
    }
    return temp;
}

template<typename _Type, size_t piec_len, size_t count>
__forceinline void MutexBlock<_Type, piec_len, count>::release(typename std::add_const<_Type>::type* p_data)
{
    if (p_data)
    {
        pos_t pos = (p_data - reinterpret_cast<_Type*>(m_data.data())) / sizeof(Block);
        if (pos < count)
        {
            m_data[pos].m_stat.store(Status::free);
        }
    }

}

template<typename _Type, size_t piec_len, size_t count>
__forceinline typename std::add_const<_Type>::type*  MutexBlock<_Type, piec_len, count>::acquire(Status stat)
{
    m_mtx.lock();
    typename std::add_const<_Type>::type* p_temp = nullptr;
    if (m_curr_pos >= count)
    {
        m_curr_pos = 0;
    }
    while (m_curr_pos < count)
    {
        if (m_data[m_curr_pos++].m_stat.load() == stat)
        {
            p_temp = m_data[m_curr_pos - 1].data;
            m_data[m_curr_pos - 1].m_stat.store(Status::locked);
            break;
        }
    };
    m_mtx.unlock();
    return p_temp;
}

template<typename _Type, size_t piec_len, size_t count>
MutexBlock<_Type, piec_len, count>::~MutexBlock()
{

}

template<typename _Type, size_t piec_len, size_t count>
MutexBlock<_Type, piec_len, count>::MutexBlock()
{
    memset(m_data.data(), 0, sizeof(Block)*count);
    m_curr_pos = 0;
}

#endif // MutexBlock_h__

