#pragma once
#ifndef RANDOMIZER_H
#define RANDOMIZER_H

#include <random>
#include <chrono>
#include <type_traits>

class Randomizer
{
public:
    Randomizer();

    ~Randomizer();

    template<typename _Set>
    typename std::add_pointer<typename _Set::value_type>::type randomize(_Set set, typename _Set::size_type _times);

    template<typename _Set>
    typename std::add_pointer<typename _Set::value_type>::type unique_randomize(_Set set, typename _Set::size_type _times);

    template<typename _SeedType>
    void set_seed(_SeedType seed);

    void change_seed(int seed);

private:
    void clear();

    void clear_distribution();

    void clear_generator();

private:
    std::default_random_engine *m_generator;
    std::uniform_int_distribution<size_t> *m_distribution;
};


#endif

//note: it just return the specific number of
//random items. if set type is char, then there would be no '\0'
// at the end of returned array!
template<typename _Set>
typename std::add_pointer<typename _Set::value_type>::type Randomizer::randomize(_Set _set, typename _Set::size_type _times)
{
    clear_distribution();
    m_distribution = new std::uniform_int_distribution<size_t>(0, _set.size() - 1);
    typename std::add_pointer<typename _Set::value_type>::type rand_array = new typename _Set::value_type[_times];
    for (int i = 0; i < _times; ++i)
    {
        rand_array[i] = _set[(*m_distribution)(*m_generator)];
    }
    return rand_array;
}

template<typename _Set>
typename std::add_pointer<typename _Set::value_type>::type Randomizer::unique_randomize(_Set _set, typename _Set::size_type _times)
{
    typename std::add_pointer<typename _Set::value_type>::type rand_array = new typename _Set::value_type[_times];
    for (int i = 0; i < _times; ++i)
    {
        clear_distribution();
        m_distribution = new std::uniform_int_distribution<size_t>(0, _set.size() - 1);
        auto rand_index = (*m_distribution)(*m_generator);
        rand_array[i] = _set[rand_index];
        _set.erase(_set.begin()+rand_index);
    }
    return rand_array;
}


template<typename _SeedType>
void Randomizer::set_seed(_SeedType _seed)
{
    clear_generator();
    m_generator = new std::default_random_engine(static_cast<unsigned int>(_seed));
}
