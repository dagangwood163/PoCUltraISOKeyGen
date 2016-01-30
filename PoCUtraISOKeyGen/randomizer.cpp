#include "randomizer.h"
Randomizer::Randomizer()
{
    m_generator = nullptr;
    m_distribution = nullptr;
}

Randomizer::~Randomizer()
{
    clear();

}



void Randomizer::change_seed(int seed)
{
    set_seed(std::chrono::system_clock::now().time_since_epoch().count()+seed);
}

void Randomizer::clear()
{
    clear_distribution();
    clear_generator();
}


void Randomizer::clear_distribution()
{
    if (m_distribution != nullptr)
        delete m_distribution;
}

void Randomizer::clear_generator()
{
    if (m_generator != nullptr)
        delete m_generator;
}


