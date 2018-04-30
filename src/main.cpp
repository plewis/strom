#include <iostream>
#include "lot.hpp"

using namespace strom;

int main(int argc, const char * argv[])
    {
    Lot lot;
    lot.setSeed(12345);
    std::cout << "Seed set to 12345" << std::endl;
    std::cout << "  Uniform(0,1) random deviate:             " << lot.uniform() << std::endl;
    std::cout << "  Uniform(0,1) random deviate (log scale): " << lot.logUniform() << std::endl;
    std::cout << "  Discrete Uniform(1,4) random deviate:    " << lot.randint(1,4) << std::endl;
    std::cout << "  Normal(0,1) random deviate:              " << lot.normal() << std::endl;
    std::cout << "  Gamma(2,1) random deviate:               " << lot.gamma(2.0,1.0) << std::endl;
    lot.setSeed(12345);
    std::cout << "\nSeed set to 12345" << std::endl;
    std::cout << "  Uniform(0,1) random deviate:             " << lot.uniform() << std::endl;
    std::cout << "  Uniform(0,1) random deviate (log scale): " << lot.logUniform() << std::endl;
    std::cout << "  Discrete Uniform(1,4) random deviate:    " << lot.randint(1,4) << std::endl;
    std::cout << "  Normal(0,1) random deviate:              " << lot.normal() << std::endl;
    std::cout << "  Gamma(2,1) random deviate:               " << lot.gamma(2.0,1.0) << std::endl;

    return 0;
    }
