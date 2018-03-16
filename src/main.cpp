#if 0

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

    std::cout << "\nMean of 10000 Gamma(2,1) deviates is expected to be 2.0:" << std::endl;
    double total = 0.0;
    for (unsigned i = 0; i < 10000; ++i)
        total += lot.gamma(2.0,1.0);
    std::cout << "  sample mean = " << (total/10000.0) << std::endl;

    return 0;
    }


#else

#include <iostream>
#include "strom.hpp"

using namespace strom;

// static data member initializations
std::string Strom::_program_name    = "strom";
unsigned    Strom::_major_version   = 1;
unsigned    Strom::_minor_version   = 0;
const double Node::_smallest_edge_length  = 1.0e-12; //POLNEW

int main(int argc, const char * argv[])
    {
    Strom strom;
    try {
        strom.processCommandLineOptions(argc, argv);
        strom.run();
    }
    catch(std::exception & x) {
        std::cerr << "Exception: " << x.what() << std::endl;
        std::cerr << "Aborted." << std::endl;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }

    return 0;
    }


#endif
