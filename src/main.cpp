#include <iostream>
#include "node.hpp"
#include "tree.hpp"
#include "tree_manip.hpp"

using namespace strom;

int main(int argc, const char * argv[])
    {
    std::cout << "Starting..." << std::endl;
    TreeManip tm;
    tm.createTestTree();
    std::cout << "\nFinished!" << std::endl;

    return 0;
    }
