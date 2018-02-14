#include <iostream>
#include "data.hpp"
#include "node.hpp"
#include "tree.hpp"
#include "tree_manip.hpp"

using namespace strom;

int main(int argc, const char * argv[])
    {
    std::cout << "Starting..." << std::endl;

    // Test Tree, Node, and TreeManip classes
    TreeManip tm;
    std::string newick = std::string("(1:0.3,2:0.3,(3:0.2,(4:0.1,5:0.1):0.1):0.1)");
    tm.buildFromNewick(newick, false, false);
    std::cout << tm.makeNewick(3) << std::endl;

    std::cout << std::endl;

    // Test Data class
    Data d;
    d.getDataFromFile("rbcL.nex");
    std::cout << "Number of taxa:     " << d.getNumTaxa() << std::endl;
    std::cout << "Sequence length:    " << d.getSeqLen() << std::endl;
    std::cout << "Number of patterns: " << d.getNumPatterns() << std::endl;

    std::cout << "\nFinished!" << std::endl;

    return 0;
    }
