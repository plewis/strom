#include <iostream>
#include "tree_summary.hpp"

using namespace strom;

int main(int argc, const char * argv[])
    {
    std::cout << "Starting..." << std::endl;

    TreeSummary sumt;
    try
        {
        sumt.readTreefile("trees.t", 1);
        }
	catch(NxsException x)
		{
        std::cerr << "Program aborting due to errors encountered reading tree file." << std::endl;
        std::cerr << x.what() << std::endl;
        std::exit(0);
		}
    sumt.showSummary();

    std::cout << "\nFinished!" << std::endl;

    return 0;
    }
