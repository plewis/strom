#include <iostream>
//#include "node.hpp"
//#include "tree.hpp"
#include "data.hpp"
#include "tree_summary.hpp"
#include "likelihood.hpp"

using namespace strom;

int main(int argc, const char * argv[])
    {
    std::cout << "Starting..." << std::endl;

    try
        {
        // Read and store data
        Data::SharedPtr d(new Data());
        try
            {
            d->getDataFromFile("rbcL.nex");
            }
        catch(NxsException x)
            {
            std::cerr << "Program aborting due to errors encountered reading data file." << std::endl;
            std::cerr << x.what() << std::endl;
            std::exit(0);
            }

        // Create a likelihood object that will compute log-likelihoods
        Likelihood::SharedPtr likelihood(new Likelihood());
        likelihood->setData(d);

        // Read in a tree
        TreeSummary::SharedPtr tree_summary(new TreeSummary());
        try
            {
            tree_summary->readTreefile("rbcLjc.tre", 0);
            }
        catch(NxsException x)
            {
            std::cerr << "Program aborting due to errors encountered reading tree file." << std::endl;
            std::cerr << x.what() << std::endl;
            std::exit(0);
            }
        Tree::SharedPtr tree = tree_summary->getTree(0);

        // Calculate the log-likelihood for the tree
        double lnL = likelihood->calcLogLikelihood(tree);
        std::cout << boost::str(boost::format("log likelihood = %.5f") % lnL) << std::endl;
        std::cout << "      (expecting -286.9238)" << std::endl;
        }
    catch (XStrom x)
        {
        std::cerr << "Strom encountered a problem:\n  " << x.what() << std::endl;
        }

    std::cout << "\nFinished!" << std::endl;

    return 0;
    }
