#include <iostream>
#include "data.hpp"
#include "tree_summary.hpp"
#include "likelihood.hpp"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char * argv[])
    {
    std::cout << "Starting..." << std::endl;

    try
        {
        // Read and store data
        Data::SharedPtr d(new Data());
        d->getDataFromFile("rbcL.nex");

        // Create a likelihood object that will compute log-likelihoods
        Likelihood::SharedPtr likelihood(new Likelihood());
        likelihood->setData(d);

        // Read in a tree
        TreeSummary::SharedPtr tree_summary(new TreeSummary());
        tree_summary->readTreefile("rbcLjc.tre", 0);
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
