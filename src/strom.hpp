#pragma once

#include <iostream>
#include "tree_summary.hpp"
#include "likelihood.hpp"
#include <boost/program_options.hpp>

namespace strom {

class Strom
    {
    public:
                            Strom();
                            ~Strom();

        void                processCommandLineOptions(int argc, const char * argv[]);
        void                run();

    private:

        std::string         _data_file_name;
        std::string         _tree_file_name;

        static std::string  _program_name;
        static unsigned     _major_version;
        static unsigned     _minor_version;

    };

inline Strom::Strom()
    {
    std::cout << "Constructing a Strom" << std::endl;
    }

inline Strom::~Strom()
    {
    std::cout << "Destroying a Strom" << std::endl;
    }

inline void Strom::processCommandLineOptions(int argc, const char * argv[])
    {
    boost::program_options::variables_map       vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile,d",  boost::program_options::value(&_data_file_name)->required(), "name of data file in NEXUS format")
        ("treefile,t",  boost::program_options::value(&_tree_file_name)->required(), "name of data file in NEXUS format")
        ;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    // If user specified --help on command line, output usage summary and quit
    if (vm.count("help") > 0)
        {
        std::cout << desc << "\n";
        std::exit(1);
        }

    // If user specified --version on command line, output version and quit
    if (vm.count("version") > 0)
        {
        std::cout << boost::str(boost::format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << std::endl;
        std::exit(1);
        }
    }

inline void Strom::run()
    {
    std::cout << "Starting..." << std::endl;

    try
        {
        // Read and store data
        Data::SharedPtr d(new Data());
        d->getDataFromFile(_data_file_name);

        // Create a likelihood object that will compute log-likelihoods
        Likelihood::SharedPtr likelihood(new Likelihood());
        likelihood->setData(d);

        // Read in a tree
        TreeSummary::SharedPtr tree_summary(new TreeSummary());
        tree_summary->readTreefile(_tree_file_name, 0);
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
    }

}

