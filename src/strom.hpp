#pragma once

#include <iostream>
#include "tree_summary.hpp"
#include "data.hpp"
#include "likelihood.hpp"
#include <boost/program_options.hpp>

namespace strom {

class Strom
    {
    public:
                            Strom();
                            ~Strom();

        void                clear();
        void                processCommandLineOptions(int argc, const char * argv[]);
        void                run();

    private:

        std::string            _data_file_name;
        std::string            _tree_file_name;
        unsigned               _tree_to_plot;

        Data::SharedPtr        _data;
        Likelihood::SharedPtr   _likelihood;
        TreeSummary::SharedPtr _tree_summary;

        static std::string     _program_name;
        static unsigned        _major_version;
        static unsigned        _minor_version;

    };

inline Strom::Strom()
    {
    //std::cout << "Constructing a Strom" << std::endl;
    clear();
    }

inline Strom::~Strom()
    {
    //std::cout << "Destroying a Strom" << std::endl;
    }

inline void Strom::clear()
    {
    _data_file_name = "";
    _tree_file_name = "";
    _tree_to_plot = 0;
    _data           = nullptr;
    _likelihood     = nullptr;
    _tree_summary   = nullptr;
    }

inline void Strom::processCommandLineOptions(int argc, const char * argv[])
    {
    boost::program_options::variables_map       vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        //("datafile,d",  boost::program_options::value(&_data_file_name)->required(), "name of data file in NEXUS format")
        ("treefile,t",  boost::program_options::value(&_tree_file_name)->required(), "name of tree file in NEXUS format")
        ("plot",  boost::program_options::value(&_tree_to_plot), "number of tree to plot")
        ;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    try
        {
        const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("strom.conf", desc, false);
        boost::program_options::store(parsed, vm);
        }
    catch(boost::program_options::reading_file & x)
        {
        std::cout << "Note: configuration file (strom.conf) not found" << std::endl;
        }
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

    // If user specified --plot on command line, check to make sure it is valid
    if (vm.count("plot") > 0)
        {
        if (_tree_to_plot == 0)
            {
            std::cout << "plot must be greater than zero" << std::endl;
            std::exit(1);
            }
        }
    }

inline void Strom::run()
    {
    std::cout << "Starting..." << std::endl;

    try
        {
        // Read and store data
        //_data = Data::SharedPtr(new Data());
        //_data->getDataFromFile(_data_file_name);

        // Create a Likelihood object that will compute log-likelihoods
        //_likelihood = Likelihood::SharedPtr(new Likelihood());
        //_likelihood->setData(_data);

        // Read in trees
        _tree_summary = TreeSummary::SharedPtr(new TreeSummary());
        _tree_summary->readTreefile(_tree_file_name, 0, _tree_to_plot);
        _tree_summary->showSummary(_tree_to_plot);

        //Tree::SharedPtr tree = _tree_summary->getTree(0);

        // Calculate the log-likelihood for the tree
        //double lnL = _likelihood->calcLogLikelihood(tree);
        //std::cout << boost::str(boost::format("log likelihood = %.5f") % lnL) << std::endl;
        //std::cout << "      (expecting -286.9238)" << std::endl;
        }
    catch (XStrom & x)
        {
        std::cerr << "Strom encountered a problem:\n  " << x.what() << std::endl;
        }

    std::cout << "\nFinished!" << std::endl;
    }

} // namespace strom
