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

        double              _expected_log_likelihood;
        double              _gamma_shape;
        unsigned            _num_categ;
        std::vector<double> _state_frequencies;
        std::vector<double> _exchangeabilities;

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
        ("expectedLnL", boost::program_options::value(&_expected_log_likelihood)->default_value(0.0), "log likelihood expected")
        ("gammashape,s", boost::program_options::value(&_gamma_shape)->default_value(0.5), "shape parameter of the Gamma among-site rate heterogeneity model")
        ("ncateg,c",     boost::program_options::value(&_num_categ)->default_value(1),     "number of categories in the discrete Gamma rate heterogeneity model")
        ("statefreq,f",  boost::program_options::value(&_state_frequencies)->multitoken()->default_value(std::vector<double> {0.25, 0.25, 0.25, 0.25}, "0.25 0.25 0.25 0.25"),  "state frequencies in the order A C G T (will be normalized to sum to 1)")
        ("rmatrix,r",    boost::program_options::value(&_exchangeabilities)->multitoken()->default_value(std::vector<double> {1, 1, 1, 1, 1, 1}, "1 1 1 1 1 1"),                "GTR exchangeabilities in the order AC AG AT CG CT GT (will be normalized to sum to 1)")
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

    // Be sure state frequencies sum to 1.0 and are all positive
    double sum_freqs = std::accumulate(_state_frequencies.begin(), _state_frequencies.end(), 0.0);
    for (auto & freq : _state_frequencies)
        {
        if (freq <= 0.0)
            throw XStrom("all statefreq entries must be positive real numbers");
        freq /= sum_freqs;
        }

    // Be sure exchangeabilities sum to 1.0 and are all positive
    double sum_xchg = std::accumulate(_exchangeabilities.begin(), _exchangeabilities.end(), 0.0);
    for (auto & xchg : _exchangeabilities)
        {
        if (xchg <= 0.0)
            throw XStrom("all rmatrix entries must be positive real numbers");
        xchg /= sum_xchg;
        }

    // Be sure gamma shape parameter is positive
    if (_gamma_shape <= 0.0)
        throw XStrom("gamma shape must be a positive real number");

    // Be sure number of gamma rate categories is greater than or equal to 1
    if (_num_categ < 1)
        throw XStrom("ncateg must be a positive integer greater than 0");
    }

inline void Strom::run()
    {
    std::cout << "Starting..." << std::endl;

    try
        {
        // Read and store data
        Data::SharedPtr d(new Data());
        d->getDataFromFile(_data_file_name);

        // Create a substitution model
        GTRModel::SharedPtr gtr = GTRModel::SharedPtr(new GTRModel());

        // old way
        //gtr->setExchangeabilitiesAndStateFreqs(
        //    {0.1307058, 0.1583282, 0.1598077, 0.2456609, 0.2202439, 0.08525337},
        //    {0.2249887, 0.2090694, 0.1375933, 0.4283486});
        //gtr->setGammaShape(1.480126);
        //gtr->setGammaNCateg(4);

        // new way
        gtr->setExchangeabilitiesAndStateFreqs(_exchangeabilities, _state_frequencies);
        gtr->setGammaShape(_gamma_shape);
        gtr->setGammaNCateg(_num_categ);

        std::cout << gtr->describeModel() << std::endl;

        // Create a likelihood object that will compute log-likelihoods
        Likelihood::SharedPtr likelihood(new Likelihood());
        likelihood->setData(d);
        likelihood->setModel(gtr);

        // Read in a tree
        TreeSummary::SharedPtr tree_summary(new TreeSummary());
        tree_summary->readTreefile(_tree_file_name, 0);
        Tree::SharedPtr tree = tree_summary->getTree(0);

        // Calculate the log-likelihood for the tree
        double lnL = likelihood->calcLogLikelihood(tree);
        std::cout << boost::str(boost::format("log likelihood = %.5f") % lnL) << std::endl;

        // old way
        //std::cout << "      (expecting -274.59466)" << std::endl;

        // new way
        std::cout << boost::str(boost::format("      (expecting %.5f)") % _expected_log_likelihood) << std::endl;
        }
    catch (XStrom x)
        {
        std::cerr << "Strom encountered a problem:\n  " << x.what() << std::endl;
        }

    std::cout << "\nFinished!" << std::endl;
    }

}

