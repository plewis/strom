#pragma once

#include <iostream>
#include "tree_summary.hpp"
#include "data.hpp"
#include "likelihood.hpp"
#include "lot.hpp"
#include "pwk.hpp"
#include "chain.hpp"
#include <boost/program_options.hpp>
#include "output_manager.hpp"

namespace strom {

class Strom
    {
    public:
                                 Strom();
                                    ~Strom();

        void                        clear();
        void                        processCommandLineOptions(int argc, const char * argv[]);
        void                        run();

    private:

        std::string                 _data_file_name;
        std::string                 _tree_file_name;

        double                      _expected_log_likelihood;
        double                      _gamma_shape;
        unsigned                    _num_categ;
        std::vector<double>         _state_frequencies;
        std::vector<double>         _exchangeabilities;

        Data::SharedPtr             _data;
        Model::SharedPtr            _model;
        Likelihood::SharedPtr       _likelihood;
        TreeSummary::SharedPtr      _tree_summary;
        Lot::SharedPtr              _lot;

        unsigned                    _random_seed;
        unsigned                    _num_iter;
        unsigned                    _num_burnin_iter;
        bool                        _using_stored_data;
        unsigned                    _sample_freq;

        unsigned                    _num_chains;
        double                      _heating_lambda;
        std::vector<Chain>          _chains;
        std::vector<double>         _heating_powers;
        std::vector<unsigned>       _swaps;

        void                        sample(unsigned iter, Chain & chain);

        void                        calcHeatingPowers();
        void                        initChains();
        void                        stopTuningChains();
        void                        stepChains(unsigned iteration, bool sampling);
        void                        swapChains();
        void                        stopChains();
        void                        swapSummary() const;
        void                        showLambdas() const;

        OutputManager::SharedPtr    _output_manager;

        static std::string          _program_name;
        static unsigned             _major_version;
        static unsigned             _minor_version;

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
    _data_file_name          = "";
    _tree_file_name          = "";
    _data                    = nullptr;
    _model                   = nullptr;
    _likelihood              = nullptr;
    _tree_summary            = nullptr;
    _lot                     = nullptr;
    _output_manager          = nullptr;
    _expected_log_likelihood = 0.0;
    _gamma_shape             = 0.5;
    _num_categ               = 1;
    _random_seed             = 1;
    _num_iter                = 1000;
    _sample_freq             = 1;
    _num_burnin_iter         = 1000;
    _heating_lambda          = 0.5;
    _num_chains              = 1;
    _using_stored_data       = true;

    _state_frequencies.resize(0);
    _exchangeabilities.resize(0);
    _chains.resize(0);
    _heating_powers.resize(0);
    _swaps.resize(0);
    }

inline void Strom::processCommandLineOptions(int argc, const char * argv[])
    {
    boost::program_options::variables_map       vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("seed,z",        boost::program_options::value(&_random_seed)->default_value(1),   "pseudorandom number seed")
        ("niter,n",       boost::program_options::value(&_num_iter)->default_value(1000),   "number of MCMC iterations")
        ("samplefreq",  boost::program_options::value(&_sample_freq)->default_value(1),   "skip this many iterations before sampling next")
        ("datafile,d",  boost::program_options::value(&_data_file_name)->required(), "name of data file in NEXUS format")
        ("treefile,t",  boost::program_options::value(&_tree_file_name)->required(), "name of data file in NEXUS format")
        ("expectedLnL", boost::program_options::value(&_expected_log_likelihood)->default_value(0.0), "log likelihood expected")
        ("gammashape,s", boost::program_options::value(&_gamma_shape)->default_value(0.5), "shape parameter of the Gamma among-site rate heterogeneity model")
        ("ncateg,c",     boost::program_options::value(&_num_categ)->default_value(1),     "number of categories in the discrete Gamma rate heterogeneity model")
        ("statefreq,f",  boost::program_options::value(&_state_frequencies)->multitoken()->default_value(std::vector<double> {0.25, 0.25, 0.25, 0.25}, "0.25 0.25 0.25 0.25"),  "state frequencies in the order A C G T (will be normalized to sum to 1)")
        ("rmatrix,r",    boost::program_options::value(&_exchangeabilities)->multitoken()->default_value(std::vector<double> {1, 1, 1, 1, 1, 1}, "1 1 1 1 1 1"),                "GTR exchangeabilities in the order AC AG AT CG CT GT (will be normalized to sum to 1)")
        ("nchains",       boost::program_options::value(&_num_chains)->default_value(1),                "number of chains")
        ("heatfactor",    boost::program_options::value(&_heating_lambda)->default_value(0.5),          "determines how hot the heated chains are")
        ("burnin",        boost::program_options::value(&_num_burnin_iter)->default_value(100),         "number of iterations used to burn in chains")
        ("usedata",       boost::program_options::value(&_using_stored_data)->default_value(true),      "use the stored data in calculating likelihoods (specify no to explore the prior)")
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

    // Be sure number of chains is greater than or equal to 1
    if (_num_chains < 1)
        throw XStrom("nchains must be a positive integer greater than 0");

    // Be sure heatfactor is between 0 and 1
    if (_heating_lambda <= 0.0 || _heating_lambda > 1.0)
        throw XStrom("heatfactor must be a real number in the interval (0.0,1.0]");

    if (!_using_stored_data)
        std::cout << "\n*** Not using stored data (posterior = prior) ***\n" << std::endl;
    }

inline void Strom::calcHeatingPowers()
    {
    // Specify chain heating power (e.g. _heating_lambda = 0.2)
    // chain_index  power
    //      0       1.000 = 1/(1 + 0.2*0)
    //      1       0.833 = 1/(1 + 0.2*1)
    //      2       0.714 = 1/(1 + 0.2*2)
    //      3       0.625 = 1/(1 + 0.2*3)
    unsigned i = 0;
    for (auto & h : _heating_powers)
        {
        h = 1.0/(1.0 + _heating_lambda*i++);
        }
    }

inline void Strom::initChains()
    {
    // Create _num_chains chains
    _chains.resize(_num_chains);

    // Create _num_chains by _num_chains swap matrix
    _swaps.assign(_num_chains*_num_chains, 0);
    std::cout << "Number of chains = " << _num_chains << std::endl;

    // Create heating power vector
    _heating_powers.assign(_num_chains, 1.0);
    calcHeatingPowers();

    // Initialize chains
    unsigned chain_index = 0;
    for (auto & c : _chains)
        {
        // Give the chain a starting tree
        std::string newick = _tree_summary->getNewick(0);
        c.setTreeFromNewick(newick);

        // Set the pseudorandom number generator
        c.setLot(_lot);

        // Create a substitution model
        Model::SharedPtr model = Model::SharedPtr(new Model());
        model->setExchangeabilitiesAndStateFreqs(_exchangeabilities, _state_frequencies);
        model->setGammaShape(_gamma_shape);
        model->setGammaNCateg(_num_categ);
        model->useStoredData(_using_stored_data);

        // Create a likelihood object that will compute log-likelihoods
        Likelihood::SharedPtr likelihood = Likelihood::SharedPtr(new Likelihood());
        likelihood->setData(_data);
        likelihood->setModel(model);
        likelihood->useStoredData(_using_stored_data);

        // Provide the chain a likelihood calculator
        c.setLikelihood(likelihood);

        // Tell the chain that it should adapt its updators (at least initially)
        c.startTuning();

        // Set heating power to precalculated value
        c.setChainIndex(chain_index);
        c.setHeatingPower(_heating_powers[chain_index]);

        // Print headers in output files and make sure each updator has its starting value
        c.start();

        if (chain_index == 0)
            {
            // Summarize model
            std::cout << model->describeModel() << std::endl;

            // Calculate the log-likelihood for the tree
            Tree::SharedPtr tree = _tree_summary->getTree(0);
            double lnL = likelihood->calcLogLikelihood(tree);
            std::cout << boost::str(boost::format("log likelihood = %.5f") % lnL) << std::endl;
            if (_expected_log_likelihood != 0.0)
                std::cout << boost::str(boost::format("      (expecting %.5f)") % _expected_log_likelihood) << std::endl;
            }

        ++chain_index;
        }
    }

inline void Strom::showLambdas() const
    {
    for (unsigned idx = 0; idx < _num_chains; ++idx)
        {
        for (auto & c : _chains)
            {
            if (c.getChainIndex() == idx)
                {
                _output_manager->outputConsole(boost::str(boost::format("Chain %d (power %.5f)") % idx % c.getHeatingPower()));
                std::vector<std::string> names = c.getUpdaterNames();
                std::vector<double> lambdas    = c.getLambdas();
                std::vector<double> acceptpcts = c.getAcceptPercentages();
                unsigned n = (unsigned)names.size();
                _output_manager->outputConsole(boost::str(boost::format("%30s %15s %15s") % "Updater" % "Tuning Param." % "Accept %"));
                for (unsigned i = 0; i < n; ++i)
                    {
                    _output_manager->outputConsole(boost::str(boost::format("%30s %15.8f %15.1f") % names[i] % lambdas[i] % acceptpcts[i]));
                    }
                }
            }
        }
    }

inline void Strom::stopTuningChains()
    {
    _swaps.assign(_num_chains*_num_chains, 0);
    for (auto & c : _chains)
        {
        c.stopTuning();
        }
    }

inline void Strom::stepChains(unsigned iteration, bool sampling)
    {
    for (auto & c : _chains)
        {
        c.nextStep(iteration);
        if (sampling)
            sample(iteration, c);
        }
    }

inline void Strom::swapChains()
    {
    if (_num_chains == 1)
        return;

    // Select two chains at random to swap
    // If _num_chains = 3...
    //  i  j  = (i + 1 + randint(0,1)) % _num_chains
    // ---------------------------------------------
    //  0  1  = (0 + 1 +      0      ) %     3
    //     2  = (0 + 1 +      1      ) %     3
    // ---------------------------------------------
    //  1  2  = (1 + 1 +      0      ) %     3
    //     0  = (1 + 1 +      1      ) %     3
    // ---------------------------------------------
    //  2  0  = (2 + 1 +      0      ) %     3
    //     1  = (2 + 1 +      1      ) %     3
    // ---------------------------------------------
    unsigned i = _lot->randint(0, _num_chains-1);
    unsigned j = i + 1 + _lot->randint(0, _num_chains-2);
    j %= _num_chains;

    assert(i != j && i >=0 && i < _num_chains && j >= 0 && j < _num_chains);

    // Determine upper and lower triangle cells in _swaps vector
    unsigned smaller = _num_chains;
    unsigned larger  = _num_chains;
    double index_i   = _chains[i].getChainIndex();
    double index_j   = _chains[j].getChainIndex();
    if (index_i < index_j)
        {
        smaller = index_i;
        larger  = index_j;
        }
    else
        {
        smaller = index_j;
        larger  = index_i;
        }
    unsigned upper = smaller*_num_chains + larger;
    unsigned lower = larger*_num_chains  + smaller;
    _swaps[upper]++;

    // Propose swap of chains i and j
	// Proposed state swap will be successful if a uniform random deviate is less than R, where
	//    R = Ri * Rj = (Pi(j) / Pi(i)) * (Pj(i) / Pj(j))
    // Chain i: power = a, kernel = pi
    // Chain j: power = b, kernel = pj
    //      pj^a         pi^b
    // Ri = ----    Rj = ----
    //      pi^a         pj^b
    // log R = (a-b) [log(pj) - log(pi)]

    double heat_i       = _chains[i].getHeatingPower();
    double log_kernel_i = _chains[i].calcLogLikelihood() + _chains[i].calcLogJointPrior();

    double heat_j       = _chains[j].getHeatingPower();
    double log_kernel_j = _chains[j].calcLogLikelihood() + _chains[j].calcLogJointPrior();

    double logR = (heat_i - heat_j)*(log_kernel_j - log_kernel_i);

    double logu = _lot->logUniform();
    if (logu < logR)
        {
        // accept swap
        _swaps[lower]++;
        _chains[j].setHeatingPower(heat_i);
        _chains[i].setHeatingPower(heat_j);
        _chains[j].setChainIndex(index_i);
        _chains[i].setChainIndex(index_j);
        std::vector<double> lambdas_i = _chains[i].getLambdas();
        std::vector<double> lambdas_j = _chains[j].getLambdas();
        _chains[i].setLambdas(lambdas_j);
        _chains[j].setLambdas(lambdas_i);
        }

    }

inline void Strom::stopChains()
    {
    for (auto & c : _chains)
        c.stop();
    }

inline void Strom::swapSummary() const
    {
    if (_num_chains > 1)
        {
        unsigned i, j;
        std::cout << "\nSwap summary (upper triangle = no. attempted swaps; lower triangle = no. successful swaps):" << std::endl;

        // column headers
        std::cout << boost::str(boost::format("%12s") % " ");
        for (i = 0; i < _num_chains; ++i)
            std::cout << boost::str(boost::format(" %12d") % i);
        std::cout << std::endl;

        // top line
        std::cout << boost::str(boost::format("%12s") % "------------");
        for (i = 0; i < _num_chains; ++i)
            std::cout << boost::str(boost::format("-%12s") % "------------");
        std::cout << std::endl;

        // table proper
        for (i = 0; i < _num_chains; ++i)
            {
            std::cout << boost::str(boost::format("%12d") % i);
            for (j = 0; j < _num_chains; ++j)
                {
                if (i == j)
                    std::cout << boost::str(boost::format(" %12s") % "---");
                else
                    std::cout << boost::str(boost::format(" %12.5f") % _swaps[i*_num_chains + j]);
                }
            std::cout << std::endl;
            }

        // bottom line
        std::cout << boost::str(boost::format("%12s") % "------------");
        for (i = 0; i < _num_chains; ++i)
            std::cout << boost::str(boost::format("-%12s") % "------------");
        std::cout << std::endl;
        }
    }

inline void Strom::sample(unsigned iteration, Chain & chain)
    {
    if (chain.getHeatingPower() == 1 && iteration % _sample_freq == 0)
        {
        double logLike = chain.calcLogLikelihood();
        double logPrior = chain.calcLogJointPrior();
        double TL = chain.getTreeManip()->calcTreeLength();
        _output_manager->outputConsole(boost::str(boost::format("%12d %12.5f %12.5f %12.5f") % iteration % logLike % logPrior % TL));
        _output_manager->outputTree(iteration, chain.getTreeManip());
        _output_manager->outputParameters(iteration, logLike, logPrior, TL, chain.getModel());
        }
    }

inline void Strom::run()
    {
    std::cout << "Starting..." << std::endl;

    try
        {
        // Read and store data
        _data = Data::SharedPtr(new Data());
        _data->getDataFromFile(_data_file_name);

        // Read in trees
        _tree_summary = TreeSummary::SharedPtr(new TreeSummary());
        _tree_summary->readTreefile(_tree_file_name, 0);

        // Create a Lot object that generates (pseudo)random numbers
        _lot = Lot::SharedPtr(new Lot);
        _lot->setSeed(_random_seed);

#if 1
        PWK pwk(_lot, _tree_summary);
        pwk.logMarginalLikelihood();
#else
        // Create  Chain objects
        initChains();

        // Create an output manager and open output files
        _output_manager.reset(new OutputManager);
        _output_manager->outputConsole(boost::str(boost::format("\n%12s %12s %12s %12s") % "iteration" % "logLike" % "logPrior" % "TL"));
        _output_manager->openTreeFile("trees.tre", _data);
        _output_manager->openParameterFile("params.txt", _chains[0].getModel());
        sample(0, _chains[0]);

        // Burn-in the chains
        std::cout << "Burning in for " << _num_burnin_iter << " iterations... " << std::endl;
        for (unsigned iteration = 1; iteration <= _num_burnin_iter; ++iteration)
            {
            stepChains(iteration, false);
            swapChains();
            }

        std::cout << "Burn-in finished, no longer tuning updaters." << std::endl;
        stopTuningChains();
        showLambdas();

        // Sample the chains
        for (unsigned iteration = 1; iteration <= _num_iter; ++iteration)
            {
            stepChains(iteration, true);
            swapChains();
            }
        showLambdas();
        stopChains();

        // Create swap summary
        swapSummary();

        // Close output files
        _output_manager->closeTreeFile();
        _output_manager->closeParameterFile();
#endif
        }
    catch (XStrom & x)
        {
        std::cerr << "Strom encountered a problem:\n  " << x.what() << std::endl;
        }

    std::cout << "\nFinished!" << std::endl;
    }

} // namespace strom
