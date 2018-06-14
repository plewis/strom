#pragma once

#include "likelihood.hpp"
#include "model.hpp"
#include "lot.hpp"
#include "tree_summary.hpp"

namespace strom
{

class PWK
	{
	public:
                                            PWK(Lot::SharedPtr lot, TreeSummary::SharedPtr sumt);
                                            ~PWK();

        double                              logMarginalLikelihood();

    private:

        double                              sortTrees();

        TreeSummary::sorted_vect_t          _sorted_trees;

        Lot::SharedPtr                      _lot;                             // pseudorandom number generator
        //Likelihood::SharedPtr               _likelihood;                      // the likelihood calculator
        //Model::SharedPtr                    _model;                           // the substitution model
        TreeSummary::SharedPtr              _tree_summary;                    // contains treeID map and newick vector stored from tree sample file
        //ParamSummary::SharedPtr             _param_summary;                 // contains vector of parameter vectors stored from parameter sample file
        //Tree::SharedPtr                     _tree;                          // used to store current tree topology
        //double                              _log_num_topol;                 // log(number of possible tree topologies) used in computing discrete uniform topology prior
        //double                              _edgelen_prior_lambda;          // rate parameter of Exponential distribution used for edge length
        //double                              _shape_prior_lambda;            // rate parameter of Exponential distribution used for gamma shape
        //double                              _exchangeability_prior_param;   // Dirichlet parameter for exchangeabilities (1 for flat prior)
        //double                              _frequency_prior_param;         // Dirichlet parameter for frequencies (1 for flat prior)
        //unsigned                            _minimum_sample_size;           // smallest sample size for which tree topology will be included
        //unsigned                            _num_shells;                    // number of concentric shells used in the pwk estimator
        double                              _log_total_sampled_trees;

	public:
		typedef boost::shared_ptr< PWK >            SharedPtr;
        typedef std::pair<double,double>            double_pair_t;

    }; // class PWK

inline PWK::PWK(Lot::SharedPtr lot, TreeSummary::SharedPtr sumt)
  : _lot(lot), _tree_summary(sumt)
    {
    }

inline PWK::~PWK()
    {
    }

inline double PWK::sortTrees()
    {
    unsigned ntrees = _tree_summary->getNumStoredTrees();

    // Show how many trees are stored in _tree_summary
    std::cout << boost::str(boost::format("Read %d trees from file") % ntrees) << std::endl;

    // Let TreeSummary do the work
    _tree_summary->sortTrees(_sorted_trees);

    return log(ntrees);
    }

inline double PWK::logMarginalLikelihood()
    {
    _log_total_sampled_trees = sortTrees();
    std::cout << boost::str(boost::format("\nTotal sampled trees: %d (log scale: %.5f)\n") % exp(_log_total_sampled_trees) % _log_total_sampled_trees);

    std::cout << "\nTopologies sorted by sample frequency:" << std::endl;
    std::cout << boost::str(boost::format("%20s %20s") % "topology" % "frequency") << std::endl;
    for (auto & ntrees_topol_pair : boost::adaptors::reverse(_sorted_trees))
        {
        unsigned n = ntrees_topol_pair.first;
        unsigned t = ntrees_topol_pair.second;
        std::cout << boost::str(boost::format("%20d %20d") % t % n) << std::endl;
        }

    return 0.0;

#if 0
    _posterior_proportion = 0.0;
    unsigned which_tree = 0;
    unsigned cum_trees_sampled = 0;
    std::vector<double> log_numerators;
    unsigned num_topol_considered = 0;
    BOOST_REVERSE_FOREACH(sorted_pair_t & p, _sortedTrees)
        {
        if (p.first < _minimum_sample_size)
            break;

        ++which_tree;
        unsigned num_trees_sampled = p.first;
        cum_trees_sampled += num_trees_sampled;

        // Get vector of indices of the p.first sampled trees having topology specified by tree ID p.second
        Split::treeid_t & tree_id = p.second;
        Split::treemap_t & tree_map = _tree_summary->getTreeIDMap();
        std::vector<unsigned> & tree_indices = tree_map[tree_id];

        storePosteriorSamples(which_tree, tree_indices);
        unsigned num_params = standardizeSamples();
        double this_numerator = matryoshka(which_tree, num_params, _num_shells);

        // debugging
        double log_c_this_topology = log(_posterior_samples.size()) - this_numerator + _log_num_topol;
        double post_prob_this_topology = (double)_posterior_samples.size()/exp(_log_total_sampled_trees);
        double info_contribution = post_prob_this_topology*log_c_this_topology;
        std::cerr << std::endl;
        std::cerr << "This tree topology only:" << std::endl;
        std::cerr << "  log(total num. topologies)    = " << _log_num_topol << std::endl;
        std::cerr << "  newick description            = " << _tree_summary->getNewick(tree_indices[0]) << std::endl;
        std::cerr << "  log(numerator)                = " << this_numerator << std::endl;
        std::cerr << "  number of samples             = " << _posterior_samples.size() << std::endl;
        std::cerr << "  log(number of samples)        = " << log(_posterior_samples.size()) << std::endl;
        std::cerr << "  log(c)                        = " << log_c_this_topology << std::endl;
        std::cerr << "  posterior probability pr(tau) = " << post_prob_this_topology << std::endl;
        std::cerr << "  pr(tau) * log(c)              = " << info_contribution << std::endl;

        log_numerators.push_back(this_numerator);
        num_topol_considered++;
        }

    double max_numerator_term = *std::max_element(log_numerators.begin(), log_numerators.end());
    double numerator_sum = 0.0;
    BOOST_FOREACH(double d, log_numerators)
        {
        numerator_sum += exp(d - max_numerator_term);
        }
    double log_numerator = max_numerator_term + log(numerator_sum);

    std::cout << std::endl;
    std::cout << "Number of topologies considered: " << num_topol_considered << std::endl;
    std::cout << "log(T) = log(number of samples): " << _log_total_sampled_trees << std::endl;
    std::cerr << "proportion posterior included  = " << _posterior_proportion << std::endl;
    std::cout << std::endl;

    //temporary!
    //std::cout << "log(num_topol_considered) = " << log(num_topol_considered) << std::endl;
    //std::cout << "_log_num_topol            = " << _log_num_topol << std::endl;
    //std::cout << "log_numerator             = " << log_numerator << std::endl;
    //std::cout << "_log_total_sampled_trees                     = " << _log_total_sampled_trees << std::endl;

    //double log_c = (log(num_topol_considered) - _log_num_topol) - (log_numerator - _log_total_sampled_trees);
    double log_c =  log(_posterior_proportion) - (log_numerator - _log_total_sampled_trees);  //temporary
    //double log_c =  _log_total_sampled_trees - log_numerator;  //temporary
#endif
    }

} // namespace strom
