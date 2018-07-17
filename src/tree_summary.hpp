#pragma once

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "split.hpp"
#include "tree_manip.hpp"
#include "xstrom.hpp"

#include "ncl/nxsmultiformat.h"

namespace strom
    {

    class TreeSummary
        {
        public:
                                        TreeSummary();
                                        ~TreeSummary();

            void                        readTreefile(const std::string filename, unsigned skip, unsigned tree_to_plot);
            void                        showSummary(unsigned tree_to_plot) const;
            typename Tree::SharedPtr    getTree(unsigned index);
            std::string                 getNewick(unsigned index);
            void                        clear();

        private:
        
            Split::treeid_t             _plotted;
            Split::treemap_t            _treeIDs;
            std::vector<std::string>    _newicks;
            std::vector<std::string>    _taxon_names;

        public:

            typedef std::shared_ptr< TreeSummary > SharedPtr;
        };

inline TreeSummary::TreeSummary()
    {
    //std::cout << "Constructing a TreeSummary" << std::endl;
    }

inline TreeSummary::~TreeSummary()
    {
    //std::cout << "Destroying a TreeSummary" << std::endl;
    }

inline Tree::SharedPtr TreeSummary::getTree(unsigned index)
    {
    if (index >= _newicks.size())
        throw XStrom("getTree called with index >= number of stored trees");

    TreeManip tm;

    // build the tree
    tm.buildFromNewick(_newicks[index], false, false);

    return tm.getTree();
    }

inline std::string TreeSummary::getNewick(unsigned index)
    {
    if (index >= _newicks.size())
        throw XStrom("getNewick called with index >= number of stored trees");

        return _newicks[index];
    }

inline void TreeSummary::clear()
    {
    _treeIDs.clear();
    _newicks.clear();
    }

inline void TreeSummary::readTreefile(const std::string filename, unsigned skip, unsigned tree_to_plot)
    {
    TreeManip tm;
    Split::treeid_t splitset;

    // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

    MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
    try {
        nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
    catch(...)
        {
        nexusReader.DeleteBlocksFromFactories();
        throw;
        }

    int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    for (int i = 0; i < numTaxaBlocks; ++i)
        {
        clear();
        NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
        std::string taxaBlockTitle = taxaBlock->GetTitle();
        
        _taxon_names.clear();
        _taxon_names = taxaBlock->GetAllLabels();

        const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
        for (unsigned j = 0; j < nTreesBlocks; ++j)
            {
            const NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, j);
            unsigned ntrees = treesBlock->GetNumTrees();
            if (skip < ntrees)
                {
                //std::cout << "Trees block contains " << ntrees << " tree descriptions.\n";
                for (unsigned t = skip; t < ntrees; ++t)
                    {
                    const NxsFullTreeDescription & d = treesBlock->GetFullTreeDescription(t);

                    // store the newick tree description
                    std::string newick = d.GetNewick();;
                    _newicks.push_back(newick);
                    unsigned tree_index = (unsigned)_newicks.size() - 1;

                    // build the tree
                    tm.buildFromNewick(newick, false, false);
                    //std::cerr << "newick = " << newick << std::endl;

                    // store set of splits
                    if (t != tree_to_plot - 1)
                        {
                        //std::cerr << "tree " << t << ", tree_to_plot = " << tree_to_plot << ": storing splits in _treeIDs" << std::endl;

                        splitset.clear();
                        tm.storeSplits(splitset);

                        // iterator iter will point to the value corresponding to key splitset
                        // or to end (if splitset is not already a key in the map)
                        Split::treemap_t::iterator iter = _treeIDs.lower_bound(splitset);

                        if (iter == _treeIDs.end() || iter->first != splitset)
                            {
                            // splitset key not found in map, need to create an entry
                            std::vector<unsigned> v(1, tree_index);  // vector of length 1 with only element set to tree_index
                            _treeIDs.insert(iter, Split::treemap_t::value_type(splitset, v));
                            }
                        else
                            {
                            // splitset key was found in map, need to add this tree's index to vector
                            iter->second.push_back(tree_index);
                            }
                        }
                    else
                        {
                        //std::cerr << "tree " << t << ", tree_to_plot = " << tree_to_plot << ": storing splits in _plotted" << std::endl;
                        _plotted.clear();
                        tm.storeSplits(_plotted);
                        }

                    } // trees loop
                } // if skip < ntrees
            } // TREES block loop
        } // TAXA block loop

    // No longer any need to store raw data from nexus file
    nexusReader.DeleteBlocksFromFactories();
    }
    
inline void TreeSummary::showSummary(unsigned tree_to_plot) const
    {
    unsigned num_input_trees = (unsigned)_newicks.size() - 1;

    // Produce some output to show that it works
    std::cout << boost::str(boost::format("\nRead %d trees from file") % _newicks.size()) << std::endl;
    std::cout << boost::str(boost::format("\nTree to plot is number %d") % tree_to_plot) << std::endl;
    std::cout << boost::str(boost::format("\nNumber of input trees %d") % num_input_trees) << std::endl;
    
    // Build tree to plot
    TreeManip tm;
    tm.buildFromNewick(_newicks[tree_to_plot-1], false, false);
    Split::treeid_t splitset;
    tm.storeSplits(splitset);

    // Show all unique topologies with a list of the trees that have that topology

    // Create a vector of (gene tree frequency, topology index) tuples sorted by increasing frequency (sorted_trees)
    typedef std::pair<unsigned, unsigned> sorted_trees_pair_t;
    std::vector< sorted_trees_pair_t > sorted_trees;

    // Also create a map (clademap) that can be used to find the gene tree frequency of any split
    typedef std::map<Split, unsigned> clademap_t;
    clademap_t clademap;

    // Here is the loop over stored trees:
    //  key_value_pair.first is a tree ID (set of splits defining the tree topology)
    //  key_value_pair.second is a vector of indices into _newicks
    unsigned topology = 0;
    for (auto & key_value_pair : _treeIDs)
        {
        topology++;
        
        // key_value_pair.first is a so-called tree ID; a set of all splits in the tree
        const Split::treeid_t & splitset = key_value_pair.first;

        // clademap is a map in which keys are splits and values are frequencies
        for (auto & s : splitset)
            {
            clademap_t::iterator lowb = clademap.lower_bound(s);
            if (lowb != clademap.end() && !(clademap.key_comp()(s, lowb->first)))
                {
                // this pattern has already been seen
                lowb->second += 1;
                }
            else
                {
                // this pattern has not yet been seen
                clademap.insert(lowb, clademap_t::value_type(s, 1));
                }
            }

        // key_value_pair.second is a vector of tree indices (of all trees having the same topology)
        unsigned ntrees = (unsigned)key_value_pair.second.size();

        // sorted_trees vector holds tuples (n,t), where t is the index of a tree topology and
        // n is the number of trees in the file with topology t
        sorted_trees.push_back(std::pair<unsigned, unsigned>(ntrees,topology));
        std::cout << "Topology " << topology << " seen in these " << ntrees << " trees:" << std::endl << "  ";
        std::copy(key_value_pair.second.begin(), key_value_pair.second.end(), std::ostream_iterator<unsigned>(std::cout, " "));
        std::cout << std::endl;
        }
        
    // Create sorted_clades vector of tuples (n,s), where n is the frequency of Split s
    typedef std::pair<unsigned, Split> sorted_clades_pair_t;
    std::vector< sorted_clades_pair_t > sorted_clades;
    for (auto & m : clademap)
        {
        sorted_clades.push_back(std::make_pair(m.second, m.first));
        }
    std::sort(sorted_clades.begin(), sorted_clades.end());

    // Show tree topologies sorted from most to least frequent (note: sorted_trees does not include tree number tree_to_plot)
    std::sort(sorted_trees.begin(), sorted_trees.end());
    //unsigned npairs = (unsigned)sorted_trees.size();
    std::cout << "\nTopologies sorted by sample frequency:" << std::endl;
    std::cout << boost::str(boost::format("%20s %20s") % "topology" % "frequency") << std::endl;
    for (auto & ntrees_topol_pair : boost::adaptors::reverse(sorted_trees))
        {
        unsigned n = ntrees_topol_pair.first;
        unsigned t = ntrees_topol_pair.second;
        std::cout << boost::str(boost::format("%20d %20d") % t % n) << std::endl;
        }

    // Show key to taxa in split representations
    std::cout << "\nTaxon names and the position of each in split representations below:" << std::endl;
    unsigned taxon_number = 0;
    for (auto & s : _taxon_names)
        {
        std::cout << boost::str(boost::format("%12d %s") % (++taxon_number) % s) << std::endl;
        }
        
    // Show support for clades in the tree_to_plot
    std::cout << boost::str(boost::format("\nGSF and IC for the %d splits in tree number %d") % _plotted.size() % tree_to_plot) << std::endl;
    std::cout << boost::str(boost::format("%20s %s") % "support" % "split") << std::endl;
    for (auto & s : _plotted)
        {
        int freq_s = -1;
        int freq_alt = -1;
        Split split_alt;
        for (auto & x : boost::adaptors::reverse(sorted_clades))
            {
            if (x.second == s)
                {
                freq_s = x.first;
                }
            if (freq_alt < 0 && s.conflictsWith(x.second))
                {
                freq_alt = x.first;
                split_alt = x.second;
                }
            }
            
        // Notes on why log(a)/log(2) equals log2(a)
        //
        // Show:
        //   log2(a) = log(a)/log(2)
        //
        // Proof:
        //   log2(a) = log(a)/log(2)
        //   log2(a) log(2) = log(a)
        //   log2(a) log(2) = log(2^{log2(a)}}) = log(a)
        //   log(a) = log(a) QED
            
        if (freq_s > 0)
            {
            // split found in tree to be plotted was also found in sorted_clades
            double focal_prop = 1.0*freq_s/num_input_trees;
            double ic = 1.0;  // if freq_alt still equals -1, then no conflicting split was found. In this case, IC = 1
            
            if (freq_alt >= 0.0)
                {
                // freq_alt is NOT still -1, so a conflicting split was found.
                double x1 = 1.0*freq_s/(freq_s + freq_alt);
                double x2 = 1.0*freq_alt/(freq_s + freq_alt);
                double logx1 = log(x1)/log(2);
                double logx2 = log(x2)/log(2);
                ic = 1.0 + x1*logx1 + x2*logx2;
                if (freq_alt > freq_s)
                    ic *= -1.0;
                std::cout << boost::str(boost::format("\n%20.3f %s") % focal_prop    % s.createPatternRepresentation()) << std::endl;
                double conflicting_prop = 1.0*freq_alt/num_input_trees;
                std::cout << boost::str(boost::format("%20.3f %s") % conflicting_prop % split_alt.createPatternRepresentation()) << std::endl;
                std::cout << boost::str(boost::format("%20.3f IC") % ic) << std::endl;
                }

            Node * nd = tm.findNodeWithSplit(s);
            assert(nd);
            nd->setSplitInfo(boost::str(boost::format("[%.5f,%.5f]") % focal_prop % ic));

            //clademap_t::iterator lowb = clademap.lower_bound(s);
            //if (lowb != clademap.end())
            //    {
            //    unsigned n = lowb->second;
            //    double pct = 100.0*n/num_input_trees;
            //    double ic = calcIC(s, );
            //    std::cout << boost::str(boost::format("%20.1f %s") % pct % s.createPatternRepresentation()) << std::endl;
            //    }
            // else
            //     {
            //     std::cout << boost::str(boost::format("%20s %s") % "null" % s.createPatternRepresentation()) << std::endl;
            //     }
            }
        else
            {
            // split found in tree to be plotted was NOT found in sorted_clades
            std::cout << boost::str(boost::format("\n%20s %s") % "null" % s.createPatternRepresentation()) << std::endl;
            }
        }
        
    std::ofstream outf("treedata.js", std::ios::out);
    outf << "var newick = \"" << tm.makeNewickNames(5, _taxon_names) << "\";" << std::endl;
    outf.close();
    }

}
