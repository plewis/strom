#pragma once

#include <climits>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/algorithm/string/replace.hpp>
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

            void                        readTreefiles(std::ostream & outf, const std::vector<std::string> filenames, unsigned skip, unsigned tree_to_plot, std::vector<std::string> & outgroup_names);
            void                        showSummary(std::ostream & outf) const;
            typename Tree::SharedPtr    getTree(unsigned index);
            std::string                 getNewick(unsigned index);
            void                        clear();

        private:
        
            unsigned                    _outgroup_node_number;
            Split::SharedPtr            _outgroup_split;
            unsigned                    _tree_to_plot;
            Split::treeid_t             _plotted;
            Split::treemap_t            _treeIDs;
            std::vector<std::string>    _newicks;
            std::vector<std::string>    _taxon_names;

            void                        readTreefile(std::ostream & outf, const std::string filename, unsigned skip, std::vector<std::string> & outgroup_names, bool tree_to_plot);

        public:

            typedef std::shared_ptr< TreeSummary > SharedPtr;
        };

inline TreeSummary::TreeSummary() : _outgroup_node_number(UINT_MAX), _tree_to_plot(0)
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

    // root the tree
    if (_outgroup_split)
        tm.rerootAtSplit(_outgroup_split);
    else if (_outgroup_node_number != UINT_MAX)
        tm.rerootAtNodeNumber(_outgroup_node_number);
        
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

inline void TreeSummary::readTreefiles(std::ostream & outf, const std::vector<std::string> filenames, unsigned skip, unsigned tree_to_plot, std::vector<std::string> & outgroup_names)
    {
    unsigned file_index = 0;
    for (auto fn : filenames)
        {
        if (++file_index == tree_to_plot) {
            readTreefile(outf, fn, skip, outgroup_names, true);
            }
        else {
            readTreefile(outf, fn, skip, outgroup_names, false);
            }
        }
    }
    
inline void TreeSummary::readTreefile(std::ostream & outf, const std::string filename, unsigned skip, std::vector<std::string> & outgroup_names, bool tree_to_plot)
    {
    TreeManip tm;
    Split::treeid_t splitset;
    
    if (tree_to_plot)
        skip = 0;

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
        //clear();
        NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
        std::string taxaBlockTitle = taxaBlock->GetTitle();
        
        //_taxon_names.clear();
        if (_taxon_names.size() == 0) {
            _taxon_names = taxaBlock->GetAllLabels();
            }
        else {
            std::vector<std::string> tmp = taxaBlock->GetAllLabels();
            for (unsigned i = 0; i < tmp.size(); ++i) {
                if (tmp[i] != _taxon_names[i])
                    throw XStrom(boost::str(boost::format("Taxon order different in file \"%s\"") % filename));
                }
        }
        
        // if one outgroup taxon was specified, store leaf node number in _outgroup_node_number
        // if more than one outgroup taxon was specified, store split of outgroup MRCA in _outgroup_split
        if (_outgroup_node_number == UINT_MAX)
            {
            if (outgroup_names.size() == 1)
                {
                // store index of taxon name in _taxon_names vector in _outgroup_node_number
                // leave _outgroup_split empty
                assert(outgroup_names[0].length() > 0);

                //unsigned t = 0;
                //for (auto nm : _taxon_names) {
                //    outf << boost::str(boost::format("%d --> %s") % (t++) % nm) << std::endl;
                //    }

                std::string outgroup_name = outgroup_names[0];
                boost::replace_all(outgroup_name, "_", " ");
                auto it = std::find(_taxon_names.begin(), _taxon_names.end(), outgroup_name);
                if (it != _taxon_names.end())
                    {
                    auto idx = std::distance(_taxon_names.begin(), it);
                    _outgroup_node_number = (unsigned)idx;
                    outf << boost::str(boost::format("Outgroup found: index of \"%s\" was %d") % outgroup_name % _outgroup_node_number)<< std::endl;
                    }
                else
                    throw XStrom(boost::str(boost::format("Outgroup NOT found: \"%s\" not located in _taxon_names") % outgroup_name));
                }
            else
                {
                // store split with bits set for all outgroup taxa in _outgroup_split
                // leave _outgroup_node_number set to UINT_MAX
                _outgroup_split = Split::SharedPtr(new Split);
                _outgroup_split->resize((unsigned)_taxon_names.size());
                for (auto outgroup_name : outgroup_names)
                    {
                    boost::replace_all(outgroup_name, "_", " ");
                    auto it = std::find(_taxon_names.begin(), _taxon_names.end(), outgroup_name);
                    if (it != _taxon_names.end())
                        {
                        auto idx = std::distance(_taxon_names.begin(), it);
                        unsigned taxon_index = (unsigned)idx;
                        _outgroup_split->setBitAt(taxon_index);
                        outf << boost::str(boost::format("Outgroup taxon found: index of \"%s\" was %d") % outgroup_name % taxon_index)<< std::endl;
                        }
                    else
                        outf << boost::str(boost::format("Outgroup taxon NOT found: \"%s\" not located in _taxon_names vector") % outgroup_name) << std::endl;
                    }
                //outf << "outgroup split: " << _outgroup_split->createPatternRepresentation() << std::endl;
                }
            }

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
                    std::string newick = d.GetNewick();
                    _newicks.push_back(newick);
                    unsigned tree_index = (unsigned)_newicks.size() - 1;

                    // build the tree
                    tm.buildFromNewick(newick, false, false);
                    if (_outgroup_split)
                        tm.rerootAtSplit(_outgroup_split);
                    else if (_outgroup_node_number != UINT_MAX)
                        tm.rerootAtNodeNumber(_outgroup_node_number);

                    // store set of splits
                    if (tree_to_plot)
                        {
                        _plotted.clear();
                        tm.storeSplits(_plotted);
                        _tree_to_plot = tree_index;
                        outf << "\n\nnewick = " << newick << "\n" << std::endl;
                        break;  // assume first tree is the one to plot
                        }
                    else
                        {
                        //outf << "tree " << t << ", tree_to_plot = " << tree_to_plot << ": storing splits in _treeIDs" << std::endl;

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

                    } // trees loop
                } // if skip < ntrees
            } // TREES block loop
        } // TAXA block loop

    // No longer any need to store raw data from nexus file
    nexusReader.DeleteBlocksFromFactories();
    }
    
inline void TreeSummary::showSummary(std::ostream & outf) const
    {
    unsigned num_input_trees = (unsigned)_newicks.size() - 1;

    // Produce some output to show that it works
    outf << boost::str(boost::format("A total of %d trees were read from files") % _newicks.size()) << std::endl;
    outf << boost::str(boost::format("%6d tree to plot") % 1) << std::endl;
    outf << boost::str(boost::format("%6d input trees") % num_input_trees) << std::endl;
    
    // Build tree to plot
    TreeManip tm;
    tm.buildFromNewick(_newicks[_tree_to_plot], false, false);
    if (_outgroup_split)
        tm.rerootAtSplit(_outgroup_split);
    else if (_outgroup_node_number != UINT_MAX)
        tm.rerootAtNodeNumber(_outgroup_node_number);
    outf << "\nTree to plot: " << _newicks[_tree_to_plot] << "\n" << std::endl;
    Split::treeid_t splitset;
    tm.storeSplits(splitset);   //TODO should be identical to _plotted

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
        unsigned topo_freq = (unsigned)key_value_pair.second.size();

        std::cerr << "Processing topology " << topology << " (frequency = " << topo_freq << ")" << std::endl; //temporary!

        // clademap is a map in which keys are splits and values are frequencies
        for (auto & s : splitset)
            {
            clademap_t::iterator lowb = clademap.lower_bound(s);
            if (lowb != clademap.end() && !(clademap.key_comp()(s, lowb->first)))
                {
                // this pattern has already been seen
                std::cerr << "  split already seen in " << lowb->second << " trees: " << s.createPatternRepresentation() << std::endl; //temporary!
                lowb->second += topo_freq;
                }
            else
                {
                // this pattern has not yet been seen
                std::cerr << "  new split: " << s.createPatternRepresentation() << std::endl; //temporary!
                clademap.insert(lowb, clademap_t::value_type(s, topo_freq));
                }
            }

        // key_value_pair.second is a vector of tree indices (of all trees having the same topology)
        unsigned ntrees = (unsigned)key_value_pair.second.size();

        // sorted_trees vector holds tuples (n,t), where t is the index of a tree topology and
        // n is the number of trees in the file with topology t
        sorted_trees.push_back(std::pair<unsigned, unsigned>(ntrees,topology));
        outf << "Topology " << topology << " seen in these " << ntrees << " trees:" << std::endl << "  ";
        std::copy(key_value_pair.second.begin(), key_value_pair.second.end(), std::ostream_iterator<unsigned>(outf, " "));
        outf << std::endl;
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
    outf << "\nTopologies sorted by sample frequency:" << std::endl;
    outf << boost::str(boost::format("%20s %20s") % "topology" % "frequency") << std::endl;
    for (auto & ntrees_topol_pair : boost::adaptors::reverse(sorted_trees))
        {
        unsigned n = ntrees_topol_pair.first;
        unsigned t = ntrees_topol_pair.second;
        outf << boost::str(boost::format("%20d %20d") % t % n) << std::endl;
        }

    // Show key to taxa in split representations
    outf << "\nTaxon names and the position of each in split representations below:" << std::endl;
    unsigned taxon_number = 0;
    for (auto & s : _taxon_names)
        {
        outf << boost::str(boost::format("%12d %s") % (++taxon_number) % s) << std::endl;
        }
        
    // Show support for clades in the tree_to_plot
    outf << boost::str(boost::format("\nGSF and IC for the %d splits in tree number %d") % _plotted.size() % _tree_to_plot) << std::endl;
    outf << boost::str(boost::format("%20s %s") % "support" % "split") << std::endl;
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
        //   a = 2^{log2(a)}               by definition of logarithm base 2
        //   log(a) = log{ 2^{log2(a)} }   take natural log of both sides
        //   log(a) = log2(a) log(2)       because log(a^b) = b log(a)
        //   log(a)/log(2) = log2(a)       divide both sides by log(2)
            
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
                outf << boost::str(boost::format("\n%20.3f %s") % focal_prop    % s.createPatternRepresentation()) << std::endl;
                double conflicting_prop = 1.0*freq_alt/num_input_trees;
                outf << boost::str(boost::format("%20.3f %s") % conflicting_prop % split_alt.createPatternRepresentation()) << std::endl;
                outf << boost::str(boost::format("%20.3f IC") % ic) << std::endl;
                }

            Node * nd = tm.findNodeWithSplit(s);
            if (!nd) {
                std::cerr << "oops" << std::endl;
                nd = tm.findNodeWithSplit(s);
                }
            assert(nd);
            nd->setSplitInfo(boost::str(boost::format("[%.5f,%.5f]") % focal_prop % ic));
            
            std::cerr << "setting split info for node " << nd->getNumber() << ": focal_prop = " << focal_prop << ", IC = " << ic << ", split_info = " << nd->getSplitInfo() << std::endl; //temporary!

            //clademap_t::iterator lowb = clademap.lower_bound(s);
            //if (lowb != clademap.end())
            //    {
            //    unsigned n = lowb->second;
            //    double pct = 100.0*n/num_input_trees;
            //    double ic = calcIC(s, );
            //    outf << boost::str(boost::format("%20.1f %s") % pct % s.createPatternRepresentation()) << std::endl;
            //    }
            // else
            //     {
            //     outf << boost::str(boost::format("%20s %s") % "null" % s.createPatternRepresentation()) << std::endl;
            //     }
            }
        else
            {
            // split found in tree to be plotted was NOT found in sorted_clades
            outf << boost::str(boost::format("\n%20s %s") % "null" % s.createPatternRepresentation()) << std::endl;
            }
        }
        
    std::ofstream jsfile("treedata.js", std::ios::out);
    jsfile << "var newick = \"" << tm.makeNewickNames(5, _taxon_names) << "\";" << std::endl;
    jsfile.close();
    }

}
