#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <map>
#include <boost/format.hpp>
#include "xstrom.hpp"

#include "ncl/nxsmultiformat.h"

namespace strom {

class Data
    {
	public:
        typedef std::vector<double>             pattern_counts_t;
        typedef std::vector<std::string>        taxon_names_t;
        typedef std::vector<int>                pattern_t;
        typedef std::map< pattern_t, unsigned > pattern_map_t;
        typedef std::vector< pattern_t >        data_matrix_t;
        typedef std::shared_ptr< Data >         SharedPtr;

                                                Data();
                                                ~Data();

        void                                    getDataFromFile(const std::string filename);

        const pattern_counts_t &                getPatternCounts() const;
        const taxon_names_t &                   getTaxonNames() const;
        const data_matrix_t &                   getDataMatrix() const;

        void                                    clear();
        unsigned                                getNumPatterns() const;
        unsigned                                getSeqLen() const;
        unsigned                                getNumTaxa() const;

    private:

        void                                    updatePatternMap(Data::pattern_t & pattern);
        void                                    compressPatterns();

        pattern_map_t                           _pattern_map;       // used as workspace
        pattern_counts_t                        _pattern_counts;
        taxon_names_t                           _taxon_names;
        data_matrix_t                           _data_matrix;
    };

inline Data::Data()
    {
    std::cout << "Creating Data object" << std::endl;
    }

inline Data::~Data()
    {
    std::cout << "Destroying Data object" << std::endl;
    }

inline const Data::pattern_counts_t & Data::getPatternCounts() const
	{
	return _pattern_counts;
	}

inline const Data::taxon_names_t & Data::getTaxonNames() const
	{
	return _taxon_names;
	}

inline const Data::data_matrix_t & Data::getDataMatrix() const
	{
	return _data_matrix;
	}

inline void Data::clear()
    {
    _pattern_map.clear();
    _pattern_counts.clear();
    _taxon_names.clear();
    _data_matrix.clear();
    }

inline unsigned Data::getNumPatterns() const
    {
    return (unsigned)_pattern_counts.size();
    }

inline unsigned Data::getNumTaxa() const
    {
    return (unsigned)_taxon_names.size();
    }

inline unsigned Data::getSeqLen() const
    {
    return (unsigned)std::accumulate(_pattern_counts.begin(), _pattern_counts.end(), 0);
    }

inline void Data::compressPatterns()
    {
    // sanity checks
    if (_data_matrix.empty())
        throw XStrom("Attempted to compress an empty data matrix");
    if (_data_matrix[0].size() < getSeqLen())
        throw XStrom("Attempted to compress an already compressed data matrix");

    // create map with keys equal to patterns and values equal to site counts
    _pattern_map.clear();

    std::vector<int> pattern;
    unsigned ntaxa = (unsigned)_data_matrix.size();
    unsigned seqlen = (unsigned)_data_matrix[0].size();
    for (unsigned i = 0; i < seqlen; ++i)
        {
        // Create vector representing pattern at site i
        pattern.clear();
        for (unsigned j = 0; j < ntaxa; ++j)
            {
            pattern.push_back(_data_matrix[j][i]);
            }

        // Add this pattern to pattern_counts
        updatePatternMap(pattern);
        }

    // resize _pattern_counts
    unsigned npatterns = (unsigned)_pattern_map.size();
    _pattern_counts.resize(npatterns);

    // resize _data_matrix so that we can use operator[] to assign values
    _data_matrix.resize(ntaxa);
    for (auto & row : _data_matrix)
        {
        row.resize(npatterns);
        }

    unsigned j = 0;
    for (auto & pc : _pattern_map)
        {
        _pattern_counts[j] = pc.second;

        unsigned i = 0;
        for (auto sc : pc.first)
            {
            _data_matrix[i][j] = sc;
            ++i;
            }

        ++j;
        }

    // Everything has been transferred to _data_matrix and _pattern_counts, so can now free this memory
    _pattern_map.clear();

    unsigned total_num_sites = std::accumulate(_pattern_counts.begin(), _pattern_counts.end(), 0);
    if (seqlen != total_num_sites)
        throw XStrom(boost::str(boost::format("Total number of sites before compaction (%d) not equal to toal number of sites after (%d)") % seqlen % total_num_sites));
    }

inline void Data::updatePatternMap(Data::pattern_t & pattern)
    {
    // If pattern is not already in pattern_map, insert it and set value to 1.
    // If it does exist, increment its current value.
    // (see item 24, p. 110, in Meyers' Efficient STL for more info on the technique used here)
    pattern_map_t::iterator lowb = _pattern_map.lower_bound(pattern);
    if (lowb != _pattern_map.end() && !(_pattern_map.key_comp()(pattern, lowb->first)))
        {
        // this pattern has already been seen
        lowb->second += 1;
        }
    else
        {
        // this pattern has not yet been seen
        _pattern_map.insert(lowb, pattern_map_t::value_type(pattern, 1));
        }
    }

inline void Data::getDataFromFile(const std::string filename)
    {
    // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for documentation
    //
    // -1 means "process all blocks found" (this is a bit field and -1 fills the bit field with 1s)
    // Here are the bits (and nexus blocks) that are defined:
    //     enum NexusBlocksToRead
    //     {
    //         NEXUS_TAXA_BLOCK_BIT = 0x01,
    //         NEXUS_TREES_BLOCK_BIT = 0x02,
    //         NEXUS_CHARACTERS_BLOCK_BIT = 0x04,
    //         NEXUS_ASSUMPTIONS_BLOCK_BIT = 0x08,
    //         NEXUS_SETS_BLOCK_BIT = 0x10,
    //         NEXUS_UNALIGNED_BLOCK_BIT = 0x20,
    //         NEXUS_DISTANCES_BLOCK_BIT = 0x40,
    //         NEXUS_UNKNOWN_BLOCK_BIT = 0x80
    //     };
    MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
    try {
        nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
    catch(...)
        {
        nexusReader.DeleteBlocksFromFactories();
        throw;
        }

    // Commit to storing new data
    clear();

    int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    for (int i = 0; i < numTaxaBlocks; ++i)
        {
        NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
        _taxon_names.clear();
        for (auto s : taxaBlock->GetAllLabels()) {
            _taxon_names.push_back(s);
            }
        unsigned ntax = (unsigned)_taxon_names.size();

        const unsigned numCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
        for (unsigned j = 0; j < numCharBlocks; ++j)
            {
            const NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, j);
            std::string charBlockTitle = taxaBlock->GetTitle();

            _data_matrix.resize(ntax);
            for (unsigned t = 0; t < ntax; ++t)
                {
                const NxsDiscreteStateRow & row = charBlock->GetDiscreteMatrixRow(t);
                unsigned seqlen = (unsigned)row.size();
                _data_matrix[t].resize(seqlen);
                unsigned k = 0;
                for (auto state_code : row) {
                    if (state_code < 0)
                        _data_matrix[t][k++] = 4;
                    else
                        _data_matrix[t][k++] = state_code;
                    }
                }

            }
        }

    // No longer any need to store raw data from nexus file
    nexusReader.DeleteBlocksFromFactories();

    // Compress _data_matrix so that it holds only unique patterns (counts stored in _pattern_counts)
    compressPatterns();
    }


}
