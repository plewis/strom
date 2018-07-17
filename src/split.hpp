#pragma once

#include <vector>
#include <memory>
#include <set>
#include <map>
#include <climits>
#include <cassert>

namespace strom
    {

    class Split
        {
        public:
                                                                Split();
                                                                Split(const Split & other);
                                                                ~Split();

            Split &                                             operator=(const Split & other);
            bool                                                operator==(const Split & other) const;
            bool                                                operator<(const Split & other) const;

            void                                                clear();
            void                                                resize(unsigned nleaves);

            typedef unsigned long                               split_unit_t;
            typedef std::vector<split_unit_t>                   split_t;
            typedef std::set<Split>                             treeid_t;
            typedef std::map< treeid_t, std::vector<unsigned> > treemap_t;
            typedef std::tuple<unsigned,unsigned,unsigned>      split_metrics_t;

            void                                                setBitAt(unsigned leaf_index);
            void                                                addSplit(const Split & other);
            
            bool                                                isCompatible(const Split & other) const;
            bool                                                conflictsWith(const Split & other) const;

            std::string                                         createPatternRepresentation() const;
            split_metrics_t                                     getSplitMetrics() const;

        private:

            split_t                                             _bits;
            unsigned                                            _bits_per_unit;
            unsigned                                            _nleaves;
            const split_unit_t                                  _unity;

        public:

            typedef std::shared_ptr< Split >                    SharedPtr;
    };

inline Split::Split() : _unity(1)
    {
    _nleaves = 0;
    //std::cout << "CHAR_BIT = " << CHAR_BIT << std::endl;
    //std::cout << "sizeof(Split::split_unit_t) = " << sizeof(Split::split_unit_t) << std::endl;
    _bits_per_unit = (CHAR_BIT)*sizeof(Split::split_unit_t);
    clear();
    //std::cout << "Constructing a Split" << std::endl;
    }

inline Split::Split(const Split & other) : _unity(1)
    {
    _nleaves = other._nleaves;
    _bits_per_unit = (CHAR_BIT)*sizeof(Split::split_unit_t);
    _bits = other._bits;
    }

inline Split::~Split()
    {
    //std::cout << "Destroying a Split" << std::endl;
    }

inline void Split::clear()
    {
    for (auto & u : _bits)
        {
        u = 0L;
        }
    }

inline Split & Split::operator=(const Split & other)
    {
    _nleaves = other._nleaves;
    _bits = other._bits;
    return *this;
    }

inline bool Split::operator==(const Split & other) const
    {
    return (_bits == other._bits);
    }

inline bool Split::operator<(const Split & other) const
    {
    assert(_bits.size() == other._bits.size());
    return (_bits < other._bits);
    }

inline void Split::resize(unsigned nleaves)
    {
    _nleaves = nleaves;
    unsigned nunits = 1 + ((nleaves - 1)/_bits_per_unit);
    _bits.resize(nunits);
    clear();
    }

inline void Split::setBitAt(unsigned leaf_index)
    {
    unsigned unit_index = leaf_index/_bits_per_unit;
    unsigned bit_index = leaf_index - unit_index*_bits_per_unit;
    split_unit_t bit_to_set = _unity << bit_index;
    _bits[unit_index] |= bit_to_set;
    }

inline void Split::addSplit(const Split & other)
    {
    unsigned nunits = (unsigned)_bits.size();
    assert(nunits == other._bits.size());
    for (unsigned i = 0; i < nunits; ++i)
        {
        _bits[i] |= other._bits[i];
        }
    }
    
// Returns true if this split and `other' are compatible. The two splits a and b are compatible if a & b is nonzero
// and also not equal to either a or b. For example, these two splits
//
//  split a: -***---*--
//  split b: ----***--*
//    a & b: ----------
//
//  are compatible, because a & b = 0. The two splits below are also compatible because a & b == b:
//
//  split a: -****-*---
//  split b: --**--*--- <-
//    a & b: --**--*--- <-
//
//  These two splits, on the other hand, are not compatible because a & b != 0 and is not equal to either a or b:
//
//  split a: -***---*--
//  split b: ---***---*
//    a & b: ---*------
//
inline bool Split::isCompatible(const Split & other) const
    {
    for (unsigned i = 0; i < _bits.size(); ++i)
        {
        split_unit_t a = _bits[i];
        split_unit_t b = other._bits[i];
        split_unit_t a_and_b = (a & b);
        bool equals_a   = (a_and_b == a);
        bool equals_b   = (a_and_b == b);
        if (a_and_b && !(equals_a || equals_b))
            {
            // A failure of any unit to be compatible makes the entire split incompatible
            return false;
            }
        }

    // None of the units were incompatible, so that means the splits are compatible
    return true;
    }

inline bool Split::conflictsWith(const Split & other) const
    {
    return !isCompatible(other);
    }

inline std::string Split::createPatternRepresentation() const
    {
    std::string s;
    unsigned ntax_added = 0;
    for (unsigned i = 0; i < _bits.size(); ++i)
        {
        for (unsigned j = 0; j < _bits_per_unit; ++j)
            {
            split_unit_t bitmask = ((split_unit_t)1 << j);
            bool bit_is_set = ((_bits[i] & bitmask) > (split_unit_t)0);
            if (bit_is_set)
                s += '*';
            else
                s += '-';
            if (++ntax_added == _nleaves)
                break;
            }
        }
    return s;
    }


    }
