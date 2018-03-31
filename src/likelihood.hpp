#pragma once

#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "libhmsbeagle/beagle.h"
#include "data.hpp"
#include "gtr_model.hpp"
#include "xstrom.hpp"

namespace strom {

class Likelihood
    {
    public:
                                    Likelihood();
                                    ~Likelihood();

        void                        useStoredData(bool using_data);
        std::string                 availableResources();

        double                      calcLogLikelihood(typename Tree::SharedPtr t);
        void                        setModel(GTRModel::SharedPtr model);
        GTRModel::SharedPtr         getModel();

        void                        setData(Data::SharedPtr d);
        Data::SharedPtr             getData();

    private:

        void                        initBeagleLib();
        void                        setTipStates();
        void                        setPatternWeights();
        void                        setDiscreteGammaShape();
        void                        setModelRateMatrix();
        void                        defineOperations(typename Tree::SharedPtr t);
        void                        updateTransitionMatrices();
        void                        calculatePartials();

        int                         _instance;
        std::map<int, std::string>  _beagle_error;
        std::vector<int>            _operations;
        std::vector<int>            _pmatrix_index;
        std::vector<double>         _edge_lengths;

        Data::SharedPtr             _data;
        unsigned                    _ntaxa;
        unsigned                    _nstates;
        unsigned                    _npatterns;
        bool                        _rooted;
        bool                        _prefer_gpu;

        bool                        _using_data;

        GTRModel::SharedPtr         _model;

    public:
        typedef std::shared_ptr< Likelihood > SharedPtr;
    };

inline Likelihood::Likelihood()
    {
    _instance   = -1;
    _ntaxa      = 0;
    _nstates    = 0;
    _npatterns  = 0;
    _rooted     = false;
    _prefer_gpu = false;
    _using_data = true;

    _model = GTRModel::SharedPtr(new GTRModel());

    // store BeagleLib error codes so that useful
    // error messages may be provided to the user
    _beagle_error[0]  = std::string("success");
    _beagle_error[-1] = std::string("unspecified error");
    _beagle_error[-2] = std::string("not enough memory could be allocated");
    _beagle_error[-3] = std::string("unspecified exception");
    _beagle_error[-4] = std::string("the instance index is out of range, or the instance has not been created");
    _beagle_error[-5] = std::string("one of the indices specified exceeded the range of the array");
    _beagle_error[-6] = std::string("no resource matches requirements");
    _beagle_error[-7] = std::string("no implementation matches requirements");
    _beagle_error[-8] = std::string("floating-point range exceeded");
    }

inline Likelihood::~Likelihood()
    {
    if (_instance >= 0)
        {
        int code = beagleFinalizeInstance(_instance);
        if (code != 0)
            std::cerr << boost::str(boost::format("Likelihood destructor failed to finalize BeagleLib instance. BeagleLib error code was %d (%s).") % code % _beagle_error[code]) << std::endl;
        }
    }

inline std::string Likelihood::availableResources()
    {
    BeagleResourceList * rsrcList = beagleGetResourceList();
	std::string s;
    for (int i = 0; i < rsrcList->length; ++i)
        {
        std::string desc = rsrcList->list[i].description;
        boost::trim(desc);
        if (desc.size() > 0)
            s += boost::str(boost::format("%d: %s (%s)\n") % i % rsrcList->list[i].name % desc);
        else
            s += boost::str(boost::format("%d: %s\n") % i % rsrcList->list[i].name);
        }
    return s;
    }

inline Data::SharedPtr Likelihood::getData()
    {
    return _data;
    }

inline void Likelihood::setData(Data::SharedPtr data)
    {
    _data = data;
    if (_instance >= 0)
        {
        // initBeagleLib function was previously called, so
        // finalize existing instance and create new one
        int code = beagleFinalizeInstance(_instance);
        if (code != 0)
            throw XStrom(boost::str(boost::format("Likelihood setData function failed to finalize BeagleLib instance. BeagleLib error code was %d (%s).") % code % _beagle_error[code]));
        _instance = -1;
        assert(_ntaxa > 0 && _nstates > 0 && _npatterns > 0);
        initBeagleLib();
        }
    }

inline void Likelihood::useStoredData(bool using_data)
    {
    _using_data = using_data;
    }

inline void Likelihood::initBeagleLib()
    {
    // a non-operation ("no-op") if a valid instance has already been created
    if (_instance >= 0)
        return;

    assert(_data);

    _ntaxa      = _data->getNumTaxa();
    _npatterns  = _data->getNumPatterns();
    _nstates    = 4;
    _rooted     = false;
    _prefer_gpu = false;

    std::cout << "Sequence length:    " << _data->getSeqLen() << std::endl;
    std::cout << "Number of taxa:     " << _ntaxa << std::endl;
    std::cout << "Number of patterns: " << _npatterns << std::endl;

    unsigned num_internals        = (_rooted ? (_ntaxa - 1) : (_ntaxa - 2));
    unsigned num_transition_probs = (_rooted ? (2*_ntaxa - 2) : (2*_ntaxa - 3));

    long requirementFlags = 0;
    requirementFlags |= BEAGLE_FLAG_PRECISION_DOUBLE;

    long preferenceFlags = 0;
    if (_prefer_gpu)
        preferenceFlags |= BEAGLE_FLAG_PROCESSOR_GPU;
    else
        preferenceFlags |= BEAGLE_FLAG_PROCESSOR_CPU;

    requirementFlags |= BEAGLE_FLAG_SCALING_MANUAL;

    BeagleInstanceDetails instance_details;
    _instance = beagleCreateInstance(
         _ntaxa,                    // tips
         num_internals,             // partials
         _ntaxa,                    // sequences
         _nstates,                  // states
         _npatterns,                // patterns
         1,                         // models
         num_transition_probs,      // transition matrices
         _model->_num_categ,        // rate categories
         num_internals + 1,         // scale buffers
         NULL,                      // resource restrictions
         0,                         // length of resource list
         preferenceFlags,           // preferred flags
         requirementFlags,          // required flags
         &instance_details);        // pointer for details

    if (_instance < 0)
        {
        // beagleCreateInstance returns one of the following:
        //   valid instance (0, 1, 2, ...)
        //   error code (negative integer)
        throw XStrom(boost::str(boost::format("Likelihood init function failed to create Likelihood instance (BeagleLib error code was %d)") % _beagle_error[_instance]));
        }

    setTipStates();
    setPatternWeights();

    //std::cout << boost::str(boost::format("BeagleLib instance (%d) created.") % _instance) << std::endl;
    }

inline void Likelihood::setTipStates()
    {
    assert(_data);

    const Data::data_matrix_t & data_matrix = _data->getDataMatrix();
    unsigned i = 0;

    typedef const std::vector<int> & ref_int_vect_t;
    for (ref_int_vect_t v : data_matrix)
        {
        int code = beagleSetTipStates(
            _instance,      // Instance number
            i,              // Index of destination compactBuffer
            &v[0]);         // Pointer to compact states vector

        if (code != 0)
            throw XStrom(boost::str(boost::format("failed to set tip state for taxon %d (\"%s\"; BeagleLib error code was %d)") % (i+1) % _data->getTaxonNames()[i] % code % _beagle_error[code]));
        ++i;
        }
    }

inline void Likelihood::setPatternWeights()
    {
    assert(_data);

    int code = 0;
    const Data::pattern_counts_t & v = _data->getPatternCounts();
    if (v.empty())
        throw XStrom(boost::str(boost::format("failed to set pattern weights because data matrix has empty pattern count vector") % code));

    code = beagleSetPatternWeights(
       _instance,     // instance number
       &v[0]);        // vector of pattern counts

    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to set pattern weights. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));
    }

inline void Likelihood::setDiscreteGammaShape()
    {
    int code = _model->setBeagleAmongSiteRateVariationRates(_instance);
    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to set category rates. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));

    code = _model->setBeagleAmongSiteRateVariationProbs(_instance);
    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to set category probabilities. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));
    }

inline void Likelihood::setModelRateMatrix()
    {
    int code = _model->setBeagleStateFrequencies(_instance);
    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to set state frequencies. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));

    code = _model->setBeagleEigenDecomposition(_instance);
    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to set eigen decomposition. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));

    code = _model->setBeagleAmongSiteRateVariationRates(_instance);
    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to set among-site rate variation rates. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));

    code = _model->setBeagleAmongSiteRateVariationProbs(_instance);
    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to set among-site rate variation weights. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));
    }

inline void Likelihood::defineOperations(typename Tree::SharedPtr t)
    {
    _operations.clear();
    _pmatrix_index.clear();
    _edge_lengths.clear();

    //BOOST_REVERSE_FOREACH(T * nd, t->_preorder)
    for (auto nd : boost::adaptors::reverse(t->_preorder))
        {
        assert(nd->_number >= 0);
        if (!nd->_left_child)
            {
            // This is a leaf
            _pmatrix_index.push_back(nd->_number);
            _edge_lengths.push_back(nd->_edge_length);
            }
        else
            {
            // This is an internal node
            _pmatrix_index.push_back(nd->_number);
            _edge_lengths.push_back(nd->_edge_length);

            // Internal nodes have partials to be calculated, so define
            // an operation to compute the partials for this node

            // 1. destination partial to be calculated
            int partial = nd->_number;
            _operations.push_back(partial);

            // 2. destination scaling buffer index to write to
            //oldway
            //_operations.push_back(BEAGLE_OP_NONE);
            //newway
            int scaler = nd->_number - _ntaxa + 1;
            _operations.push_back(scaler);

            // 3. destination scaling buffer index to read from
            _operations.push_back(BEAGLE_OP_NONE);

            // 4. left child partial index
            partial = nd->_left_child->_number;
            _operations.push_back(partial);

            // 5. left child transition matrix index
            int tmatrix = nd->_left_child->_number;
            _operations.push_back(tmatrix);

            // 6. right child partial index
            assert(nd->_left_child);
            assert(nd->_left_child->_right_sib);
            partial = nd->_left_child->_right_sib->_number;
            _operations.push_back(partial); // assumes binary tree

            // 7. right child transition matrix index
            tmatrix = nd->_left_child->_right_sib->_number;
            _operations.push_back(tmatrix);
        }
    }

    // if tree is unrooted and thus "rooted" at a leaf, need to
    // use the transition matrix associated with the leaf
    if (!t->_is_rooted)
        _pmatrix_index[_pmatrix_index.size()-1] = t->_root->_number;
    }

inline void Likelihood::updateTransitionMatrices()
    {
    int code = beagleUpdateTransitionMatrices(
        _instance,                      // Instance number
        0,                              // Index of eigen-decomposition buffer
        &_pmatrix_index[0],             // transition probability matrices to update
        NULL,                           // first derivative matrices to update
        NULL,                           // second derivative matrices to update
        &_edge_lengths[0],              // List of edge lengths
        (int)_pmatrix_index.size());    // Length of lists

    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to update transition matrices. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));
    }

inline void Likelihood::calculatePartials()
    {
    int code = beagleResetScaleFactors(_instance, 0);
    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to reset scale factors in calculatePartials. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));

    // Calculate or queue for calculation partials using a list of operations
    int totalOperations = (int)(_operations.size()/7);
    code = beagleUpdatePartials(
        _instance,                              // Instance number
        (BeagleOperation *) &_operations[0],    // BeagleOperation list specifying operations
        totalOperations,                        // Number of operations
        0);                                     // Index number of scaleBuffer to store accumulated factors

    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to update partials. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));
    }

inline double Likelihood::calcLogLikelihood(typename Tree::SharedPtr t)
    {
    if (t->_is_rooted)
        throw XStrom("can only compute likelihoods for unrooted trees currently");

    if (!_data)
        throw XStrom("must call setData before calcLogLikelihood");

    initBeagleLib(); // this is a no-op if a valid instance already exists

    // Assuming "root" is leaf 0
    assert(t->_root->_number == 0 && t->_root->_left_child == t->_preorder[0] && !t->_preorder[0]->_right_sib);

    // Assuming there are as many transition matrices as there are edge lengths
    assert(_pmatrix_index.size() == _edge_lengths.size());

    setModelRateMatrix();
    setDiscreteGammaShape();
    defineOperations(t);
    updateTransitionMatrices();
    calculatePartials();

    // The beagleCalculateEdgeLogLikelihoods function integrates a list of partials
    // at a parent and child node with respect to a set of partials-weights and
    // state frequencies to return the log likelihood and first and second derivative sums
    int stateFrequencyIndex  = 0;
    int categoryWeightsIndex = 0;
    int cumulativeScalingIndex = 0;
    double log_likelihood = 0.0;

    // index_focal_child is the root node
    int index_focal_child  = t->_root->_number;

    // index_focal_parent is the only child of root node
    int index_focal_parent = t->_preorder[0]->_number;

    int code = beagleCalculateEdgeLogLikelihoods(
        _instance,                  // instance number
        &index_focal_parent,        // indices of parent partialsBuffers
        &index_focal_child,         // indices of child partialsBuffers
        &index_focal_child,         // transition probability matrices for this edge
        NULL,                       // first derivative matrices
        NULL,                       // second derivative matrices
        &categoryWeightsIndex,      // weights to apply to each partialsBuffer
        &stateFrequencyIndex,       // state frequencies for each partialsBuffer
        &cumulativeScalingIndex,    // scaleBuffers containing accumulated factors
        1,                          // Number of partialsBuffer
        &log_likelihood,            // destination for log likelihood
        NULL,                       // destination for first derivative
        NULL);                      // destination for second derivative

    if (code != 0)
        throw XStrom(boost::str(boost::format("failed to calculate edge logLikelihoods in CalcLogLikelihood. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));

    return log_likelihood;
    }

inline void Likelihood::setModel(GTRModel::SharedPtr model)
    {
    _model = model;
    if (_instance >= 0)
        {
        // init function was previously called, so set the model and create new BeagleLib instance
        int code = beagleFinalizeInstance(_instance);
        if (code != 0)
            throw XStrom(boost::str(boost::format("Likelihood setModel function failed to finalize BeagleLib instance. BeagleLib error code was %d (%s).") % code % _beagle_error[code]));
        _instance = -1;
        assert(_ntaxa > 0 && _nstates > 0 && _npatterns > 0);
        initBeagleLib();
        }
    }

inline GTRModel::SharedPtr Likelihood::getModel()
    {
    return _model;
    }

}
