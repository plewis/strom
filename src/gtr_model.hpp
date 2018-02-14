#pragma once

#include <algorithm>
#include <vector>
#include "libhmsbeagle/beagle.h"
#include <boost/math/distributions/gamma.hpp>
#include <Eigen/Dense>

namespace strom
    {

    class GTRModel
        {
        friend class Likelihood;

        public:
            typedef Eigen::Matrix<double, 4, 4, Eigen::RowMajor>    EigenMatrix4d;
            typedef Eigen::Vector4d                                 EigenVector4d;


                                        GTRModel();
                                        ~GTRModel();

            std::string                 describeModel() const;

            void                        setGammaShape(double shape);
            void                        setGammaNCateg(unsigned ncateg);
            void                        setExchangeabilities(const std::vector<double> & exchangeabilities);
            void                        setStateFreqs(const std::vector<double> & state_frequencies);
            void                        setExchangeabilitiesAndStateFreqs(const std::vector<double> & exchangeabilities, const std::vector<double> & state_frequencies);

            int                         setBeagleEigenDecomposition(int beagle_instance);
            int                         setBeagleStateFrequencies(int beagle_instance);
            int                         setBeagleAmongSiteRateVariationRates(int beagle_instance);
            int                         setBeagleAmongSiteRateVariationProbs(int beagle_instance);

                                        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            typedef boost::shared_ptr< GTRModel > SharedPtr;

        private:

            void                        clear();
            void                        recalcRateMatrix();
            void                        recalcGammaRates();

            // substitution model specification
            std::vector<double>         _state_freqs;
            std::vector<double>         _exchangeabilities;

            // for computing eigenvectors/eigenvalues
            EigenMatrix4d               _sqrtPi;
            EigenMatrix4d               _sqrtPiInv;
            EigenMatrix4d               _qmatrix;
            EigenMatrix4d               _eigenvectors;
            EigenMatrix4d               _inverse_eigenvectors;
            EigenVector4d               _eigenvalues;

            // among-site rate heterogeneity specification
            unsigned                    _num_categ;
            double                      _gamma_shape;
            std::vector<double>         _relative_rates;
            std::vector<double>         _categ_boundaries;
            std::vector<double>         _rate_probs;

        };

inline GTRModel::GTRModel()
    {
    clear();
    }

inline GTRModel::~GTRModel()
    {
    }

inline void GTRModel::clear()
    {
    _num_categ = 1;
    _gamma_shape = 0.5;

    // Set up GTR rate matrix representing the JC69 model by default
    setExchangeabilitiesAndStateFreqs({1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, {0.25, 0.25, 0.25, 0.25});

    _relative_rates.assign(1, 1.0);
    _categ_boundaries.assign(1, 0.0);
    _rate_probs.assign(1, 1.0);

    recalcRateMatrix();
    }

inline std::string GTRModel::describeModel() const
    {
    std::string s;
    s += "\n----------------- GTR Model Info ------------------\n";
    s += boost::str(boost::format("State frequencies:   \n  piA = %g\n  piC = %g\n  piG = %g\n  piT = %g") % _state_freqs[0] % _state_freqs[1] % _state_freqs[2] % _state_freqs[3]);
    s += boost::str(boost::format("\nRelative rates:    \n  rAC = %g\n  rAG = %g\n  rAT = %g\n  rCG = %g\n  rCT = %g\n  rGT = %g") % _exchangeabilities[0] % _exchangeabilities[1] % _exchangeabilities[2] % _exchangeabilities[3] % _exchangeabilities[4] % _exchangeabilities[5]);
    s += boost::str(boost::format("\nRate categories:   \n  %d") % _num_categ);
    s += boost::str(boost::format("\nGamma shape:       \n  %g") % _gamma_shape);
    if (_num_categ > 1)
        {
        s += "\nCategory boundaries and relative rate means:";
        s += boost::str(boost::format("\n%12s %12s %12s %12s") % "category" %  "lower" % "upper" % "mean rate");
        for (unsigned i = 1; i < _num_categ; ++i)
            s += boost::str(boost::format("\n%12d %12.5f %12.5f %12.5f") % i % _categ_boundaries[i-1] % _categ_boundaries[i] % _relative_rates[i-1]);
        s += boost::str(boost::format("\n%12d %12.5f %12s %12.5f") % _num_categ % _categ_boundaries[_num_categ-1] % "infinity" % _relative_rates[_num_categ-1]);
        }
    s += "\n---------------------------------------------------\n";
    return s;
    }

inline void GTRModel::setGammaNCateg(unsigned ncateg)
    {
    if (ncateg < 1)
        throw XStrom(boost::str(boost::format("number of categories used for among-site rate variation must be greater than zero but the value %d was supplied") % ncateg));
    _num_categ = ncateg;
    recalcGammaRates();
    }

inline void GTRModel::setGammaShape(double shape)
    {
    if (shape <= 0.0)
        throw XStrom(boost::str(boost::format("gamma shape must be greater than zero but the value %.5f was supplied") % shape));
    _gamma_shape = shape;
    recalcGammaRates();
    }

inline void GTRModel::setExchangeabilities(const std::vector<double> & exchangeabilities)
    {
    _exchangeabilities.assign(exchangeabilities.begin(), exchangeabilities.end());
    recalcRateMatrix();
    }

inline void GTRModel::setStateFreqs(const std::vector<double> & state_frequencies)
    {
    _state_freqs.assign(state_frequencies.begin(), state_frequencies.end());

    Eigen::Map<const Eigen::Array4d> tmp(&_state_freqs[0]);
    _sqrtPi = tmp.sqrt().matrix().asDiagonal();
    _sqrtPiInv = _sqrtPi.inverse();

    recalcRateMatrix();
    }

inline void GTRModel::setExchangeabilitiesAndStateFreqs(const std::vector<double> & exchangeabilities, const std::vector<double> & state_frequencies)
    {
    _exchangeabilities.assign(exchangeabilities.begin(), exchangeabilities.end());

    _state_freqs.assign(state_frequencies.begin(), state_frequencies.end());

    Eigen::Map<const Eigen::Array4d> tmp(&_state_freqs[0]);
    _sqrtPi = tmp.sqrt().matrix().asDiagonal();
    _sqrtPiInv = _sqrtPi.inverse();

    recalcRateMatrix();
    }

inline void GTRModel::recalcRateMatrix()
    {
    double piA = _state_freqs[0];
    double piC = _state_freqs[1];
    double piG = _state_freqs[2];
    double piT = _state_freqs[3];

    double rAC = _exchangeabilities[0];
    double rAG = _exchangeabilities[1];
    double rAT = _exchangeabilities[2];
    double rCG = _exchangeabilities[3];
    double rCT = _exchangeabilities[4];
    double rGT = _exchangeabilities[5];

    double inverse_scaling_factor = piA*(rAC*piC + rAG*piG + rAT*piT) + piC*(rAC*piA + rCG*piG + rCT*piT) + piG*(rAG*piA + rCG*piC + rGT*piT) + piT*(rAT*piA + rCT*piC + rGT*piG);
    double scaling_factor = 1.0/inverse_scaling_factor;

    _qmatrix(0,0) = -scaling_factor*(rAC*piC + rAG*piG + rAT*piT);
    _qmatrix(0,1) = scaling_factor*rAC*piC;
    _qmatrix(0,2) = scaling_factor*rAG*piG;
    _qmatrix(0,3) = scaling_factor*rAT*piT;

    _qmatrix(1,0) = scaling_factor*rAC*piA;
    _qmatrix(1,1) = -scaling_factor*(rAC*piA + rCG*piG + rCT*piT);
    _qmatrix(1,2) = scaling_factor*rCG*piG;
    _qmatrix(1,3) = scaling_factor*rCT*piT;

    _qmatrix(2,0) = scaling_factor*rAG*piA;
    _qmatrix(2,1) = scaling_factor*rCG*piC;
    _qmatrix(2,2) = -scaling_factor*(rAG*piA + rCG*piC + rGT*piT);
    _qmatrix(2,3) = scaling_factor*rGT*piT;

    _qmatrix(3,0) = scaling_factor*rAT*piA;
    _qmatrix(3,1) = scaling_factor*rCT*piC;
    _qmatrix(3,2) = scaling_factor*rGT*piG;
    _qmatrix(3,3) = -scaling_factor*(rAT*piA + rCT*piC + rGT*piG);

    // S is a symmetric matrix
    EigenMatrix4d S = EigenMatrix4d(_sqrtPi*_qmatrix*_sqrtPiInv);

    // Can use efficient eigensystem solver because S is symmetric
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(S);
    if (solver.info() != Eigen::Success)
        {
        throw XStrom("Error in the calculation of eigenvectors and eigenvalues of the GTR rate matrix");
        }

    _eigenvectors           = _sqrtPiInv*solver.eigenvectors();
    _inverse_eigenvectors   = solver.eigenvectors().transpose()*_sqrtPi;
    _eigenvalues            = solver.eigenvalues();
    }

inline void GTRModel::recalcGammaRates()
    {
    assert(_num_categ > 0);
    _relative_rates.assign(_num_categ, 1.0);
    _categ_boundaries.assign(_num_categ, 0.0);
    _rate_probs.assign(_num_categ, 1.0/_num_categ);

    if (_num_categ == 1)
        return;

    assert(_gamma_shape > 0.0);
    double alpha = _gamma_shape;
    double beta = 1.0/_gamma_shape;

    boost::math::gamma_distribution<> my_gamma(_gamma_shape, 1.0/_gamma_shape);
    boost::math::gamma_distribution<> my_gamma_plus(_gamma_shape + 1.0, 1.0/_gamma_shape);

    double cum_upper        = 0.0;
    double cum_upper_plus   = 0.0;
    double upper            = 0.0;
    double cum_prob         = 0.0;
    for (unsigned i = 1; i <= _num_categ; ++i)
        {
        double lower                = upper;
        double cum_lower_plus       = cum_upper_plus;
        double cum_lower            = cum_upper;
        cum_prob                    += _rate_probs[i-1];

        if (i < _num_categ)
            {
            upper                   = boost::math::quantile(my_gamma, cum_prob);
            cum_upper_plus          = boost::math::cdf(my_gamma_plus, upper);
            cum_upper               = boost::math::cdf(my_gamma, upper);
            }
        else
            {
            cum_upper_plus          = 1.0;
            cum_upper               = 1.0;
            }

        double numer                = cum_upper_plus - cum_lower_plus;
        double denom                = cum_upper - cum_lower;
        double r_mean               = (denom > 0.0 ? (alpha*beta*numer/denom) : 0.0);
        _relative_rates[i-1]        = r_mean;
        _categ_boundaries[i-1]      = lower;
        }
    }

inline int GTRModel::setBeagleEigenDecomposition(int beagle_instance)
    {
    int code = beagleSetEigenDecomposition(
        beagle_instance,
        0,
        &_eigenvectors.data()[0],
        &_inverse_eigenvectors.data()[0],
        &_eigenvalues.data()[0]);

    return code;
    }

inline int GTRModel::setBeagleStateFrequencies(int beagle_instance)
    {
    int code = beagleSetStateFrequencies(
         beagle_instance,
         0,
         &_state_freqs[0]);

    return code;
    }

inline int GTRModel::setBeagleAmongSiteRateVariationRates(int beagle_instance)
    {
    int code = beagleSetCategoryRates(
        beagle_instance,
        &_relative_rates[0]);

    return code;
    }

inline int GTRModel::setBeagleAmongSiteRateVariationProbs(int beagle_instance)
    {
    int code = beagleSetCategoryWeights(
        beagle_instance,
        0,
        &_rate_probs[0]);

    return code;
    }

}
