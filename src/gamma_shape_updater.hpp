#pragma once

#include "updater.hpp"

namespace strom
{

    class GammaShapeUpdater : public Updater
    {
        public:

                                        GammaShapeUpdater();
                                        ~GammaShapeUpdater();

            virtual void                clear();
            virtual double              calcLogPrior() const;
            virtual void                pullCurrentStateFromModel();
            virtual void                pushCurrentStateToModel() const;
            virtual void                proposeNewState();
            virtual void                revert();

            double                      getCurrentPoint() const;

            double                      calcGammaShapePrior() const;

        private:

            double                      _prev_point;
            double                      _curr_point;

        public:
            typedef std::shared_ptr< GammaShapeUpdater > SharedPtr;
    };

inline GammaShapeUpdater::GammaShapeUpdater()
    {
    //std::cout << "GammaShapeUpdater being created" << std::endl;
    _name = "shape";
    clear();
    }

inline GammaShapeUpdater::~GammaShapeUpdater()
    {
    //std::cout << "GammaShapeUpdater being destroyed" << std::endl;
    }

inline double GammaShapeUpdater::getCurrentPoint() const
    {
    return _curr_point;
    }

inline void GammaShapeUpdater::clear()
    {
    Updater::clear();
    _prev_point             = 0.0;
    _curr_point             = 0.0;
    reset();
    }

inline double GammaShapeUpdater::calcLogPrior() const
    {
    assert(_prior_parameters.size() == 2);
    double prior_a = _prior_parameters[0];
    double prior_b = _prior_parameters[1];

    double log_prior = 0.0;
    log_prior += (prior_a - 1.0)*log(_curr_point);
    log_prior -= _curr_point/prior_b;
    log_prior -= prior_b*log(prior_a);
    log_prior -= std::lgamma(prior_a);
    return log_prior;
    }

inline void GammaShapeUpdater::pullCurrentStateFromModel()
    {
    GTRModel::SharedPtr gtr = _likelihood->getModel();
    _curr_point = gtr->getGammaShape();
    }

inline void GammaShapeUpdater::pushCurrentStateToModel() const
    {
    GTRModel::SharedPtr gtr = _likelihood->getModel();
    gtr->setGammaShape(_curr_point);
    }

inline void GammaShapeUpdater::proposeNewState()
    {
    // Save copy of _curr_point in case revert is necessary.
    _prev_point = _curr_point;

    // Let _curr_point be proposed value
    double m = exp(_lambda*(_lot->uniform() - 0.5));
    _curr_point = m*_prev_point;

    //std::cerr << boost::str(boost::format("_curr_point    = %12.5f") % _curr_point) << std::endl;

    // calculate log of Hastings ratio
    _log_hastings_ratio = log(m);
    }

inline void GammaShapeUpdater::revert()
    {
    _curr_point = _prev_point;
    }

inline double GammaShapeUpdater::calcGammaShapePrior() const
    {
    return calcLogPrior();
    }

}
