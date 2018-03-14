#pragma once

#include <memory>
#include <boost/format.hpp>
#include "lot.hpp"
#include "data.hpp"
#include "tree.hpp"
#include "likelihood.hpp"
#include "tree_manip.hpp"
#include "gamma_shape_updater.hpp"
#include "statefreq_updater.hpp"
#include "exchangeability_updater.hpp"

namespace strom
    {

    class Likelihood;

    class Chain
        {
        friend class Likelihood;

        public:
                                                    Chain();
                                                    ~Chain();

            void                                    clear();

            void                                    startTuning();
            void                                    stopTuning();

            void                                    start();
            void                                    stop();
            void                                    nextStep(int iteration, unsigned sampling_freq);

            void                                    setTreeManip(TreeManip::SharedPtr tm);
            void                                    setLikelihood(typename Likelihood::SharedPtr likelihood);
            void                                    setLot(typename Lot::SharedPtr lot);
            void                                    setSliceLambda(double slicelambda);

            void                                    setHeatingPower(double p);
            double                                  getHeatingPower() const;

            std::vector<std::string>                getLambdaNames() const;
            std::vector<double>                     getLambdas() const;
            void                                    setLambdas(std::vector<double> & v);

            double                                  calcLogLikelihood() const;
            double                                  calcLogJointPrior() const;

            typedef std::shared_ptr< Chain >        SharedPtr;

        private:

            Likelihood::SharedPtr               _likelihood;
            TreeManip::SharedPtr                _tree_manipulator;

            GammaShapeUpdater::SharedPtr        _shape_updater;
            StateFreqUpdater::SharedPtr         _statefreq_updater;
            ExchangeabilityUpdater::SharedPtr   _exchangeability_updater;

            double                              _heating_power;
            double                              _log_likelihood;
            int                                 _iter;

        };

inline Chain::Chain()
    {
    //std::cout << "Chain being created" << std::endl;
    clear();
    }

inline Chain::~Chain()
    {
    //std::cout << "Chain being destroyed" << std::endl;
    }

inline void Chain::setHeatingPower(double p)
    {
    _heating_power = p;
    _shape_updater->setHeatingPower(p);
    _statefreq_updater->setHeatingPower(p);
    _exchangeability_updater->setHeatingPower(p);
    }

inline double Chain::getHeatingPower() const
    {
    return _heating_power;
    }

//TODO shouldn't this just be the updater names?
inline std::vector<std::string> Chain::getLambdaNames() const
    {
    std::vector<std::string> v;
    v.push_back("Gamma shape lambda");
    return v;
    }

inline std::vector<double> Chain::getLambdas() const
    {
    std::vector<double> v;
    v.push_back(_shape_updater->getLambda());
    v.push_back(_statefreq_updater->getLambda());
    v.push_back(_exchangeability_updater->getLambda());
    return v;
    }

inline void Chain::setLambdas(std::vector<double> & v)
    {
    assert(v.size() >= 3);
    _shape_updater->setLambda(v[0]);
    _statefreq_updater->setLambda(v[1]);
    _exchangeability_updater->setLambda(v[2]);
    }

inline void Chain::startTuning()
    {
    _shape_updater->setTuning(true);
    _statefreq_updater->setTuning(true);
    _exchangeability_updater->setTuning(true);
    }

inline void Chain::stopTuning()
    {
    _shape_updater->setTuning(false);
    _statefreq_updater->setTuning(false);
    _exchangeability_updater->setTuning(false);
    }

inline void Chain::setTreeManip(TreeManip::SharedPtr tm)
    {
    _tree_manipulator = tm;
    _shape_updater->setTreeManip(_tree_manipulator);
    _statefreq_updater->setTreeManip(_tree_manipulator);
    _exchangeability_updater->setTreeManip(_tree_manipulator);
    }

inline void Chain::setLikelihood(typename Likelihood::SharedPtr likelihood)
    {
    _likelihood = likelihood;
    _shape_updater->setLikelihood(likelihood);
    _statefreq_updater->setLikelihood(likelihood);
    _exchangeability_updater->setLikelihood(likelihood);
    }

inline void Chain::setLot(typename Lot::SharedPtr lot)
    {
    _shape_updater->setLot(lot);
    _statefreq_updater->setLot(lot);
    _exchangeability_updater->setLot(lot);
    }

inline void Chain::clear()
    {
    _log_likelihood = 0.0;

    std::vector<double> edge_length_prior_parameters(2, 1.0);
    edge_length_prior_parameters[1] = 10.0;

    _shape_updater.reset(new GammaShapeUpdater);
    _shape_updater->setLambda(1.0);
    _shape_updater->setTargetAcceptanceRate(0.3);
    _shape_updater->setPriorParameters(std::vector<double>(2, 1.0));

    _statefreq_updater.reset(new StateFreqUpdater);
    _statefreq_updater->setLambda(0.001);
    _statefreq_updater->setTargetAcceptanceRate(0.3);
    _statefreq_updater->setPriorParameters(std::vector<double>(4, 1.0));

    _exchangeability_updater.reset(new ExchangeabilityUpdater);
    _exchangeability_updater->setLambda(0.001);
    _exchangeability_updater->setTargetAcceptanceRate(0.3);
    _exchangeability_updater->setPriorParameters(std::vector<double>(6, 1.0));

    setHeatingPower(1.0);
    startTuning();
    }

inline void Chain::start()
    {
    _shape_updater->pullCurrentStateFromModel();
    _statefreq_updater->pullCurrentStateFromModel();
    _exchangeability_updater->pullCurrentStateFromModel();
    _log_likelihood = calcLogLikelihood();
    }

inline void Chain::stop()
    {
    }

inline double Chain::calcLogLikelihood() const
    {
    return _shape_updater->calcLogLikelihood();
    }

inline double Chain::calcLogJointPrior() const
    {
    double lnP = 0.0;
    lnP += _shape_updater->calcLogPrior();
    lnP += _statefreq_updater->calcLogPrior();
    lnP += _exchangeability_updater->calcLogPrior();
    return lnP;
    }

inline void Chain::nextStep(int iteration, unsigned sampling_freq)
    {
    GTRModel::SharedPtr gtr = _likelihood->getModel();
    if (gtr->getGammaNCateg() > 1)
        _log_likelihood = _shape_updater->update(_log_likelihood);
    _log_likelihood = _statefreq_updater->update(_log_likelihood);
    _log_likelihood = _exchangeability_updater->update(_log_likelihood);
    }

}
