#pragma once

#include <memory>
#include <boost/format.hpp>
#include "lot.hpp"
#include "data.hpp"
#include "tree.hpp"
#include "likelihood.hpp"
#include "tree_manip.hpp"
#include "gamma_shape_updater.hpp"
//#include "statefreq_updater.hpp"
//#include "exchangeability_updater.hpp"
//#include "tree_updater.hpp"
//#include "tree_length_updater.hpp"
//#include "output_manager.hpp"

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
            void                                    nextStep(int iteration);

            void                                    setTreeFromNewick(std::string & newick);
            void                                    setTreeManip(TreeManip::SharedPtr tm);
            void                                    setLikelihood(typename Likelihood::SharedPtr likelihood);
            void                                    setLot(typename Lot::SharedPtr lot);

            TreeManip::SharedPtr                    getTreeManip();
            Model::SharedPtr                        getModel();

            void                                    setHeatingPower(double p);
            double                                  getHeatingPower() const;

            void                                    setChainIndex(unsigned idx);
            double                                  getChainIndex() const;

            std::vector<std::string>                getUpdaterNames() const;
            std::vector<double>                     getAcceptPercentages() const;
            std::vector<double>                     getLambdas() const;
            void                                    setLambdas(std::vector<double> & v);

            double                                  calcLogLikelihood() const;
            double                                  calcLogJointPrior() const;

            typedef std::shared_ptr< Chain >        SharedPtr;

        private:

            Likelihood::SharedPtr               _likelihood;
            TreeManip::SharedPtr                _tree_manipulator;

            GammaShapeUpdater::SharedPtr        _shape_updater;
            //StateFreqUpdater::SharedPtr         _statefreq_updater;
            //ExchangeabilityUpdater::SharedPtr   _exchangeability_updater;
            //TreeUpdater::SharedPtr              _tree_updater;
            //TreeLengthUpdater::SharedPtr         _tree_length_updater;

            unsigned                            _chain_index;
            double                              _heating_power;
            double                              _log_likelihood;
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

inline TreeManip::SharedPtr Chain::getTreeManip()
    {
    return _tree_manipulator;
    }

inline Model::SharedPtr Chain::getModel()
    {
    return _likelihood->getModel();
    }

inline void Chain::setHeatingPower(double p)
    {
    _heating_power = p;
    _shape_updater->setHeatingPower(p);
    //_statefreq_updater->setHeatingPower(p);
    //_exchangeability_updater->setHeatingPower(p);
    //_tree_updater->setHeatingPower(p);
    //_tree_length_updater->setHeatingPower(p);
    }

inline double Chain::getHeatingPower() const
    {
    return _heating_power;
    }

inline void Chain::setChainIndex(unsigned idx)
    {
    _chain_index = idx;
    }

inline double Chain::getChainIndex() const
    {
    return _chain_index;
    }

inline std::vector<std::string> Chain::getUpdaterNames() const
    {
    std::vector<std::string> v;
    v.push_back(_shape_updater->getUpdaterName());
    //v.push_back(_statefreq_updater->getUpdaterName());
    //v.push_back(_exchangeability_updater->getUpdaterName());
    //v.push_back(_tree_updater->getUpdaterName());
    //v.push_back(_tree_length_updater->getUpdaterName());
    return v;
    }

inline std::vector<double> Chain::getAcceptPercentages() const
    {
    std::vector<double> v;
    v.push_back(_shape_updater->getAcceptPct());
    //v.push_back(_statefreq_updater->getAcceptPct());
    //v.push_back(_exchangeability_updater->getAcceptPct());
    //v.push_back(_tree_updater->getAcceptPct());
    //v.push_back(_tree_length_updater->getAcceptPct());
    return v;
    }

inline std::vector<double> Chain::getLambdas() const
    {
    std::vector<double> v;
    v.push_back(_shape_updater->getLambda());
    //v.push_back(_statefreq_updater->getLambda());
    //v.push_back(_exchangeability_updater->getLambda());
    //v.push_back(_tree_updater->getLambda());
    //v.push_back(_tree_length_updater->getLambda());
    return v;
    }

inline void Chain::setLambdas(std::vector<double> & v)
    {
    assert(v.size() >= 4);
    _shape_updater->setLambda(v[0]);
    //_statefreq_updater->setLambda(v[1]);
    //_exchangeability_updater->setLambda(v[2]);
    //_tree_updater->setLambda(v[3]);
    //_tree_length_updater->setLambda(v[3]);
    }

inline void Chain::startTuning()
    {
    _shape_updater->setTuning(true);
    //_statefreq_updater->setTuning(true);
    //_exchangeability_updater->setTuning(true);
    //_tree_updater->setTuning(true);
    //_tree_length_updater->setTuning(true);
    }

inline void Chain::stopTuning()
    {
    _shape_updater->setTuning(false);
    //_statefreq_updater->setTuning(false);
    //_exchangeability_updater->setTuning(false);
    //_tree_updater->setTuning(false);
    //_tree_length_updater->setTuning(false);
    }

inline void Chain::setTreeFromNewick(std::string & newick)
    {
    if (!_tree_manipulator)
        _tree_manipulator.reset(new TreeManip);
    _tree_manipulator->buildFromNewick(newick, false, false);

    _shape_updater->setTreeManip(_tree_manipulator);
    //_statefreq_updater->setTreeManip(_tree_manipulator);
    //_exchangeability_updater->setTreeManip(_tree_manipulator);
    //_tree_updater->setTreeManip(_tree_manipulator);
    //_tree_length_updater->setTreeManip(_tree_manipulator);
    }

inline void Chain::setLikelihood(typename Likelihood::SharedPtr likelihood)
    {
    _likelihood = likelihood;
    _shape_updater->setLikelihood(likelihood);
    //_statefreq_updater->setLikelihood(likelihood);
    //_exchangeability_updater->setLikelihood(likelihood);
    //_tree_updater->setLikelihood(likelihood);
    //_tree_length_updater->setLikelihood(likelihood);
    }

inline void Chain::setLot(typename Lot::SharedPtr lot)
    {
    _shape_updater->setLot(lot);
    //_statefreq_updater->setLot(lot);
    //_exchangeability_updater->setLot(lot);
    //_tree_updater->setLot(lot);
    //_tree_length_updater->setLot(lot);
    }

inline void Chain::clear()
    {
    _log_likelihood = 0.0;

    _shape_updater.reset(new GammaShapeUpdater);
    _shape_updater->setLambda(1.0);
    _shape_updater->setTargetAcceptanceRate(0.3);
    _shape_updater->setPriorParameters({1.0, 1.0});

    //_statefreq_updater.reset(new StateFreqUpdater);
    //_statefreq_updater->setLambda(0.001);
    //_statefreq_updater->setTargetAcceptanceRate(0.3);
    //_statefreq_updater->setPriorParameters({1.0, 1.0, 1.0, 1.0});

    //_exchangeability_updater.reset(new ExchangeabilityUpdater);
    //_exchangeability_updater->setLambda(0.001);
    //_exchangeability_updater->setTargetAcceptanceRate(0.3);
    //_exchangeability_updater->setPriorParameters({1.0, 1.0, 1.0, 1.0, 1.0, 1.0});

    //double tree_length_shape = 1.0;
    //double tree_length_scale = 10.0;
    //double dirichlet_param   = 1.0;

    //_tree_updater.reset(new TreeUpdater);
    //_tree_updater->setLambda(0.2);
    //_tree_updater->setTargetAcceptanceRate(0.3);
    //_tree_updater->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});

    //_tree_length_updater.reset(new TreeLengthUpdater);
    //_tree_length_updater->setLambda(0.2);
    //_tree_length_updater->setTargetAcceptanceRate(0.3);
    //_tree_length_updater->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});

    _chain_index = 0;
    setHeatingPower(1.0);
    startTuning();
    }

inline void Chain::start()
    {
    _shape_updater->pullCurrentStateFromModel();
    //_statefreq_updater->pullCurrentStateFromModel();
    //_exchangeability_updater->pullCurrentStateFromModel();
    //_tree_updater->pullCurrentStateFromModel();
    //_tree_length_updater->pullCurrentStateFromModel();
    _log_likelihood = calcLogLikelihood();

    // Output column headers and first line of output showing starting state (iteration 0)
    std::cout << boost::str(boost::format("%12s %12s %12s %12s %12s\n")
        % "iteration"
        % "lnLike"
        % "lnPrior"
        % "shape"
        % "accept");
    double lnP = calcLogJointPrior();
    std::cout << boost::str(boost::format("%12d %12.5f %12.5f %12.5f %12s\n")
        % 0
        % _log_likelihood
        % lnP
        % _shape_updater->getCurrentPoint()
        % "---");
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
    //lnP += _statefreq_updater->calcLogPrior();
    //lnP += _exchangeability_updater->calcLogPrior();
    //lnP += _tree_updater->calcLogPrior();
    return lnP;
    }

inline void Chain::nextStep(int iteration)
    {
    Model::SharedPtr model = _likelihood->getModel();
    if (model->getGammaNCateg() > 1)
        _log_likelihood = _shape_updater->update(_log_likelihood);
    //_log_likelihood = _statefreq_updater->update(_log_likelihood);
    //_log_likelihood = _exchangeability_updater->update(_log_likelihood);
    //_log_likelihood = _tree_updater->update(_log_likelihood);
    //_log_likelihood = _tree_length_updater->update(_log_likelihood);

    double log_prior = calcLogJointPrior();

    std::cout << boost::str(boost::format("%12d %12.5f %12.5f %12.5f %12.1f\n")
        % iteration
        % _log_likelihood
        % log_prior
        % _shape_updater->getCurrentPoint()
        % _shape_updater->getAcceptPct());
    }

}
