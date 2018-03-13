#pragma once

#include <memory>
#include <boost/format.hpp>
#include "lot.hpp"
#include "data.hpp"
#include "tree.hpp"
#include "likelihood.hpp"
#include "tree_manip.hpp"
#include "gamma_shape_updater.hpp"

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

            double                              _tmp;

        private:

            Likelihood::SharedPtr               _likelihood;
            TreeManip::SharedPtr                _tree_manipulator;

            GammaShapeUpdater::SharedPtr        _shape_updater;

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
    return v;
    }

inline void Chain::setLambdas(std::vector<double> & v)
    {
    _shape_updater->setLambda(v[0]);
    }

inline void Chain::startTuning()
    {
    _shape_updater->setTuning(true);
    }

inline void Chain::stopTuning()
    {
    _shape_updater->setTuning(false);
    }

inline void Chain::setTreeManip(TreeManip::SharedPtr tm)
    {
    _tree_manipulator = tm;
    _shape_updater->setTreeManip(_tree_manipulator);
    }

inline void Chain::setLikelihood(typename Likelihood::SharedPtr likelihood)
    {
    _likelihood = likelihood;
    _shape_updater->setLikelihood(likelihood);
    }

inline void Chain::setLot(typename Lot::SharedPtr lot)
    {
    _shape_updater->setLot(lot);
    }

inline void Chain::clear()
    {
    _tmp = 0.0;
    _log_likelihood = 0.0;

    std::vector<double> edge_length_prior_parameters(2, 1.0);
    edge_length_prior_parameters[1] = 10.0;

    _shape_updater.reset(new GammaShapeUpdater);
    _shape_updater->setLambda(1.0);
    _shape_updater->setTargetAcceptanceRate(0.3);
    _shape_updater->setPriorParameters(std::vector<double>(2, 1.0));
    //std::vector<double> v = {2.0,10.0};
    //_shape_updater->setPriorParameters(v);

    setHeatingPower(1.0);
    startTuning();
    }

inline void Chain::start()
    {
    _shape_updater->pullCurrentStateFromModel();
    _log_likelihood = calcLogLikelihood();
    if (_heating_power == 1.0)
        {
        std::cout << boost::str(boost::format("%12s %12s %12s %12s %12s")
            % "iteration"
            % "lnLike"
            % "lnPrior"
            % "shape"
            % "accept"
            ) << std::endl;

        double log_prior = calcLogJointPrior();
        std::cout << boost::str(boost::format("%12d %12.5f %12.5f %12.5f %12s")
                % 0
                % _log_likelihood
                % log_prior
                % _shape_updater->getCurrentPoint()
                % "---") << std::endl;
        }
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
    lnP += _shape_updater->calcGammaShapePrior();
    return lnP;
    }

inline void Chain::nextStep(int iteration, unsigned sampling_freq)
    {
    GTRModel::SharedPtr gtr = _likelihood->getModel();
    if (gtr->getGammaNCateg() > 1)
        _log_likelihood = _shape_updater->update(_log_likelihood);
    double log_prior = calcLogJointPrior();
    if (_heating_power == 1.0)
        {
        if (sampling_freq > 0 && iteration % sampling_freq == 0)
            {
            if (iteration > 0)
                _tmp += _shape_updater->getCurrentPoint();
            std::cout << boost::str(boost::format("%12d %12.5f %12.5f %12.5f %12.1f")
                    % iteration
                    % _log_likelihood
                    % log_prior
                    % _shape_updater->getCurrentPoint()
                    % _shape_updater->getAcceptPct()) << std::endl;
            }
        }
    }

}
