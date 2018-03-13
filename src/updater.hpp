#pragma once

#include "tree.hpp"
#include "tree_manip.hpp"
#include "lot.hpp"
#include "xstrom.hpp"
#include "likelihood.hpp"

namespace strom
{
    class Chain;

    class Updater
        {
        friend class Chain;

        public:

                                    Updater();
            virtual                 ~Updater();

            void                    setLikelihood(typename Likelihood::SharedPtr likelihood);
            void                    setTreeManip(typename TreeManip::SharedPtr treemanip);
            void                    setLot(Lot::SharedPtr lot);
            void                    setLambda(double lambda);
            void                    setHeatingPower(double p);
            void                    setTuning(bool on);
            void                    setTargetAcceptanceRate(double target);
            void                    setPriorParameters(const std::vector<double> & c);

            TreeManip::SharedPtr    getTreeManip() const;
            double                  getLambda() const;
            double                  getAcceptPct() const;

            virtual void            clear();

            virtual double          calcLogPrior() const = 0;
            double                  calcLogLikelihood() const;
            virtual double          update(double prev_lnL);

        protected:

            virtual void            reset();
            virtual void            tune(bool accepted);

            virtual void            revert() = 0;
            virtual void            proposeNewState() = 0;
            virtual void            pullCurrentStateFromModel() = 0;
            virtual void            pushCurrentStateToModel() const = 0;

            Lot::SharedPtr          _lot;
            Likelihood::SharedPtr   _likelihood;
            TreeManip::SharedPtr    _tree_manipulator;
            std::string             _name;
            double                  _lambda;
            double                  _log_hastings_ratio;
            double                  _target_acceptance;
            unsigned                _naccepts;
            unsigned                _nattempts;
            bool                    _tuning;
            //double                  _prev_log_likelihood;
            //double                  _prev_log_prior;
            std::vector<double>     _prior_parameters;

            double                  _heating_power;

        };

inline Updater::Updater()
    {
    //std::cout << "Updater constructor called" << std::endl;
    clear();
    }

inline Updater::~Updater()
    {
    //std::cout << "Updater destructor called" << std::endl;
    }

inline void Updater::clear()
    {
    _name                   = "updater";
    _tuning                 = true;
    _lambda                 = 0.0001;
    _target_acceptance      = 0.3;
    _naccepts               = 0;
    _nattempts              = 0;
    _heating_power          = 1.0;
    _prior_parameters.clear();
    reset();
    }

inline void Updater::reset()
    {
    _log_hastings_ratio = 0.0;
    }

inline void Updater::setLikelihood(typename Likelihood::SharedPtr likelihood)
    {
    _likelihood = likelihood;
    }

inline void Updater::setTreeManip(typename TreeManip::SharedPtr treemanip)
    {
    _tree_manipulator = treemanip;
    }

inline TreeManip::SharedPtr Updater::getTreeManip() const
    {
    return _tree_manipulator;
    }

inline void Updater::setLot(Lot::SharedPtr lot)
    {
    _lot = lot;
    }

inline void Updater::setHeatingPower(double p)
    {
    _heating_power = p;
    }

inline void Updater::setLambda(double lambda)
    {
    _lambda = lambda;
    }

void Updater::setTuning(bool do_tune)
    {
    _tuning = do_tune;
    _naccepts = 0;
    _nattempts = 0;
    }

inline void Updater::tune(bool accepted)
    {
    _nattempts++;
    if (_tuning)
        {
        double gamma_n = 10.0/(100.0 + (double)_nattempts);
        if (accepted)
            _lambda *= 1.0 + gamma_n*(1.0 - _target_acceptance)/(2.0*_target_acceptance);
        else
            _lambda *= 1.0 - gamma_n*0.5;
        }
    }

inline void Updater::setTargetAcceptanceRate(double target)
    {
    _target_acceptance = target;
    }

inline void Updater::setPriorParameters(const std::vector<double> & c)
    {
    _prior_parameters.clear();
    _prior_parameters.assign(c.begin(), c.end());
    }

inline double Updater::getLambda() const
    {
    return _lambda;
    }

inline double Updater::getAcceptPct() const
    {
    return (_nattempts == 0 ? 0.0 : (100.0*_naccepts/_nattempts));
    }

inline double Updater::calcLogLikelihood() const
    {
    Tree::SharedPtr t = _tree_manipulator->getTree();
    return _likelihood->calcLogLikelihood(t);
    //return _likelihood->calcLogLikelihood(_tree_manipulator->getTree());
    }

inline double Updater::update(double prev_lnL)
    {
    // Copy current state from model into _curr_point.
    pullCurrentStateFromModel();

    double prev_log_prior      = calcLogPrior();

    // Set model to proposed state and calculate _log_hastings_ratio
    proposeNewState();
    pushCurrentStateToModel();

    double log_likelihood = calcLogLikelihood();
    double log_prior = calcLogPrior();
    double log_diff = _log_hastings_ratio;
    log_diff += _heating_power*((log_likelihood + log_prior) - (prev_lnL + prev_log_prior));

    //std::cerr << boost::str(boost::format("prev_lnL       = %12.5f") % prev_lnL) << std::endl;
    //std::cerr << boost::str(boost::format("log_likelihood = %12.5f") % log_likelihood) << std::endl;
    //std::cerr << boost::str(boost::format("prev_log_prior = %12.5f") % prev_log_prior) << std::endl;
    //std::cerr << boost::str(boost::format("log_prior      = %12.5f") % log_prior) << std::endl;
    //std::cerr << std::endl;

    bool accept = true;
    double logu = _lot->logUniform();
    if (logu > log_diff)
        accept = false;

    if (accept)
        _naccepts++;
    else
        {
        revert();
        pushCurrentStateToModel();
        log_likelihood = prev_lnL;
        }

    tune(accept);
    reset();

    return log_likelihood;
    }

}
