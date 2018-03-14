#pragma once

#include "dirichlet_updater.hpp"

namespace strom
{

    class StateFreqUpdater : public DirichletUpdater
        {
        public:

                                        StateFreqUpdater();
                                        ~StateFreqUpdater();

            virtual void                pullCurrentStateFromModel();
            virtual void                pushCurrentStateToModel() const;

            std::vector<double>         getCurrentPoint() const;

        public:
            typedef std::shared_ptr< StateFreqUpdater > SharedPtr;
        };

inline StateFreqUpdater::StateFreqUpdater()
    {
    //std::cout << "Creating a StateFreqUpdater" << std::endl;
    DirichletUpdater::clear();
    _name = "freq";
    }

inline StateFreqUpdater::~StateFreqUpdater()
    {
    //std::cout << "Destroying a StateFreqUpdater" << std::endl;
    }

inline std::vector<double> StateFreqUpdater::getCurrentPoint() const
    {
    return _curr_point;
    }

inline void StateFreqUpdater::pullCurrentStateFromModel()
    {
    GTRModel::SharedPtr gtr = _likelihood->getModel();
    const std::vector<double> & freqs = gtr->getStateFreqs();
    _curr_point.assign(freqs.begin(), freqs.end());
    }

inline void StateFreqUpdater::pushCurrentStateToModel() const
    {
    // sanity checks
    assert(_curr_point.size() == 4);
    assert(fabs(std::accumulate(_curr_point.begin(), _curr_point.end(), 0.0) - 1.0) < 1.e-8);

    GTRModel::SharedPtr gtr = _likelihood->getModel();
    gtr->setStateFreqs(_curr_point);
    }

}
