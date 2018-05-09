#pragma once

#include "updater.hpp"

namespace strom
{

    class TreeLengthUpdater : public Updater
    {
        public:

                                        TreeLengthUpdater();
                                        ~TreeLengthUpdater();

            virtual void                clear();
            virtual void                pullCurrentStateFromModel();
            virtual void                pushCurrentStateToModel() const;
            virtual void                proposeNewState();
            virtual void                revert();

            virtual double              calcLogPrior() const;

        private:

            double                      _prev_point;
            double                      _curr_point;

        public:

            typedef std::shared_ptr< TreeLengthUpdater > SharedPtr;
    };

inline TreeLengthUpdater::TreeLengthUpdater()
    {
    // std::cout << "Creating a TreeLengthUpdater..." << std::endl;
    clear();
    _name = "Tree Length"; //POLNEW
    }

inline TreeLengthUpdater::~TreeLengthUpdater()
    {
    // std::cout << "Destroying a TreeLengthUpdater..." << std::endl;
    }

inline void TreeLengthUpdater::clear()
    {
    Updater::clear();
    _prev_point     = 0.0;
    _curr_point     = 0.0;
    reset();
    }

inline double TreeLengthUpdater::calcLogPrior() const
    {
    return Updater::calcEdgeLengthPrior();
    }

inline void TreeLengthUpdater::pullCurrentStateFromModel()
    {
    _curr_point = _tree_manipulator->calcTreeLength();
    }

inline void TreeLengthUpdater::pushCurrentStateToModel() const
    {
    double scaler = _curr_point/_prev_point;
    _tree_manipulator->scaleAllEdgeLengths(scaler);
    }

inline void TreeLengthUpdater::proposeNewState()
    {
    // Save copy of _curr_point in case revert is necessary.
    _prev_point = _curr_point;

    // Let _curr_point be proposed value
    double m = exp(_lambda*(_lot->uniform() - 0.5));
    _curr_point = m*_prev_point;

    // calculate log of Hastings ratio under GammaDir parameterization
    _log_hastings_ratio = log(m);
    }

inline void TreeLengthUpdater::revert()
    {
    // swap _curr_point and _prev_point so that edge length scaler
    // in pushCurrentStateToModel will be correctly calculated
    double tmp = _curr_point;
    _curr_point = _prev_point;
    _prev_point = tmp;
    }

}
