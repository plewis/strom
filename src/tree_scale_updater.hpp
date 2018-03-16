#pragma once

#include "updater.hpp"

namespace strom
{

    class TreeScaleUpdater : public Updater
    {
        public:

                                        TreeScaleUpdater();
                                        ~TreeScaleUpdater();

            virtual void                clear();
            virtual void                pullCurrentStateFromModel();
            virtual void                pushCurrentStateToModel() const;
            virtual void                proposeNewState();
            virtual void                revert();

            virtual double              calcLogPrior() const;

        private:

            double                      _prev_point;
            double                      _curr_point;
            unsigned                    _nedges;

            mutable double              _tree_length;

        public:

            typedef std::shared_ptr< TreeScaleUpdater > SharedPtr;
    };

    inline TreeScaleUpdater::TreeScaleUpdater()
        {
        _name = "scale";
        clear();
        }

    inline TreeScaleUpdater::~TreeScaleUpdater()
        {
        }

    inline void TreeScaleUpdater::clear()
        {
        Updater::clear();
        _prev_point     = 0.0;
        _curr_point     = 0.0;
        _tree_length    = 0.0;
        _nedges         = 1;
        reset();
        }

    inline double TreeScaleUpdater::calcLogPrior() const
        {
        return Updater::calcEdgeLengthPrior();
        }

    inline void TreeScaleUpdater::pullCurrentStateFromModel()
        {
        _nedges = _tree_manipulator->getTree()->numNodes() - 1;
        _curr_point = _tree_manipulator->calcTreeLength();
        _tree_length = _curr_point;
        }

    inline void TreeScaleUpdater::pushCurrentStateToModel() const
        {
        double scaler = _curr_point/_tree_length;
        _tree_length = _curr_point;
        _tree_manipulator->scaleAllEdgeLengths(scaler);
        }

    inline void TreeScaleUpdater::proposeNewState()
        {
        // Save copy of _curr_point in case revert is necessary.
        _prev_point = _curr_point;

#if 1
        // Let _curr_point be proposed value
        double m = exp(_lambda*(_lot->uniform() - 0.5));
        _curr_point = m*_prev_point;

        // calculate log of Hastings ratio
        _log_hastings_ratio = _nedges*log(m);
#else
        _log_hastings_ratio = 0.0;
#endif
        }

    inline void TreeScaleUpdater::revert()
        {
        _curr_point = _prev_point;
        }

}

