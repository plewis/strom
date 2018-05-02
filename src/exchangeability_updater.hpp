#pragma once

#include "dirichlet_updater.hpp"

namespace strom
    {

    class ExchangeabilityUpdater : public DirichletUpdater
        {
        public:

                                        ExchangeabilityUpdater();
                                        ~ExchangeabilityUpdater();

            virtual void                pullCurrentStateFromModel();
            virtual void                pushCurrentStateToModel() const;

            std::vector<double>         getCurrentPoint() const;

        public:
            typedef std::shared_ptr< ExchangeabilityUpdater > SharedPtr;
        };

inline ExchangeabilityUpdater::ExchangeabilityUpdater()
    {
    // std::cout << "Creating an ExchangeabilityUpdater" << std::endl;
    DirichletUpdater::clear();
    _name = "Exchangeabilities";
    }

inline ExchangeabilityUpdater::~ExchangeabilityUpdater()
    {
    // std::cout << "Destroying an ExchangeabilityUpdater" << std::endl;
    }

inline std::vector<double> ExchangeabilityUpdater::getCurrentPoint() const
    {
    return _curr_point;
    }

inline void ExchangeabilityUpdater::pullCurrentStateFromModel()
    {
    Model::SharedPtr model = _likelihood->getModel();
    const std::vector<double> & xchg = model->getExchangeabilities();
    _curr_point.assign(xchg.begin(), xchg.end());
    }

inline void ExchangeabilityUpdater::pushCurrentStateToModel() const
    {
    // sanity checks
    assert(_curr_point.size() == 6);
    assert(fabs(std::accumulate(_curr_point.begin(), _curr_point.end(), 0.0) - 1.0) < 1.e-8);

    Model::SharedPtr model = _likelihood->getModel();
    model->setExchangeabilities(_curr_point);
    }
}
