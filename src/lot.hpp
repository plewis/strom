#pragma once

#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace strom
    {

    class Lot
        {
        public:
                                    Lot();
                                    ~Lot();

            void                    setSeed(unsigned seed);
            double                  uniform();
            int                     randint(int low, int high);
            double                  normal();
            double                  gamma(double shape, double scale);
            double                  logUniform();

            typedef boost::shared_ptr<Lot> SharedPtr;

        private:

            typedef boost::variate_generator<boost::mt19937 &, boost::uniform_real<> > uniform_variate_generator_t;
            typedef boost::variate_generator<boost::mt19937 &, boost::normal_distribution<> > normal_variate_generator_t;
            typedef boost::variate_generator<boost::mt19937 &, boost::gamma_distribution<> > gamma_variate_generator_t;

            unsigned                        _seed;
            boost::mt19937                  _generator;
            uniform_variate_generator_t *   _uniform_variate_generator;
            normal_variate_generator_t  *   _normal_variate_generator;
        };

    inline Lot::Lot() : _seed(0), _generator(1), _uniform_variate_generator(0), _normal_variate_generator(0)
        {
        _generator.seed(static_cast<unsigned int>(std::time(0)));
        _uniform_variate_generator = new uniform_variate_generator_t(_generator, boost::uniform_real<>(0,1));
        _normal_variate_generator = new normal_variate_generator_t(_generator, boost::normal_distribution<>(0,1));
        }

    inline Lot::~Lot()
        {
        }

    inline void Lot::setSeed(unsigned seed)
        {
        _seed = seed;
        _generator.seed(_seed > 0 ? _seed : static_cast<unsigned int>(std::time(0)));
        }

    inline double Lot::uniform()
        {
        return (*_uniform_variate_generator)();
        }

    inline double Lot::normal()
        {
        return (*_normal_variate_generator)();
        }

    inline double Lot::gamma(double shape, double scale)
        {
        return gamma_variate_generator_t(_generator, boost::gamma_distribution<>(shape,scale))();
        }

    inline double Lot::logUniform()
        {
        double u = (*_uniform_variate_generator)();
        assert(u > 0.0);
        return std::log(u);
        }

    inline int Lot::randint(int low, int high)
        {
        // Example: return random int k between low and high (inclusive)
        // low = -2, high = 3, 3-(-2)+1 = 6
        // Draw u ~ Uniform(0,1)
        //   k  u between  details of calculation
        //  -2  0/6 - 1/6  i = 0 retval = 0 + -2 = -2 = low
        //  -1  1/6 - 2/6  i = 1 retval = 1 + -2 = -1
        //   0  2/6 - 3/6  i = 2 retval = 2 + -2 =  0
        //   1  3/6 - 4/6  i = 3 retval = 3 + -2 =  1
        //   2  4/6 - 5/6  i = 4 retval = 4 + -2 =  2
        //   3  5/6 - 6/6  i = 5 retval = 5 + -2 =  3 = high
        double u = (*_uniform_variate_generator)();
        double n = (double)(high - low + 1);
        int retval = n;
        double cumprob = 0.0;
        for (unsigned i = 0; i < n; ++i)
            {
            cumprob += 1.0/n;
            if (u < cumprob)
                {
                retval = i;
                break;
                }
            }
        retval += low;
        return retval;
        }
    }
