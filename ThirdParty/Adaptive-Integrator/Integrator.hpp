#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_
#include <functional>

template <class T>
class Integrator {

public:
    virtual double integrate(std::function<double(const double)>, const double a, const double b, const double tol) = 0;

    virtual ~Integrator(){}
};

#endif // INTEGRATOR_HPP_

