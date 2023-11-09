Gauss-Lobatto Adaptive Integrator
===================

A C++ port of the adaptive quadrature algorithm by [Gander and Gautschi](http://link.springer.com/article/10.1023%2FA%3A1022318402393). The same algorithm is implemented in [Numerical Recipes](http://www.nr.com/), but their license is very restrictive. I translated the code from the original paper so it can be used freely, as outlined in the [MIT License](http://opensource.org/licenses/MIT).

Usage
-------------------
A sample example file,`example.cpp`, is provided.
Simply include the header `AdaptiveIntegrator.hpp` in your code:
```cpp
    #include "AdaptiveIntegrator.hpp"
```
Then create an integrator object specialized for the type of functional 
you wish to integrate. In the simplest case this is just a normal function 
with signature `double(const double)`, so the declaration looks like
```cpp
    AdaptiveIntegrator<double(const double)> test
```
Finally to integrate your function, `func`, invoke the `integrate` method 
with integral bounds `a` and `b` and desired tolerance `tol`:
```cpp
    test.integrate(func, a, b, tol)
```
