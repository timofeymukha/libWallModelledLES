#ifndef ADAPTIVEINTEGRATOR_HPP_
#define ADAPTIVEINTEGRATOR_HPP_

#include "Integrator.hpp"
#include <cmath>
using std::sqrt;
using std::fabs;
#include <limits>
using std::numeric_limits;
#include <iostream>
using std::cerr;
#include <functional>

typedef std::function<double(const double)> IntegrandFunc;

template <class T>
class AdaptiveIntegrator : Integrator<T> {
private:
    const static double alpha, beta, x1, x2, x3;
    const static double x[12];

    bool terminated;
public:
    AdaptiveIntegrator();

    void set_tol(double tol);

    double integrate(IntegrandFunc func, const double a, const double b, const double tol);

    double adaptlobstp(IntegrandFunc func, const double a, const double b, const double fa,
         const double fb, const double is);
};
template <class T>
const double AdaptiveIntegrator<T>::alpha = sqrt(2./3.);
template <class T>
const double AdaptiveIntegrator<T>::beta = 1./sqrt(5.);
template <class T>
const double AdaptiveIntegrator<T>::x1 = .94288241569547971905635175843185720232;
template <class T>
const double AdaptiveIntegrator<T>::x2 = .64185334234578130578123554132903188354;
template <class T>
const double AdaptiveIntegrator<T>::x3 = .23638319966214988028222377349205292599;

template <class T>
AdaptiveIntegrator<T>::AdaptiveIntegrator() : terminated(false) {}

template <class T>
double AdaptiveIntegrator<T>::integrate(IntegrandFunc f, const double a, const double b, const double tol_)
{
    double tol, eps;
    eps = numeric_limits<double>::epsilon();
    tol = (tol_ < eps) ?  eps : tol_;

    double m, h;
    m = (a+b)/2.; h = (b-a)/2.;

    double y[13] = {f(a),f(m-x1*h),f(m-alpha*h),f(m-x2*h),f(m-beta*h),f(m-x3*h),f(m),
                    f(m+x3*h),f(m+beta*h),f(m+x2*h),f(m+alpha*h),f(m+x1*h) ,f(b)};

    double fa, fb;
    fa = y[0]; fb = y[12];

    double i1, i2, is;

    i2 = (h/6.)*(y[0] + y[12] + 5.*(y[4] + y[8]));
    i1 = (h/1470.)*(77.*(y[0]+y[12]) + 432.*(y[2]+y[10]) + 625.*(y[4]+y[8]) + 672.*y[6]);
    is = h*(.0158271919734802*(y[0]+y[12]) + .0942738402188500*(y[1]+y[11])
           + .155071987336585*(y[2]+y[10]) + .188821573960182*(y[3]+y[9])
           + .199773405226859*(y[4]+y[8]) + .224926465333340*(y[5]+y[7])
           + .242611071901408*y[6]);

    double erri1, erri2, R;
    erri1 = fabs(i1 - is); erri2 = fabs(i2 - is);
    R = (erri2 != 0.) ? erri1/erri2 : 1.;

    tol = (R > 0. and R < 1.) ? tol_/R : tol_;
    is = fabs(is)*tol/eps;
    if (is == 0.) is = b-a;

    return adaptlobstp(f, a, b, fa, fb, is);

}

template <class T>
double AdaptiveIntegrator<T>::adaptlobstp(IntegrandFunc f, const double a, const double b,
                                          const double fa, const double fb, const double is )
{
    double m, h;
    m = (a+b)/2.; h = (b-a)/2.;

    double mll, ml, mr, mrr;
    mll = m - alpha*h; ml = m - beta*h; mr = m + beta*h; mrr = m + alpha*h;

    double fmll, fml, fm, fmr, fmrr;
    fmll = f(mll); fml = f(ml); fm = f(m); fmr = f(mr); fmrr = f(mrr);

    double i1, i2;
    i2 = (h/6.)*(fa + fb + 5.*(fml+fmr));
    i1 = (h/1470.)*(77*(fa+fb) + 432.*(fmll + fmrr) + 625.*(fml + fmr) + 672.*fm);

    if (is + (i1-i2) == is or mll <= a or b <= mrr)
    {
        if ( (m <= a or b <= m) and !terminated)
        {
            cerr << "No machine number in the interval. Requested tolerance may not be met.\n";
            terminated = true;
        }
        return i1;
    }
    else
    {
        return adaptlobstp(f,a,mll,fa,fmll,is)
             + adaptlobstp(f,mll,ml,fmll,fml,is)
             + adaptlobstp(f,ml,m,fml,fm,is)
             + adaptlobstp(f,m,mr,fm,fmr,is)
             + adaptlobstp(f,mr,mrr,fmr,fmrr,is)
             + adaptlobstp(f,mrr,b,fmrr,fb,is);
    }
}



#endif // ADAPTIVEINTEGRATOR_HPP_
