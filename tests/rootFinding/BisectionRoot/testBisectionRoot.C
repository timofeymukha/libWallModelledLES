#include "fvCFD.H"
#include "BisectionRootFinder.H"
#include "RootFinder.H"
#include <functional>
#include "gtest.h"
#undef Log

using std::placeholders::_1;

namespace
{
    const scalar cubicRootOfFour = Foam::cbrt(scalar(4));
    const scalar rootTolerance = 1e-6;

    struct Cubic
    {
        scalar val(scalar x, scalar a)
        {
            return x*x*x - a;
        }

        scalar deriv(scalar x)
        {
            return 3*x*x;
        }
    };
}


TEST(BisectionRootFinder, FullConstructor)
{
    Cubic cubic;
    label maxIter = 100;
    std::function<scalar(scalar)> f =
        std::bind(&Cubic::val, &cubic, _1, 3);
    std::function<scalar(scalar)> d =
        std::bind(&Cubic::deriv, &cubic, _1);

    BisectionRootFinder bisection =
        BisectionRootFinder("Bisection", f, d, maxIter);

    ASSERT_EQ(bisection.maxIter(), maxIter);
    ASSERT_DOUBLE_EQ(bisection.f(1), -2);
    ASSERT_DOUBLE_EQ(bisection.d(1), 3);
}


TEST(BisectionRootFinder, DictDefaultValues)
{
    dictionary dict;
    BisectionRootFinder bisection(dict);

    ASSERT_EQ(bisection.maxIter(), 30);
}


TEST(BisectionRootFinder, Root)
{
    Cubic cubic;
    label maxIter = 100;
    std::function<scalar(scalar)> value =
        std::bind(&Cubic::val, &cubic, _1, 4);
    std::function<scalar(scalar)> deriv =
        std::bind(&Cubic::deriv, &cubic, _1);

    BisectionRootFinder rootFinder =
        BisectionRootFinder("Bisection", value, deriv, maxIter);

    std::pair<scalar, label> result = rootFinder.root(1, 1, 2);

    ASSERT_NEAR(result.first, cubicRootOfFour, rootTolerance);
    ASSERT_GT(result.second, 0);
    ASSERT_LE(result.second, maxIter);
}


TEST(BisectionRootFinder, ReversedBounds)
{
    Cubic cubic;
    std::function<scalar(scalar)> value =
        std::bind(&Cubic::val, &cubic, _1, 4);
    std::function<scalar(scalar)> deriv =
        std::bind(&Cubic::deriv, &cubic, _1);

    BisectionRootFinder rootFinder =
        BisectionRootFinder("Bisection", value, deriv, 100);

    std::pair<scalar, label> result = rootFinder.root(1, 2, 1);

    ASSERT_NEAR(result.first, cubicRootOfFour, rootTolerance);
}


TEST(BisectionRootFinder, RootAtBound)
{
    Cubic cubic;
    std::function<scalar(scalar)> value =
        std::bind(&Cubic::val, &cubic, _1, 1);
    std::function<scalar(scalar)> deriv =
        std::bind(&Cubic::deriv, &cubic, _1);

    BisectionRootFinder rootFinder =
        BisectionRootFinder("Bisection", value, deriv, 100);

    std::pair<scalar, label> result = rootFinder.root(0, 1, 2);

    ASSERT_DOUBLE_EQ(result.first, 1);
    ASSERT_EQ(result.second, 0);
}


TEST(BisectionRootFinder, RuntimeSelectionTable)
{
    Cubic cubic;
    std::function<scalar(scalar)> value =
        std::bind(&Cubic::val, &cubic, _1, 4);
    std::function<scalar(scalar)> deriv =
        std::bind(&Cubic::deriv, &cubic, _1);

    autoPtr<RootFinder> rootFinder =
        RootFinder::New("Bisection", value, deriv, 100);

    std::pair<scalar, label> result = rootFinder->root(1, 1, 2);

    ASSERT_NEAR(result.first, cubicRootOfFour, rootTolerance);
}
