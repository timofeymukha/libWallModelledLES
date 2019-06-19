#include "fvCFD.H"
#include "NewtonRootFinder.H"
#include <functional>
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"

using std::placeholders::_1;

struct Foo
{
    scalar val(scalar f, scalar a)
    {
        return f*f*f - a;
    }

    scalar deriv(scalar f)
    {
        return 3*f*f;
    }
};

TEST(NetwonRootFinder, FullConstructor)
{
    Foo foo;
    label maxIter = 100;
    scalar eps = 1e-5;
    std::function<scalar(scalar)> f = std::bind(&Foo::val, &foo, _1, 3);
    std::function<scalar(scalar)> d = std::bind(&Foo::deriv, &foo, _1);

    NewtonRootFinder newton = NewtonRootFinder("Newton", f, d, eps, maxIter);
    ASSERT_EQ(newton.maxIter(), maxIter);
    ASSERT_DOUBLE_EQ(newton.eps(), eps);
    ASSERT_DOUBLE_EQ(newton.f(1), -2.);
    ASSERT_DOUBLE_EQ(newton.d(1), 3.);
}

TEST(NetwonRootFinder, FuncDerivDictConstructor)
{
    Foo foo;
    std::function<scalar(scalar)> f = std::bind(&Foo::val, &foo, _1, 3);
    std::function<scalar(scalar)> d = std::bind(&Foo::deriv, &foo, _1);
    dictionary dict = dictionary();
    dict.add("eps", 1e-5);
    dict.add("maxIter", 100);

    NewtonRootFinder newton = NewtonRootFinder(f, d, dict);
    ASSERT_EQ(newton.maxIter(), 100);
    ASSERT_DOUBLE_EQ(newton.eps(), 1e-5);
    ASSERT_DOUBLE_EQ(newton.f(1), -2.);
    ASSERT_DOUBLE_EQ(newton.d(1), 3.);
}

TEST(NetwonRootFinder, DictConstructor)
{
    dictionary dict = dictionary();
    dict.add("eps", 1e-5);
    dict.add("maxIter", 100);

    NewtonRootFinder newton = NewtonRootFinder(dict);
    ASSERT_EQ(newton.maxIter(), 100);
    ASSERT_DOUBLE_EQ(newton.eps(), 1e-5);
    ASSERT_DOUBLE_EQ(newton.f(10), 0);
    ASSERT_DOUBLE_EQ(newton.d(10), 0);
}

TEST(NetwonRootFinder, DictDefaultValues)
{
    dictionary dict = dictionary();

    NewtonRootFinder newton = NewtonRootFinder(dict);
    ASSERT_EQ(newton.maxIter(), 30);
    ASSERT_DOUBLE_EQ(newton.eps(), 1e-3);
}


TEST(NewtonRootFinder, Root)
{

    Foo foo;
    label maxIter = 100;

    std::function<scalar(scalar)> value = std::bind(&Foo::val, &foo, _1, 0);
    std::function<scalar(scalar)> deriv = std::bind(&Foo::deriv, &foo, _1);

    // Construct normally
    NewtonRootFinder rootFinder = 
        NewtonRootFinder("Newton", value, deriv, 1e-14, maxIter);

    ASSERT_NEAR(rootFinder.root(2.), 0.0, 1e-10);


    // Now through the RTS
    Foam::autoPtr<RootFinder> rootFinder2 =
        RootFinder::New("Newton", value, deriv, 1e-10, maxIter);

    ASSERT_NEAR(rootFinder.root(2.), 0.0, 1e-10);
}
