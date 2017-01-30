#include "fvCFD.H"
#include "NewtonRoot.H"
#include <functional>

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

int main(int argc, char *argv[])
{
    argList args(argc, argv);

    Foo foo;
    label maxIter = 20;

    std::function<scalar(scalar)> value = std::bind(&Foo::val, &foo, _1, 1);
    std::function<scalar(scalar)> deriv = std::bind(&Foo::deriv, &foo, _1);

    NewtonRoot rootFinder = NewtonRoot(value, deriv, 1e-10, maxIter);

    Info<< rootFinder.root(2.);

    return 0;
}
