#include "fvCFD.H"
#include "NewtonRoot.H"

scalar func(scalar f)
{
    return f*f*f - 10;
}

scalar deriv(scalar f)
{
    return 3*f*f;
}

int main(int argc, char *argv[])
{
    argList args(argc, argv);

//    xSquaredDeriv deriv = xSquaredDeriv();

    //NewtonSecantRoot<xSquared,xSquaredDeriv> rootFinder = NewtonSecantRoot<xSquared,xSquaredDeriv>(func, deriv, 1e-10);

    NewtonRoot rootFinder = NewtonRoot(func, deriv, 1e-10);
    Info<< rootFinder.root(2.);
    

    return 0;
}
