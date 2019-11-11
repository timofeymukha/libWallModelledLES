#include "codeRules.H"
#include "fvCFD.H"
#undef Log
#include "gtest.h"


Foam::argList * mainArgs;

int main(int argc, char **argv)
{
    Foam::argList::noBanner();

    ::testing::InitGoogleTest(&argc, argv);

    mainArgs = new Foam::argList(argc, argv);

    return RUN_ALL_TESTS();
    delete mainArgs;
}

