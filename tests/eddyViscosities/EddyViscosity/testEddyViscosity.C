#include "codeRules.H"
#include "fvCFD.H"
#include "EddyViscosity.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
class DummyEddyViscosity : public EddyViscosity
{
public:
    TypeName("DummyEddyViscosity");

    //- Construct from dictionary
    DummyEddyViscosity(const dictionary & dict)
    :
    EddyViscosity(dict)
    {}

    //- Construct from TypeName and dictionary
    DummyEddyViscosity
    (
        const word& lawName,
        const dictionary & dict
    )
    :
    EddyViscosity(lawName, dict)
    {}

    virtual autoPtr<EddyViscosity> clone() const override
    {
        return autoPtr<EddyViscosity>
        (
            new DummyEddyViscosity(*this) 
        );
    }

    virtual void printCoeffs() const override
    {}

    //- Return the value of the implicit function defining the law
    virtual scalarList
    value
    (
        const SingleCellSampler & sampler,
        label index,
        const scalarList & y,
        scalar uTau,
        scalar nu
    ) const override
    {
        return scalarList(1, 0.0);
    }
};

    defineTypeNameAndDebug(DummyEddyViscosity, 0);
    addToRunTimeSelectionTable
    (
        EddyViscosity,
        DummyEddyViscosity,
        Dictionary
    );
    addToRunTimeSelectionTable
    (
        EddyViscosity,
        DummyEddyViscosity,
        TypeAndDictionary
    );
}

class EddyViscosityTest : public ChannelFlow
{};

TEST_F(EddyViscosityTest, ConstructFromDict)
{
    dictionary dict = dictionary();
    dict.add("Test1", 0.395);
    dict.add("Test2", 4);
    DummyEddyViscosity law = DummyEddyViscosity(dict);


    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test1", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test2", 0.0), 4);
}

TEST_F(EddyViscosityTest, CopyConstructor)
{
    dictionary dict = dictionary();
    dict.add("Test1", 0.395);
    dict.add("Test2", 4);
    DummyEddyViscosity law = DummyEddyViscosity(dict);
    DummyEddyViscosity law2(law);    

    dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test1", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test2", 0.0), 4);
}

TEST_F(EddyViscosityTest, NewTypeNameDict)
{
    dictionary dict = dictionary();
    dict.add("Test1", 0.395);
    dict.add("Test2", 4);
    autoPtr<EddyViscosity> law = DummyEddyViscosity::New("DummyEddyViscosity", dict);

    dict = law().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test1", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test2", 0.0), 4);
    ASSERT_EQ(law->type(), word("DummyEddyViscosity"));
}

TEST_F(EddyViscosityTest, NewDict)
{
    dictionary dict = dictionary();
    dict.add("Test1", 0.395);
    dict.add("Test2", 4);
    dict.add("type", "DummyEddyViscosity");
    autoPtr<EddyViscosity> law = DummyEddyViscosity::New(dict);

    dict = law().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test1", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test2", 0.0), 4);
    ASSERT_EQ(law->type(), word("DummyEddyViscosity"));
}