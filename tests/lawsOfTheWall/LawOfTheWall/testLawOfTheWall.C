#include "codeRules.H"
#include "fvCFD.H"
#include "LawOfTheWall.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
class DummyLawOfTheWall : public LawOfTheWall
{
public:
    TypeName("DummyLawOfTheWall");

    //- Construct from dictionary
    DummyLawOfTheWall(const dictionary & dict)
    :
    LawOfTheWall(dict)
    {}

    //- Construct from TypeName and dictionary
    DummyLawOfTheWall
    (
        const word& lawName,
        const dictionary & dict
    )
    :
    LawOfTheWall(lawName, dict)
    {}

    virtual autoPtr<LawOfTheWall> clone() const override
    {
        return autoPtr<LawOfTheWall>
        (
            new DummyLawOfTheWall(*this) 
        );
    }

    virtual void printCoeffs() const override
    {}

    //- Return the value of the implicit function defining the law
    virtual scalar value(const SingleCellSampler & sampler, label index,
                            scalar uTau, scalar nu) const override
    {
        return 0;
    }
    
    virtual scalar 
    derivative(const SingleCellSampler & sampler, label index,
                scalar uTau, scalar nu) const override
    {
        return 0;
    }
};

    defineTypeNameAndDebug(DummyLawOfTheWall, 0);
    addToRunTimeSelectionTable
    (
        LawOfTheWall,
        DummyLawOfTheWall,
        Dictionary
    );
    addToRunTimeSelectionTable
    (
        LawOfTheWall,
        DummyLawOfTheWall,
        TypeAndDictionary
    );
}

class LawOfTheWallTest : public ChannelFlow
{};

TEST_F(LawOfTheWallTest, ConstructFromDict)
{
    dictionary dict = dictionary();
    dict.add("Test1", 0.395);
    dict.add("Test2", 4);
    DummyLawOfTheWall law = DummyLawOfTheWall(dict);


    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test1", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test2", 0.0), 4);
}

TEST_F(LawOfTheWallTest, CopyConstructor)
{
    dictionary dict = dictionary();
    dict.add("Test1", 0.395);
    dict.add("Test2", 4);
    DummyLawOfTheWall law = DummyLawOfTheWall(dict);
    DummyLawOfTheWall law2(law);    

    dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test1", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test2", 0.0), 4);
}

TEST_F(LawOfTheWallTest, NewTypeNameDict)
{
    dictionary dict = dictionary();
    dict.add("Test1", 0.395);
    dict.add("Test2", 4);
    autoPtr<LawOfTheWall> law = DummyLawOfTheWall::New("DummyLawOfTheWall", dict);

    dict = law().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test1", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test2", 0.0), 4);
    ASSERT_EQ(law->type(), word("DummyLawOfTheWall"));
}

TEST_F(LawOfTheWallTest, NewDict)
{
    dictionary dict = dictionary();
    dict.add("Test1", 0.395);
    dict.add("Test2", 4);
    dict.add("type", "DummyLawOfTheWall");
    autoPtr<LawOfTheWall> law = DummyLawOfTheWall::New(dict);

    dict = law().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test1", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("Test2", 0.0), 4);
    ASSERT_EQ(law->type(), word("DummyLawOfTheWall"));
}