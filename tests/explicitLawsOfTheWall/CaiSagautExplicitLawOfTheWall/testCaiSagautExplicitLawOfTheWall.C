#include "codeRules.H"
#include "fvCFD.H"
#include "CaiSagautExplicitLawOfTheWall.H"
#undef Log
#include "gtest.h"

class CaiSagautExplicitLawOfTheWallTest : public ::testing::Test
{};

TEST_F(CaiSagautExplicitLawOfTheWallTest, ConstructFromConstants)
{
    CaiSagautExplicitLawOfTheWall law(0.395, 4.0, 1.2, 220.0);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4.0);
    ASSERT_DOUBLE_EQ(law.p(), 1.2);
    ASSERT_DOUBLE_EQ(law.s(), 220.0);

    dictionary dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4.0);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("p", 0.0), 1.2);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("s", 0.0), 220.0);
}


TEST_F(CaiSagautExplicitLawOfTheWallTest, ConstructFromDictionary)
{
    dictionary dict;
    dict.add("kappa", 0.395);
    dict.add("B", 4.0);
    dict.add("p", 1.2);
    dict.add("s", 220.0);

    CaiSagautExplicitLawOfTheWall law(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4.0);
    ASSERT_DOUBLE_EQ(law.p(), 1.2);
    ASSERT_DOUBLE_EQ(law.s(), 220.0);
}


TEST_F(CaiSagautExplicitLawOfTheWallTest, ConstructFromDictionaryDefaultValues)
{
    dictionary dict;
    CaiSagautExplicitLawOfTheWall law(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.4);
    ASSERT_DOUBLE_EQ(law.B(), 5.5);
    ASSERT_DOUBLE_EQ(law.p(), 1.138);
    ASSERT_DOUBLE_EQ(law.s(), 217.8);
}


TEST_F(CaiSagautExplicitLawOfTheWallTest, Clone)
{
    CaiSagautExplicitLawOfTheWall law(0.395, 4.0, 1.2, 220.0);
    autoPtr<ExplicitLawOfTheWall> law2 = law.clone();

    const CaiSagautExplicitLawOfTheWall& clonedLaw =
        dynamic_cast<const CaiSagautExplicitLawOfTheWall&>(law2());

    ASSERT_DOUBLE_EQ(clonedLaw.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(clonedLaw.B(), 4.0);
    ASSERT_DOUBLE_EQ(clonedLaw.p(), 1.2);
    ASSERT_DOUBLE_EQ(clonedLaw.s(), 220.0);
}
