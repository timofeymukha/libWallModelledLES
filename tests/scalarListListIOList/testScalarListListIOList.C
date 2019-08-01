#include "codeRules.H"
#include "fvCFD.H"
#include "scalarListListIOList.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"

class ScalarListListIOList : public ::testing::Test
{};


TEST(ScalarListListIOList, Constructors)
{
    scalarListListList list(10);
    ASSERT_EQ(list.size(), 10);

    scalarListListList list2(2, scalarListList(3, scalarList(2, 3.0)));
    ASSERT_EQ(list2.size(), 2);

    forAll(list2, i)
    {
        ASSERT_EQ(list2[i].size(), 3);

        forAll(list2[i], j)
        {
            ASSERT_EQ(list2[i][j].size(), 2);

            forAll(list2[i][j], k)
            {
                ASSERT_DOUBLE_EQ(list2[i][j][k], 3.0);
            }
        }
    }


}