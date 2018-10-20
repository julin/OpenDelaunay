#if defined(_ENABLE_UNITEST)
#include <gtest/gtest.h>

#define EXPECT_PERCENT_NEAR(val1, val2, percent) \
  EXPECT_NEAR(val1,val2,(double)val1*percent)

#define EXPECT_PERCENT0(val1, val2) \
  EXPECT_PERCENT_NEAR(val1,val2,0.25)

#endif
