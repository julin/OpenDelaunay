#if defined(_ENABLE_UNITEST)
#include <gtest/gtest.h>

#define EXPECT_PERCENT_NEAR(val1, val2, percent) \
  EXPECT_NEAR(val1,val2,(double)val1*percent)

#define EXPECT_PERCENT0(val1, val2) \
  EXPECT_PERCENT_NEAR(val1,val2,0.25)

#if defined (_WINDOWS)
#if defined (_DEBUG)
static int gRunAllTest = 0;
#else
static int gRunAllTest = 1;
#endif
#define MEXT_ENABLE_ASSERT _set_abort_behavior(1, _WRITE_ABORT_MSG)
#else
#define MEXT_ENABLE_ASSERT 
static int gRunAllTest = 1;
#endif
#endif
