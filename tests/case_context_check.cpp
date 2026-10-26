/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "tests/test_harness.h"

static void lambda_case() {
  [] {
    HS_CONTEXT("sample", 7);
    HS_EXPECT_EQ(1, 2);
  }();
}

static void argument_case(int &calls) {
  ++calls;
  HS_EXPECT_EQ(hs_test::context_depth(), 1);
  HS_EXPECT_STREQ(hs_test::context_frames()[0].label, "argument_case(calls)");
}

template <int A, int B> static void template_case() { HS_EXPECT_EQ(A + B, 3); }

int main() {
  HS_RUN(lambda_case);
  HS_EXPECT_EQ(hs_test::context_depth(), 0);
  int calls = 0;
  HS_RUN(argument_case, calls);
  HS_EXPECT_EQ(calls, 1);
  HS_RUN((template_case<1, 2>));
  HS_EXPECT_EQ(hs_test::context_depth(), 0);
  return hs_test::stats().failed == 1 ? 0 : 1;
}
