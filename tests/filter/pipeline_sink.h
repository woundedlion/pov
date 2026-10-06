/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Pipeline sink + get<T>()
// ============================================================================

struct ReplacingTerminalWithoutHistory : Filter::Is2D {
  static constexpr bool is_terminal = true;
  static constexpr bool terminal_replaces = true;
  int flushes = 0;
  float last_alpha = 0.0f;
  void flush(Canvas &, float alpha) {
    ++flushes;
    last_alpha = alpha;
  }
};

inline void test_replacing_terminal_without_history_flushes() {
  StubEffect effect(32, 16);
  Canvas canvas(effect);
  Pipeline<32, 16, ReplacingTerminalWithoutHistory> pipeline;
  const auto &terminal = pipeline.get<ReplacingTerminalWithoutHistory>();
  static_assert(!decltype(pipeline)::has_history);
  (void)pipeline.begin_frame(canvas, 0.25f);
  HS_EXPECT_EQ(terminal.flushes, 1);
  HS_EXPECT_EQ(terminal.last_alpha, 0.25f);
  (void)pipeline.begin_frame(canvas, 0.75f);
  HS_EXPECT_EQ(terminal.flushes, 2);
  HS_EXPECT_EQ(terminal.last_alpha, 0.75f);
}

/**
 * @brief Verifies the bare (filter-free) pipeline sink is 2D.
 */
inline void test_pipeline_sink_is_2d() {
  HS_EXPECT_TRUE((Pipeline<32, 32>::is_2d));
  HS_EXPECT_FALSE((Pipeline<32, 32>::is_terminal));
  using Terminal = Pipeline<32, 32, Filter::Pixel::Feedback<32, 32>>;
  using NonTerminal = Pipeline<32, 32, Filter::Screen::AntiAlias<32, 32>>;
  using Prepared =
      std::remove_reference_t<decltype(std::declval<Terminal &>().begin_frame(
          std::declval<Canvas &>(), 1.0f))>;
  static_assert(Terminal::is_terminal);
  static_assert(!NonTerminal::is_terminal);
  static_assert(!RawFramePlotter<Terminal>);
  static_assert(!TerminalFlusher<Terminal>);
  static_assert(!TerminalFlusher<Filter::Pixel::Feedback<32, 32>>);
  using RetrievedFeedback = std::remove_reference_t<
      decltype(std::declval<Terminal &>()
                   .template get<Filter::Pixel::Feedback<32, 32>>())>;
  static_assert(!TerminalFlusher<RetrievedFeedback>);
  static_assert(ReplacementFrameStarter<Terminal>);
  static_assert(std::is_empty_v<Prepared>);
  static_assert(Filter::PipelineFoldSurface<Prepared>);
  static_assert(RawFramePlotter<Prepared>);
  static_assert(!TerminalFlusher<Prepared>);
  static_assert(!ReplacementFrameStarter<Prepared>);
  static_assert(RawFramePlotter<NonTerminal>);
  static_assert(!ReplacementFrameStarter<NonTerminal>);
}

/**
 * @brief Verifies get<T>() resolves each composed filter to the correctly-typed
 *        base subobject, for both the head node and tail nodes, in const and
 *        non-const pipelines.
 */
inline void test_pipeline_get_returns_correct_filter() {
  constexpr int W = 32, H = 32;
  using AA = Filter::Screen::AntiAlias<W, H>;
  using CS = Filter::Pixel::ChromaticShift<W>;
  using Blur = Filter::Screen::Blur<W, H>;

  Pipeline<W, H, AA, Blur, CS> pipe(AA{}, Blur{1.0f}, CS{});

  Blur &bl = pipe.get<Blur>();
  CS &cs = pipe.get<CS>();

  static_assert(std::is_same_v<decltype(pipe.get<AA>()), AA &>,
                "get<AA>() returns AA&");
  static_assert(std::is_same_v<decltype(pipe.get<Blur>()), Blur &>,
                "get<Blur>() returns Blur&");
  static_assert(std::is_same_v<decltype(pipe.get<CS>()), CS &>,
                "get<CS>() returns CS&");

  static_assert(!std::is_convertible_v<Pipeline<W, H, AA, Blur, CS> &, AA &>);
  static_assert(!std::is_convertible_v<Pipeline<W, H, Blur, CS> &, Blur &>);
  static_assert(!std::is_convertible_v<Pipeline<W, H, CS> &, CS &>);
  HS_EXPECT_TRUE(static_cast<const void *>(&bl) !=
                 static_cast<const void *>(&cs));

  const Pipeline<W, H, AA, Blur, CS> &cpipe = pipe;
  const Blur &cbl = cpipe.get<Blur>();
  HS_EXPECT_TRUE(&cbl == &bl);

  // get<T>() on an absent stage is a hard compile error, not SFINAE-detectable.
  using Absent = Filter::Pixel::Feedback<W, H>;

  // A duplicated stage type makes get<T>() ambiguous; stage_count is what the
  // guard in get<T>() reads. The ambiguous case is a hard compile error, so
  // only the count itself can be asserted here.
  static_assert(Pipeline<W, H, AA, Blur, CS>::stage_count<Blur> == 1);
  static_assert(Pipeline<W, H, AA, Blur, CS>::stage_count<Absent> == 0);
  static_assert(Pipeline<W, H, Blur, Blur, CS>::stage_count<Blur> == 2);
}
