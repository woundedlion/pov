/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Pixel parity of every filter-pipeline effect against recorded frame hashes.
 */
#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include "core/render/spacetime.h"
#include "hardware/pov_segment_map.h"
#include "targets/effects.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace plot_parity_tests {

/** @brief Render width of every parity scenario. */
constexpr int PARITY_W = 288;
/** @brief Render height of every parity scenario. */
constexpr int PARITY_H = 144;
/** @brief Frames folded into each scenario's hash. */
constexpr int PARITY_FRAMES = 48;
/** @brief Segments the clipped scenario's display is divided into. */
constexpr int PARITY_SEGMENTS = 4;
/** @brief Segment the clipped scenario renders. */
constexpr int PARITY_SEGMENT_ID = 1;

/** @brief Recorded render sequences for one effect across the parity scenarios. */
struct ParityHashes {
  const char *effect;       /**< Effect name. */
  RecordedGolden temporal;  /**< Full frames, temporal sampling enabled. */
  RecordedGolden frame_end; /**< Full frames, temporal sampling disabled. */
  RecordedGolden segment_clipped; /**< One display segment's clip per frame. */
};

/** @brief Computed sequences for one effect's parity scenarios. */
struct ParityRenders {
  const char *effect;
  PixelGolden temporal;
  PixelGolden frame_end;
  PixelGolden segment_clipped;
};

/** @brief Recorded temporal, frame-end and segment-clipped render sequences. */
inline constexpr ParityHashes PARITY_HASHES[] = {
    {"Comets",
     {{0x8b8536ade45aba11ull,
       {{172892154ull, 169225874ull, 291812097ull},
        46754718597ull,
        44923120391ull}},
      {0x9a2b810fd287b024ull,
       {{172870233ull, 169201730ull, 291771823ull},
        46749147387ull,
        44917245103ull}}},
     {{0x7d05d3cd68d1cbb9ull,
       {{204934466ull, 190186630ull, 320628647ull},
        51402881614ull,
        49449636374ull}},
      {0xaac5ff81ea1847c2ull,
       {{204908186ull, 190159673ull, 320585057ull},
        51396855544ull,
        49443198522ull}}},
     {{0x7c59ab4460ae24e3ull,
       {{52911474ull, 38924709ull, 70432523ull},
        15867213715ull,
        13013449707ull}},
      {0x702011234b2a480cull,
       {{52906405ull, 38920834ull, 70425601ull},
        15865693364ull,
        13012176838ull}}}},
    {"DisplacementField",
     {{0xc4bc19266261803dull,
       {{2787083657ull, 2695513000ull, 836806405ull},
        450897683462ull,
        906526832411ull}},
      {0xfc3c76f9f754aa8cull,
       {{2787083568ull, 2695513336ull, 836806778ull},
        450897730540ull,
        906526868524ull}}},
     {{0x47fa576d9a94eb76ull,
       {{2814892854ull, 2721147015ull, 845090131ull},
        455301753460ull,
        914315475807ull}},
      {0xe52d3552a02a5bacull,
       {{2814892936ull, 2721146911ull, 845090372ull},
        455301750293ull,
        914315475888ull}}},
     {{0xc0cc816e5cbabb92ull,
       {{818129378ull, 781269291ull, 243072894ull},
        191221248945ull,
        265291098058ull}},
      {0xd3ff0c99319f8491ull,
       {{818129667ull, 781269796ull, 243073120ull},
        191221331806ull,
        265291231249ull}}}},
    {"DreamBalls",
     {{0x9c69e67167f317a2ull,
       {{3223945041ull, 0ull, 1249196611ull},
        319905595071ull,
        642543300551ull}},
      {0x0f86f70ebec4f9a3ull,
       {{3223973493ull, 0ull, 1249201930ull},
        319909006098ull,
        642547303863ull}}},
     {{0x448c45bc7689b9e1ull,
       {{3430548532ull, 0ull, 1293235207ull},
        338345398609ull,
        676160025527ull}},
      {0x56bfc11b83edc88full,
       {{3430522125ull, 0ull, 1293232404ull},
        338343993751ull,
        676154836463ull}}},
     {{0xad6de67b304e0a6aull,
       {{877116409ull, 0ull, 339827076ull}, 123296319362ull, 176424641669ull}},
      {0x1532ecf59ca7ecd9ull,
       {{877116553ull, 0ull, 339823893ull}, 123295910512ull, 176423850253ull}},
      0.001,
      0.025}},
    {"Dynamo",
     {{0x3abd50e13f110816ull,
       {{401965812ull, 422139093ull, 125840104ull},
        23865250665ull,
        122170733415ull}},
      {0x4dee96fe6b829209ull,
       {{401959805ull, 422131756ull, 125837968ull},
        23865197885ull,
        122169032786ull}}},
     {{0x3abd50e13f110816ull,
       {{401965812ull, 422139093ull, 125840104ull},
        23865250665ull,
        122170733415ull}},
      {0x4dee96fe6b829209ull,
       {{401959805ull, 422131756ull, 125837968ull},
        23865197885ull,
        122169032786ull}}},
     {{0x3abd50e13f110816ull,
       {{401965812ull, 422139093ull, 125840104ull},
        23865250665ull,
        122170733415ull}},
      {0x4dee96fe6b829209ull,
       {{401959805ull, 422131756ull, 125837968ull},
        23865197885ull,
        122169032786ull}}}},
    {"Fishbowl",
     {{0xe7715a41d1d9b309ull,
       {{436224441ull, 155153696ull, 63283327ull},
        50313714350ull,
        98343222973ull}},
      {0x31ed1af1437d1a70ull,
       {{436248087ull, 155155081ull, 63280457ull},
        50315350845ull,
        98349420264ull}}},
     {{0xe470fb4fe3b1eed2ull,
       {{496608704ull, 176173768ull, 71952184ull},
        57068880631ull,
        110121445273ull}},
      {0x38d60878c8b029cfull,
       {{496611792ull, 176173916ull, 71952064ull},
        57069270621ull,
        110121680429ull}}},
     {{0x4ea8cf3a5cefad6full,
       {{141170590ull, 50838477ull, 21074010ull},
        22481057347ull,
        28395846726ull}},
      {0xd499c268880adf66ull,
       {{141167077ull, 50823154ull, 21031193ull},
        22477083243ull,
        28383868176ull}},
      0.005,
      0.03}},
    {"Galaxies",
     {{0xe5a998917a1e8e25ull,
       {{765242538ull, 750982966ull, 758812506ull},
        161201577515ull,
        322654514559ull}},
      {0x3a878926e03a80a7ull,
       {{765200649ull, 750941047ull, 758770590ull},
        161192518433ull,
        322637149432ull}}},
     {{0x4637663771062f8full,
       {{728771023ull, 714497557ull, 722269521ull},
        153407654202ull,
        309685769331ull}},
      {0x1a32fb7059d74036ull,
       {{728729745ull, 714456293ull, 722228198ull},
        153398744211ull,
        309668578051ull}}},
     {{0xfeb6f3baf57701b6ull,
       {{263432933ull, 257854735ull, 259531567ull},
        79566008860ull,
        100162522764ull}},
      {0x8c73049a65315255ull,
       {{263419848ull, 257841618ull, 259518433ull},
        79561541490ull,
        100158045109ull}}}},
    {"GnomonicStars",
     {{0x81bb46212d3052a6ull,
       {{4149279533ull, 2287397446ull, 665172744ull},
        384043592310ull,
        1021650624965ull}},
      {0xfc766ac66d5b2db6ull,
       {{4147939011ull, 2286656164ull, 664972493ull},
        383919970319ull,
        1021325076058ull}}},
     {{0x16a3021ea493af52ull,
       {{4136992719ull, 2282720521ull, 666248722ull},
        384261503398ull,
        1019145656674ull}},
      {0x2a766799b7c4de47ull,
       {{4135656400ull, 2281984609ull, 666047974ull},
        384138456679ull,
        1018822763895ull}}},
     {{0x19e77fedacd5f287ull,
       {{129932712ull, 738207572ull, 316694037ull},
        131064937704ull,
        144280265357ull}},
      {0x2485b886e9f6467aull,
       {{129890587ull, 737980178ull, 316603826ull},
        131025822186ull,
        144239089811ull}}}},
    {"HankinSolids",
     {{0x9cb8fc1dd116511dull,
       {{51096339711ull, 20556697841ull, 22164777386ull},
        6578967900036ull,
        13507150281694ull}},
      {0xbce80fc0fb38b526ull,
       {{51096337729ull, 20556700999ull, 22164774983ull},
        6578968338896ull,
        13507150748490ull}}},
     {{0x530957c44028e83bull,
       {{52045641609ull, 20313181878ull, 21906244358ull},
        6604299715281ull,
        13574336121858ull}},
      {0x95d1e5ae44d09b08ull,
       {{52045635842ull, 20313182622ull, 21906242134ull},
        6604299329264ull,
        13574335018959ull}}},
     {{0x0cdcd78b958faf19ull,
       {{12976237619ull, 5133423594ull, 5709152697ull},
        2496729450233ull,
        3348262188445ull}},
      {0x077ffbfc653a5f87ull,
       {{12976236385ull, 5133423720ull, 5709152015ull},
        2496729294276ull,
        3348261813228ull}}}},
    {"HopfFibration",
     {{0xcd1f5106e2ab7670ull,
       {{6846649521ull, 2031499668ull, 102927963ull},
        589312491240ull,
        1621896930474ull}},
      {0x5ece142cc16585cdull,
       {{6846854380ull, 2031528781ull, 102933213ull},
        589328552345ull,
        1621957783668ull}}},
     {{0xde39c32a03287cc3ull,
       {{7889386020ull, 2372297080ull, 117391407ull},
        681170912161ull,
        1881784411794ull}},
      {0xdb9955d075930144ull,
       {{7889973195ull, 2372461609ull, 117395564ull},
        681202397228ull,
        1881905077709ull}}},
     {{0xba7b1f8954e81e32ull,
       {{1676000839ull, 482079390ull, 29970668ull},
        220481116842ull,
        413403740546ull}},
      {0x6c98130728c8bfa8ull,
       {{1676032262ull, 482079560ull, 29972372ull},
        220485202201ull,
        413410568154ull}}}},
    {"IslamicStars",
     {{0x96b8c46b2f41db54ull,
       {{20514060628ull, 13061337707ull, 11724899181ull},
        3297521524496ull,
        6542034908525ull}},
      {0x620dbf50829d4cdbull,
       {{20514056975ull, 13061337027ull, 11724893603ull},
        3297520759497ull,
        6542033494061ull}}},
     {{0x66b52924497c1972ull,
       {{21284175998ull, 13233937867ull, 11795349947ull},
        3369947135083ull,
        6666789759476ull}},
      {0x721a158e6e619894ull,
       {{21284151814ull, 13233953733ull, 11795353253ull},
        3369946715432ull,
        6666789058733ull}}},
     {{0x8ef51aa75425b232ull,
       {{5270716509ull, 3501368213ull, 3171629901ull},
        1304734204421ull,
        1791789268435ull}},
      {0x3c500e3a11a0e6deull,
       {{5270715524ull, 3501367982ull, 3171628342ull},
        1304733834374ull,
        1791788915235ull}}}},
    {"MeshFeedback",
     {{0xc562f2896902ca10ull,
       {{6842427473ull, 4877140903ull, 850861680ull},
        905602387243ull,
        1784300196852ull}},
      {0xc6430066319f568dull,
       {{6842432059ull, 4877112135ull, 850848297ull},
        905601599220ull,
        1784293666687ull}}},
     {{0x20674c39f14b693cull,
       {{7426028320ull, 5282426023ull, 911152059ull},
        973538725766ull,
        1947118629591ull}},
      {0x67671acf08daf21bull,
       {{7426015675ull, 5282421190ull, 911164592ull},
        973539211619ull,
        1947117308053ull}}},
     {{0xc562f2896902ca10ull,
       {{6842427473ull, 4877140903ull, 850861680ull},
        905602387243ull,
        1784300196852ull}},
      {0xc6430066319f568dull,
       {{6842432059ull, 4877112135ull, 850848297ull},
        905601599220ull,
        1784293666687ull}}}},
    {"MindSplatter",
     {{0x2584dda7eedf202eull,
       {{321487231ull, 189659287ull, 294994098ull},
        58558422538ull,
        113031486146ull}},
      {0x2e68ad4891689869ull,
       {{321517439ull, 189680026ull, 295013816ull},
        58570505235ull,
        113049885455ull}},
      0.001,
      0.03},
     {{0x2d2ed377b25b0f51ull,
       {{383602584ull, 228490737ull, 344437828ull},
        69593929666ull,
        133221574453ull}},
      {0xa4ff267717354839ull,
       {{383631591ull, 228502409ull, 344475106ull},
        69598133215ull,
        133234979122ull}}},
     {{0xddec7ce73b346ca2ull,
       {{89764594ull, 51158528ull, 80489940ull},
        22873164291ull,
        34217516832ull}},
      {0x7406edb900344f70ull,
       {{89768814ull, 51161411ull, 80500249ull},
        22875467042ull,
        34222399672ull}},
      0.001,
      0.03}},
    {"MobiusRings",
     {{0x67d2b1e0df5fd6feull,
       {{178117751ull, 251522000ull, 593368873ull},
        74428538493ull,
        128702953045ull}},
      {0x3a9df2b742fa7236ull,
       {{178121027ull, 251526115ull, 593373725ull},
        74429143546ull,
        128703602706ull}}},
     {{0x00f0e15d48db5d71ull,
       {{177317518ull, 250333328ull, 590824222ull},
        74088283202ull,
        128042154053ull}},
      {0x0da342e762e8d5c0ull,
       {{177320777ull, 250337377ull, 590829405ull},
        74088851107ull,
        128042842837ull}}},
     {{0x00b3fed29e618219ull,
       {{58854347ull, 85295677ull, 205936344ull},
        35100668099ull,
        44672992878ull}},
      {0xc393b331421358a8ull,
       {{58855275ull, 85296324ull, 205938581ull},
        35101222112ull,
        44674219510ull}}}},
    {"PetalFlow",
     {{0x83741882ab777c5cull,
       {{222054493ull, 358384195ull, 355218898ull},
        67119597655ull,
        124484309217ull}},
      {0xb8c1e9d6f281faddull,
       {{222060491ull, 358395221ull, 355226661ull},
        67121177532ull,
        124488075472ull}}},
     {{0x1312633d0088b40full,
       {{251857489ull, 424850991ull, 423627755ull},
        78770040171ull,
        147514659413ull}},
      {0x80b59f436a2b0fe4ull,
       {{251866952ull, 424875588ull, 423644641ull},
        78773986929ull,
        147520649727ull}}},
     {{0x0bb0443594b408d8ull,
       {{72707811ull, 114596899ull, 111493616ull},
        27529925389ull,
        40145681957ull}},
      {0xf6bf48e270538696ull,
       {{72706890ull, 114595780ull, 111490893ull},
        27529270837ull,
        40145507354ull}}}},
    {"RingShower",
     {{0x7624bf028d7c940aull,
       {{53346736ull, 78042168ull, 88473434ull},
        11277465447ull,
        24698293606ull}},
      {0xd54068f6624d8a17ull,
       {{53356836ull, 78053811ull, 88494956ull},
        11279329338ull,
        24700934656ull}},
      0.001,
      0.04},
     {{0x7624bf028d7c940aull,
       {{53346736ull, 78042168ull, 88473434ull},
        11277465447ull,
        24698293606ull}},
      {0xd54068f6624d8a17ull,
       {{53356836ull, 78053811ull, 88494956ull},
        11279329338ull,
        24700934656ull}},
      0.001,
      0.04},
     {{0x685b21e41b9424d0ull,
       {{14053870ull, 20640578ull, 18929111ull}, 4507106461ull, 4884882489ull}},
      {0x3181142ae4d30cccull,
       {{14053787ull, 20640635ull, 18929348ull},
        4507122601ull,
        4884905246ull}}}},
    {"RingSpin",
     {{0x6c0a0a7a0d9b5f37ull,
       {{10285788237ull, 5132252692ull, 2160504337ull},
        1256866073320ull,
        2607641516863ull}},
      {0x161b15be5129365dull,
       {{10286029879ull, 5132277292ull, 2160424559ull},
        1256879422362ull,
        2607668217060ull}}},
     {{0x67e9653ae68abed2ull,
       {{10480462304ull, 5207345772ull, 2193893865ull},
        1278541711006ull,
        2642914455422ull}},
      {0xa6413c609d7fb5e6ull,
       {{10480634488ull, 5207328956ull, 2193831181ull},
        1278548395970ull,
        2642927734486ull}}},
     {{0xbaff0e1259f62c26ull,
       {{2682370673ull, 1325021899ull, 564697615ull},
        492585843081ull,
        819439663828ull}},
      {0x591ea86458fb35fcull,
       {{2682431296ull, 1325028348ull, 564676325ull},
        492593519574ull,
        819440051825ull}}}},
    {"ShapeShifter",
     {{0x776048d2aefa8d00ull,
       {{15314518083ull, 3096011298ull, 727963878ull},
        1367589396922ull,
        2354887172707ull}},
      {0x3be17111ca0edf4cull,
       {{15314516902ull, 3096021654ull, 727963231ull},
        1367589934462ull,
        2354887630331ull}}},
     {{0x47d2347a023fe8dfull,
       {{15314395345ull, 3095913116ull, 727964304ull},
        1367573636585ull,
        2354868049363ull}},
      {0xca76a3c14c761399ull,
       {{15314394185ull, 3095923502ull, 727963658ull},
        1367574177663ull,
        2354868511799ull}}},
     {{0xdb96e15259c022f0ull,
       {{4034868311ull, 815070839ull, 195160576ull},
        456477063309ull,
        545799058758ull}},
      {0xc636851155147df7ull,
       {{4034866669ull, 815073387ull, 195159315ull},
        456476878718ull,
        545798747312ull}}}},
    {"Thrusters",
     {{0x04886484f4367d8full,
       {{164906368ull, 20489923ull, 73700538ull},
        18866796612ull,
        39749124979ull}},
      {0xac65d889d60adb09ull,
       {{164907014ull, 20490142ull, 73700336ull},
        18866842244ull,
        39749250493ull}}},
     {{0x04886484f4367d8full,
       {{164906368ull, 20489923ull, 73700538ull},
        18866796612ull,
        39749124979ull}},
      {0xac65d889d60adb09ull,
       {{164907014ull, 20490142ull, 73700336ull},
        18866842244ull,
        39749250493ull}}},
     {{0xae12ff3819824aafull,
       {{56718367ull, 8304818ull, 18133530ull}, 7946569184ull, 12075564744ull}},
      {0x574e2c900d48b210ull,
       {{56718338ull, 8304875ull, 18133519ull},
        7946566522ull,
        12075578579ull}}}},
};

/**
 * @brief Renders PARITY_FRAMES frames of one effect and folds every displayed
 *        frame into one fingerprint.
 * @tparam E Effect class template.
 * @param temporal Temporal sampling selection for the whole run.
 * @param segment_clipped Render one display segment's clip each frame, where
 *        the segmented driver clips the effect.
 * @param lit Cleared when no displayed frame has a lit pixel.
 * @return The folded sequence fingerprint.
 */
template <template <int, int> class E>
inline PixelGolden render_fold(bool temporal, bool segment_clipped, bool &lit) {
  reset_globals();
  Spacetime::set_flat_surface(!temporal);
  pin_frame_clock(0);
  E<PARITY_W, PARITY_H> effect;
  effect.init();
  constexpr int S = PARITY_H * 2;
  const pov::SegmentMap map =
      pov::segment_map(PARITY_SEGMENT_ID, S, PARITY_SEGMENTS);
  const bool clip_applies = pov::segment_clip_applies(effect.needs_full_frame(),
                                                      effect.persists_pixels());
  PixelGolden fold;
  bool any_lit = false;
  for (int f = 0; f < PARITY_FRAMES; ++f) {
    if (segment_clipped && clip_applies) {
      const pov::SegmentClip clip =
          pov::segment_clip(map, (f & 1) == 0, S, PARITY_SEGMENTS, PARITY_W);
      effect.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
    }
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    fold_pixel_golden<PARITY_W, PARITY_H>(fold, effect);
    any_lit = any_lit || count_lit_region<PARITY_W, PARITY_H>(effect) != 0;
  }
  lit = lit && any_lit;
  Spacetime::set_flat_surface(false);
  hs::clear_mock_time();
  return fold;
}

/**
 * @brief Renders every parity scenario of one effect.
 * @tparam E Effect class template.
 * @param name Effect name.
 * @return The effect's render fingerprints.
 */
template <template <int, int> class E>
inline ParityRenders render_hashes(const char *name) {
  bool lit = true;
  const ParityRenders hashes{name, render_fold<E>(true, false, lit),
                             render_fold<E>(false, false, lit),
                             render_fold<E>(true, true, lit)};
  HS_CONTEXT(name);
  HS_EXPECT(lit, "every parity scenario must display a lit pixel");
  return hashes;
}

/**
 * @brief Every filter-pipeline effect matches its recorded frame hashes.
 */
inline void test_filter_pipeline_effects_match_recorded_frames() {
  const ParityRenders computed[] = {
      render_hashes<Comets>("Comets"),
      render_hashes<DisplacementField>("DisplacementField"),
      render_hashes<DreamBalls>("DreamBalls"),
      render_hashes<Dynamo>("Dynamo"),
      render_hashes<Fishbowl>("Fishbowl"),
      render_hashes<Galaxies>("Galaxies"),
      render_hashes<GnomonicStars>("GnomonicStars"),
      render_hashes<HankinSolids>("HankinSolids"),
      render_hashes<HopfFibration>("HopfFibration"),
      render_hashes<IslamicStars>("IslamicStars"),
      render_hashes<MeshFeedback>("MeshFeedback"),
      render_hashes<MindSplatter>("MindSplatter"),
      render_hashes<MobiusRings>("MobiusRings"),
      render_hashes<PetalFlow>("PetalFlow"),
      render_hashes<RingShower>("RingShower"),
      render_hashes<RingSpin>("RingSpin"),
      render_hashes<ShapeShifter>("ShapeShifter"),
      render_hashes<Thrusters>("Thrusters"),
  };
  std::printf("  float environment: 0x%016llxull\n",
              static_cast<unsigned long long>(float_environment_fingerprint()));
  HS_EXPECT_EQ(std::size(computed), std::size(PARITY_HASHES));
  for (size_t i = 0; i < std::size(computed) && i < std::size(PARITY_HASHES);
       ++i) {
    const ParityHashes &want = PARITY_HASHES[i];
    const ParityRenders &got = computed[i];
    HS_CONTEXT(got.effect);
    HS_EXPECT(std::strcmp(want.effect, got.effect) == 0,
              "recorded effect order");
    expect_pixel_golden("temporal", got.temporal, want.temporal);
    expect_pixel_golden("frame_end", got.frame_end, want.frame_end);
    expect_pixel_golden("segment_clipped", got.segment_clipped,
                        want.segment_clipped);
  }
}

/** @brief Spatial fingerprints reject energy-preserving image changes. */
inline void test_render_fingerprint_rejects_channel_and_spatial_changes() {
  struct Frame {
    std::array<Pixel, 35> pixels{};
    const Pixel &get_pixel(int x, int y) const { return pixels[y * 7 + x]; }
  };
  Frame original;
  original.pixels[1] = Pixel(1000, 4000, 7000);
  original.pixels[2 * 7 + 4] = Pixel(9000, 500, 3000);
  original.pixels[3 * 7 + 5] = Pixel(300, 11000, 600);
  auto capture = [](const Frame &frame) {
    PixelGolden golden;
    fold_pixel_golden<7, 5>(golden, frame);
    return golden;
  };
  const PixelGolden expected = capture(original);
  HS_EXPECT_EQ(expected.fingerprint.channels[0], uint64_t{10300});
  HS_EXPECT_EQ(expected.fingerprint.channels[1], uint64_t{15500});
  HS_EXPECT_EQ(expected.fingerprint.channels[2], uint64_t{10600});
  HS_EXPECT_EQ(expected.fingerprint.row, uint64_t{60700});
  HS_EXPECT_EQ(expected.fingerprint.column, uint64_t{121500});
  HS_EXPECT_TRUE(render_fingerprint_matches(expected.fingerprint,
                                            expected.fingerprint, 1e-3, 0.01));
  Frame rounded = original;
  ++rounded.pixels[1].r;
  HS_EXPECT_TRUE(render_fingerprint_matches(capture(rounded).fingerprint,
                                            expected.fingerprint, 1e-3, 0.01));
  for (int mutation = 0; mutation < 5; ++mutation) {
    HS_CONTEXT("mutation", mutation);
    Frame changed;
    for (int y = 0; y < 5; ++y)
      for (int x = 0; x < 7; ++x) {
        Pixel pixel = original.get_pixel(x, y);
        if (mutation == 0)
          std::swap(pixel.r, pixel.b);
        const int xx = mutation == 1 ? 6 - x : mutation == 3 ? (x + 1) % 7 : x;
        const int yy = mutation == 2 ? 4 - y : mutation == 4 ? (y + 1) % 5 : y;
        changed.pixels[yy * 7 + xx] = pixel;
      }
    const PixelGolden actual = capture(changed);
    HS_EXPECT_EQ(actual.fingerprint.energy(), expected.fingerprint.energy());
    HS_EXPECT_FALSE(render_fingerprint_matches(
        actual.fingerprint, expected.fingerprint, 1e-3, 0.01));
  }
  HS_EXPECT_FALSE(
      render_fingerprint_matches({}, expected.fingerprint, 1e-3, 0.01));
}

/**
 * @brief Runs every plot parity test case.
 * @return The module's failure count.
 */
inline int run_plot_parity_tests() {
  hs_test::ModuleFixture fixture("plot_parity");
  HS_RUN(test_render_fingerprint_rejects_channel_and_spatial_changes);
  HS_RUN(test_filter_pipeline_effects_match_recorded_frames);
  return fixture.result();
}

} // namespace plot_parity_tests
} // namespace hs_test
