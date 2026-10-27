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
     {{0xeef0796b35fc207bull,
       {{3387524892ull, 0ull, 1291610249ull},
        334904181064ull,
        671031083453ull}},
      {0xd3ed8534a3a3c0c3ull,
       {{3387533729ull, 0ull, 1291610880ull},
        334905494734ull,
        671032545999ull}}},
     {{0x448c45bc7689b9e1ull,
       {{3430548532ull, 0ull, 1293235207ull},
        338345398609ull,
        676160025527ull}},
      {0x56bfc11b83edc88full,
       {{3430522125ull, 0ull, 1293232404ull},
        338343993751ull,
        676154836463ull}}},
     {{0xab2f61b06de4f335ull,
       {{921364614ull, 0ull, 350420563ull}, 130974334749ull, 184708235540ull}},
      {0x3c615f788aa9e8b9ull,
       {{921364407ull, 0ull, 350419297ull}, 130974151946ull, 184708259492ull}},
      0.001,
      0.025}},
    {"Dynamo",
     {{0x156b14f07cb00b21ull,
       {{396035167ull, 416799849ull, 124183317ull},
        23866900258ull,
        120265239757ull}},
      {0xfe53c3ba8b3285c2ull,
       {{396036715ull, 416800826ull, 124183479ull},
        23866916924ull,
        120265778614ull}}},
     {{0x3abd50e13f110816ull,
       {{401965812ull, 422139093ull, 125840104ull},
        23865250665ull,
        122170733415ull}},
      {0x4dee96fe6b829209ull,
       {{401959805ull, 422131756ull, 125837968ull},
        23865197885ull,
        122169032786ull}}},
     {{0x156b14f07cb00b21ull,
       {{396035167ull, 416799849ull, 124183317ull},
        23866900258ull,
        120265239757ull}},
      {0xfe53c3ba8b3285c2ull,
       {{396036715ull, 416800826ull, 124183479ull},
        23866916924ull,
        120265778614ull}}}},
    {"Fishbowl",
     {{0xbe8c99ad2b9e2658ull,
       {{436583740ull, 155182791ull, 63283427ull},
        50368287197ull,
        98347198502ull}},
      {0x521f7c057ce42447ull,
       {{436607419ull, 155184176ull, 63280558ull},
        50369928705ull,
        98353397011ull}}},
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
     {{0x17757f081dfd369dull,
       {{4148753139ull, 2287128305ull, 665179201ull},
        384008576219ull,
        1021508100845ull}},
      {0x53396c7069537a90ull,
       {{4147410252ull, 2286385617ull, 664978601ull},
        383884670155ull,
        1021181917759ull}}},
     {{0x16a3021ea493af52ull,
       {{4136992719ull, 2282720521ull, 666248722ull},
        384261503398ull,
        1019145656674ull}},
      {0x2a766799b7c4de47ull,
       {{4135656400ull, 2281984609ull, 666047974ull},
        384138456679ull,
        1018822763895ull}}},
     {{0x85863196bc09e538ull,
       {{129938039ull, 738145260ull, 316699748ull},
        131059483775ull,
        144255352687ull}},
      {0x6070d94724e78258ull,
       {{129895609ull, 737917199ull, 316609199ull},
        131020231105ull,
        144214067662ull}}}},
    {"HankinSolids",
     {{0x93c556c0cd01c6cbull,
       {{51096326162ull, 20556689001ull, 22164776162ull},
        6578967127937ull,
        13507147750863ull}},
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
     {{0x2edfcdc6789d0201ull,
       {{12976236670ull, 5133422967ull, 5709152603ull},
        2496729290498ull,
        3348262118094ull}},
      {0x077ffbfc653a5f87ull,
       {{12976236385ull, 5133423720ull, 5709152015ull},
        2496729294276ull,
        3348261813228ull}}}},
    {"HopfFibration",
     {{0x38d350fba24597a0ull,
       {{6851508713ull, 2032836043ull, 103070720ull},
        589808483302ull,
        1622681908231ull}},
      {0x439806d1b0cc4484ull,
       {{6851676885ull, 2032852904ull, 103075815ull},
        589824526983ull,
        1622736675523ull}}},
     {{0xde39c32a03287cc3ull,
       {{7889386020ull, 2372297080ull, 117391407ull},
        681170912161ull,
        1881784411794ull}},
      {0xdb9955d075930144ull,
       {{7889973195ull, 2372461609ull, 117395564ull},
        681202397228ull,
        1881905077709ull}}},
     {{0x122e27fdb1caf405ull,
       {{1677505253ull, 482438350ull, 30039266ull},
        220756174080ull,
        413591719164ull}},
      {0x87927c756f894badull,
       {{1677536578ull, 482438480ull, 30040966ull},
        220760239466ull,
        413598526378ull}}}},
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
     {{0x90611df1ee0811ddull,
       {{7428838334ull, 5263851795ull, 901779940ull},
        984400573817ull,
        1924639324174ull}},
      {0xfbe5f5b5eb81f5e8ull,
       {{7428780966ull, 5263797008ull, 901769641ull},
        984384384876ull,
        1924627968936ull}}},
     {{0x20674c39f14b693cull,
       {{7426028320ull, 5282426023ull, 911152059ull},
        973538725766ull,
        1947118629591ull}},
      {0x67671acf08daf21bull,
       {{7426015675ull, 5282421190ull, 911164592ull},
        973539211619ull,
        1947117308053ull}}},
     {{0x90611df1ee0811ddull,
       {{7428838334ull, 5263851795ull, 901779940ull},
        984400573817ull,
        1924639324174ull}},
      {0xfbe5f5b5eb81f5e8ull,
       {{7428780966ull, 5263797008ull, 901769641ull},
        984384384876ull,
        1924627968936ull}}}},
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
     {{0xd8321f5a30b713e5ull,
       {{177224536ull, 248277625ull, 591242785ull},
        73981424996ull,
        127674892612ull}},
      {0x932666bcc2b400fdull,
       {{177226734ull, 248278565ull, 591246984ull},
        73982244904ull,
        127675874362ull}}},
     {{0x00f0e15d48db5d71ull,
       {{177317518ull, 250333328ull, 590824222ull},
        74088283202ull,
        128042154053ull}},
      {0x0da342e762e8d5c0ull,
       {{177320777ull, 250337377ull, 590829405ull},
        74088851107ull,
        128042842837ull}}},
     {{0x8d43a967946b1300ull,
       {{58843772ull, 84428535ull, 204991443ull},
        34940961349ull,
        44478309016ull}},
      {0x938504a4c6e2ca8eull,
       {{58845561ull, 84429107ull, 204994667ull},
        34941695585ull,
        44479209493ull}}}},
    {"PetalFlow",
     {{0x2e1f446b4a53b805ull,
       {{222130307ull, 358485999ull, 355339207ull},
        67142681398ull,
        124519478675ull}},
      {0xf559e2e0bc787d5eull,
       {{222136298ull, 358497234ull, 355347473ull},
        67144363032ull,
        124523303203ull}}},
     {{0x1312633d0088b40full,
       {{251857489ull, 424850991ull, 423627755ull},
        78770040171ull,
        147514659413ull}},
      {0x80b59f436a2b0fe4ull,
       {{251866952ull, 424875588ull, 423644641ull},
        78773986929ull,
        147520649727ull}}},
     {{0x722bf08bee491fa6ull,
       {{72722988ull, 114629303ull, 111572386ull},
        27547862049ull,
        40166531878ull}},
      {0x10623e86e501088cull,
       {{72722066ull, 114628397ull, 111570155ull},
        27547307447ull,
        40166415210ull}}}},
    {"RingShower",
     {{0xe6b468b39c66e1e0ull,
       {{55682989ull, 79374227ull, 91275599ull},
        11731109338ull,
        25481017760ull}},
      {0x6be800863cce7b63ull,
       {{55685410ull, 79382264ull, 91280735ull},
        11731550332ull,
        25482511912ull}},
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
     {{0x82622108b88e4fd5ull,
       {{16188710ull, 22040007ull, 20422628ull}, 4892461452ull, 5365971910ull}},
      {0xfb044c3bf293025dull,
       {{16188643ull, 22039848ull, 20422476ull},
        4892429774ull,
        5365921207ull}}}},
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
     {{0x4d4c8a02a82cbf2aull,
       {{10096038003ull, 2188976309ull, 352284387ull},
        902431436622ull,
        1514810896650ull}},
      {0x77560a39dda92621ull,
       {{10096032609ull, 2188976446ull, 352304114ull},
        902432513143ull,
        1514812852172ull}}},
     {{0xe3bb5e344a6165b2ull,
       {{10196840911ull, 2210889935ull, 352427582ull},
        911267015550ull,
        1528857863272ull}},
      {0x9a97936632102acdull,
       {{10196835284ull, 2210889801ull, 352446924ull},
        911267994166ull,
        1528859796596ull}}},
     {{0x006b9339fb8223baull,
       {{2640496704ull, 573917087ull, 91862010ull},
        296414285012ull,
        338263043550ull}},
      {0x9c9fa1971e08fa70ull,
       {{2640495858ull, 573917256ull, 91867125ull},
        296414736250ull,
        338263670027ull}}}},
    {"Thrusters",
     {{0x4142ef3789f0b3d9ull,
       {{164917678ull, 20727560ull, 73689934ull},
        18840263761ull,
        39899806495ull}},
      {0xaaa42142f01d2c66ull,
       {{164917620ull, 20727564ull, 73689784ull},
        18840247712ull,
        39899777339ull}}},
     {{0x04886484f4367d8full,
       {{164906368ull, 20489923ull, 73700538ull},
        18866796612ull,
        39749124979ull}},
      {0xac65d889d60adb09ull,
       {{164907014ull, 20490142ull, 73700336ull},
        18866842244ull,
        39749250493ull}}},
     {{0x26a7491160496413ull,
       {{56364288ull, 8313585ull, 17273634ull}, 7842625095ull, 11948154129ull}},
      {0x3278ca5199f76df7ull,
       {{56364205ull, 8313561ull, 17273632ull},
        7842611443ull,
        11948125884ull}}}},
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

/** @brief Optional platform fast-math references retain the native fallback. */
inline void test_platform_fastmath_golden_reference() {
  RecordedGolden recorded;
  recorded.fastmath.fingerprint.channels = {10, 20, 30};
  recorded.linux_fastmath.hash = 17;
  HS_EXPECT_EQ(&fastmath_golden_reference(recorded, recorded.linux_fastmath),
               &recorded.fastmath);
  HS_EXPECT_EQ(&fastmath_golden_reference(recorded, recorded.wasm_fastmath),
               &recorded.fastmath);
  recorded.linux_fastmath.fingerprint.channels = {40, 50, 60};
  HS_EXPECT_EQ(&fastmath_golden_reference(recorded, recorded.linux_fastmath),
               &recorded.linux_fastmath);
  HS_EXPECT_EQ(&fastmath_golden_reference(recorded, recorded.wasm_fastmath),
               &recorded.fastmath);
  recorded.wasm_fastmath.fingerprint.channels = {70, 80, 90};
  HS_EXPECT_EQ(&fastmath_golden_reference(recorded, recorded.wasm_fastmath),
               &recorded.wasm_fastmath);
  HS_EXPECT_EQ(&fastmath_golden_reference(recorded, recorded.linux_fastmath),
               &recorded.linux_fastmath);
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
  HS_RUN(test_platform_fastmath_golden_reference);
  HS_RUN(test_render_fingerprint_rejects_channel_and_spatial_changes);
  HS_RUN(test_filter_pipeline_effects_match_recorded_frames);
  return fixture.result();
}

} // namespace plot_parity_tests
} // namespace hs_test
