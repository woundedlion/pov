#include "pullback_manifest.generated.h"

#include <cstddef>
#include <cstdint>

int main() {
  static_assert(PullbackManifest::PROGRAMS.size() == 11);
  static_assert(PullbackManifest::ORACLES.size() == 2);
  static_assert(PullbackManifest::BASE_SHA.size() == 40);
  static_assert(PullbackManifest::MANIFEST_SHA256.size() == 64);
  uint16_t preset_mask = 0;
  for (const PullbackManifest::ProgramEntry &program :
       PullbackManifest::PROGRAMS)
    preset_mask |= program.preset_mask;
  return preset_mask == 0x0fff ? 0 : 1;
}
