/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <cstring>
#include <string_view>
#include "core/control/params.h"
#include "core/render/pullback/composed_descriptors.h"

/** @file composed_admission.h
 * @brief `ComposedAdmission`, the GUI write admission of a composed effect.
 */

namespace Pullback {

/**
 * @brief Admits GUI edits of a composed effect's parameter families and
 *        modulation documents; empty unless the GUI bridge is built.
 * @details CRTP base of `Effect`, which befriends it.
 * @tparam Effect The composed effect.
 */
template <typename Effect> class ComposedAdmission {
#if HS_ENABLE_PARAM_GUI_BRIDGE
public:
  /** @brief Visits current typed modulation resources without retaining them.
   * @tparam Visitor Receives the exact resource key and immutable recipe.
   * @param visitor Synchronous resource visitor. */
  template <typename Visitor>
  void visit_document_modulations(Visitor &&visitor) const {
    effect().params.visit([&]<typename Resource>(const auto &family) {
      if constexpr (Resource::KIND == ResourceKind::MODULATION)
        visitor(Resource::KEY.view(), family.modulation);
    });
  }

  /** @brief Admits an edit to one current modulation resource.
   * @tparam Editor Edits a candidate and returns its write result.
   * @param resource Exact resource key.
   * @param editor Synchronous candidate editor.
   * @return The editor's result, or a resource/admission refusal. */
  template <typename Editor>
  ParamSetResult edit_document_modulation(std::string_view resource,
                                          Editor &&editor) {
    ParamSetResult result = ParamSetResult::UNKNOWN_PARAM;
    effect().params.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == ResourceKind::MODULATION) {
        if (Resource::KEY.view() != resource)
          return;
        auto candidate = family.modulation;
        result = editor(candidate);
        if (result != ParamSetResult::APPLIED)
          return;
        if ((!Resource::Family::DYNAMIC &&
             !Modulation::same_topology(family.modulation, candidate)) ||
            !Modulation::admitted(candidate)) {
          result = ParamSetResult::INADMISSIBLE;
          return;
        }
        family.modulation = candidate;
        effect().setAnimationsPaused(true);
        effect().parameter_written();
      }
    });
    return result;
  }

protected:
  /**
   * @brief Warning text for a slider whose last write was refused.
   * @param name Parameter name.
   * @return The family warning for the last refused parameter, else null.
   */
  const char *refused_parameter_warning(const char *name) const {
    return refused_name != nullptr && std::strcmp(name, refused_name) == 0
               ? refused_warning
               : nullptr;
  }

  /**
   * @brief Refuses a write that would make its parameter family invalid.
   * @param parameter Descriptor being written.
   * @param value Proposed value.
   * @return False when the write is refused.
   */
  bool admit_parameter_write(const ParamDef &parameter, float value) {
    bool admitted = true;
    effect().params.visit([&]<typename Resource>(const auto &family) {
      if constexpr (Resource::KIND == ResourceKind::LENS ||
                    Resource::KIND == ResourceKind::SURFACE ||
                    Resource::KIND == ResourceKind::COLOR ||
                    Resource::KIND == ResourceKind::MODULATION) {
        const uintptr_t begin = reinterpret_cast<uintptr_t>(&family);
        const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
        if (target >= begin && target - begin < sizeof(family)) {
          auto candidate = family;
          ParamDef proposed = parameter;
          proposed.target =
              reinterpret_cast<unsigned char *>(&candidate) + (target - begin);
          Effect::write_parameter_unchecked(proposed, value);
          refused_warning = Pullback::refusal(candidate);
          admitted = refused_warning == nullptr;
        }
      }
    });
    if (!admitted) {
      refused_name = parameter.name;
      return false;
    }
    refused_name = nullptr;
    return true;
  }

private:
  Effect &effect() { return static_cast<Effect &>(*this); }
  const Effect &effect() const { return static_cast<const Effect &>(*this); }

  /// Name of the last refused parameter, or null after an admitted write.
  const char *refused_name = nullptr;
  /// Warning associated with the last refused family write.
  const char *refused_warning = nullptr;
#endif
};

} // namespace Pullback
