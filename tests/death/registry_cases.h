/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Registry death cases.

/** @brief Builds a registration carrying only a name and stable ID. */
constexpr EffectRegistration named_registration(std::string_view name,
                                                std::string_view stable_id) {
  EffectRegistration entry{};
  entry.name = name;
  entry.stable_id = stable_id;
  return entry;
}

static_assert(!registration_names_unique(std::array{
    named_registration("DeathEmptyId", "")}));
static_assert(!registration_names_unique(std::array{
    named_registration("", "death-empty-name")}));

/**
 * @brief Death case: a class name equal to another effect's stable ID must
 *        trap.
 */
inline void case_effect_registry_name_matches_stable_id() {
  EffectRegistration first{};
  first.name = "DeathFirst";
  first.stable_id = "DeathPersistedAlias";

  EffectRegistration second{};
  second.name = "DeathPersistedAlias";
  second.stable_id = "death-second";
  validate_effect_registrations(std::array{first, second});
}

/**
 * @brief Death case: registering two effects under one name must trap.
 */
inline void case_effect_registry_duplicate_name() {
  EffectRegistration first{};
  first.name = "DeathDuplicate";
  first.stable_id = "death-dup-a";

  EffectRegistration second{};
  second.name = "DeathDuplicate";
  second.stable_id = "death-dup-b";
  validate_effect_registrations(std::array{first, second});
}

/** @brief Death case: two effects declaring the same stable ID must trap. */
inline void case_effect_registry_duplicate_stable_id() {
  EffectRegistration first{};
  first.name = "DeathStableA";
  first.stable_id = "death-stable";

  EffectRegistration second{};
  second.name = "DeathStableB";
  second.stable_id = "death-stable";
  validate_effect_registrations(std::array{first, second});
}

/** @brief Death case: a stable ID equal to another effect's name must trap. */
inline void case_effect_registry_stable_id_matches_name() {
  EffectRegistration first{};
  first.name = "DeathClassAlias";
  first.stable_id = "death-first";

  EffectRegistration second{};
  second.name = "DeathOther";
  second.stable_id = "DeathClassAlias";
  validate_effect_registrations(std::array{first, second});
}
