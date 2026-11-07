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
  entry.identity.display_name = name;
  entry.identity.machine_id = stable_id;
  return entry;
}

// registration_names_unique rejects an empty name or stable ID.
static_assert(!registration_names_unique(std::array{
    named_registration("DeathEmptyId", "")}));
static_assert(!registration_names_unique(std::array{
    named_registration("", "death-empty-name")}));

/**
 * @brief Death case: a class name equal to another effect's stable ID must
 *        trap.
 */
inline void case_effect_registry_name_matches_stable_id() {
  validate_effect_registrations(
      std::array{named_registration("DeathFirst", "DeathPersistedAlias"),
                 named_registration("DeathPersistedAlias", "death-second")});
}

/**
 * @brief Death case: registering two effects under one name must trap.
 */
inline void case_effect_registry_duplicate_name() {
  validate_effect_registrations(
      std::array{named_registration("DeathDuplicate", "death-dup-a"),
                 named_registration("DeathDuplicate", "death-dup-b")});
}

/** @brief Death case: two effects declaring the same stable ID must trap. */
inline void case_effect_registry_duplicate_stable_id() {
  validate_effect_registrations(
      std::array{named_registration("DeathStableA", "death-stable"),
                 named_registration("DeathStableB", "death-stable")});
}

/** @brief Death case: a stable ID equal to another effect's name must trap. */
inline void case_effect_registry_stable_id_matches_name() {
  validate_effect_registrations(
      std::array{named_registration("DeathClassAlias", "death-first"),
                 named_registration("DeathOther", "DeathClassAlias")});
}
