/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file choreography_hooks.h
 * @brief Explicit compile-time bindings for choreography customizations. */

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Control {

/** @brief Choreography customization points. */
enum class ChoreographyHook {
  PARAMETER_FIELDS,     ///< Static parameter descriptions.
  VALID_PARAMS,         ///< Static parameter validation.
  PRESET_ROW,           ///< Static authored preset row.
  INITIAL_PARAMS,       ///< Static startup parameters.
  PRESET_PARAMS,        ///< Runtime preset target parameters.
  ADOPT_PARAMS,         ///< Snap target adoption.
  BLEND_PARAMS,         ///< In-flight interpolation.
  TRANSITION_ARMED,     ///< Lerp endpoint-state capture.
  FINISH_BLEND,         ///< Automatic target adoption.
  PRESET_OPACITY,       ///< Fade opacity update.
  AFTER_PROFILE_PRESET, ///< Profiling preset adjustment.
  AFTER_COMPOSED_INIT   ///< Composed lifecycle initialization.
};

/** @brief Ownership of a composed effect's hook declarations. */
enum class HookOwnership {
  AUTHORED, ///< Derived declares its owned hook list.
  INHERITED ///< Composed base supplies its default hook list.
};

namespace HookDetail {
/** @brief Function-pointer signature traits. @tparam Function Pointer type. */
template <typename Function> struct Signature;
/** @brief Static function signature.
 * @tparam ReturnT Return type.
 * @tparam Args Argument types. */
template <typename ReturnT, typename... Args>
struct Signature<ReturnT (*)(Args...)> {
  using Return = ReturnT;                ///< Return type.
  using Class = void;                    ///< No member owner.
  using Arguments = std::tuple<Args...>; ///< Exact argument types.
  static constexpr bool MEMBER = false;  ///< Static function.
};
/** @brief Member function signature.
 * @tparam ReturnT Return type.
 * @tparam Owner Member owner.
 * @tparam Args Argument types. */
template <typename ReturnT, typename Owner, typename... Args>
struct Signature<ReturnT (Owner::*)(Args...)>
    : Signature<ReturnT (*)(Args...)> {
  using Class = Owner;                 ///< Member owner.
  static constexpr bool MEMBER = true; ///< Member function.
};
template <typename ReturnT, typename Owner, typename... Args>
struct Signature<ReturnT (Owner::*)(Args...) const>
    : Signature<ReturnT (Owner::*)(Args...)> {};
template <typename ReturnT, typename... Args>
struct Signature<ReturnT (*)(Args...) noexcept>
    : Signature<ReturnT (*)(Args...)> {};
template <typename ReturnT, typename Owner, typename... Args>
struct Signature<ReturnT (Owner::*)(Args...) noexcept>
    : Signature<ReturnT (Owner::*)(Args...)> {};
template <typename ReturnT, typename Owner, typename... Args>
struct Signature<ReturnT (Owner::*)(Args...) const noexcept>
    : Signature<ReturnT (Owner::*)(Args...)> {};
} // namespace HookDetail

/** @brief Binds one customization point to a function.
 * @tparam Kind Customization point.
 * @tparam Function Static or member function pointer. */
template <ChoreographyHook Kind, auto Function> struct HookBinding {
  static constexpr ChoreographyHook KIND = Kind; ///< Customization point.
  static constexpr auto FUNCTION = Function;     ///< Bound function pointer.
  using Signature =
      HookDetail::Signature<decltype(Function)>; ///< Function signature.
};

/** @brief Type-only list of explicitly bound hooks.
 * @tparam OwnerT Type declaring this list.
 * @tparam Bindings HookBinding entries. */
template <typename OwnerT, typename... Bindings> struct HookSet {
  using Owner = OwnerT; ///< Declaration owner.
  static_assert(
      [] {
        constexpr ChoreographyHook KINDS[] = {
            Bindings::KIND..., ChoreographyHook::PARAMETER_FIELDS};
        for (size_t i = 0; i < sizeof...(Bindings); ++i)
          for (size_t j = i + 1; j < sizeof...(Bindings); ++j)
            if (KINDS[i] == KINDS[j])
              return false;
        return true;
      }(),
      "choreography hooks: duplicate binding");

  /** @brief Whether a customization point is bound.
   * @tparam Kind Customization point. */
  template <ChoreographyHook Kind>
  static constexpr bool HAS = ((Bindings::KIND == Kind) || ...);

  /** @brief Adds bindings and assigns their declaration owner.
   * @tparam NewOwner New declaration owner.
   * @tparam More Additional HookBinding entries.
   * @return Extended type-only list. */
  template <typename NewOwner, typename... More>
  static consteval auto extend() {
    return HookSet<NewOwner, Bindings..., More...>{};
  }

  /** @brief The function bound to a customization point.
   * @tparam Kind Customization point.
   * @return Bound function pointer. */
  template <ChoreographyHook Kind> static consteval auto function() {
    static_assert(HAS<Kind>, "choreography hooks: unbound hook");
    return find<Kind, Bindings...>();
  }

  /** @brief Invokes a bound hook on an effect.
   * @tparam Kind Customization point.
   * @tparam Effect Effect type.
   * @tparam Args Argument types.
   * @param effect Effect receiving a member call.
   * @param args Hook arguments.
   * @return Hook result. */
  template <ChoreographyHook Kind, typename Effect, typename... Args>
  static constexpr decltype(auto) invoke(Effect &effect, Args &&...args) {
    constexpr auto FUNCTION = function<Kind>();
    if constexpr (std::is_member_function_pointer_v<decltype(FUNCTION)>)
      return (effect.*FUNCTION)(std::forward<Args>(args)...);
    else
      return FUNCTION(std::forward<Args>(args)...);
  }

  /** @brief Invokes a static hook.
   * @tparam Kind Customization point.
   * @tparam Args Argument types.
   * @param args Hook arguments.
   * @return Hook result. */
  template <ChoreographyHook Kind, typename... Args>
  static constexpr decltype(auto) invoke_static(Args &&...args) {
    return function<Kind>()(std::forward<Args>(args)...);
  }

  /** @brief Checks exact argument and return types of the bound hooks.
   * @tparam Params Effect parameter aggregate.
   * @tparam Row Preset row type.
   * @return True when every binding has its declared signature. */
  template <typename Params, typename Row> static consteval bool valid() {
    return (valid_binding<Bindings, Params, Row>() && ...);
  }

private:
  template <ChoreographyHook Kind, typename First, typename... Rest>
  static consteval auto find() {
    if constexpr (First::KIND == Kind)
      return First::FUNCTION;
    else
      return find<Kind, Rest...>();
  }

  template <typename Binding, typename Params, typename Row>
  static consteval bool valid_binding() {
    using Signature = typename Binding::Signature;
    using Args = typename Signature::Arguments;
    using Return = typename Signature::Return;
    constexpr auto KIND = Binding::KIND;
    if constexpr (Signature::MEMBER &&
                  !std::is_base_of_v<typename Signature::Class, Owner>)
      return false;
    else if constexpr (KIND == ChoreographyHook::PARAMETER_FIELDS)
      return !Signature::MEMBER && std::is_same_v<Args, std::tuple<>>;
    else if constexpr (KIND == ChoreographyHook::VALID_PARAMS)
      return !Signature::MEMBER &&
             std::is_same_v<Args, std::tuple<const Params &>> &&
             std::is_same_v<Return, bool>;
    else if constexpr (KIND == ChoreographyHook::PRESET_ROW)
      return !Signature::MEMBER && std::is_same_v<Args, std::tuple<size_t>> &&
             std::is_same_v<Return, Row>;
    else if constexpr (KIND == ChoreographyHook::INITIAL_PARAMS)
      return !Signature::MEMBER && std::is_same_v<Args, std::tuple<>> &&
             std::is_same_v<Return, Params>;
    else if constexpr (KIND == ChoreographyHook::PRESET_PARAMS)
      return Signature::MEMBER && std::is_same_v<Args, std::tuple<size_t>> &&
             std::is_same_v<Return, Params>;
    else if constexpr (KIND == ChoreographyHook::ADOPT_PARAMS ||
                       KIND == ChoreographyHook::TRANSITION_ARMED ||
                       KIND == ChoreographyHook::FINISH_BLEND)
      return Signature::MEMBER &&
             std::is_same_v<Args, std::tuple<const Params &>> &&
             std::is_same_v<Return, void>;
    else if constexpr (KIND == ChoreographyHook::BLEND_PARAMS ||
                       KIND == ChoreographyHook::PRESET_OPACITY)
      return Signature::MEMBER && std::is_same_v<Args, std::tuple<float>> &&
             std::is_same_v<Return, void>;
    else if constexpr (KIND == ChoreographyHook::AFTER_PROFILE_PRESET ||
                       KIND == ChoreographyHook::AFTER_COMPOSED_INIT)
      return Signature::MEMBER && std::is_same_v<Args, std::tuple<>> &&
             std::is_same_v<Return, void>;
    else
      return false;
  }
};
} // namespace Control
