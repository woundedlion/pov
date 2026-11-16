"""Read and validate firmware memory-budget declarations."""

from __future__ import annotations

import json
from pathlib import Path

from teensy_memory import MEMORY_MAP


def _strip_jsonc_comments(text: str) -> str:
    """Remove line comments from JSON text, leaving strings untouched."""
    out = []
    i, n = 0, len(text)
    in_string = False
    while i < n:
        c = text[i]
        if in_string:
            out.append(c)
            if c == "\\" and i + 1 < n:
                # Escaped char: copy verbatim; it can't terminate the string.
                out.append(text[i + 1])
                i += 2
                continue
            if c == '"':
                in_string = False
            i += 1
        elif c == '"':
            in_string = True
            out.append(c)
            i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "/":
            i += 2  # // line comment: drop through end of line (keep the \n)
            while i < n and text[i] != "\n":
                i += 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def declares_components(budget: dict) -> bool:
    """True if any region of `budget` carries per-component ceilings.

    The `size -A` fallback synthesizes region totals with no component
    breakdown, so such a budget cannot be evaluated from it at all.
    """
    return any(spec.get("components")
               for spec in budget.get("regions", {}).values())


class BudgetSchemaError(ValueError):
    """A budgets entry carries a key the gate would never read, or omits a
    required region / layout symbol."""


# Every key the gate reads, per nesting level. Region names, component names and
# layout-symbol keys stay free-form: a typo there is already a loud
# region-missing / component-missing / symbol-not-found violation.
_BUDGET_KEYS = frozenset({"regions", "symbols"})
_REGION_KEYS = frozenset({"max_bytes", "free_min_bytes", "components"})
_COMPONENT_KEYS = frozenset({"max_bytes", "max_banks_from_stack_floor"})
# Every key is also required: the derivation reads all three unconditionally.
_DERIVED_KEYS = frozenset(
    {"bank_bytes", "total_banks", "min_headroom_bytes"})
_SYMBOL_KEYS = frozenset({"name", "region", "min_bytes", "max_bytes"})

# Keys whose absence leaves an object schema-valid but inert (every ceiling is
# read with `.get()`).
_REGION_REQUIRED_KEYS = frozenset({"max_bytes"})
_SYMBOL_REQUIRED_KEYS = frozenset({"name", "region"})

# A component budget declares at least one ceiling.
_COMPONENT_ONE_OF_KEYS = frozenset({"max_bytes", "max_banks_from_stack_floor"})

# Regions whose 'max_bytes' is the hardware size, which the linker hits first,
# so the free floor is required.
_REGION_REQUIRED_KEYS_BY_REGION: dict[str, frozenset[str]] = {
    "ram2": _REGION_REQUIRED_KEYS | {"free_min_bytes"},
}

# Layout symbols whose invariant is magnitude: a shrunken arena frees RAM1 and
# trips no other check.
_SYMBOL_REQUIRED_KEYS_BY_SYMBOL: dict[str, frozenset[str]] = {
    "arena": _SYMBOL_REQUIRED_KEYS | {"min_bytes", "max_bytes"},
}

# Region objects and layout symbols every target budget must declare; evaluate()
# enforces only what the budget carries.
_REQUIRED_REGIONS = frozenset({"flash", "ram1", "ram2"})
_REQUIRED_SYMBOLS = frozenset({"arena", "framebuffer_a", "framebuffer_b"})

# Layout symbols required of specific targets only.
_REQUIRED_SYMBOLS_BY_ENV: dict[str, frozenset[str]] = {
    "phantasm": frozenset({"reaction_graph", "reaction_graph_index",
                           "reaction_graph_count", "dma_tx_buffer"}),
}
_REQUIRED_COMPONENTS_BY_ENV: dict[str, dict[str, frozenset[str]]] = {
    "phantasm": {"ram1": frozenset({"code"})},
}
_COMPONENT_REQUIRED_KEYS_BY_COMPONENT: dict[str, frozenset[str]] = {
    "code": frozenset({"max_banks_from_stack_floor"}),
}


def _check_keys(spec: object, allowed: frozenset[str], where: str,
                required: frozenset[str] = frozenset(),
                one_of: frozenset[str] = frozenset()) -> dict:
    """Reject unknown / missing keys in one budgets object; return it."""
    if not isinstance(spec, dict):
        raise BudgetSchemaError(
            f"{where}: expected an object, got {type(spec).__name__}.")
    unknown = sorted(set(spec) - allowed)
    if unknown:
        raise BudgetSchemaError(
            f"{where}: unknown key(s) {', '.join(repr(k) for k in unknown)} - "
            f"the gate reads only {', '.join(sorted(allowed))}. Every ceiling is "
            f"optional, so an unrecognized key silently disables the check it "
            f"was meant to configure.")
    missing = sorted(required - set(spec))
    if missing:
        raise BudgetSchemaError(
            f"{where}: missing required key(s) "
            f"{', '.join(repr(k) for k in missing)}.")
    null_required = sorted(key for key in required if spec[key] is None)
    if null_required:
        raise BudgetSchemaError(
            f"{where}: required key(s) cannot be null: "
            f"{', '.join(repr(k) for k in null_required)}.")
    if one_of and not any(spec.get(key) is not None for key in one_of):
        raise BudgetSchemaError(
            f"{where}: must declare at least one of "
            f"{', '.join(repr(k) for k in sorted(one_of))} - an object "
            f"carrying neither ceiling passes the schema and enforces nothing.")
    return spec


def _check_byte_limits(spec: dict, keys: tuple[str, ...], where: str,
                       required: frozenset[str] = frozenset()) -> None:
    for key in keys:
        value = spec.get(key)
        if value is None:
            continue
        minimum = 1 if key in required else 0
        if type(value) is not int or value < minimum:
            domain = "positive" if minimum else "non-negative"
            raise BudgetSchemaError(
                f"{where}: '{key}' must be a {domain} integer.")
    lo, hi = spec.get("min_bytes"), spec.get("max_bytes")
    if lo is not None and hi is not None and lo > hi:
        raise BudgetSchemaError(
            f"{where}: 'min_bytes' must not exceed 'max_bytes'.")


def _child_map(parent: dict, key: str, where: str) -> dict:
    """Fetch an optional nested object, rejecting a non-object value."""
    value = parent.get(key)
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise BudgetSchemaError(
            f"{where}: '{key}' must be an object, got {type(value).__name__}.")
    return value


def _require_present(present: dict, required: frozenset[str], where: str,
                     kind: str) -> None:
    """Reject a budgets mapping that omits a required entry."""
    missing = sorted(required - set(present))
    if missing:
        raise BudgetSchemaError(
            f"{where}: missing required {kind}(s) "
            f"{', '.join(repr(k) for k in missing)} - the gate iterates only "
            f"the entries the budget declares, so a deleted one removes its "
            f"ceiling / invariant and the gate reports PASS.")


def validate_budgets(budgets: object) -> dict:
    """Validate a parsed budgets mapping against the gate's schema."""
    if not isinstance(budgets, dict):
        raise BudgetSchemaError(
            f"budgets: expected an object of env -> budget, got "
            f"{type(budgets).__name__}.")
    for env, budget in budgets.items():
        _check_keys(budget, _BUDGET_KEYS, f"env '{env}'")
        regions = _child_map(budget, "regions", f"env '{env}'")
        _require_present(regions, _REQUIRED_REGIONS, f"env '{env}'", "region")
        for region, spec in regions.items():
            rwhere = f"env '{env}' region '{region}'"
            required = _REGION_REQUIRED_KEYS_BY_REGION.get(
                region, _REGION_REQUIRED_KEYS)
            _check_keys(spec, _REGION_KEYS, rwhere,
                        required=required)
            _check_byte_limits(spec, ("max_bytes", "free_min_bytes"),
                               rwhere, required)
            components = _child_map(spec, "components", rwhere)
            _require_present(
                components,
                _REQUIRED_COMPONENTS_BY_ENV.get(env, {}).get(
                    region, frozenset()),
                rwhere, "component")
            for cname, cspec in components.items():
                cwhere = f"{rwhere} component '{cname}'"
                _check_keys(cspec, _COMPONENT_KEYS, cwhere,
                            required=_COMPONENT_REQUIRED_KEYS_BY_COMPONENT.get(
                                cname, frozenset()),
                            one_of=_COMPONENT_ONE_OF_KEYS)
                _check_byte_limits(cspec, ("max_bytes",), cwhere)
                derived = cspec.get("max_banks_from_stack_floor")
                if derived is not None:
                    _check_keys(derived, _DERIVED_KEYS,
                                f"{cwhere} max_banks_from_stack_floor",
                                required=_DERIVED_KEYS)
                    min_headroom = derived["min_headroom_bytes"]
                    if (not isinstance(min_headroom, int) or
                            isinstance(min_headroom, bool) or
                            min_headroom < 0):
                        raise BudgetSchemaError(
                            f"{cwhere} max_banks_from_stack_floor: "
                            f"'min_headroom_bytes' must be a non-negative "
                            f"integer.")
                    # The bank geometry divides the derivation; a bad value
                    # would raise out of evaluate() with the violation code.
                    for key in ("bank_bytes", "total_banks"):
                        count = derived[key]
                        if (not isinstance(count, int) or
                                isinstance(count, bool) or count <= 0):
                            raise BudgetSchemaError(
                                f"{cwhere} max_banks_from_stack_floor: "
                                f"'{key}' must be a positive integer.")
        syms = _child_map(budget, "symbols", f"env '{env}'")
        _require_present(
            syms,
            _REQUIRED_SYMBOLS | _REQUIRED_SYMBOLS_BY_ENV.get(env, frozenset()),
            f"env '{env}'", "layout symbol")
        for key, spec in syms.items():
            swhere = f"env '{env}' symbol '{key}'"
            required = _SYMBOL_REQUIRED_KEYS_BY_SYMBOL.get(
                key, _SYMBOL_REQUIRED_KEYS)
            _check_keys(spec, _SYMBOL_KEYS, swhere, required=required)
            _check_byte_limits(spec, ("min_bytes", "max_bytes"), swhere, required)
            if not isinstance(spec["name"], str) or not spec["name"].strip():
                raise BudgetSchemaError(
                    f"{swhere}: 'name' must be a non-empty string.")
            if (not isinstance(spec["region"], str) or
                    spec["region"] not in {name for name, _, _ in MEMORY_MAP} | {"OTHER"}):
                raise BudgetSchemaError(
                    f"{swhere}: 'region' must name a known memory region.")
    return budgets


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise BudgetSchemaError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_budgets(path: str | Path) -> dict:
    """Load tools/teensy_budgets.json with optional line comments."""
    raw = Path(path).read_text(encoding="utf-8")
    return validate_budgets(json.loads(_strip_jsonc_comments(raw),
                                       object_pairs_hook=_unique_object))
