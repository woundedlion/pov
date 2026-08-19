/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <cstring>
#include <new>

#include "render/pullback/operator_table.h"

/**
 * @file interpreter.h
 * @brief ChainProgram: the transactional two-arena stage-program interpreter.
 */

namespace Pullback {

namespace Interp {

inline constexpr size_t MAX_CHAIN_OPS = 32;
inline constexpr size_t MAX_INSTANCE_ID = 48;
/** One arena's capacity — the exported budget; the engine holds two. */
inline constexpr size_t CHAIN_ARENA_BYTES = 96 * 1024;
inline constexpr size_t MAX_CHAIN_PARAMS = 224;
/** Fixed per-entry arena cost beyond the cataloged blocks: the instance-id
    copy reservation. */
inline constexpr size_t PER_OP_OVERHEAD_BYTES = MAX_INSTANCE_ID + 1;
/** Longest schema field id, in bytes, excluding the NUL; pinned so the
    per-field name reservation below is a table-independent constant. */
inline constexpr size_t MAX_FIELD_ID = 31;
/** Fixed arena reservation per schema field for its registered
    "{instance}.{field-id}" parameter name: worst-case instance id, the dot,
    worst-case field id, and the NUL. */
inline constexpr size_t PER_PARAM_NAME_BYTES =
    MAX_INSTANCE_ID + 1 + MAX_FIELD_ID + 1;

namespace Detail {

template <typename List> struct SlotTraits;
template <typename... Carriers>
struct SlotTraits<Pullback::Detail::TypeList<Carriers...>> {
  static_assert((std::is_trivially_destructible_v<Carriers> && ...),
                "chain interpreter: every canonical carrier must be trivially "
                "destructible");
  static constexpr size_t SIZE = std::max({sizeof(Carriers)...});
  static constexpr size_t ALIGN = std::max({alignof(Carriers)...});
};

} // namespace Detail

/** @brief Ping-pong slot ABI, derived from the closed carrier set. */
inline constexpr size_t SLOT_SIZE = Detail::SlotTraits<CarrierList>::SIZE;
inline constexpr size_t SLOT_ALIGN = Detail::SlotTraits<CarrierList>::ALIGN;

enum class ChainStatus : uint8_t {
  OK,
  NOT_CHAIN_EFFECT,
  MALFORMED_PAYLOAD,
  EMPTY,
  TOO_LONG,
  UNKNOWN_OPERATOR,
  DUPLICATE_INSTANCE,
  MALFORMED_INSTANCE,
  ENTRY_FAMILY,
  EXIT_FAMILY,
  CARRIER_MISMATCH,
  ARENA_OVERFLOW,
  PARAM_OVERFLOW,
  MIGRATE_FAILED
};

/** @brief Wire spelling of @p status; the embind boundary renames a committed
    OK to "APPLIED". */
inline const char *chain_status_name(ChainStatus status) {
  switch (status) {
  case ChainStatus::OK:
    return "OK";
  case ChainStatus::NOT_CHAIN_EFFECT:
    return "NOT_CHAIN_EFFECT";
  case ChainStatus::MALFORMED_PAYLOAD:
    return "MALFORMED_PAYLOAD";
  case ChainStatus::EMPTY:
    return "EMPTY";
  case ChainStatus::TOO_LONG:
    return "TOO_LONG";
  case ChainStatus::UNKNOWN_OPERATOR:
    return "UNKNOWN_OPERATOR";
  case ChainStatus::DUPLICATE_INSTANCE:
    return "DUPLICATE_INSTANCE";
  case ChainStatus::MALFORMED_INSTANCE:
    return "MALFORMED_INSTANCE";
  case ChainStatus::ENTRY_FAMILY:
    return "ENTRY_FAMILY";
  case ChainStatus::EXIT_FAMILY:
    return "EXIT_FAMILY";
  case ChainStatus::CARRIER_MISMATCH:
    return "CARRIER_MISMATCH";
  case ChainStatus::ARENA_OVERFLOW:
    return "ARENA_OVERFLOW";
  case ChainStatus::PARAM_OVERFLOW:
    return "PARAM_OVERFLOW";
  case ChainStatus::MIGRATE_FAILED:
    return "MIGRATE_FAILED";
  }
  __builtin_unreachable();
}

/** @brief Structured refusal: the code plus the offending entry index
    (-1 refuses the chain as a whole). */
struct ChainRefusal {
  ChainStatus code;
  int16_t entry_index;
};

/** Every shipped schema field id must fit the fixed name reservation. */
consteval bool operator_schema_ids_fit_names() {
  for (const OperatorDescriptor &op : OPERATOR_TABLE)
    for (uint16_t index = 0; index < op.schema_count; ++index)
      if (std::string_view(op.schema[index].id).size() > MAX_FIELD_ID)
        return false;
  return true;
}
static_assert(operator_schema_ids_fit_names(),
              "chain operator table: a schema field id exceeds MAX_FIELD_ID");

/** @brief One wire chain entry; the views are caller-owned and only read
    during compile(). */
struct ChainEntryRequest {
  std::string_view instance_id;
  std::string_view operator_id;
};

/** @brief Whether @p id matches [a-z][a-z0-9]*(-[a-z0-9]+)* within the length
    cap: one kebab segment, no dots. */
inline bool instance_id_wellformed(std::string_view id) {
  if (id.empty() || id.size() > MAX_INSTANCE_ID)
    return false;
  if (id.front() < 'a' || id.front() > 'z')
    return false;
  bool previous_dash = false;
  for (size_t index = 1; index < id.size(); ++index) {
    const char c = id[index];
    if (c == '-') {
      if (previous_dash)
        return false;
      previous_dash = true;
      continue;
    }
    if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9')))
      return false;
    previous_dash = false;
  }
  return !previous_dash;
}

/**
 * @brief An arena-backed chain program: transactional compile into the
 *        inactive arena, per-frame advance/prepare, ping-pong evaluate.
 * @details Storage is two caller-owned fixed blocks bound once; steady-state
 * recompiles never allocate. Compilation is transactional: any refusal leaves
 * the previous program and all live instance state untouched, tearing down
 * only the candidate. No heap anywhere.
 */
class ChainProgram {
public:
  /** @brief One compiled entry; offsets index the owning side's block. */
  struct ChainOp {
    const OperatorDescriptor *op;
    const char *instance; /**< Copy in the owning arena. */
    uint32_t stable_hash;
    uint32_t param_offset;
    uint32_t prepared_offset;
    uint32_t state_offset;
    uint32_t name_offset; /**< "{instance}.{field-id}" slots, PER_PARAM_NAME_
                               BYTES apart, one per schema field. */
  };

  /**
   * @brief Binds the two program arenas and the operator table.
   * @param block_a First arena block, at least @p block_capacity bytes.
   * @param block_b Second arena block, same size.
   * @param block_capacity One arena's capacity; the budget compile() checks.
   * @param operator_table Resolution table; defaults to OPERATOR_TABLE.
   * @details Blocks must be max_align-aligned and outlive the program. Bind
   * once, before the first compile().
   */
  void bind_storage(uint8_t *block_a, uint8_t *block_b,
                    size_t block_capacity = CHAIN_ARENA_BYTES,
                    std::span<const OperatorDescriptor> operator_table =
                        std::span<const OperatorDescriptor>(OPERATOR_TABLE)) {
    HS_CHECK(block_a != nullptr && block_b != nullptr,
             "ChainProgram::bind_storage: null block");
    HS_CHECK(reinterpret_cast<uintptr_t>(block_a) % alignof(std::max_align_t) ==
                 0,
             "ChainProgram::bind_storage: block_a misaligned");
    HS_CHECK(reinterpret_cast<uintptr_t>(block_b) % alignof(std::max_align_t) ==
                 0,
             "ChainProgram::bind_storage: block_b misaligned");
    HS_CHECK(blocks[0] == nullptr, "ChainProgram::bind_storage: rebinding");
    blocks[0] = block_a;
    blocks[1] = block_b;
    capacity = block_capacity;
    table = operator_table;
  }

  /**
   * @brief Compiles a program shape transactionally.
   * @param request Ordered {instance_id, operator_id} entries.
   * @return {OK, -1} on commit; otherwise the refusal, with the previous
   *         program, its parameter blocks and all live instance state intact.
   * @details Order: shape checks against the request alone, then budgets
   * against one arena's capacity, then layout into the inactive arena, then
   * state migration/init (a failing migrate tears the candidate down as a
   * unit), then the active-index flip and loser teardown. Parameter blocks
   * are default-constructed; values re-apply through the value channel after
   * commit.
   */
  ChainRefusal compile(std::span<const ChainEntryRequest> request) {
    HS_CHECK(blocks[0] != nullptr, "ChainProgram::compile before bind_storage");
    if (request.empty())
      return {ChainStatus::EMPTY, -1};
    if (request.size() > MAX_CHAIN_OPS)
      return {ChainStatus::TOO_LONG, -1};
    const OperatorDescriptor *resolved[MAX_CHAIN_OPS];
    for (size_t index = 0; index < request.size(); ++index) {
      if (!instance_id_wellformed(request[index].instance_id))
        return {ChainStatus::MALFORMED_INSTANCE, static_cast<int16_t>(index)};
      for (size_t previous = 0; previous < index; ++previous)
        if (request[previous].instance_id == request[index].instance_id)
          return {ChainStatus::DUPLICATE_INSTANCE, static_cast<int16_t>(index)};
      resolved[index] = nullptr;
      for (const OperatorDescriptor &candidate : table)
        if (request[index].operator_id == candidate.operator_id) {
          resolved[index] = &candidate;
          break;
        }
      if (resolved[index] == nullptr)
        return {ChainStatus::UNKNOWN_OPERATOR, static_cast<int16_t>(index)};
    }
    if (resolved[0]->input != CarrierId::SPHERE)
      return {ChainStatus::ENTRY_FAMILY, 0};
    if (resolved[request.size() - 1]->output != CarrierId::COLOR)
      return {ChainStatus::EXIT_FAMILY,
              static_cast<int16_t>(request.size() - 1)};
    for (size_t index = 1; index < request.size(); ++index)
      if (resolved[index]->input != resolved[index - 1]->output)
        return {ChainStatus::CARRIER_MISMATCH, static_cast<int16_t>(index)};

    size_t field_total = 0;
    for (size_t index = 0; index < request.size(); ++index)
      field_total += resolved[index]->schema_count;
    if (field_total > MAX_CHAIN_PARAMS)
      return {ChainStatus::PARAM_OVERFLOW, -1};
    if (plan_layout(resolved, request, nullptr) > capacity)
      return {ChainStatus::ARENA_OVERFLOW, -1};

    Side &candidate = sides[active ^ 1];
    candidate.count = 0;
    candidate.used_bytes = plan_layout(resolved, request, &candidate);

    const Side &current = sides[active];
    for (size_t index = 0; index < candidate.count; ++index) {
      ChainOp &op = candidate.ops[index];
      const InstanceId id{op.instance, op.op->operator_id, op.stable_hash};
      const ChainOp *survivor = nullptr;
      if (has_program)
        for (size_t existing = 0; existing < current.count; ++existing) {
          const ChainOp &live = current.ops[existing];
          if (std::string_view(live.instance) == op.instance &&
              std::string_view(live.op->operator_id) == op.op->operator_id) {
            survivor = &live;
            break;
          }
        }
      if (survivor != nullptr) {
        const Status migrated = op.op->runtime.migrate(
            block_ptr(active ^ 1, op.state_offset),
            block_ptr(active, survivor->state_offset), id);
        if (migrated != Status::OK) {
          // The candidate arena tears down as a unit: every state compile
          // already constructed, inits and completed migrates alike.
          for (size_t constructed = 0; constructed < index; ++constructed)
            candidate.ops[constructed].op->runtime.destroy(
                block_ptr(active ^ 1, candidate.ops[constructed].state_offset));
          candidate.count = 0;
          return {ChainStatus::MIGRATE_FAILED, static_cast<int16_t>(index)};
        }
      } else {
        op.op->runtime.init(block_ptr(active ^ 1, op.state_offset), id);
      }
    }

    Side &loser = sides[active];
    if (has_program)
      for (size_t index = 0; index < loser.count; ++index)
        loser.ops[index].op->runtime.destroy(
            block_ptr(active, loser.ops[index].state_offset));
    loser.count = 0;
    loser.used_bytes = 0;
    active ^= 1;
    has_program = true;
    needs_prepare = true;
    return {ChainStatus::OK, -1};
  }

  /** @brief Steps every op's per-frame clocks; once per frame. */
  void advance() {
    Side &side = sides[active];
    for (size_t index = 0; index < side.count; ++index) {
      ChainOp &op = side.ops[index];
      op.op->runtime.advance(block_ptr(active, op.state_offset),
                             const_block_ptr(active, op.param_offset));
    }
  }

  /** @brief Resolves every op's prepared block from @p ctx. */
  void prepare(const FrameContext &ctx) {
    Side &side = sides[active];
    for (size_t index = 0; index < side.count; ++index) {
      ChainOp &op = side.ops[index];
      op.op->runtime.prepare(
          ctx, const_block_ptr(active, op.param_offset),
          block_ptr(active, op.state_offset),
          static_cast<uint8_t *>(block_ptr(active, op.prepared_offset)));
    }
    needs_prepare = false;
  }

  /** @brief Seeds the entry carrier and runs the chain over @p view. */
  Color4 evaluate(const Vector &view, const FrameContext &ctx) const {
    HS_CHECK(has_program, "ChainProgram::evaluate without a program");
    HS_CHECK(!needs_prepare, "ChainProgram::evaluate before prepare()");
    alignas(SLOT_ALIGN) uint8_t slot_a[SLOT_SIZE];
    alignas(SLOT_ALIGN) uint8_t slot_b[SLOT_SIZE];
    void *in = slot_a;
    void *out = slot_b;
    ::new (in) SphereSample{view, 0.0f};
    const Side &side = sides[active];
    for (size_t index = 0; index < side.count; ++index) {
      const ChainOp &op = side.ops[index];
      op.op->runtime.run(in, out, ctx, const_block_ptr(active, op.param_offset),
                         const_block_ptr(active, op.prepared_offset));
      std::swap(in, out);
    }
    return *std::launder(reinterpret_cast<const Color4 *>(in));
  }

  /** @brief The committed program's entries, in chain order. */
  std::span<const ChainOp> ops() const {
    return {sides[active].ops, sides[active].count};
  }

  bool compiled() const { return has_program; }

  /** @brief Bytes the committed program occupies in its arena. */
  size_t used_bytes() const { return sides[active].used_bytes; }

  /** @brief Mutable param block of entry @p index; the value channel's
      storage and the registration target. */
  uint8_t *param_block(size_t index) {
    HS_CHECK(index < sides[active].count, "ChainProgram::param_block range");
    return static_cast<uint8_t *>(
        block_ptr(active, sides[active].ops[index].param_offset));
  }

  /** @brief Registered parameter name "{instance}.{field-id}" of schema entry
      @p schema_index of entry @p index.
      @details Storage lives in the winning arena, so the pointer dies at the
      next successful compile — the caller must re-register immediately after
      every commit, before anything reads a previously registered name. */
  const char *param_name(size_t index, uint16_t schema_index) const {
    HS_CHECK(index < sides[active].count, "ChainProgram::param_name range");
    const ChainOp &op = sides[active].ops[index];
    HS_CHECK(schema_index < op.op->schema_count,
             "ChainProgram::param_name schema range");
    return reinterpret_cast<const char *>(const_block_ptr(
        active, op.name_offset + schema_index * PER_PARAM_NAME_BYTES));
  }

  /** @brief Instance state of entry @p index (LUT bake inputs, inspection). */
  const void *state_block(size_t index) const {
    HS_CHECK(index < sides[active].count, "ChainProgram::state_block range");
    return const_block_ptr(active, sides[active].ops[index].state_offset);
  }

  /** @brief Prepared block of entry @p index; valid after prepare(). */
  const uint8_t *prepared_block(size_t index) const {
    HS_CHECK(index < sides[active].count, "ChainProgram::prepared_block range");
    return const_block_ptr(active, sides[active].ops[index].prepared_offset);
  }

  /** @brief Destroys the committed program's instance state and empties the
      program; storage stays bound. */
  void clear() {
    if (!has_program)
      return;
    Side &side = sides[active];
    for (size_t index = 0; index < side.count; ++index)
      side.ops[index].op->runtime.destroy(
          block_ptr(active, side.ops[index].state_offset));
    side.count = 0;
    side.used_bytes = 0;
    has_program = false;
  }

private:
  struct Side {
    ChainOp ops[MAX_CHAIN_OPS];
    size_t count = 0;
    size_t used_bytes = 0;
  };

  static size_t align_up(size_t offset, size_t alignment) {
    return (offset + alignment - 1) & ~(alignment - 1);
  }

  void *block_ptr(uint8_t side, uint32_t offset) {
    return blocks[side] + offset;
  }
  const uint8_t *const_block_ptr(uint8_t side, uint32_t offset) const {
    return blocks[side] + offset;
  }

  /**
   * @brief Computes the exact arena layout of a resolved chain.
   * @param out Side to fill with offsets, instance-id copies, and parameter
   *        name slots, or null for the budget dry run — both paths share this
   *        arithmetic, so the budget check equals the layout.
   * @return Total bytes consumed.
   */
  size_t plan_layout(const OperatorDescriptor *const *resolved,
                     std::span<const ChainEntryRequest> request, Side *out) {
    size_t cursor = 0;
    uint8_t *base = out != nullptr ? blocks[active ^ 1] : nullptr;
    for (size_t index = 0; index < request.size(); ++index) {
      const OperatorRuntime &runtime = resolved[index]->runtime;
      cursor = align_up(cursor, runtime.param.align);
      const size_t param_offset = cursor;
      cursor += runtime.param.size;
      cursor = align_up(cursor, runtime.prepared.align);
      const size_t prepared_offset = cursor;
      cursor += runtime.prepared.size;
      cursor = align_up(cursor, runtime.state.align);
      const size_t state_offset = cursor;
      cursor += runtime.state.size;
      const size_t id_offset = cursor;
      cursor += PER_OP_OVERHEAD_BYTES;
      const size_t name_offset = cursor;
      cursor += static_cast<size_t>(resolved[index]->schema_count) *
                PER_PARAM_NAME_BYTES;
      if (out != nullptr) {
        char *id_copy = reinterpret_cast<char *>(base + id_offset);
        const std::string_view instance = request[index].instance_id;
        std::memcpy(id_copy, instance.data(), instance.size());
        id_copy[instance.size()] = '\0';
        char *names = reinterpret_cast<char *>(base + name_offset);
        for (uint16_t field = 0; field < resolved[index]->schema_count;
             ++field) {
          const std::string_view field_id{resolved[index]->schema[field].id};
          HS_CHECK(field_id.size() <= MAX_FIELD_ID,
                   "chain compile: schema field id exceeds MAX_FIELD_ID");
          char *slot = names + field * PER_PARAM_NAME_BYTES;
          std::memcpy(slot, instance.data(), instance.size());
          slot[instance.size()] = '.';
          std::memcpy(slot + instance.size() + 1, field_id.data(),
                      field_id.size());
          slot[instance.size() + 1 + field_id.size()] = '\0';
        }
        runtime.construct_params(base + param_offset);
        out->ops[index] = ChainOp{
            resolved[index],
            id_copy,
            instance_hash(instance, resolved[index]->operator_id),
            static_cast<uint32_t>(param_offset),
            static_cast<uint32_t>(prepared_offset),
            static_cast<uint32_t>(state_offset),
            static_cast<uint32_t>(name_offset),
        };
        out->count = index + 1;
      }
    }
    return cursor;
  }

  Side sides[2];
  uint8_t *blocks[2] = {nullptr, nullptr};
  size_t capacity = 0;
  std::span<const OperatorDescriptor> table;
  uint8_t active = 0;
  bool has_program = false;
  bool needs_prepare = false;
};

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
