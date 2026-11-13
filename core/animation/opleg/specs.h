/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file specs.h
 * @brief Conway morph leg specifications and phase types. */

#include "core/mesh/conway_graph.h"
#include "core/color/palettes.h"

namespace Animation::detail {

/** Leg kind selected by the constructor. */
enum class LegKind : uint8_t {
  CONWAY_SWEEP,
  HANKIN_SWEEP,
  RELAX_SLERP,
  MEDIAL_SLERP,
  GATED_SWAP
};

/** Partition operator a GATED_SWAP leg swaps to. */
enum class SwapOp : uint8_t { KIS, DUAL };

/**
 * @brief Graph-edge Conway sweep: the edge's single op swept between the
 * edge's endpoint parameters.
 */
struct EdgeSweepSpec {
  const ConwayGraph::EdgeSpec *edge =
      nullptr;           /**< Graph edge being traversed. */
  bool reverse = false;  /**< True when traversing to_node -> from_node. */
  int sweep_frames = 0;  /**< Operator-sweep frames (N). */
  int settle_frames = 0; /**< Relax-slerp frames (S); positive on a settling
                            edge, 0 otherwise. */
};

/** @brief Hankin contact-angle sweep on a fixed seed. */
struct HankinSweepSpec {
  float theta_start = 0.0f; /**< Contact angle at frame 0, radians;
                               negatives clamp to 0, and it enters only as
                               the opening slerp fraction theta_start /
                               theta_end, floored at K_EPS. */
  float theta_end = 0.0f;   /**< Arrival contact angle, radians; floored at
                               THETA_EPS. */
  int sweep_frames = 0;     /**< Angle-sweep frames (N). */
};

/**
 * @brief Relax leg: every vertex slerps from its seed position to its
 * relaxed one.
 */
struct RelaxSpec {
  int iterations = 0; /**< Spring-relaxation passes of the arrival form. */
  const MeshOps::RelaxBake *bake =
      nullptr; /**< Shipped converged relaxation to land on; null instead
                  runs `iterations` live steps, and one of the two is
                  required. */
  int sweep_frames = 0; /**< Slerp frames (N). */
};

/** @brief Medial-slerp leg (the Conway-dual bridge's middle leg). */
struct MedialSpec {
  int sweep_frames = 0; /**< Slerp frames (N). */
};

/** @brief Reconcile leg closing a smooth kis/needle path. */
struct ReconcileSpec {
  const math::Vector *to_positions =
      nullptr;          /**< Authored vertex positions, one per seed vertex
                  (index-corresponded through the caller's bijection); slerp
                  endpoints. */
  size_t to_count = 0;  /**< Number of authored endpoint positions. */
  int sweep_frames = 0; /**< Slerp frames (N). */
};

/** @brief Gated-swap leg: a partition op with no sweep. */
struct GatedSwapSpec {
  SwapOp op = SwapOp::KIS; /**< Partition operator applied at the swap
                              frame. */
  int gate_frames = 0;     /**< Frames on each side of the swap (F_gate); the
                          leg runs 2 * gate_frames + 1 frames. */
};

/** @brief Face-order relationship between departed and swept meshes. */
enum class FaceCorrespondence : uint8_t {
  GEOMETRIC,    /**< Matches departed faces by centroid geometry. */
  IDENTITY,     /**< Carries palettes by unchanged face emission order. */
  DUAL_CLOSING, /**< Maps the dual bridge closing faces to the departed mesh. */
};

/**
 * @brief Palette provenance of the departed node.
 * @details prev_face_palette describes the node base mesh the leg departs
 * from, in emission order; read only during construction. Supplying
 * prev_face_centroid makes the mapping geometric.
 */
struct PaletteHandoff {
  const BakedPaletteBank *bank =
      nullptr; /**< The effect's baked source LUTs. */
  const uint8_t *prev_face_palette =
      nullptr;           /**< Per-face palette of the departed base mesh. */
  size_t prev_faces = 0; /**< Face count of the departed base mesh. */
  const math::Vector *prev_face_centroid =
      nullptr; /**< Unit centroid per departed base face; enables the
                  geometric provenance mapping. */
  FaceCorrespondence correspondence =
      FaceCorrespondence::GEOMETRIC;    /**< Palette provenance mapping. */
  bool preserve_without_births = false; /**< Hold colors without births. */
};

/**
 * @brief Bookend grouping of the arrival node.
 * @details topology[f] is the class the effect displays arrival face f
 * with at the closing bookend. Color targets key on it, so faces it merges
 * converge to one color by w = 1. A null topology keys targets on the swept
 * classification instead. Palette-preserving legs retain departed colors.
 */
struct BookendClasses {
  const uint16_t *topology = nullptr; /**< Bookend class per arrival face. */
  size_t faces = 0;                   /**< Arrival base-mesh face count. */
};

/** @brief Split of a swept leg's eased progress into its two stages. */
struct SweepPhase {
  float k;            /**< Operator-sweep fraction in [0, 1]. */
  float settle_alpha; /**< Relax-slerp fraction in [0, 1]. */
};

/** @brief Per-leg switches for init_conway. */
struct ConwayLegFlags {
  /** Whether the leg relax-slerps: at arrival on a forward leg, at
   * departure on a reverse leg. */
  bool settle = false;
  /** Whether the leg is the jitterbug bridge (vertex-orbit faces survive
   * into the node mesh). */
  bool jitterbug = false;
  /** Dual-bridge provenance: start centroids come from the leg's own start
   * mesh when closing, or arrival mesh when opening. */
  bool bridge_provenance = false;
};

} // namespace Animation::detail
