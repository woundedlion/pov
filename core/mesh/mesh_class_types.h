/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstddef>
#include <cstdint>

#include "engine/memory.h"
#include "render/sdf/face.h"

/**
 * @file mesh_class_types.h
 * @brief Records a congruence-class bake produces and the rasterizer reads.
 *
 * The runtime half of the congruence-class LUT feature: the class id space and
 * the three record types Scan::Mesh binds per frame. The clustering and LUT
 * bake that fill them live in mesh_classes.h, which the rasterizer does not
 * need — this header keeps the mesh-building machinery out of every rasterizer
 * translation unit.
 */
namespace MeshOps {

/** Sentinel class id: face keeps the per-face exact path. */
inline constexpr uint8_t NO_CLASS = 0xFF;
/** Congruence-class capacity per mesh (census max over the registry is 24;
 *  overflow degrades the excess faces to NO_CLASS, it never traps). */
inline constexpr int MAX_CONGRUENCE_CLASSES = 32;
static_assert(MAX_CONGRUENCE_CLASSES < NO_CLASS,
              "class ids must fit uint8_t without aliasing NO_CLASS");

/**
 * @brief Per-face congruence record, baked once per spawned mesh.
 */
struct FaceClassRec {
  uint8_t class_id; /**< Congruence class, or NO_CLASS. */
  uint8_t
      vert_offset; /**< Cyclic offset aligning mesh vertex order to canonical. */
  uint8_t reflected; /**< Non-zero for the mirror family. */
};

static_assert(SDF::FaceScratchBuffer::MAX_VERTS <= UINT8_MAX,
              "FaceClassRec::vert_offset must hold any face vertex index");

/**
 * @brief One congruence class: canonical shape + optional distance LUT.
 */
struct CongruenceClass {
  SDF::ClassLut
      lut; /**< Distance LUT; data == nullptr when the class has none. */
  const float *canon_xy =
      nullptr;           /**< Canonical centered 2D polygon, x/y pairs. */
  float canon_sq = 0.0f; /**< Sum of squared canonical vertex magnitudes. */
  int n_verts = 0;       /**< Vertex count. */
  int topo_id = -1;      /**< Seeding topology class. */
  uint16_t members = 0;  /**< Faces assigned to this class. */
  bool concave = false;  /**< Canonical shape is concave (LUT-eligible). */
};

/**
 * @brief Congruence-class bake for one spawned mesh (persistent-arena owned).
 * @details Lives in the effect's persistent arena and dies at every arena
 * compaction — rebake unconditionally after compact_keep_front or the fading
 * mesh silently degrades to NO_CLASS (correct but slower).
 */
struct MeshClassBake {
  ArenaVector<CongruenceClass> classes; /**< Congruence classes, dense ids. */
  ArenaVector<FaceClassRec>
      face_recs;                  /**< Per-face records, mesh face order. */
  float worst_residual_px = 0.0f; /**< Worst accepted RMS residual (pixels). */
  float predicted_hit_share =
      0.0f; /**< Safe-cell share on LUT-bound faces (quality; gate >= MIN_CLASS_HIT_SHARE). */
  float lut_face_share =
      0.0f; /**< Faces bound to a LUT / all faces (coverage). */
  uint16_t shared_faces = 0;  /**< Faces in classes with >= 2 members. */
  uint16_t concave_faces = 0; /**< Faces in concave (LUT-eligible) classes. */
  uint16_t lut_faces = 0;     /**< Faces whose class received a LUT. */
  uint16_t luts_built = 0;    /**< Classes that received a LUT. */
  uint16_t degraded_classes =
      0; /**< Eligible classes whose grid was shrunk to fit the budget (some are then discarded by the hit-share gate). */
  uint16_t dropped_classes =
      0; /**< Eligible classes denied a LUT once the budget could not fit even the min grid. */
  uint16_t dropped_faces =
      0; /**< Faces in dropped classes (kept the exact path). */
  uint16_t lowq_classes =
      0; /**< Classes whose built LUT missed MIN_CLASS_HIT_SHARE and was discarded. */
  uint16_t overflow_faces =
      0; /**< Faces left NO_CLASS because the class table was already full. */
  size_t lut_bytes =
      0; /**< Persistent bytes charged to kept LUTs (<= the byte budget). */
  size_t aux_bytes =
      0; /**< Persistent bytes claimed outside the LUT budget: the class and per-face record vectors plus every founded class's canonical polygon. */
};

} // namespace MeshOps
