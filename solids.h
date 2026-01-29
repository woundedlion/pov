/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "geometry.h"

// 1. TETRAHEDRON (4 Verts, 4 Faces)
struct Tetrahedron {
  std::vector<Vector> vertices = {
    Vector(0.57735, 0.57735, 0.57735),
    Vector(0.57735, -0.57735, -0.57735),
    Vector(-0.57735, 0.57735, -0.57735),
    Vector(-0.57735, -0.57735, 0.57735)
  };

  std::vector<std::vector<int>> faces = {
    {0, 3, 1}, {0, 2, 3}, {0, 1, 2}, {1, 3, 2}
  };

  Tetrahedron() {
    for (auto& v : vertices) v.normalize();
  }
};

// 2. CUBE (8 Verts, 6 Faces)
struct Cube {
  const VertexList vertices = {
    Vector(1, 1, 1),   // 0
    Vector(1, 1, -1),  // 1
    Vector(1, -1, 1),  // 2
    Vector(1, -1, -1), // 3
    Vector(-1, 1, 1),  // 4
    Vector(-1, 1, -1), // 5
    Vector(-1, -1, 1), // 6
    Vector(-1, -1, -1) // 7
  };

  const std::vector<std::vector<int>> faces = {
    {0, 3, 2, 1}, // Bottom
    {0, 1, 5, 4}, // Front
    {0, 4, 7, 3}, // Left
    {6, 5, 1, 2}, // Right
    {6, 2, 3, 7}, // Back
    {6, 7, 4, 5}  // Top
  };

  const AdjacencyList eulerPath = {
    {1, 2, 4}, // 0
    {3, 5},    // 1
    {3, 6},    // 2
    {7},       // 3
    {5, 6},    // 4
    {7},       // 5
    {7},       // 6
    {}         // 7
  };
};

// 3. OCTAHEDRON (6 Verts, 8 Faces)
struct Octahedron {
  std::vector<Vector> vertices = {
    Vector(1, 0, 0),
    Vector(-1, 0, 0),
    Vector(0, 1, 0),
    Vector(0, -1, 0),
    Vector(0, 0, 1),
    Vector(0, 0, -1)
  };

  std::vector<std::vector<int>> faces = {
    {4, 0, 2}, {4, 2, 1}, {4, 1, 3}, {4, 3, 0},
    {5, 2, 0}, {5, 1, 2}, {5, 3, 1}, {5, 0, 3}
  };

  Octahedron() {
      for (auto& v : vertices) v.normalize();
  }
};

// 4. ICOSAHEDRON (12 Verts, 20 Faces)
struct Icosahedron {
  std::vector<Vector> vertices;
  std::vector<std::vector<int>> faces = {
    {0, 1, 4}, {0, 4, 9}, {9, 4, 5}, {4, 8, 5}, {4, 1, 8},
    {8, 1, 10}, {8, 10, 3}, {5, 8, 3}, {5, 3, 2}, {2, 3, 7},
    {7, 3, 10}, {7, 10, 6}, {7, 6, 11}, {11, 6, 0}, {0, 6, 1},
    {6, 10, 1}, {9, 11, 0}, {9, 2, 11}, {9, 5, 2}, {7, 11, 2}
  };

  Icosahedron() {
      const float X = 0.525731112119f;
      const float Z = 0.850650808352f;
      vertices = {
        Vector(-X, 0.0, Z), Vector(X, 0.0, Z), Vector(-X, 0.0, -Z), Vector(X, 0.0, -Z),
        Vector(0.0, Z, X), Vector(0.0, Z, -X), Vector(0.0, -Z, X), Vector(0.0, -Z, -X),
        Vector(Z, X, 0.0), Vector(-Z, X, 0.0), Vector(Z, -X, 0.0), Vector(-Z, -X, 0.0)
      };
      for (auto& v : vertices) v.normalize();
  }
};

// 5. DODECAHEDRON (20 Verts, 12 Faces)
struct Dodecahedron {
  VertexList vertices = {
    {1, 1, 1},       // 0
    {1, -1, 1},      // 1
    {1, 1, -1},      // 2
    {1, -1, -1},     // 3
    {-1, 1, 1},      // 4
    {-1, -1, 1},     // 5
    {-1, 1, -1},     // 6
    {-1, -1, -1},    // 7

    {0, PHI, 1 / PHI},   //8
    {0, -PHI, 1 / PHI},  // 9
    {0, PHI, -1 / PHI},  // 10
    {0, -PHI, -1 / PHI}, // 11

    {1 / PHI, 0, PHI},   // 12
    {1 / PHI, 0, -PHI},  // 13
    {-1 / PHI, 0, PHI},  // 14
    {-1 / PHI, 0, -PHI}, // 15

    {PHI, 1 / PHI, 0},   // 16
    {PHI, -1 / PHI, 0},  // 17
    {-PHI, 1 / PHI, 0},  // 18
    {-PHI, -1 / PHI, 0}, // 19
  };

  AdjacencyList edges = {
    {8, 12, 16},  // 0
    {9, 12, 17},  // 1
    {10, 13, 16}, // 2
    {11, 13, 17}, // 3
    {8, 14, 18},  // 4
    {9, 14, 19},  // 5
    {10, 15, 18}, // 6
    {11, 15, 19}, // 7
    {0, 4, 10},   // 8
    {1, 5, 11},   // 9
    {2, 6, 8},    // 10
    {3, 7, 9},    // 11
    {0, 1, 14},   // 12
    {2, 3, 15},   // 13
    {4, 5, 12},   // 14
    {6, 7, 13},   // 15
    {0, 2, 17},   // 16
    {1, 3, 16},   // 17
    {4, 6, 19},   // 18
    {5, 7, 18},   // 19
  };

  std::vector<std::vector<int>> faces = {
        {0, 8, 9, 4, 16},
        {0, 12, 13, 1, 8},
        {0, 16, 17, 2, 12},
        {8, 1, 18, 5, 9},
        {12, 2, 10, 3, 13},
        {16, 4, 14, 6, 17},
        {9, 5, 15, 14, 4},
        {6, 11, 10, 2, 17},
        {3, 19, 18, 1, 13},
        {7, 15, 5, 18, 19},
        {7, 11, 6, 14, 15},
        {7, 19, 3, 10, 11}
  };

  AdjacencyList euler_path = {
    {8, 12, 16},  // 0
    {9, 12, 17},  // 1
    {10, 13, 16}, // 2
    {11, 13, 17}, // 3
    {8, 14, 18},  // 4
    {9, 14, 19},  // 5
    {10, 15, 18}, // 6
    {11, 15, 19}, // 7
    {10},   // 8
    {11},   // 9
    {8},    // 10
    {9},    // 11
    {14},   // 12
    {15},   // 13
    {12},   // 14
    {13},   // 15
    {17},   // 16
    {16},   // 17
    {19},   // 18
    {18},   // 19
  };

  /**
   * @brief Constructor that normalizes all vertices to the unit sphere.
   */
  Dodecahedron() {
    for (auto& v : vertices) {
      v.normalize();
    }
  }

  /**
   * @brief Prints the polyhedron's data to the serial output for debugging.
   */
  void dump() const {
    Serial.printf("Dodecahedron\n");
    Serial.printf("vertices [%d] = \n", vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
      Serial.printf("%d: [%f, %f, %f]\n", i, vertices[i].i, vertices[i].j, vertices[i].k);
    }
    Serial.printf("euler path [%d] = \n", euler_path.size());
    for (size_t i = 0; i < euler_path.size(); ++i) {
      Serial.printf("%d: [", i);
      for (size_t j = 0; j < euler_path[i].size(); ++j) {
        Serial.printf("%d, ", euler_path[i][j]);
      }
      Serial.printf("]\n", i);
    }
  }
};
