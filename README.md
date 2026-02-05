# Holosphere Firmware (POV)

**Holosphere** is a high-performance graphics engine for Teensy 4.x microcontrollers, designed to drive Spherical Persistence of Vision (POV) displays. While traditional POV displays often rely on simple 2D polar mappings, Holosphere treats the display as a **3D virtual unit sphere**, enabling complex vector graphics, volumetric solids, and relativistic distortions (Möbius transformations) with 60 FPS fluidity.

This project is the firmware implementation of the **[Daydream](../daydream)** Digital Twin, maintaining strict nomenclature and math parity with its WebGL counterpart.

---

## 📜 License

This project operates under a **Dual License** model to protect the creative works (Effects) while keeping the engine open for non-commercial educational use.

### 1. The Core Engine
The firmware core, including all mathematics, rendering pipelines, and driver logic, is licensed under the **[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0)**.
*   **Allowed**: Personal use, education, hobbyist projects.
*   **Prohibited**: Commercial distribution or selling hardware running this engine.
*   **Files**: `3dmath.h`, `geometry.h`, `led.h`, `plot.h`, `scan.h`, `filter.h`, `color.h`, `animation.h`, `effects_engine.h`.

### 2. The Effects Library
The specific visual sketches found in the `effects/` directory are **Proprietary** and **All Rights Reserved**.
*   **Restricted**: These files may not be copied, modified, or redistributed without explicit written permission from the author.
*   **Files**: All contents of the `effects/` directory.

---
*See [Daydream](../daydream/README.md) for the WebGL simulator used to develop these effects.*

---

## 🛰️ Graphics Engine Architecture

The Holosphere engine is built on four architectural pillars:

1.  **3D Spherical Coordinate Space**: All geometry, physics, and lighting calculations are performed on a unit sphere. This eliminates the "seam" and "pole" distortion issues common in 2D POV implementations.
2.  **Orientation History (Temporal Smoothing)**: To overcome the low persistence of mechanical POV displays, the engine maintains a high-resolution history buffer of rotations for every moving object. By "tweening" through this history within a single frame, the engine renders continuous **Great Circle arcs** instead of discrete points, resulting in hardware-accelerated motion blur.
3.  **Adaptive Rasterization**: The engine dynamically calculates the required sampling density for every line or shape based on the current display resolution (`MAX_W` x `H`), ensuring that graphics are neither "dashed" (under-sampled) nor "aliased" (over-sampled).
4.  **The Register Bank (`v0-v3`)**: Shaders in Holosphere communicate through a standardized fragment register bank. This ensures that effects developed in the JavaScript Digital Twin can be ported to C++ by simply matching the register-to-color mapping logic.

---

## ⚙️ Rendering Pipeline

Data flows through the engine in four distinct stages:

### 1. The Kinetic Step (`draw_frame()`)
The effect class updates its internal state. This involves stepping the `Timeline`, resolving physics (velocities/positions), and pushing new rotation quaternions into `Orientation` history buffers.

### 2. Geometry Emission
An effect calls a drawing primitive from either the `Plot` or `Scan` namespaces:
*   **Vector (`Plot`)**: Primitives like `Line` or `Ring` sample their paths into a stream of **Fragments**.
*   **Volumetric (`Scan`)**: The engine iterates through scanline intervals, using **Signed Distance Fields (SDF)** to determine the coverage of solids like `Polygon` or `HarmonicBlob`.

### 3. The Filter Pipeline
Atomic Fragments enter a template-resolved `Pipeline`. Filters are applied in a strict, zero-allocation sequence:
*   **`FilterMobius`**: Warps the 3D position using complex-plane transformations.
*   **`FilterOrient`**: Projects the world through the object's `Orientation` history, generating motion blur.
*   **`TemporalFilter`**: Smoothly accumulates energy over multiple frames to eliminate mechanical flicker.

### 4. Hardware Rasterization
Processed colors are plotted into a `Canvas`. The `POVDisplay` class uses a hardware `IntervalTimer` to pulse these columns to the physical LED strip via the `FastLED` driver, synchronized with the motor's RPM.

---

## 📁 File Manifest

### Core Engine & Mathematics
*   **`3dmath.h`**: The foundational algebra library. Implements high-performance `Vector` (Cartesian), `Quaternion` (Rotation), and `Complex` (Mobius) types specialized for the Teensy FPU.
*   **`geometry.h`**: Defines the engine's atomic units: the `Fragment` (position + registers) and the `Orientation` history buffer. It also provides `MeshState` for handling 3D geometry data.
*   **`led.h` & `rotate.h`**: Hardware abstraction layer. Defines physical LED pins, resolution constants (`MAX_W`, `H`), and provides inverse stereographic projection and basis-creation helpers.
*   **`spatial.h`**: Spatial indexing and acceleration structures (`BVH`, `KDTree`, `SpatialHash`) for optimizing mesh intersections and neighbor searches.
*   **`static_circular_buffer.h`**: A zero-heap, fixed-capacity container used for orientation history and fragment pools to prevent memory fragmentation in long-running firmware.
*   **`util.h`**: Generic utility functions for random number generation and modulo-wrapping arithmetic.
*   **`FastNoiseLite.h`**: A lightweight, header-only noise generation library used for procedural textures and movement (Simplex/Perlin noise).

### Rendering & Shading
*   **`plot.h`**: The vector graphics engine. Contains implementations for Geodesic Lines, Rings, Stars, Flowers, and the high-fidelity `ParticleSystem` (v0: trail progress, v1: ID).
*   **`scan.h`**: The volumetric rasterizer. Defines the `SDF` namespace (signed distance math for solids) and the scanline `rasterize` engine used for filled shapes and meshes.
*   **`filter.h`**: The post-processing suite. Contains the `Pipeline` class and all fragment-level filters (Mobius, Orientation, Temporal Interpolation).
*   **`color.h`**: The core color type definitions (`Color4`) and blending logic.
*   **`palettes.h`**: Defines `ProceduralPalette`, `GenerativePalette`, and decorators like `FalloffPalette` to manage color gradients and transitions.
*   **`solids.h`**: A library of pre-defined 3D mesh data (Platonic solids like Icosahedrons, Dodecahedrons) used for testing and rigid body effects.

### Animation & Logic
*   **`animation.h`**: The kinetic engine. Contains the `Timeline`, `Sprite` lifecycles, and physics-based movement classes (`Mutation`, `Transition`, `Motion`, `ParticleSystem`).
*   **`easing.h`**: A collection of standard easing functions (e.g., `ease_in_out_cubic`, `ease_mid`) used for smooth animation transitions.
*   **`effects_engine.h`**: The interface bridge. Defines engine-wide concepts, template aliases, and the standardized shader signatures required for parity with the Digital Twin.

### Application Entry & Management
*   **`Holosphere.ino`**: The main Arduino/Teensy entry point. Initializes the hardware, sets up the `POVDisplay`, and manages effect switching.
*   **`effects.h`**: A registry file that includes all active effect headers, making them available to the main application.
*   **`effects_legacy.h`**: Contains legacy or deprecated effects that are preserved for reference or regression testing.

---

## 📦 Included Effects

*   **`MindSplatter`**: Complex particle simulation with relativistic Mobius warping and alpha-masking.
*   **`DreamBalls`**: A demonstration of the volumetric SDF engine rendering moving metaballs.
*   **`MobiusGrid`**: A visualization of conformal mappings with dynamic axis labeling and antipode stabilization.
*   **`FlowField`**: 4D OpenSimplex noise driving thousands of particles across the sphere surface.

---

## 🔧 Setup

1.  **Hardware**: Requires Teensy 4.0/4.1.
2.  **Toolchain**: PlatformIO or Arduino IDE with Teensyduino.
3.  **Dependencies**: FastLED.
4.  **Calibration**: Adjust `RPM` and `NUM_PIXELS` in `led.h` to match your physical hardware.
