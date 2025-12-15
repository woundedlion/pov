# Holosphere Firmware (POV)

**Holosphere** is a graphics engine built for Teensy 4.x microcontrollers to drive high-speed Spherical Persistence of Vision (POV) displays. Unlike traditional POV renderers that work on a 2D polar grid, Holosphere abstracts the display as a **virtual 3D unit sphere**, allowing for complex vector graphics, particle physics, and 3D transformations.

This project is the hardware counterpart to **[Daydream]([../daydream](https://github.com/woundedlion/daydream))**, a WebGL simulator used for prototyping these effects.

---

## 🚀 System Architecture & Component Guide

The engine is architected as a modular pipeline. Data flows from high-level **Animations** $\to$ **Drawing Primitives** $\to$ **Abstract Dots** $\to$ **Filter Chain** $\to$ **Physical LEDs**.

### 1. Core Mathematics (`3dmath.h`)
The foundation of the engine, providing optimized floating-point primitives for the Teensy's FPU.
* **`Vector`**: A custom 3D Cartesian vector struct $(x, y, z)$. It supports operator overloading for intuitive vector math (addition, cross products, normalization) required for 3D physics and geometry.
* **`Quaternion`**: Implements 4D rotation math ($w, x, y, z$). This is crucial for the engine's camera system, avoiding Gimbal lock and allowing for smooth spherical interpolation (SLERP) of rotations.
* **`Spherical`**: A helper struct for converting Cartesian vectors into Polar coordinates $(\theta, \phi)$ just before rasterization.
* **`MobiusParams`**: A structure holding complex number coefficients ($a, b, c, d$) used for conformal mapping transformations.

### 2. Scene Graph & Geometry (`geometry.h`)
Defines how objects exist and move in the virtual world.
* **`Orientation`**: **(Key Feature)** This class solves the "strobing" issue inherent in POV displays. Instead of storing a single rotation state, it maintains a **history buffer** of quaternions for the current frame. When an object moves, the engine "tweens" through this history to draw it as a continuous arc rather than a single point, generating high-quality temporal motion blur.
* **`Dot`**: The atomic unit of rendering. It bundles a `Vector` (position) with a `Color4` (RGBA color).
* **`Projection`**: Contains the logic to map normalized 3D vectors onto the 2D "unwrapped" map of the LED strip (`pixel_to_vector` and `vector_to_pixel`).

### 3. The Render Pipeline (`filter.h`)
The engine uses a recursive template-based pipeline to process graphics. This allows for zero-cost abstractions where filter chains are resolved at compile time.
* **`Pipeline<...>`**: A variadic template that chains filters together. Data passed into the pipeline traverses every filter before reaching the canvas.
* **`FilterOrient`**: The "Camera" of the engine. It intercepts every 3D point and rotates it by the global `Orientation` quaternion, allowing the entire world to spin or tumble.
* **`FilterDecay`**: Implements visual trails. Instead of drawing to the screen, it records points into a `DecayBuffer`. In subsequent frames, it replays these points with diminishing brightness based on their age.
* **`FilterMobius`**: Performs a complex-plane Möbius transformation $f(z) = \frac{az+b}{cz+d}$ on the 3D vectors, creating hyperbolic and elliptic geometry distortions.
* **`FilterAntiAlias`**: The "Rasterizer". It takes the final floating-point coordinates and distributes the pixel's energy to the four nearest physical LEDs using a **Quintic Kernel** (SmootherStep) for ultra-smooth sub-pixel rendering.

### 4. Animation Timeline (`animation.h`)
A scripting system that replaces complex state machines with a linear, scheduled timeline.
* **`Timeline`**: The master scheduler. You add events to it with a specific start time.
* **`Sprite`**: Manages the lifecycle of a visual object. It handles the `draw()` callback, and automatically manages fade-in and fade-out opacity curves over the object's lifespan.
* **`Transition`**: Smoothly interpolates a scalar value (like a radius or color index) from $A$ to $B$ over time using a specific easing function.
* **`Mutation`**: Similar to Transition, but drives a value using a continuous function (like a Sine wave) rather than a start/end point.
* **`Motion`**: Moves an `Orientation` object along a defined 3D `Path`.

### 5. Drawing Primitives (`draw.h`)
Helper functions that generate `Dots` to be fed into the pipeline.
* **`draw_line`**: Draws a Great Circle arc between two vectors on the sphere surface. It uses adaptive sampling to ensure the line resolution matches the screen resolution.
* **`draw_ring`**: Generates a circle on the sphere surface defined by a normal vector and a radius.
* **`draw_polyhedron`**: Renders wireframe meshes (like Dodecahedrons) defined in `geometry.h` by iterating their edge lists.

### 6. Color Theory (`color.h`)
* **`GenerativePalette`**: A procedural color system. Instead of hardcoded colors, it generates harmonious palettes on the fly using rules (e.g., "Split-Complementary", "Triadic") and shape distributions (e.g., "Vignette", "Circular").
* **`ProceduralPalette`**: Defines gradients using cosine-wave coefficients, allowing for smooth, mathematical color evolution without lookup tables.

### 7. Hardware Abstraction (`led.h`)
* **`POVDisplay`**: Manages the `FastLED` integration and the high-speed hardware timer (`IntervalTimer`). It synchronizes the column-display rate to the motor's RPM.
* **`Effect` / `Canvas`**: Implements **Double Buffering**.
    * Buffer A is being written to by the `draw()` loop.
    * Buffer B is being read by the high-priority Interrupt Service Routine (ISR) to light the LEDs.
    * The `Canvas` class uses RAII to automatically swap and lock these buffers at the start/end of a frame.

---

## 📦 Included Effects

The `effects/` directory contains specific visual sketches demonstrating the engine's capabilities:

* **`GSReactionDiffusion` / `BZReactionDiffusion`**: Simulates biological pattern formation (Gray-Scott and Belousov-Zhabotinsky models) on a spherical Fibonacci lattice.
* **`FlowField`**: Uses 4D OpenSimplex noise to drive thousands of particles across the sphere surface, creating fluid-like motion.
* **`MobiusGrid`**: Visualizes conformal mappings by projecting a grid through a dynamic Möbius transformation.
* **`Thrusters`**: A physics simulation of particles emitted from a moving source, affected by drag and momentum.
* **`Comets`**: Agents that follow complex Lissajous curves on the sphere surface, leaving long decay trails.

## 🔧 Setup & Configuration

1.  **Hardware**: Requires a **Teensy 4.0 or 4.1**.
2.  **Libraries**: Install **FastLED**.
3.  **Config (`led.h`)**:
    ```cpp
    static constexpr int NUM_PIXELS = 40;   // Number of physical LEDs
    static constexpr unsigned int RPM = 480; // Rotation speed of the motor
    ```

---
*See [Daydream](../daydream/README.md) for the WebGL simulator used to develop these effects.*
