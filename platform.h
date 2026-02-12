/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifdef ARDUINO
    #include <Arduino.h>
    #include <FastLED.h>
    
    // Alias Pixel to CRGB for consistency with non-Arduino builds if needed, 
    // though the codebase seems to use CRGB directly often. 
    // We will ensure Pixel type alias is consistent.
    using Pixel = CRGB;

    namespace hs {
       inline void log(const char* msg) { Serial.println(msg); }
       inline unsigned long millis() { return ::millis(); }
       inline void disable_interrupts() { noInterrupts(); }
       inline void enable_interrupts() { interrupts(); }
    }

#else
    // Non-Arduino / PC Simulation Platform
    #include <cstdint>
    #include <cmath>
    #include <algorithm>
    #include <iostream>
    #include <vector>
    #include <cstring>
    #include <chrono>

    // --- Mock Arduino Constants/Types ---
    using  byte = uint8_t;
    using  boolean = bool;

    // --- Mock FastLED Types ---
    struct CRGB {
        uint8_t r, g, b;
        constexpr CRGB() : r(0), g(0), b(0) {}
        constexpr CRGB(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
        constexpr CRGB(uint8_t gray) : r(gray), g(gray), b(gray) {} // explicit?
        
        // predefined colors
        static const CRGB Black;
        static const CRGB White;
        static const CRGB Red;
        static const CRGB Green;
        static const CRGB Blue;

        // Operators
        bool operator==(const CRGB& rhs) const { return r == rhs.r && g == rhs.g && b == rhs.b; }
        bool operator!=(const CRGB& rhs) const { return !(*this == rhs); }
        
        CRGB& operator+=(const CRGB& rhs) {
            // Saturated add
            r = (r + rhs.r > 255) ? 255 : r + rhs.r;
            g = (g + rhs.g > 255) ? 255 : g + rhs.g;
            b = (b + rhs.b > 255) ? 255 : b + rhs.b;
            return *this;
        }

        // Methods matching FastLED
        CRGB& lerp16(const CRGB& other, uint16_t frac) {
            // frac is 0..65535
            r = static_cast<uint8_t>((static_cast<uint32_t>(r) * (65535 - frac) + static_cast<uint32_t>(other.r) * frac) >> 16);
            g = static_cast<uint8_t>((static_cast<uint32_t>(g) * (65535 - frac) + static_cast<uint32_t>(other.g) * frac) >> 16);
            b = static_cast<uint8_t>((static_cast<uint32_t>(b) * (65535 - frac) + static_cast<uint32_t>(other.b) * frac) >> 16);
            return *this;
        }
    };
    
    inline const CRGB CRGB::Black = CRGB(0,0,0);
    inline const CRGB CRGB::White = CRGB(255,255,255);
    inline const CRGB CRGB::Red = CRGB(255,0,0);
    inline const CRGB CRGB::Green = CRGB(0,255,0);
    inline const CRGB CRGB::Blue = CRGB(0,0,255);

    using Pixel = CRGB;

    struct CHSV {
        uint8_t h, s, v;
        constexpr CHSV() : h(0), s(0), v(0) {}
        constexpr CHSV(uint8_t h, uint8_t s, uint8_t v) : h(h), s(s), v(v) {}
    };
    
    // --- Mock FastLED Functions ---
    // qadd8: saturated addition for 8-bit integers
    inline uint8_t qadd8(uint8_t i, uint8_t j) {
        int t = i + j;
        if (t > 255) t = 255;
        return t;
    }
    
    // qsub8: saturated subtraction
    inline uint8_t qsub8(uint8_t i, uint8_t j) {
        int t = i - j;
        if (t < 0) t = 0;
        return t;
    }

    // --- Mock Arduino Functions ---
    inline int random(int max) { return rand() % max; }
    inline int random(int min, int max) { return min + (rand() % (max - min)); }
    inline long map(long x, long in_min, long in_max, long out_min, long out_max) {
      return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    }

    // --- System Mock ---
    struct SerialMock {
        void println(const char* msg) { std::cout << msg << std::endl; }
        void println(int val) { std::cout << val << std::endl; }
        void printf(const char* fmt, ...) {
             // simplified
             std::cout << fmt << std::endl; 
        }
    };
    static SerialMock Serial;

    // --- FastLED Mocks ---
    inline uint8_t random8() { return rand() % 256; }
    inline uint8_t random8(uint8_t top) { return rand() % top; }
    inline uint16_t random16() { return rand() % 65536; }
    inline void random16_add_entropy(uint16_t) {}

    inline uint8_t beatsin8(uint16_t bpm, uint8_t low, uint8_t high, uint16_t phase = 0, uint16_t offset = 0) {
        // Simplified implementation using system clock
        // Real implementation would need a unified timebase
        return low; 
    }
    
    inline uint16_t beatsin16(uint16_t bpm, uint16_t low, uint16_t high, uint16_t phase = 0, uint16_t offset = 0) {
        return low;
    }
    
    inline uint8_t addmod8(uint8_t a, uint8_t b, uint8_t m) {
        return (a + b) % m;
    }
    
    inline uint8_t map8(uint8_t x, uint8_t in_min, uint8_t in_max) {
        return (x - in_min) * 255 / (in_max - in_min); // rough approximation
    }
    
    inline uint8_t triwave8(uint8_t in) {
        if (in & 0x80) {
            in = 255 - in;
        }
        return in << 1;
    }

    #define FASTRUN

    // Mock EVERY_N_MILLIS using a simple static checker
    #define EVERY_N_MILLIS(N) \
      static unsigned long __last_##__LINE__ = 0; \
      unsigned long __now_##__LINE__ = hs::millis(); \
      if (__now_##__LINE__ - __last_##__LINE__ >= (N) && (__last_##__LINE__ = __now_##__LINE__))

    #define EVERY_N_SECONDS(N) EVERY_N_MILLIS((N) * 1000)
    #define EVERY_N_MILLISECONDS(N) EVERY_N_MILLIS(N)

    namespace hs {
       inline unsigned long millis() { 
           using namespace std::chrono;
           return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
       }
       inline void disable_interrupts() {}
       inline void enable_interrupts() {}
    }
    
    // Global millis/micros if needed, though prefer namespaced
    inline unsigned long millis() { return hs::millis(); }
    inline unsigned long micros() { return hs::millis() * 1000; } // approx

#endif
