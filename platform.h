/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifdef ARDUINO
    #include <Arduino.h>
    #include <FastLED.h>

    // using Pixel = CRGB; // Moved to color.h for Linear 16-bit support

    namespace hs {
       inline void log(const char* fmt, ...) {
            char buf[128];
            va_list args;
            va_start(args, fmt);
            vsnprintf(buf, sizeof(buf), fmt, args);
            va_end(args);
            Serial.println(buf);
       }
       inline unsigned long millis() { return ::millis(); }
       inline void disable_interrupts() { noInterrupts(); }
       inline void enable_interrupts() { interrupts(); }

       /**
       * @brief Generates a pseudo-random floating-point number between 0.0 and 1.0.
       * @return A random float in the range [0.0, 1.0].
       */
       float rand_f() {
           return static_cast<float>(::random(0, std::numeric_limits<int32_t>::max()))
           / std::numeric_limits<int32_t>::max();
       }

       float rand_f(float min, float max) {
           return min + rand_f() * (max - min);
       }

       /**
       * @brief Generates a pseudo-random integer within a specified range.
       * @param min The minimum value (inclusive).
       * @param max The maximum value (exclusive).
       * @return A random integer in the range [min, max).
       */
       int rand_int(int min, int max) {
           return ::random(min, max);
       }
       
       static constexpr int H_OFFSET = 3;
    }

#else 

    // Non-Arduino / PC Simulation Platform
    #ifndef DMAMEM
    #define DMAMEM
    #endif

    #ifdef __EMSCRIPTEN__
        #include <emscripten.h>
        #include <emscripten/bind.h>
    #endif

    #include <cstdint>
    #include <cmath>
    #include <algorithm>
    #include <iostream>
    #include <vector>
    #include <cstring>
    #include <chrono>
    #include <iostream>

    template <typename T, typename U>
    constexpr T lerp(T a, T b, U t) {
        return (T)(a + (b - a) * t);
    }

    struct CHSV {
        uint8_t h, s, v;
        constexpr CHSV() : h(0), s(0), v(0) {}
        constexpr CHSV(uint8_t h, uint8_t s, uint8_t v) : h(h), s(s), v(v) {}
    };

    // --- Mock FastLED Types ---
    struct CRGB {
        uint8_t r, g, b;
        constexpr CRGB() : r(0), g(0), b(0) {}
        constexpr CRGB(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
        constexpr CRGB(uint8_t gray) : r(gray), g(gray), b(gray) {}
        
        // Convert HSV to RGB (Basic implementation)
        constexpr CRGB(const CHSV& hsv) {
            unsigned char region, remainder, p, q, t;

            if (hsv.s == 0) {
                r = hsv.v;
                g = hsv.v;
                b = hsv.v;
                return;
            }

            region = hsv.h / 43;
            remainder = (hsv.h - (region * 43)) * 6; 

            p = (hsv.v * (255 - hsv.s)) >> 8;
            q = (hsv.v * (255 - ((hsv.s * remainder) >> 8))) >> 8;
            t = (hsv.v * (255 - ((hsv.s * (255 - remainder)) >> 8))) >> 8;

            switch (region) {
                case 0: r = hsv.v; g = t; b = p; break;
                case 1: r = q; g = hsv.v; b = p; break;
                case 2: r = p; g = hsv.v; b = t; break;
                case 3: r = p; g = q; b = hsv.v; break;
                case 4: r = t; g = p; b = hsv.v; break;
                default: r = hsv.v; g = p; b = q; break;
            }
        }

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
        CRGB lerp16(const CRGB& other, uint16_t frac) const {
            CRGB ret;
            // frac is 0..65535
            ret.r = static_cast<uint8_t>((static_cast<uint32_t>(r) * (65535 - frac) + static_cast<uint32_t>(other.r) * frac) >> 16);
            ret.g = static_cast<uint8_t>((static_cast<uint32_t>(g) * (65535 - frac) + static_cast<uint32_t>(other.g) * frac) >> 16);
            ret.b = static_cast<uint8_t>((static_cast<uint32_t>(b) * (65535 - frac) + static_cast<uint32_t>(other.b) * frac) >> 16);
            return ret;
        }
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

    #define PI 3.1415926535897932384626433832795

    namespace hs {
       inline void log(const char* fmt, ...) {
           va_list args;
           va_start(args, fmt);
           vprintf(fmt, args);
           va_end(args);
           printf("\n");
       }
    }

    // --- Mock Arduino Constants/Types ---
    using  byte = uint8_t;
    using  boolean = bool;

    // --- Mock FastLED Constants ---
    enum FastLEDCheck {
        UncorrectedColor,
        TypicalLEDStrip,
        UncorrectedTemperature,
        Candle
    };

    struct FastLEDMock {
        void setCorrection(int) {}
        void setTemperature(int) {}
        template<typename T, int P1, int P2, int P3, int P4>
        void addLeds(CRGB* data, int nLeds) {}
        void show() {}
        void showColor(const CRGB&) {}
    };
    static FastLEDMock FastLED;

    // Helper for addLeds template args
    enum LEDType { WS2801 };
    enum ColorOrder { RGB };
    #define DATA_RATE_MHZ(x) x

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

       /**
       * @brief Generates a pseudo-random floating-point number between 0.0 and 1.0.
       * @return A random float in the range [0.0, 1.0].
       */
       float rand_f() {
           return static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
       }

       float rand_f(float min, float max) {
           return min + rand_f() * (max - min);
       }

       /**
       * @brief Generates a pseudo-random integer within a specified range.
       * @param min The minimum value (inclusive).
       * @param max The maximum value (exclusive).
       * @return A random integer in the range [min, max).
       */
       int rand_int(int min, int max) {
           return std::rand() % (max - min) + min;
       }
       
       static constexpr int H_OFFSET = 0;
    }
    
    // Global millis/micros if needed, though prefer namespaced
    inline unsigned long millis() { return hs::millis(); }
    inline unsigned long micros() { return hs::millis() * 1000; } // approx

#endif
