#pragma once

#include "types.h"

#define ASSERT(exp) \
    if (exp) { } \
    else \
    { \
        assert_handler(#exp, __FILE__, __LINE__); \
        (__debugbreak(), 1); \
    }

#ifndef FORCEINLINE
#if defined(_MSC_VER) // MSVC
    #define FORCEINLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) // GCC/Clang
    #define FORCEINLINE inline __attribute__((always_inline))
#else
    #define FORCEINLINE inline
#endif
#endif
