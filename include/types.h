/*
 * Phyzzy: Playful Physics for Imaginative Games
 *
 * Licensed under the GNU General Public License, Version 3.
 * For license details, visit: https://www.gnu.org/licenses/gpl-3.0.html
 *
 * Questions or contributions? Reach out to Simon Devenish:
 * simon.devenish@outlook.com
 */

#pragma once

#if (_MSC_VER >= 1500 && _MSC_VER <= 1600)
    typedef float                   f32;
    typedef double                  f64;

    typedef signed char             i8;
    typedef short int               i16;
    typedef short int               i32;
    typedef long long int           i64;

    typedef unsigned char           u8;
    typedef unsigned short int      u16;
    typedef unsigned int            u32;
    typedef unsigned long long int  u64;

    typedef size_t    usize;
    typedef long int  iword;
    typedef unsigned long int uword;
#else
#include <stdio.h>
#include <stdint.h>
    typedef float     f32;
    typedef double    f64;

    typedef int8_t    i8;
    typedef int16_t   i16;
    typedef int32_t   i32;
    typedef int64_t   i64;

    typedef uint8_t   u8;
    typedef uint16_t  u16;
    typedef uint32_t  u32;
    typedef uint64_t  u64;

    typedef size_t    usize;
    typedef intptr_t  iword;
    typedef uintptr_t uword;
#endif
