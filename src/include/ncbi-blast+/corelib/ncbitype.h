#ifndef CORELIB___NCBITYPE__H
#define CORELIB___NCBITYPE__H

/*  $Id: ncbitype.h 375760 2012-09-24 14:15:47Z ucko $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Denis Vakatov
 *
 *
 */

/**
 * @file ncbitype.h
 *
 * Defines NCBI C/C++ fixed-size types.
 *
 *   -  Char, Uchar
 *   -  Int1, Uint1
 *   -  Int2, Uint2
 *   -  Int4, Uint4
 *   -  Int8, Uint8
 *   -  Ncbi_BigScalar
 *   -  Macros for constant values definition.
 *
 */

/** @addtogroup Portability
 *
 * @{
 */

#include <../include/ncbi-blast+/ncbiconf.h>

#ifdef HAVE_INTTYPES_H
#  ifndef   __STDC_CONSTANT_MACROS
#    define __STDC_CONSTANT_MACROS
#  endif /*!__STDC_CONSTANT_MACROS*/
#  ifndef   __STDC_FORMAT_MACROS
#    define __STDC_FORMAT_MACROS
#  endif /*!__STDC_FORMAT_MACROS*/
#  include <inttypes.h>
#endif


/* Char, Uchar, Int[1,2,4], Uint[1,2,4]
 */

#if (SIZEOF_CHAR != 1)
#  error "Unsupported size of char(must be 1 byte)"
#endif
#if (SIZEOF_SHORT != 2)
#  error "Unsupported size of short int(must be 2 bytes)"
#endif
#if (SIZEOF_INT != 4)
#  error "Unsupported size of int(must be 4 bytes)"
#endif

#ifdef Int1
#  undef Int1
#  undef Uint1
#  undef Int2
#  undef Uint2
#  undef Int4
#  undef Uint4
#endif
#ifdef Int8
#  undef Int8
#  undef Uint8
#endif

typedef          char  Char;    /**< Alias for char */
typedef signed   char  Schar;   /**< Alias for signed char */
typedef unsigned char  Uchar;   /**< Alias for unsigned char */

#ifdef HAVE_INTTYPES_H
typedef int8_t   Int1;  /**< 1-byte  (8-bit)   signed integer */
typedef uint8_t  Uint1; /**< 1-byte  (8-bit) unsigned integer */
typedef int16_t  Int2;  /**< 2-byte (16-bit)   signed integer */
typedef uint16_t Uint2; /**< 2-byte (16-bit) unsigned integer */
typedef int32_t  Int4;  /**< 4-byte (32-bit)   signed integer */
typedef uint32_t Uint4; /**< 4-byte (32-bit) unsigned integer */
typedef int64_t  Int8;  /**< 8-byte (64-bit)   signed integer */
typedef uint64_t Uint8; /**< 8-byte (64-bit) unsigned integer */
/* We still need to know (u)int64_t's ultimate type when defining functions
 * or operators with one variant per distinct integral type. :-/ */
#  if SIZEOF_LONG == 8  &&  !defined(NCBI_OS_DARWIN)
#    define NCBI_INT8_IS_LONG      1
#  elif SIZEOF_LONG_LONG == 8
#    define NCBI_INT8_IS_LONG_LONG 1
#  elif SIZEOF___INT64 == 8
#    define NCBI_INT8_IS_INT64     1
#  endif
#else
typedef signed   char  Int1;    /**< Alias for signed char */
typedef unsigned char  Uint1;   /**< Alias for unsigned char */
typedef signed   short Int2;    /**< Alias for signed short */
typedef unsigned short Uint2;   /**< Alias for unsigned short */
typedef signed   int   Int4;    /**< Alias for signed int */
typedef unsigned int   Uint4;   /**< Alias for unsigned int */


/* Int8, Uint8
 */

#  if   (SIZEOF_LONG == 8)
#    define NCBI_INT8_TYPE         long
#    define NCBI_INT8_IS_LONG      1
#  elif (SIZEOF_LONG_LONG == 8)
#    define NCBI_INT8_TYPE         long long
#    define NCBI_INT8_IS_LONG_LONG 1
#  elif (SIZEOF___INT64 == 8)
#    define NCBI_INT8_TYPE         __int64
#    define NCBI_INT8_IS_INT64     1
/** @deprecated  Use NCBI_INT8_IS_INT64 instead */
#    define NCBI_USE_INT64         1
#  else
#    error "This platform does not support 8-byte integer"
#  endif

/** Signed 8 byte sized integer */
typedef signed   NCBI_INT8_TYPE Int8;    

/** Unsigned 8 byte sized integer */
typedef unsigned NCBI_INT8_TYPE Uint8;

#endif /* HAVE_INTTYPES_H */

/* BigScalar
 */

#define NCBI_BIG_TYPE Int8
#define SIZEOF_NCBI_BIG 8
#if (SIZEOF_LONG_DOUBLE > SIZEOF_NCBI_BIG)
#  undef  NCBI_BIG_TYPE
#  undef  SIZEOF_NCBI_BIG
#  define NCBI_BIG_TYPE   long double 
#  define SIZEOF_NCBI_BIG SIZEOF_LONG_DOUBLE
#endif
#if (SIZEOF_DOUBLE > SIZEOF_NCBI_BIG)
#  undef  NCBI_BIG_TYPE
#  undef  SIZEOF_NCBI_BIG
#  define NCBI_BIG_TYPE   double
#  define SIZEOF_NCBI_BIG SIZEOF_DOUBLE
#endif
#if (SIZEOF_VOIDP > SIZEOF_NCBI_BIG)
#  undef  NCBI_BIG_TYPE
#  undef  SIZEOF_NCBI_BIG
#  define NCBI_BIG_TYPE   void*
#  define SIZEOF_NCBI_BIG SIZEOF_VOIDP
#endif

/**
 * Define large scalar type.
 *
 * This is platform dependent. It could be an Int8, long double, double
 * or void*.
 */
typedef NCBI_BIG_TYPE Ncbi_BigScalar;


#ifndef HAVE_INTPTR_T
#  if SIZEOF_INT == SIZEOF_VOIDP
typedef int intptr_t;
#  elif SIZEOF_LONG == SIZEOF_VOIDP
typedef long intptr_t;
#  elif SIZEOF_LONG_LONG == SIZEOF_VOIDP
typedef long long intptr_t;
#  else
#    error No integer type is the same size as a pointer!
#  endif
#endif

#ifndef HAVE_UINTPTR_T
#  if SIZEOF_INT == SIZEOF_VOIDP
typedef unsigned int uintptr_t;
#  elif SIZEOF_LONG == SIZEOF_VOIDP
typedef unsigned long uintptr_t;
#  elif SIZEOF_LONG_LONG == SIZEOF_VOIDP
typedef unsigned long long uintptr_t;
#  else
#    error No integer type is the same size as a pointer!
#  endif
#endif


/* Macros for constant values definition 
 */

#ifdef HAVE_INTTYPES_H
#  define NCBI_CONST_INT8(v)     INT64_C(v)
#  define NCBI_CONST_UINT8(v)   UINT64_C(v)
#  define NCBI_INT8_FORMAT_SPEC  PRId64
#  define NCBI_UINT8_FORMAT_SPEC PRIu64
#elif (SIZEOF_LONG == 8)
#  define NCBI_CONST_INT8(v)   v##L
#  define NCBI_CONST_UINT8(v)  v##UL
#  define NCBI_INT8_FORMAT_SPEC   "ld"
#  define NCBI_UINT8_FORMAT_SPEC  "lu"
#elif (SIZEOF_LONG_LONG == 8)
#  define NCBI_CONST_INT8(v)   v##LL
#  define NCBI_CONST_UINT8(v)  v##ULL
#  if defined(__MINGW32__)  ||  defined(__MINGW64__)
#    define NCBI_INT8_FORMAT_SPEC   "I64d"
#    define NCBI_UINT8_FORMAT_SPEC  "I64u"
#  else
#    define NCBI_INT8_FORMAT_SPEC   "lld"
#    define NCBI_UINT8_FORMAT_SPEC  "llu"
#  endif
#elif defined(NCBI_USE_INT64)
#  define NCBI_CONST_INT8(v)   v##i64
#  define NCBI_CONST_UINT8(v)  v##ui64
#  define NCBI_INT8_FORMAT_SPEC   "I64d"
#  define NCBI_UINT8_FORMAT_SPEC  "I64u"
#else
#  define NCBI_CONST_INT8(v)   v
#  define NCBI_CONST_UINT8(v)  v
#  define NCBI_INT8_FORMAT_SPEC   "d"
#  define NCBI_UINT8_FORMAT_SPEC  "u"
#endif
#if (SIZEOF_LONG_DOUBLE > SIZEOF_DOUBLE)
#  define NCBI_CONST_LONGDOUBLE(v)   v##L
#else
#  define NCBI_CONST_LONGDOUBLE(v)   v
#endif


/* Undef auxiliaries
 */

#undef SIZEOF_NCBI_BIG
#undef NCBI_BIG_TYPE
#undef NCBI_INT8_TYPE


#endif  /* CORELIB___NCBITYPE__H */


/* @} */
