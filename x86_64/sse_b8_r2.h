
/**
 * @file sse_b8_r2.h
 *
 * @brief a header for macros of packed 8-bit SSE4 instructions.
 *
 * @detail
 * This is a collection of wrapper macros of SSE4.1 SIMD intrinsics.
 * Each macro corresponds to a set of several intrinsics defined in
 * smmintrin.h. The details of the intrinsics are found in the intel's
 * website: https://software.intel.com/sites/landingpage/IntrinsicsGuide/
 * The required set of macros are documented in porting section of 
 * README.md in the top directory of the library.
 *
 * @sa sse.h
 */
#ifndef _SSE_B8_R2_H_INCLUDED
#define _SSE_B8_R2_H_INCLUDED

#include <smmintrin.h>

/**
 * register declarations. 
 */
#define DECLARE_VEC_CELL(v)			__m128i v##1, v##2
#define DECLARE_VEC_CELL_REG(v)		__m128i register v##1, v##2
#define DECLARE_VEC_CHAR_REG(v)		__m128i register v##1, v##2

/**
 * substitution to cell vectors
 */
#define VEC_ASSIGN(a, b) { \
	(a##1) = (b##1); (a##2) = (b##2); \
}

#define VEC_SET(v, i) { \
	(v##1) = _mm_set1_epi8(i); \
	(v##2) = _mm_set1_epi8(i); \
}

#define VEC_SETZERO(v) { \
	(v##1) = _mm_setzero_si128(); \
	(v##2) = _mm_setzero_si128(); \
}

#define VEC_SETONES(v) { \
	(v##1) = _mm_set1_epi8(0xff); \
	(v##2) = _mm_set1_epi8(0xff); \
}

/**
 * substitution to char vectors
 */
#define VEC_CHAR_SETZERO(v) { \
	VEC_SETZERO(v); \
}

#define VEC_CHAR_SETONES(v) { \
	VEC_SETONES(v); \
}

/**
 * special substitution macros
 */
#define VEC_SET_LHALF(v, i) { \
	(v##1) = _mm_set1_epi8(i); \
	(v##2) = _mm_setzero_si128(); \
}

#define VEC_SET_UHALF(v, i) { \
	(v##1) = _mm_setzero_si128(); \
	(v##2) = _mm_set1_epi8(i); \
}

#define VEC_SETF_MSB(v) { \
	VEC_SETZERO(v); VEC_INSERT_MSB(v, 0xf0); \
}

#define VEC_SETF_LSB(v) { \
	VEC_SETZERO(v); VEC_INSERT_LSB(v, 0x0f); \
}

/**
 * insertion and extraction macros
 */
#define VEC_INSERT_MSB(v, a) { \
	(v##2) = _mm_insert_epi8((v##2), (a), sizeof(__m128i)-1); \
}

#define VEC_INSERT_LSB(v, a) { \
	(v##1) = _mm_insert_epi8((v##1), (a), 0); \
}

#define VEC_MSB(v)		( (signed char)_mm_extract_epi8((v##2), sizeof(__m128i)-1) )
#define VEC_LSB(v)		( (signed char)_mm_extract_epi8((v##1), 0) )
#define VEC_CENTER(v)	( (signed char)_mm_extract_epi8((v##2), 0) )

/**
 * arithmetic and logic operations
 */
#define VEC_OR(a, b, c) { \
	(a##1) = _mm_or_si128((b##1), (c##1)); \
	(a##2) = _mm_or_si128((b##2), (c##2)); \
}

#define VEC_ADD(a, b, c) { \
	(a##1) = _mm_adds_epi8((b##1), (c##1)); \
	(a##2) = _mm_adds_epi8((b##2), (c##2)); \
}

#define VEC_ADDS(a, b, c) { \
	(a##1) = _mm_adds_epu8((b##1), (c##1)); \
	(a##2) = _mm_adds_epu8((b##2), (c##2)); \
}

#define VEC_SUB(a, b, c) { \
	(a##1) = _mm_subs_epi8((b##1), (c##1)); \
	(a##2) = _mm_subs_epi8((b##2), (c##2)); \
}

#define VEC_SUBS(a, b, c) { \
	(a##1) = _mm_subs_epu8((b##1), (c##1)); \
	(a##2) = _mm_subs_epu8((b##2), (c##2)); \
}

#define VEC_MAX(a, b, c) { \
	(a##1) = _mm_max_epi8((b##1), (c##1)); \
	(a##2) = _mm_max_epi8((b##2), (c##2)); \
}

#define VEC_MIN(a, b, c) { \
	(a##1) = _mm_min_epi8((b##1), (c##1)); \
	(a##2) = _mm_min_epi8((b##2), (c##2)); \
}

/**
 * shift operations
 */
#define VEC_SHIFT_R(a) { \
	(a##1) = _mm_or_si128( \
		_mm_srli_si128((a##1), 1), \
		_mm_slli_si128((a##2), sizeof(__m128i)-1)); \
	(a##2) = _mm_srli_si128((a##2), 1); \
}

#define VEC_SHIFT_L(a) { \
	(a##2) = _mm_or_si128( \
		_mm_slli_si128((a##2), 1), \
		_mm_srli_si128((a##1), sizeof(__m128i)-1)); \
	(a##1) = _mm_slli_si128((a##1), 1); \
}

/**
 * compare and select
 */
#define VEC_COMPARE(a, b, c) { \
	(a##1) = _mm_cmpeq_epi8((b##1), (c##1)); \
	(a##2) = _mm_cmpeq_epi8((b##2), (c##2)); \
}

#define VEC_SELECT(a, b, c, d) { \
	(a##1) = _mm_blendv_epi8((b##1), (c##1), (d##1)); \
	(a##2) = _mm_blendv_epi8((b##2), (c##2), (d##2)); \
}

/**
 * load and store operations
 */
#define VEC_STORE(p, v) { \
	_mm_store_si128((__m128i *)(p), v##1); p += sizeof(__m128i); \
	_mm_store_si128((__m128i *)(p), v##2); p += sizeof(__m128i); \
}

#define VEC_STORE_PACKED(p, dv, dh) { \
	_mm_store_si128( \
		(__m128i *)(p), \
		_mm_or_si128(_mm_slli_epi64(dh##1, 4), dv##1)); \
	p += sizeof(__m128i); \
	_mm_store_si128( \
		(__m128i *)(p), \
		_mm_or_si128(_mm_slli_epi64(dh##2, 4), dv##2)); \
	p += sizeof(__m128i); \
}

/**
 * char vector operations
 */
#define VEC_CHAR_SHIFT_R(a) { \
	VEC_SHIFT_R(a); \
}

#define VEC_CHAR_SHIFT_L(a) { \
	VEC_SHIFT_L(a); \
}

#define VEC_CHAR_INSERT_MSB(x, y) { \
	VEC_INSERT_MSB(x, y); \
}

#define VEC_CHAR_INSERT_LSB(x, y) { \
	VEC_INSERT_LSB(x, y); \
}

#endif /* #ifndef _SSE_B8_R2_H_INCLUDED */
/**
 * end of sse_b8_r2.h
 */
