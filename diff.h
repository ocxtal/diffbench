
/**
 * @file diff.h
 *
 * @brief a header for a variant management for diff algorithms.
 */
#ifndef _DIFF_H_INCLUDED
#define _DIFF_H_INCLUDED

#include <limits.h>
#include "sea.h"

#ifndef _SIMD_H_INCLUDED
#define _SIMD_H_INCLUDED

/**
 * Constants representing algorithms
 *
 * Notice: This constants must be consistent with the sea_flags_alg in sea.h.
 */
#define SW 								( 1 )
#define SEA 							( 2 )
#define XSEA 							( 3 )
#define NW 								( 6 )

/**
 * max and min
 */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MAX3(x,y,z) 	( MAX2(x, MAX2(y, z)) )
#define MAX4(w,x,y,z) 	( MAX2(MAX2(w, x), MAX2(y, z)) )

#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )
#define MIN3(x,y,z) 	( MIN2(x, MIN2(y, z)) )
#define MIN4(w,x,y,z) 	( MIN2(MIN2(w, x), MIN2(y, z)) )

/**
 * architecture flag definitions
 */
#define SSE 		( 1 )
#define AVX 		( 2 )

/**
 * @struct pos
 * @brief a struct containing a position
 */
struct pos {
	sea_int_t i, j, p, q;
};

/**
 * @struct mpos
 * @biref contains multiple position, max pos and end pos
 */
struct mpos {
	struct pos m;			/** max score pos */
	struct pos e;			/** end pos */
};

/**
 * bitwidth selection
 * The BIT_WIDTH constant represents the number of bits per cell on the memory.
 * The 8-bit wide diff algorithm uses the upper 4-bit of a byte to store DV
 * (the vertical difference) and the lower 4-bit of a byte to store DH.
 * The 16-bit wide diff algorithm uses a pair of bytes to store a pair of DV and DH.
 */
#if BIT_WIDTH == 8
 	/**
 	 * the bit width of the cell on memory: 8
 	 * the bit width of each DV or DH: 4
 	 * the bit width of the SIMD packed type: 8
 	 */
 	#define SIMD_BIT_WIDTH 				( 8 )
 	#define SIMD_BAND_WIDTH				( BAND_WIDTH )
	#define CELL_TYPE					char
	#define BYTES_PER_CELL				( sizeof(CELL_TYPE) )
	#define BYTES_PER_LINE 				( sizeof(CELL_TYPE)*BAND_WIDTH )
	#define DH(c, g)					( (*((unsigned CELL_TYPE *)(c))>>4) + g )
	#define DV(c, g)					( (*((unsigned CELL_TYPE *)(c)) & 0x0f) + g )
	#define DE(c, g)					( (*(((unsigned CELL_TYPE *)(c) + BAND_WIDTH))>>4) + g )
	#define DF(c, g)					( (*(((unsigned CELL_TYPE *)(c) + BAND_WIDTH)) & 0x0f) + g )
	#define VEC_STORE_DVDH(p, dv, dh)	{ VEC_STORE_PACKED(p, dv, dh); }
	#define CELL_MAX					( 0x0f )
	#define CELL_MIN					( 0 )
	#define MSB_MARKER					( 0xf0 )
	#define LSB_MARKER					( 0x0f )
 	#define SCORE_MIN					( INT_MIN )

#elif BIT_WIDTH == 16
 	/**
 	 * the bit width of the cell on memory: 16
 	 * the bit width of each DV or DH: 8
 	 * the bit width of the SIMD packed type: 8
 	 */
 	#define SIMD_BIT_WIDTH 				( 8 )
 	#define SIMD_BAND_WIDTH				( BAND_WIDTH )
 	#define CELL_TYPE					char
	#define BYTES_PER_CELL				( 2*sizeof(CELL_TYPE) )
	#define BYTES_PER_LINE 				( 2*sizeof(CELL_TYPE)*BAND_WIDTH )
	#define DH(c, g)					( (*((unsigned CELL_TYPE *)(c) + BAND_WIDTH)) + g )
	#define DV(c, g)					( (*(unsigned CELL_TYPE *)(c)) + g )
	#define DE(c, g)					( (*((unsigned CELL_TYPE *)(c) + 3*BAND_WIDTH)) + g )
	#define DF(c, g)					( (*((unsigned CELL_TYPE *)(c) + 2*BAND_WIDTH)) + g )
	#define VEC_STORE_DVDH(p, dv, dh)	{ VEC_STORE(p, dv); VEC_STORE(p, dh); }
	#define CELL_MAX					( UCHAR_MAX )
	#define CELL_MIN					( 0 )
	#define MSB_MARKER					( 0xf0 )
	#define LSB_MARKER					( 0x0f )
 	#define SCORE_MIN					( INT_MIN )

#else
 	#error "the BIT_WIDTH must be 8 or 16 in diff algorithms."
#endif

/**
 * p-direction length of search area.
 */
#define SEARCH_LEN				( BAND_WIDTH * 2 )

/**
 * char vector shift operations
 */
#define PUSHQ(x, y)					{ VEC_CHAR_SHIFT_L(y); VEC_CHAR_INSERT_LSB(y, x); }
#define PUSHT(x, y)					{ VEC_CHAR_SHIFT_R(y); VEC_CHAR_INSERT_MSB(y, x); }

/**
 * coordinate conversion macros (common for all algorithms)
 */
#define COX(p, q)				( ((p)>>1) - (q) )
#define COY(p, q)				( (((p)+1)>>1) + (q) )
#define COP(x, y)				( (x) + (y) )
#define COQ(x, y) 				( ((y)-(x))>>1 )

/**
 * direction macros. represents which direction the band came from. (from upward or leftward)
 */
#define DIR_V 					( 0x01 )
#define DIR_H 					( 0 )
#define DIR(p)					( ((p) & 0x01) ? DIR_V : DIR_H )

/**
 * address calculation macros for the linear-gap cost algorithms
 */
#define	ADDR(p, q)				( (BYTES_PER_LINE)*(p)+(sizeof(CELL_TYPE) * ((q)+BAND_WIDTH/2)) )
#define TOPQ(p, q) 				( - !((p)&0x01) * sizeof(CELL_TYPE) )
#define LEFTQ(p, q) 			( ((p)&0x01) * sizeof(CELL_TYPE) )
#define TOP(p, q)				( -(BYTES_PER_LINE) + TOPQ(p, q) )
#define LEFT(p, q)				( -(BYTES_PER_LINE) + LEFTQ(p, q) )
#define TOPLEFT(p, q) 			( -2*(BYTES_PER_LINE) )

#define DTOPQ(dir) 				( - !((dir)&0x01) * sizeof(CELL_TYPE) )
#define DLEFTQ(dir) 			( ((dir)&0x01) * sizeof(CELL_TYPE) )
#define DTOP(dir)				( -(BYTES_PER_LINE) + DTOPQ(dir) )
#define DLEFT(dir)				( -(BYTES_PER_LINE) + DLEFTQ(dir) )
#define DTOPLEFT(dir) 			( DTOP(dir) + DLEFT((dir)>>1) )

/**
 * address calculation macros for the affine-gap cost algorithms
 */
#define	AADDR(p, q)				( 2*(BYTES_PER_LINE)*(p)+(sizeof(CELL_TYPE) * ((q)+BAND_WIDTH/2)) )
#define ATOPQ(p, q) 			( - !((p)&0x01) * sizeof(CELL_TYPE) )
#define ALEFTQ(p, q) 			( ((p)&0x01) * sizeof(CELL_TYPE) )
#define ATOP(p, q)				( -(2 * BYTES_PER_LINE) + ATOPQ(p, q) )
#define ALEFT(p, q)				( -(2 * BYTES_PER_LINE) + ALEFTQ(p, q) )
#define ATOPLEFT(p, q) 			( -2*(2 * BYTES_PER_LINE) )

#define DATOPQ(dir) 			( - !((dir)&0x01) * sizeof(CELL_TYPE) )
#define DALEFTQ(dir) 			( ((dir)&0x01) * sizeof(CELL_TYPE) )
#define DATOP(dir)				( -(2 * BYTES_PER_LINE) + DATOPQ(dir) )
#define DALEFT(dir)				( -(2 * BYTES_PER_LINE) + DALEFTQ(dir) )
#define DATOPLEFT(dir) 			( DATOP(dir) + DALEFT((dir)>>1) )

/**
 * include SIMD intrinsic macros, depending on the value of BIT_WIDTH.
 */
#if defined(__AVX2__)
	#include "x86_64/avx.h"
#elif defined(__SSE4_1__)
	#include "x86_64/sse.h"
#else
 	#error "unsupported architecture. check definition of the 'ARCH' constant."
#endif

#endif /* _SIMD_H_INCLUDED */

struct mpos
diff_linear_dynamic_banded_fill(
	struct sea_result *aln,
	struct sea_params param,
	char *mat);

struct mpos
diff_linear_dynamic_banded_search(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
diff_linear_dynamic_banded_trace(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
diff_linear_dynamic_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth);


struct mpos
diff_affine_dynamic_banded_fill(
	struct sea_result *aln,
	struct sea_params param,
	char *mat);

struct mpos
diff_affine_dynamic_banded_search(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
diff_affine_dynamic_banded_trace(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
diff_affine_dynamic_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth);


#endif /* #ifndef _DIFF_H_INCLUDED */

/**
 * end of diff.h
 */
