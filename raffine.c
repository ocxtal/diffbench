
/**
 * @file diff_affine_dynamic_banded.c
 *
 * @brief a affine-gap cost, banded difference dynamic programming algorithm.
 *
 * @detail
 * This is a difference-variant of the diag_affine_banded.c. The Smith-Waterman
 * algorithm is not supported in the diff algorithm.
 */

#include <stdint.h>
#include <stdlib.h>						/* for definition of the NULL */

#define BIT_WIDTH 			8
#define BAND_WIDTH 			32

#include "seqreader.h"
#include "diff.h"

#define ALG 				XSEA

/**
 * function declarations
 */
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

/**
 * @fn diff_affine_dynamic_banded
 *
 * @brief a affine-gap cost banded SIMD implementation.
 *
 * @param[ref] aln : a pointer to a sea_result structure. aln must NOT be NULL. (a structure which contains an alignment result.)
 * @param[in] aln->a : a pointer to a query sequence a.
 * @param[inout] aln->apos : (input) the start position of the search section on sequence a. (0 <= apos < length(sequence a)) (output) the start position of the alignment.
 * @param[inout] aln->alen : (input) the length of the search section on sequence a. (0 < alen) (output) the length of the alignment.
 * @param[in] aln->b : a pointer to a query sequence b.
 * @param[inout] aln->bpos : same as apos.
 * @param[inout] aln->blen : same as alen.
 *
 * @param[out] aln->aln : a pointer to the alignment string.
 * @param[inout] aln->len : (input) the reserved length of the aln->aln. (output) the length of the alignment string.
 * @param[out] aln->score : an alignment score.
 *
 * @param[in] scmat : unused parameter.
 * @param[in] xdrop : unused parameter. 
 * @param[in] bandwidth : unused parameter.
 * @param[in] mat : a pointer to the DP matrix.
 * @param[in] matsize : the size of the DP matrix.
 *
 * @return status. zero means success. see sea.h for the other error code.
 */
sea_int_t
diff_affine_dynamic_banded(
	struct sea_result *aln,
	struct sea_params param,
	void *mat)
{

	/**
	 * the Smith-Waterman algorithm is not supported in the diff algorithms.
	 */
	if(ALG == SW) {
		return 0;
	} else if(ALG == NW || ALG == SEA || ALG == XSEA) {
		sea_int_t retval = SEA_ERROR;
		struct mpos o;
		/**
		 * fill in a matrix
		 */
		o = diff_affine_dynamic_banded_fill(
			aln, param, (char *)mat);
		/**
		 * search maximum score position
		 */
		o = diff_affine_dynamic_banded_search(
			aln, param,
			(char *)mat,
			o);
		/**
		 * if aln->aln is not NULL, do traceback.
		 */
		if(o.m.p < 0) { return o.m.p; }
		if(aln->aln == NULL) { return SEA_SUCCESS; }
		retval = diff_affine_dynamic_banded_trace(
			aln, param,
			(char *)mat,
			o);
		return(retval);
	}
}


/**
 * @fn diff_affine_dynamic_banded_fill
 *
 * @brief a matrix fill-in function.
 *
 * @param[ref] aln : a pointer to a sea_result struct.
 * @param[in] xdrop : currently unused.
 * @param[in] mat : a pointer to a matrix.
 * @param[in] matlen : the length of the p-direction of the matrix.
 *
 * @return an end position of the extension.
 */
struct mpos
diff_affine_dynamic_banded_fill(
	struct sea_result *aln,
	struct sea_params param,
	char *mat)
{
	sea_int_t const bw = BAND_WIDTH;
	sea_int_t i, j;
	sea_int_t m = param.m,
			  x = param.x,
			  gi = param.gi,
			  ge = param.ge,
			  xdrop = param.xdrop;
	sea_int_t apos = 0,
			  bpos = 0;
	sea_int_t alen = aln->alen,
			  blen = aln->blen;
	sea_int_t alim = alen + bw/2,
			  blim = blen + bw/2;
	sea_int_t len = COP(alim, blim);
	sea_int_t score = 0,
			  max = 0,
			  scu, scl;
	sea_int_t dir;
	char *pdir = mat + AADDR(len, -bw/2);
	struct mpos o = {{0, 0, 0, 0}, {0, 0, 0, 0}};

	/**
	 * SIMD vector registers.
	 */
	DECLARE_VEC_CELL(mggv);
	DECLARE_VEC_CELL(xggv);
	DECLARE_VEC_CELL(geiv);
	DECLARE_VEC_CELL_REG(wq);
	DECLARE_VEC_CELL_REG(wt);
	DECLARE_VEC_CELL_REG(dv);
	DECLARE_VEC_CELL_REG(dh);
	DECLARE_VEC_CELL_REG(df);
	DECLARE_VEC_CELL_REG(de);
	DECLARE_VEC_CELL_REG(dv_);
	DECLARE_VEC_CELL_REG(tmp);
	DECLARE_VEC_CELL_REG(sum);

	/**
	 * seqreader declarations.
	 */
	DECLARE_SEQ(a);
	DECLARE_SEQ(b);

	CLEAR_SEQ(a, aln->a, aln->apos, aln->alen);
	CLEAR_SEQ(b, aln->b, aln->bpos, aln->blen);

	VEC_SET(mggv, m - 2*gi);		/** (m-2gi, m-2gi, ..., m-2gi) : constant vector */
	VEC_SET(xggv, x - 2*gi);		/** (x-2gi, x-2gi, ..., x-2gi) : constant vector */
	VEC_SET(geiv, ge - gi);			/** (ge-gi, ge-gi, ..., ge-gi) : constant vector */

	VEC_SET_LHALF(dv, x - 2*gi);	/** init vertical vector with (0, 0, ..., 0, m-2gi, ..., m-2gi) */
	VEC_SET_UHALF(dh, x - 2*gi);	/** init vertical vector with (m-2gi, ..., m-2gi, 0, 0, ..., 0) */
	VEC_SHIFT_L(dv);				/** wind back dv vector from p = 1 to p = 0 */
	VEC_SETZERO(df);				/** initialize df and de with 0 */
	VEC_SETZERO(de);
	VEC_STORE_DVDH(mat, dv, dh); VEC_STORE_DVDH(mat, df, de);
	*pdir++ = DIR_V;
	score = 0;						/** score at (0, 0) */
	scl = scu = (2*gi - m) * bw/2;	/** score at (0, -bw/2), (0, bw/2) */
	xdrop -= gi;					/** compensation */

	/**
	 * prefeed of sequence buffers
	 */
	VEC_SETZERO(wq);				/** a buffer for seq.a */
	VEC_SETONES(wt);				/** a buffer for seq.b */
	for(apos = 0; apos < bw/2; apos++) {
		FETCH(a, apos); PUSHQ((apos < alen) ? DECODE(a) : 0,    wq);
	}
	for(bpos = 0; bpos < bw/2-1; bpos++) {
		FETCH(b, bpos); PUSHT((bpos < blen) ? DECODE(b) : 0xff, wt);
	}

	/**
	 * @macro UPDATE_FIRST()
	 * @brief a set of the former half of the update operations.
	 */
	#define UPDATE_FIRST() { \
		VEC_ADD(de, de, geiv); \
		VEC_SUBS(de, de, dh); \
		VEC_ADD(df, df, geiv); \
		VEC_SUBS(df, df, dv); \
	}
	/**
	 * @macro UPDATE_MIDDLE()
	 * @brief a set of the latter half of the update operations.
	 */
	#define UPDATE_MIDDLE() { \
		VEC_COMPARE(tmp, wq, wt); \
		VEC_SELECT(tmp, xggv, mggv, tmp); \
		VEC_ADD(sum, dv, de); \
		VEC_MAX(tmp, tmp, sum); \
		VEC_ADD(sum, dh, df); \
		VEC_MAX(tmp, tmp, sum); \
		VEC_SUB(dv_, tmp, dh); \
		VEC_SUB(dh, tmp, dv); \
		VEC_ASSIGN(dv, dv_); \
	}
	/**
	 * @macro UPDATE_LAST
	 * @brief update score on the center, the upper edge and the lower edge.
	 */
	#define UPDATE_LAST() { \
		if(dir == DIR_V) { VEC_ASSIGN(tmp, dv); } else { VEC_ASSIGN(tmp, dh); } \
		score += (VEC_CENTER(tmp) + gi); \
		scl += (VEC_MSB(tmp) + gi); scu += (VEC_LSB(tmp) + gi); \
	}

	i = 0; j = 0;					/** the center cell of the init vector */		
	while(i < alim && j < blim) {
		UPDATE_FIRST();
		if(scl > scu) {
			VEC_SHIFT_R(dv); VEC_SHIFT_R(de);
			j++; dir = DIR_V;		/** move downward */
			FETCH(b, bpos); PUSHT((bpos < blen) ? DECODE(b) : 0xff, wt); bpos++;
		} else {
			VEC_SHIFT_L(dh); VEC_SHIFT_L(df);
			i++; dir = DIR_H;		/** move left */
			FETCH(a, apos); PUSHQ((apos < alen) ? DECODE(a) : 0, wq); apos++;
		}
		UPDATE_MIDDLE();
		VEC_STORE_DVDH(mat, dv, dh); VEC_STORE_DVDH(mat, df, de);
		*pdir++ = dir;				/** save direction */
		UPDATE_LAST();

		if((ALG == SEA || ALG == XSEA) && score >= max) {
			max = score;
			o.m.i = i; o.m.j = j;	/* examine cells on main diagonal and save maximum */
			o.m.p = COP(i, j); o.m.q = 0;
		}
		if(ALG == NW && COP(i, j) == COP(alen, blen)) { break; }
		if(ALG == XSEA && score + xdrop - max < 0) { break; }
	}
	#undef UPDATE_FIRST
	#undef UPDATE_MIDDLE
	#undef UPDATE_LAST

	aln->len = len;
	o.e.i = i; o.e.j = j; o.e.p = COP(i, j); o.e.q = COQ(i, j);
	if(ALG == NW) {
		o.m.i = i; o.m.j = j; o.m.p = COP(i, j); o.m.q = 0;
		aln->score = score;
	} else if(ALG == SEA || ALG == XSEA) {
		aln->score = max;
	}
	return(o);
}

/**
 * @fn diff_affine_dynamic_banded_search
 *
 * @brief search a cell with maximal score.
 *
 * @return a struct mpos which contains maximal score position.
 */
struct mpos
diff_affine_dynamic_banded_search(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o)
{
	sea_int_t const bw = BAND_WIDTH,
					bl = BYTES_PER_LINE;
	sea_int_t m = param.m,
			  x = param.x,
			  gi = param.gi,
			  ge = param.ge;
	sea_int_t sp, sq, ep;			/** position to wind back */
	sea_int_t mi = o.m.i,
			  mj = o.m.j,
			  mp = o.m.p,
			  mq = o.m.q;			/** temporary variables */
	sea_int_t p, q;					/** iteration variables */
	sea_int_t s, cs,				/** score temporary */
			  max = SCORE_MIN,		/** maximum */
			  score = aln->score;	/** score found in the fill-in step */
	char *smat = mat + AADDR(mp, mq), *cptr,
		 *pdir = mat + AADDR(aln->len, -bw/2) + mp;

	/**
	 * Maximum score search strategy:
	 * The maximum score search step consists of three step. At the end of
	 * the fill-in step, we know the score at the last-calculated cell on the 
	 * center line. We use the score of the cell, winding back the score toward
	 * the initial cell of the search, then restore the scores inside the 
	 * search area.
	 */
	#define EXACT(param) 		( (param.flags & SEA_FLAGS_MASK_POLICY) == SEA_EXACT )

	/**
	 * first determine the position to which the score is wound back.
	 */
	if(ALG == NW) {
		sp = o.e.p; sq = COQ(aln->alen, aln->blen) - o.e.q;
	} else if(ALG == SEA || ALG == XSEA) {
		if(!EXACT(param)) {
			sp = MAX2(mp + bw * (m+ge-3*gi) / (2*x), 0); sq = -bw/2;
			ep = MIN2(mp + bw * (m+ge-3*gi) / (2*m), o.e.p);
		} else {
			sp = 0; sq = -bw/2; ep = o.e.p;
			mi = 0; mj = 0; mp = 0; mq = 0; score = 0; max = aln->score;
			smat = mat + ADDR(mp, mq);			
		}
	}

	/**
	 * windback along with p-coordinate, then q-coordinate.
	 */
	while(mp > sp) {
		if(*pdir-- == DIR_V) {
			score -= DV(smat, gi); mj--;
		} else {
			score -= DH(smat, gi); mi--;
		}
		smat -= 2*bl; mp--;
	}
	while(mq < sq) { score += (DV(smat+1, gi) - DH(smat, gi)); smat++; mq++; }
	while(mq > sq) { score += (DH(smat-1, gi) - DV(smat, gi)); smat--; mq--; }

	/**
	 * search the score (in the seed-and-extend algorithm)
	 */
	if(ALG == NW) {
		o.m.i = aln->alen; o.m.j = aln->blen;
		o.m.p = mp; o.m.q = mq; aln->score = score;
	} else if(ALG == SEA || ALG == XSEA) {
		smat = mat + AADDR(mp, mq);
		cptr = mat + AADDR(mp, 0);	/** mp == 0 */
		pdir = mat + AADDR(aln->len, -bw/2) + mp;
		cs = 0;						/** score on the center line */
		for(p = mp; p <= ep; p++) {
			if(!EXACT(param) || cs > max - bw*(m+ge-3*gi)/2) {
				for(q = -bw/2, s = score; q < bw/2; q++) {
					if(s >= max) {
						max = s; o.m.i = mi-q; o.m.j = mj+q; o.m.p = p; o.m.q = q;
					}
					s += (DV(smat+1, gi) - DH(smat, gi)); smat++;		/** sizeof(CELL_TYPE) == 1 */
				}
				smat += (2*bl - bw); cptr += 2*bl;
			} else {
				smat += 2*bl; cptr += 2*bl;
			}
			if(*(++pdir) == DIR_V) {
				cs += DV(cptr, gi); score += DV(smat, gi); mj++;
			} else {
				cs += DH(cptr, gi); score += DH(smat, gi); mi++;
			}
		}
		aln->alen = o.m.i;
		aln->blen = o.m.j;
		aln->score = max;
	}
	return(o);
}

/**
 * @fn diff_affine_dynamic_banded_trace
 *
 * @brief traceback function.
 *
 * @param aln : a pointer to struct sea_result.
 * @param mat : a pointer to dynamic programming matrix.
 * @param matlen : the length of the p-coordinate of the matrix.
 * @param mpos : the start position of the trace.
 */
sea_int_t
diff_affine_dynamic_banded_trace(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o)
{
	sea_int_t const bw = BAND_WIDTH;
	sea_int_t mi = o.m.i,
			  mj = o.m.j,
			  mp = o.m.p,
			  mq = o.m.q;
	sea_int_t ge = param.ge,
			  gi = param.gi;
	sea_int_t dir;								/** 2-bit direction indicator (VV, VH, HV, or HH) */
	char *tmat = (char *)mat + AADDR(mp, mq),
		 *pdir = mat + AADDR(aln->len, -bw/2) + mp;

	char *p = aln->aln + aln->len - 1;
	char type = '\0';
	int64_t len = 0;

	#define PUSH(_p, _type, _len, _c) { \
		if((_type) != (_c)) { \
			*(_p)-- = (_type); \
			while((_len) > 0) { \
				int64_t _r = (_len) % 10; \
				(_len) = (_len) / 10; \
				*(_p)-- = _r + '0'; \
			} \
			(_type) = (_c); \
			(_len) = 1; \
		} else { \
			(_len)++; \
		} \
	}

	#define DET_DIR(dir, pdir) { \
		(dir) = ((dir)>>1) | ((*pdir--)<<1); \
	}

	dir = 0; DET_DIR(dir, pdir);
	while(mi > 0 && mj > 0) {
		DET_DIR(dir, pdir);
		if(DV(tmat, gi) == 0) {
			while(mj > 0 && DF(tmat, gi) == gi-ge) {
				tmat += DATOP(dir); DET_DIR(dir, pdir);
				mj--;
				PUSH(p, type, len, 'I');
			}
			tmat += DATOP(dir);
			mj--;
			PUSH(p, type, len, 'I');
		} else if(DH(tmat, gi) == 0) {
			while(mi > 0 && DE(tmat, gi) == gi-ge) {
				tmat += DALEFT(dir); DET_DIR(dir, pdir);
				mi--;
				PUSH(p, type, len, 'D');
			}
			tmat += DALEFT(dir);
			mi--;
			PUSH(p, type, len, 'D');
		} else {
			tmat += DTOPLEFT(dir); DET_DIR(dir, pdir);
			mi--;
			mj--;
			PUSH(p, type, len, 'M');
		}
	}
	while(mj > 0) { mj--; PUSH(p, type, len, 'I'); }
	while(mi > 0) { mi--; PUSH(p, type, len, 'D'); }

	char *b = aln->aln;
	int64_t l = 0;
	while(p < (char *)aln->aln + aln->len - 1) {
		*b++ = *++p;
		l++;
	}
	aln->len = l;

	#undef DET_DIR

	return SEA_SUCCESS;
}

/**
 * @fn diff_affine_dynamic_banded_matsize
 *
 * @brief returns the size of matrix for diff_affine_banded, in bytes.
 *
 * @param[in] alen : the length of the input sequence a (in base pairs).
 * @param[in] blen : the length of the input sequence b (in base pairs).
 * @param[in] bandwidth : unused. give zero for practice.
 *
 * @return the size of a matrix in bytes.
 */
sea_int_t
diff_affine_dynamic_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth)
{
	return(4 * (alen + blen + bandwidth + 4) * BYTES_PER_LINE);
}

#if 0
/**
 * @fn main
 *
 * @brief gets two ascii strings from stdin, align strings with naive_affine_full, and print the result.
 */
#ifdef MAIN

int main(int argc, char *argv[])
{
	int matlen, alnlen;
	void *mat;
	struct sea_result aln;
	struct sea_params param;

	param.flags = 0;
	param.m = 2;
	param.x = -3;
	param.gi = -4;
	param.ge = -1;
	param.xdrop = 12;
	param.bandwidth = 32;

	aln.a = argv[1];
	aln.alen = strlen(aln.a);
	aln.apos = 0;
	aln.b = argv[2];
	aln.blen = strlen(aln.b);
	aln.bpos = 0;
	alnlen = aln.alen + aln.blen;
	aln.len = alnlen;

	alnlen = aln.alen + aln.blen;
	matlen = CALL_FUNC(diff_affine_dynamic_banded_matsize, SUFFIX)(
		aln.alen, aln.blen, param.bandwidth);

	aln.aln = (void *)malloc(alnlen);
	mat = (void *)malloc(matlen);
	CALL_FUNC(diff_affine_dynamic_banded, SUFFIX)(&aln, param, mat);

	printf("%d, %d, %s\n", aln.score, aln.len, aln.aln);

	free(mat);
	free(aln.aln);
	return 0;
}

#endif /* #ifdef MAIN */

/**
 * unittest functions.
 *
 * give a definition of TEST to a compiler, e.g. -DTEST=1, to make a library for tests.
 */
#ifdef TEST
#if SEQ == ascii && ALN == ascii

extern sea_int_t
DECLARE_FUNC_GLOBAL(naive_affine_banded_matsize, DEFAULT_SUFFIX)(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth);

extern sea_int_t
DECLARE_FUNC_GLOBAL(naive_affine_banded, DEFAULT_SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	void *mat);

/**
 * @fn test_2_diff_affine_dynamic_banded
 *
 * @brief a unittest function of diff_affine_dynamic_banded.
 *
 * @detail
 * This function is an aggregation of simple fixed ascii queries.
 * In this function, a sea_assert_align macro in tester.h is called. This macro
 * calls sea_align function with given context, checks if the score and the alignment
 * string is the same as the given score and string. If both or one of the results
 * are different, the macro prints the failure status (filename, lines, input sequences,
 * and the content of sea_result) and dumps a content of dynamic programming matrix
 * to dumps.log.
 */
void
DECLARE_FUNC_GLOBAL(test_2_diff_affine_dynamic_banded, SUFFIX)(
	void)
{
	sea_int_t m = 2,
			  x = -3,
			  gi = -4,
			  ge = -1;
	struct sea_context *ctx;

	ctx = sea_init_fp(
		SEA_BANDWIDTH_64,
		CALL_FUNC(diff_affine_dynamic_banded, SUFFIX),
		CALL_FUNC(diff_affine_dynamic_banded_matsize, SUFFIX),
		m, x, gi, ge,		/** the default blast scoring scheme */
		12);				/** xdrop threshold */

	/**
	 * when both sequences are empty
	 */
	sea_assert_align(ctx, "", 				"", 			0,			"");
	/**
	 * when one sequence is empty
	 */
	sea_assert_align(ctx, "AAA", 			"", 			0,			"");
	sea_assert_align(ctx, "TTTTTTT", 		"", 			0,			"");
	sea_assert_align(ctx, "", 				"AAA", 			0,			"");
	sea_assert_align(ctx, "", 				"TTTTTTT", 		0,			"");
	sea_assert_align(ctx, "CTAG",			"", 			0,			"");

	/**
	 * a match
	 */
	sea_assert_align(ctx, "A",				"A", 			m,			"M");
	sea_assert_align(ctx, "C", 				"C", 			m,			"M");
	sea_assert_align(ctx, "G", 				"G", 			m,			"M");
	sea_assert_align(ctx, "T", 				"T", 			m,			"M");

	/**
	 * a mismatch
	 */
	if(ALG == NW) {
		/**
		 * the Needleman-Wunsch algorithm
		 */
		sea_assert_align(ctx, "A", 				"C", 			x,			"X");
		sea_assert_align(ctx, "C", 				"A", 			x,			"X");
		sea_assert_align(ctx, "G", 				"T", 			x,			"X");
		sea_assert_align(ctx, "T", 				"A", 			x,			"X");
	} else if(ALG == SEA || ALG == SW || ALG == XSEA) {
		/**
		 * the Smith-Waterman algorithm and the seed-and-extend algorithm
		 */
		sea_assert_align(ctx, "A", 				"C", 			0,			"");
		sea_assert_align(ctx, "C", 				"A", 			0,			"");
		sea_assert_align(ctx, "G", 				"T", 			0,			"");
		sea_assert_align(ctx, "T", 				"A", 			0,			"");
	}

	/**
	 * homopolymers with different lengths.
	 */
	if(ALG == NW) {
		/**
		 * the Needleman-Wunsch algorithm
		 */
		sea_assert_align(ctx, "A", 				"AA", 			gi+m,			"IM");
		sea_assert_align(ctx, "A", 				"AAA", 			gi+ge+m,		"IIM");
		sea_assert_align(ctx, "AAAA", 			"AA", 			gi+ge+m+m,		"DDMM");
		sea_assert_align(ctx, "TTTT", 			"TTTTTTTT", 	gi+3*ge+4*m,	"IIIIMMMM");
	} else if(ALG == SEA || ALG == SW || ALG == XSEA) {
		/**
		 * the Smith-Waterman algorithm and the seed-and-extend algorithm
		 */
		sea_assert_align(ctx, "A", 				"AA", 			m,				"M");
		sea_assert_align(ctx, "A", 				"AAA", 			m,				"M");
		sea_assert_align(ctx, "AAAA", 			"AA", 			m+m,			"MM");
		sea_assert_align(ctx, "TTTT", 			"TTTTTTTT", 	4*m,			"MMMM");
	}

	/**
	 * when mismatches occurs.
	 */
	sea_assert_align(ctx, "AAAAAAAA", 		"AAAATAAA", 	7*m+x,			"MMMMXMMM");
	sea_assert_align(ctx, "TTTTTTTT", 		"TTTCTTTT", 	7*m+x,			"MMMXMMMM");
	sea_assert_align(ctx, "CCCCCCCC", 		"CCTCCCCC", 	7*m+x,			"MMXMMMMM");
	sea_assert_align(ctx, "GGGGGGGG", 		"GGCGGTGG", 	6*m+2*x,		"MMXMMXMM");

	/**
	 * when gaps with a base occurs on seq a (insertion).
	 */
	sea_assert_align(ctx, "AAAAATTTT", 		"AAAAAGTTTT", 	9*m+gi,			"MMMMMIMMMM");
	sea_assert_align(ctx, "TTTTCCCCC", 		"TTTTACCCCC", 	9*m+gi,			"MMMMIMMMMM");
	sea_assert_align(ctx, "CCCGGGGGG", 		"CCCTGGGGGG", 	9*m+gi,			"MMMIMMMMMM");
	sea_assert_align(ctx, "GGGAATTT", 		"GGGCAAGTTT", 	8*m+2*gi,		"MMMIMMIMMM");

	/**
	 * when gaps with a base occurs on seq b (deletion).
	 */
	sea_assert_align(ctx, "AAAAAGTTTT", 	"AAAAATTTT", 	9*m+gi,			"MMMMMDMMMM");
	sea_assert_align(ctx, "TTTTACCCCC", 	"TTTTCCCCC", 	9*m+gi,			"MMMMDMMMMM");
	sea_assert_align(ctx, "CCCTGGGGGG", 	"CCCGGGGGG", 	9*m+gi,			"MMMDMMMMMM");
	sea_assert_align(ctx, "GGGCAAGTTT", 	"GGGAATTT", 	8*m+2*gi,		"MMMDMMDMMM");

	/**
	 * when a gap longer than two bases occurs on seq a.
	 */
	sea_assert_align(ctx, "AAAATTTT", 		"AAAAGGTTTT", 	8*m+gi+ge,		"MMMMIIMMMM");
	sea_assert_align(ctx, "GGGGCCCC", 		"GGGGTTTCCCC", 	8*m+gi+2*ge,	"MMMMIIIMMMM");
	sea_assert_align(ctx, "GGGGGCCCCC", 	"GGGGGTTTTCCCCC",10*m+gi+3*ge,	"MMMMMIIIIMMMMM");
	sea_assert_align(ctx, "TTTTAAGGGG", 	"TTTTCCAACCGGGG",10*m+2*gi+2*ge,"MMMMIIMMIIMMMM");

	/**
	 * when a gap longer than two bases occurs on seq b.
	 */
	sea_assert_align(ctx, "AAAAGGTTTT",	 	"AAAATTTT", 	8*m+gi+ge,		"MMMMDDMMMM");
	sea_assert_align(ctx, "GGGGTTTCCCC",	"GGGGCCCC", 	8*m+gi+2*ge,	"MMMMDDDMMMM");
	sea_assert_align(ctx, "GGGGGTTTTCCCCC",	"GGGGGCCCCC", 	10*m+gi+3*ge,	"MMMMMDDDDMMMMM");
	sea_assert_align(ctx, "TTTTCCAACCGGGG", "TTTTAAGGGG", 	10*m+2*gi+2*ge, "MMMMDDMMDDMMMM");

	/**
	 * X-drop test (here xdrop threshold is 12)
	 */
	if(ALG == XSEA) {
		sea_assert_align(ctx, "TTTTTTAAAAAAAA",	"GGGGGGAAAAAAAA",	0, 			"");
		sea_assert_align(ctx, "TTTTTAAAAAAAAA",	"GGGGGAAAAAAAAA",	0, 			"");
		sea_assert_align(ctx, "TTTTAAAAAAAA",	"GGGGAAAAAAAA",		8*m+4*x, 	"XXXXMMMMMMMM");
		sea_assert_align(ctx, "TTTAAAAAAAAA",	"GGGAAAAAAAAA",		9*m+3*x, 	"XXXMMMMMMMMM");

		sea_assert_align(ctx, "AAAATTTTTAAAAAAA","AAAAGGGGGAAAAAAA",	4*m, 		"MMMM");
		sea_assert_align(ctx, "AAAATTTAAAAAAAAA","AAAAGGGAAAAAAAAA",	13*m+3*x, 	"MMMMXXXMMMMMMMMM");
	}

	sea_clean(ctx);
	return;
}

/**
 * @fn test_8_cross_diff_affine_dynamic_banded
 *
 * @brief cross test between naive_affine_banded and diff_affine_dynamic_banded
 */
#if HAVE_NAIVE_BANDED

void
DECLARE_FUNC_GLOBAL(test_8_cross_diff_affine_dynamic_banded, SUFFIX)(
	void)
{
	int i;
	int const cnt = 5;
	sea_int_t m = 2,
			  x = -3,
			  gi = -4,
			  ge = -1;
	char *a, *b;
	struct sea_context *full, *band;

	#if BIT_WIDTH == 8
		int const len = 50;
	#elif BIT_WIDTH == 16
		int const len = 1000;
	#else
		#error "bit width must be 8 or 16."
	#endif

	full = sea_init_fp(
		SEA_BANDWIDTH_64,
		CALL_FUNC(naive_affine_banded, DEFAULT_SUFFIX),
		CALL_FUNC(naive_affine_banded_matsize, DEFAULT_SUFFIX),
		m, x, gi, ge,
		10000);
	band = sea_init_fp(
		SEA_BANDWIDTH_64,
		CALL_FUNC(diff_affine_dynamic_banded, SUFFIX),
		CALL_FUNC(diff_affine_dynamic_banded_matsize, SUFFIX),
		m, x, gi, ge,
		10000);


	for(i = 0; i < cnt; i++) {
		a = rseq(len);
		b = mseq(a, 20, 100, 100);
		sea_assert_cross(
			full, band, a, b);
		free(a); free(b);
	}

	sea_clean(full);
	sea_clean(band);
	return;
}

#endif /* #if HAVE_NAIVE_BANDED */

#endif /* #if SEQ == ascii && ALN == ascii */
#endif /* #ifdef TEST */

#endif
/**
 * end of diff_affine_dynamic_banded.c
 */
