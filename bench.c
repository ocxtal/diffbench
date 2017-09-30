
//#define DEBUG

/**
 * @file bench.c
 * @brief speed benchmark of libsea
 */
#include <alloca.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "bench.h"
#include "aed.h"
#include "gaba/gaba.h"
#include "edlib.h"
#include "ksw.h"
#include "kvec.h"
#include "align.h"

int enlarge_vector(Work_Data *work, int newmax);
int enlarge_points(Work_Data *work, int newmax);
int forward_wave(Work_Data *work, Align_Spec *spec, Alignment *align, Path *bpath,
                        int *mind, int maxd, int mida, int minp, int maxp);
int reverse_wave(Work_Data *work, Align_Spec *spec, Alignment *align, Path *bpath,
                        int mind, int maxd, int mida, int minp, int maxp);


#ifdef EDLIB_COUNT_BLOCKS
extern uint64_t edlib_blocks;
#endif

int seqan_editdist(uint8_t const *a, uint8_t const *b);


#define BIT_WIDTH			8
#define BAND_WIDTH 			32

#define M 					( 1 )
#define X 					( 1 )
#define GI 					( 1 )
#define GE 					( 1 )
#define XDROP				( 40 )


#include "ddiag.h"
#include "diff.h"

/**
 * @fn print_usage
 */
void print_usage(void)
{
	fprintf(stderr, "usage: bench -l <len> -c <cnt> -x <mismatch rate> -d <indel rate>\n");
}

static inline
void *aligned_malloc(
	size_t size)
{
	void *ptr = NULL;
	if(posix_memalign(&ptr, 32, size) != 0) {
		return(NULL);
	}
	return(ptr);
}

static inline
char random_base(void)
{
	char const table[4] = {'A', 'C', 'G', 'T'};
	// char const table[4] = { 0x01, 0x02, 0x04, 0x08 };
	return(table[rand() % 4]);
}

static inline
char *generate_random_sequence(int len)
{
	int i;
	char *seq;		/** a pointer to sequence */
	seq = (char *)malloc(sizeof(char) * (len + 32 + 1));
	if(seq == NULL) { return NULL; }
	for(i = 0; i < len; i++) {
		seq[i] = random_base();
	}
	seq[len] = '\0';
	return seq;
}

static inline
char *generate_mutated_sequence(char *seq, int len, double x, double d, int bw)
{
	int i, j, wave = 0;			/** wave is q-coordinate of the alignment path */
	char *mutated_seq;

	if(seq == NULL) { return NULL; }
	mutated_seq = (char *)malloc(sizeof(char) * (len + 32 + 1));
	if(mutated_seq == NULL) { return NULL; }
	for(i = 0, j = 0; i < len; i++) {
		if(((double)rand() / (double)RAND_MAX) < x) {
			mutated_seq[i] = random_base();	j++;	/** mismatch */
		} else if(((double)rand() / (double)RAND_MAX) < d) {
			if(rand() & 0x01 && wave > -bw+1) {
				mutated_seq[i] = (j < len) ? seq[j++] : random_base();
				j++; wave--;						/** deletion */
			} else if(wave < bw-2) {
				mutated_seq[i] = random_base();
				wave++;								/** insertion */
			} else {
				mutated_seq[i] = (j < len) ? seq[j++] : random_base();
			}
		} else {
			mutated_seq[i] = (j < len) ? seq[j++] : random_base();
		}
	}
	mutated_seq[len] = '\0';
	return mutated_seq;
}

static inline
char *add_tail(
	char *seq,
	char c,
	int64_t tail_len)
{
	int64_t len = strlen(seq);
	seq = realloc(seq, len + tail_len + 1);

	for(int64_t i = 0; i < tail_len; i++) {
		seq[len + i] = (c == 0) ? random_base() : c;
	}
	seq[len + tail_len] = '\0';
	return(seq);
}

static inline
char *add_margin(
	uint8_t *seq,
	int64_t len,
	int16_t head_margin,
	int16_t tail_margin)
{
	char *p = malloc(len + head_margin + tail_margin + 1);
	memset(p, 0, head_margin);
	memcpy(p + head_margin, seq, len);
	memset(p + head_margin + len, 0, tail_margin + 1);
	free(seq);
	return(p);
}

static inline
uint8_t encode_base(
	char c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases { A = 0x01, C = 0x02, G = 0x04, T = 0x08 };
	uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		[_b('N')] = A,		/* treat 'N' as 'A' */
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b((uint8_t)c)]);

	#undef _b
}

static inline
void encode(char *ptr, int64_t len)
{
	for(int64_t i = 0; i < len; i++) {
		ptr[i] = encode_base(ptr[i]);
	}
	return;
}

static inline
void alignment_to_cigar(char *cig, char const *aln, int64_t len)
{
	if(len <= 0) {
		return;
	}
	char state = aln[0], *p = cig;
	for(int64_t i = 1, cnt = 1; i < len; i++) {
		if(aln[i] == state) {
			cnt++;
			continue;
		}
		char buf[16] = { '0' };
		int64_t j = 0;
		for(j = 0; cnt; j++) {
			buf[j] = '0' + (cnt % 10); cnt /= 10;
		}
		for(int64_t k = 0; k <= j; k++) {
			*p++ = buf[j - k];
		}
		*p++ = state;
		// p += sprintf(p, "%ld%c", cnt, state);
		state = aln[i]; cnt = 1;
	}
	return;
}

struct params {
	char const *file;
	uint64_t alen, blen, mlen;
	uint64_t cnt, max;
	double frac;
	kvec_t(char) buf;
	kvec_t(char *) seq;
	kvec_t(uint64_t) len;
	int table;
};

int parse_args(struct params *p, int c, char *arg)
{
	switch(c) {
		case 'i': p->file = arg; break;
		case 'l': p->max = atoi(arg); break;
		case 'c': p->cnt = atoi(arg); break;
		case 'f': p->frac = atof(arg); break;
		case 'a': printf("%s\n", arg); return 0;
		case 't': p->table = 1; return 0;
		/**
		 * the others: print help message
		 */
		case 'h':
		default: print_usage(); return -1;
	}
	return 0;
}


struct naive_result_s {
	int32_t score;
	uint32_t path_length;
	int64_t apos, bpos;
	int64_t alen, blen;
	char *path;
};

static inline
struct naive_result_s naive_affine(
	char const *a,
	char const *b)
{
	/* utils */
	#define _a(p, q, plen)	( (q) * ((plen) + 1) + (p) )
	#define s(p, q)			_a(p, 3*(q), alen)
	#define e(p, q)			_a(p, 3*(q)+1, alen)
	#define f(p, q)			_a(p, 3*(q)+2, alen)
	#define m(p, q)			( a[(p) - 1] == b[(q) - 1] ? m : x )

	/* load gap penalties */
	int8_t m = M;
	int8_t x = -X;
	int8_t gi = -GI;
	int8_t ge = -GE;

	/* calc lengths */
	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	/* calc min */
	int64_t min = INT16_MIN - x - 2*gi;

	/* malloc matrix */
	int16_t *mat = (int16_t *)malloc(
		3 * (alen + 1) * (blen + 1) * sizeof(int16_t));

	/* init */
	struct naive_maxpos_s {
		int16_t score;
		int64_t apos;
		int64_t bpos;
	};

	struct naive_maxpos_s max = { 0, 0, 0 };

	mat[s(0, 0)] = mat[e(0, 0)] = mat[f(0, 0)] = 0;
	for(int64_t i = 1; i < alen+1; i++) {
		mat[s(i, 0)] = mat[e(i, 0)] = MAX2(min, gi + i * ge);
		mat[f(i, 0)] = MAX2(min, gi + i * ge + gi - 1);
	}
	for(int64_t j = 1; j < blen+1; j++) {
		mat[s(0, j)] = mat[f(0, j)] = MAX2(min, gi + j * ge);
		mat[e(0, j)] = MAX2(min, gi + j * ge + gi - 1);
	}

	for(int64_t j = 1; j < blen+1; j++) {
		for(int64_t i = 1; i < alen+1; i++) {
			int16_t score_e = mat[e(i, j)] = MAX2(
				mat[s(i - 1, j)] + gi + ge,
				mat[e(i - 1, j)] + ge);
			int16_t score_f = mat[f(i, j)] = MAX2(
				mat[s(i, j - 1)] + gi + ge,
				mat[f(i, j - 1)] + ge);
			int16_t score = mat[s(i, j)] = MAX4(min,
				mat[s(i - 1, j - 1)] + m(i, j),
				score_e, score_f);
			if(score > max.score
			|| (score == max.score && (i + j) < (max.apos + max.bpos))) {
				max = (struct naive_maxpos_s){
					score, i, j
				};
			}
		}
	}
	if(max.score == 0) {
		max = (struct naive_maxpos_s){ 0, 0, 0 };
	}

	struct naive_result_s result = {
		.score = max.score,
		.apos = max.apos,
		.bpos = max.bpos,
		.path_length = max.apos + max.bpos + 1,
		.path = (char *)malloc(max.apos + max.bpos + 64)
	};
	int64_t path_index = max.apos + max.bpos + 1;
	while(max.apos > 0 || max.bpos > 0) {
		/* M > I > D > X */
		if(mat[s(max.apos, max.bpos)] == mat[f(max.apos, max.bpos)]) {
			while(mat[f(max.apos, max.bpos)] == mat[f(max.apos, max.bpos - 1)] + ge) {
				max.bpos--;
				result.path[--path_index] = 'I';
			}
			max.bpos--;
			result.path[--path_index] = 'I';
		} else if(mat[s(max.apos, max.bpos)] == mat[e(max.apos, max.bpos)]) {
			while(mat[e(max.apos, max.bpos)] == mat[e(max.apos - 1, max.bpos)] + ge) {
				max.apos--;
				result.path[--path_index] = 'D';
			}
			max.apos--;
			result.path[--path_index] = 'D';
		} else {
			result.path[--path_index] = a[max.apos - 1] == b[max.bpos - 1] ? 'M' : 'X';
			// result.path[--path_index] = 'D';
			max.apos--;
			max.bpos--;
		}
	}

	result.alen = result.apos - max.apos;
	result.blen = result.bpos - max.bpos;
	result.apos = max.apos;
	result.bpos = max.bpos;

	result.path_length -= path_index;
	for(uint64_t i = 0; i < result.path_length; i++) {
		result.path[i] = result.path[path_index++];
	}
	result.path[result.path_length] = '\0';
	free(mat);

	#undef _a
	#undef s
	#undef e
	#undef f
	#undef m
	return(result);
}


struct bench_pair_s {
	bench_t fill, trace, conv;
	int64_t score, fail;
};

static inline
struct bench_pair_s bench_adaptive_editdist(
	struct params p)
{
	void *base = aligned_malloc(2 * 4 * sizeof(uint64_t) * (p.alen + p.blen + 65));
	memset(base, 0, 4 * sizeof(uint64_t));
	void *ptr = base + 4 * sizeof(uint64_t);

	char *buf = (char *)aligned_malloc(p.alen + p.blen);
	char *cigar = malloc(sizeof(char) * (p.alen + p.blen + 1));

	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		int64_t klim = MAX2((int64_t)(1.5 * (double)(kv_at(p.len, i * 2) + kv_at(p.len, i * 2 + 1)) * 0.2 + 0.5), 10);

		bench_start(fill);
		struct aed_fill_s f = aed_fill(ptr, (uint8_t const *)kv_at(p.seq, i * 2), kv_at(p.len, i * 2), (uint8_t const *)kv_at(p.seq, i * 2 + 1), kv_at(p.len, i * 2 + 1), klim);
		score += f.score;
		bench_end(fill);

		bench_start(trace);
		int64_t len = aed_trace(buf, kv_at(p.len, i * 2) + kv_at(p.len, i * 2 + 1), ptr, f);
		bench_end(trace);
		// fprintf(stderr, "len(%ld)\n", len);

		bench_start(conv);
		alignment_to_cigar(cigar, buf, len);
		bench_end(conv);
	}
	free(base);
	free(buf);
	free(cigar);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.score = score
	});
}

static inline
struct bench_pair_s bench_ddiag_linear(
	struct params p)
{
	struct sea_params param = { 0, M, -X, -(GI + GE), 0, 0, 0, 50, 32 };


	struct sea_result *aln = (struct sea_result *)aligned_malloc(sizeof(char) * (p.alen + p.blen + 32 + 1) + sizeof(struct sea_result));
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diag_linear_dynamic_banded_matsize(p.alen, p.blen, 32);
	void *base = aligned_malloc(alnsize);
	memset(base, 0, 32 * 2 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 2 * sizeof(int16_t);

	char *cigar = malloc(sizeof(char) * (p.alen + p.blen + 1));

	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	int64_t score = 0, fail = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		aln->score = 0;
		aln->a = kv_at(p.seq, i * 2);     aln->apos = 0; aln->alen = kv_at(p.len, i * 2);
		aln->b = kv_at(p.seq, i * 2 + 1); aln->bpos = 0; aln->blen = kv_at(p.len, i * 2 + 1);
		aln->len = alnsize;

		bench_start(fill);
		struct mpos o = diag_linear_dynamic_banded_fill(aln, param, mat);
		o = diag_linear_dynamic_banded_search(aln, param, mat, o);
		score += aln->score;
		bench_end(fill);


		bench_start(trace);
		diag_linear_dynamic_banded_trace(aln, param, mat, o);
		bench_end(trace);
		if(0.8 * (kv_at(p.len, i * 2) - p.mlen) > aln->alen || 0.8 * (kv_at(p.len, i * 2 + 1) - p.mlen) > aln->blen) {
			fail++;
		}


		bench_start(conv);
		alignment_to_cigar(cigar, aln->aln, aln->len);
		bench_end(conv);
	}
	free(aln);
	free(base);
	free(cigar);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.score = score,
		.fail = fail
	});
}

static inline
struct bench_pair_s bench_ddiag_affine(
	struct params p)
{
	struct sea_params param = { 0, M, -X, -(GI + GE), -GE, 0, 0, 50, 32 };


	struct sea_result *aln = (struct sea_result *)aligned_malloc(sizeof(char) * (p.alen + p.blen + 32 + 1) + sizeof(struct sea_result));
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diag_affine_dynamic_banded_matsize(p.alen, p.blen, 32);
	void *base = aligned_malloc(alnsize);
	memset(base, 0, 32 * 6 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 6 * sizeof(int16_t);

	char *cigar = malloc(sizeof(char) * (p.alen + p.blen + 1));

	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	int64_t score = 0, fail = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		aln->score = 0;
		aln->a = kv_at(p.seq, i * 2);     aln->apos = 0; aln->alen = kv_at(p.len, i * 2);
		aln->b = kv_at(p.seq, i * 2 + 1); aln->bpos = 0; aln->blen = kv_at(p.len, i * 2 + 1);
		aln->len = alnsize;


		bench_start(fill);
		struct mpos o = diag_affine_dynamic_banded_fill(aln, param, mat);
		o = diag_affine_dynamic_banded_search(aln, param, mat, o);
		score += aln->score;
		bench_end(fill);


		bench_start(trace);
		diag_affine_dynamic_banded_trace(aln, param, mat, o);
		bench_end(trace);
		if(0.8 * (kv_at(p.len, i * 2) - p.mlen) > aln->alen || 0.8 * (kv_at(p.len, i * 2 + 1) - p.mlen) > aln->blen) {
			fail++;
		}
		// fprintf(stderr, "(%ld, %ld), (%d, %d), %d\n", kv_at(p.len, i * 2), kv_at(p.len, i * 2 + 1), aln->alen, aln->blen, aln->score);


		bench_start(conv);
		alignment_to_cigar(cigar, aln->aln, aln->len);
		bench_end(conv);
	}
	free(aln);
	free(base);
	free(cigar);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.score = score,
		.fail = fail
	});
}

static inline
struct bench_pair_s bench_diff_linear(
	struct params p)
{
	struct sea_params param = { 0, M, -X, -(GI + GE), 0, 0, 0, 50, 32 };


	struct sea_result *aln = (struct sea_result *)aligned_malloc(sizeof(char) * (p.alen + p.blen + 32 + 1) + sizeof(struct sea_result) + 1000000);
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diff_linear_dynamic_banded_matsize(p.alen, p.blen, 32);
	void *base = aligned_malloc(alnsize);
	memset(base, 0, 32 * 2 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 2 * sizeof(int16_t);

	char *cigar = malloc(sizeof(char) * (p.alen + p.blen + 1));

	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	int64_t score = 0, fail = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		aln->score = 0;
		aln->a = kv_at(p.seq, i * 2);     aln->apos = 0; aln->alen = kv_at(p.len, i * 2);
		aln->b = kv_at(p.seq, i * 2 + 1); aln->bpos = 0; aln->blen = kv_at(p.len, i * 2 + 1);
		aln->len = alnsize;


		bench_start(fill);
		struct mpos o = diff_linear_dynamic_banded_fill(aln, param, mat);
		bench_end(fill);


		bench_start(trace);
		o = diff_linear_dynamic_banded_search(aln, param, mat, o);
		score += aln->score;
		diff_linear_dynamic_banded_trace(aln, param, mat, o);
		bench_end(trace);
		if(0.8 * (kv_at(p.len, i * 2) - p.mlen) > aln->alen || 0.8 * (kv_at(p.len, i * 2 + 1) - p.mlen) > aln->blen) {
			fail++;
		}


		bench_start(conv);
		alignment_to_cigar(cigar, aln->aln, aln->len);
		bench_end(conv);
	}
	free(aln);
	free(base);
	free(cigar);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.score = score,
		.fail = fail
	});
}

static inline
struct bench_pair_s bench_diff_affine(
	struct params p)
{
	struct sea_params param = { 0, M, -X, -(GI + GE), -GE, 0, 0, 50, 32 };


	struct sea_result *aln = (struct sea_result *)aligned_malloc(sizeof(char) * (p.alen + p.blen + 32 + 1) + sizeof(struct sea_result));
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diff_affine_dynamic_banded_matsize(p.alen, p.blen, 32);
	void *base = aligned_malloc(alnsize);
	memset(base, 0, 32 * 6 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 6 * sizeof(int16_t);

	char *cigar = malloc(sizeof(char) * (p.alen + p.blen + 1));

	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	int64_t score = 0, fail = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		aln->score = 0;
		aln->a = kv_at(p.seq, i * 2);     aln->apos = 0; aln->alen = kv_at(p.len, i * 2);
		aln->b = kv_at(p.seq, i * 2 + 1); aln->bpos = 0; aln->blen = kv_at(p.len, i * 2 + 1);
		aln->len = alnsize;


		bench_start(fill);
		struct mpos o = diff_affine_dynamic_banded_fill(aln, param, mat);
		bench_end(fill);


		bench_start(trace);
		o = diff_affine_dynamic_banded_search(aln, param, mat, o);
		score += aln->score;
		diff_affine_dynamic_banded_trace(aln, param, mat, o);
		bench_end(trace);
		if(0.8 * (kv_at(p.len, i * 2) - p.mlen) > aln->alen || 0.8 * (kv_at(p.len, i * 2 + 1) - p.mlen) > aln->blen) {
			fail++;
		}
		// fprintf(stderr, "(%ld, %ld), (%d, %d), %d\n", kv_at(p.len, i * 2), kv_at(p.len, i * 2 + 1), aln->alen, aln->blen, aln->score);


		bench_start(conv);
		alignment_to_cigar(cigar, aln->aln, aln->len);
		bench_end(conv);
	}
	free(aln);
	free(base);
	free(cigar);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.score = score,
		.fail = fail
	});
}

static inline
struct bench_pair_s bench_gaba_linear(
	struct params p)
{
	char *c = (char *)malloc(p.alen + p.blen + 10);

	/** init context */
	gaba_t *ctx = gaba_init(GABA_PARAMS(
		.xdrop = XDROP,
		GABA_SCORE_SIMPLE(M, X, 0, GI + GE)));

	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	void const *lim = (void const *)0x800000000000;

	gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);
	int64_t score = 0, fail = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		gaba_dp_flush(dp, lim, lim);
		struct gaba_section_s asec = gaba_build_section(0, (uint8_t const *)kv_at(p.seq, i * 2),     kv_at(p.len, i * 2));
		struct gaba_section_s bsec = gaba_build_section(2, (uint8_t const *)kv_at(p.seq, i * 2 + 1), kv_at(p.len, i * 2 + 1));

		bench_start(fill);
		struct gaba_fill_s *f = gaba_dp_fill_root(dp, &asec, 0, &bsec, 0);
		score += f->max;
		bench_end(fill);
		
		bench_start(trace);
		struct gaba_alignment_s *r = gaba_dp_trace(dp, f, NULL, NULL);
		bench_end(trace);
		if(0.8 * (kv_at(p.len, i * 2) - p.mlen) > r->sec->alen || 0.8 * (kv_at(p.len, i * 2 + 1) - p.mlen) > r->sec->blen) {
			fail++;
		}

		bench_start(conv);
		gaba_dp_dump_cigar_forward(c, p.alen + p.blen, r->path->array, 0, r->path->len);
		bench_end(conv);
	}
	gaba_dp_clean(dp);
	gaba_clean(ctx);
	free(c);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.score = score,
		.fail = fail
	});
}

static inline
void print_path(
	struct gaba_alignment_s *r)
{
	uint64_t plen = r->path->len, cnt = 0;
	uint32_t const *path = r->path->array;
	uint64_t path_array = (uint64_t)path[0] | (((uint64_t)path[1])<<32);
	char *ptr = alloca(plen);
	char *p = ptr;

	path += 2;
	while(plen > 0) {
		if((path_array & 0x03) == 0x01) {
			*p++ = 'M'; plen -= 2;
			path_array >>= 1;
			if(++cnt == 32) {
				path_array |= ((uint64_t)*path++)<<32;
				cnt = 0;
			}
			path_array >>= 1;
			if(++cnt == 32) {
				path_array |= ((uint64_t)*path++)<<32;
				cnt = 0;
			}
		} else {
			*p++ = (path_array & 0x01) ? 'I' : 'D'; plen--;
			path_array >>= 1;
			if(++cnt == 32) {
				path_array |= ((uint64_t)*path++)<<32;
				cnt = 0;
			}
		}
	}
	*p = '\0';
	fprintf(stderr, "%s\n", ptr);
	return;
}

static inline
struct bench_pair_s bench_gaba_affine(
	struct params p)
{
	char *c = (char *)malloc(p.alen + p.blen + 10);

	/** init context */
	gaba_t *ctx = gaba_init(GABA_PARAMS(
		.xdrop = 50,
		GABA_SCORE_SIMPLE(M, X, GI, GE)));

	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	void const *lim = (void const *)0x800000000000;

	gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);
	int64_t score = 0, fail = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		gaba_dp_flush(dp, lim, lim);
		struct gaba_section_s asec = gaba_build_section(0, (uint8_t const *)kv_at(p.seq, i * 2),     kv_at(p.len, i * 2));
		struct gaba_section_s bsec = gaba_build_section(2, (uint8_t const *)kv_at(p.seq, i * 2 + 1), kv_at(p.len, i * 2 + 1));

		bench_start(fill);
		struct gaba_fill_s *f = gaba_dp_fill_root(dp, &asec, 0, &bsec, 0);
		score += f->max;
		bench_end(fill);
		
		bench_start(trace);
		struct gaba_alignment_s *r = gaba_dp_trace(dp, f, NULL, NULL);
		bench_end(trace);


		if(0.8 * (kv_at(p.len, i * 2) - p.mlen) > r->sec->alen || 0.8 * (kv_at(p.len, i * 2 + 1) - p.mlen) > r->sec->blen) {
			fail++;

#ifdef DEBUG
			struct naive_result_s nr = naive_affine(kv_at(p.seq, i * 2), kv_at(p.seq, i * 2 + 1));

			fprintf(stderr, "(%ld, %ld), (%u, %u), %ld, (%ld, %ld)\n", kv_at(p.len, i * 2), kv_at(p.len, i * 2 + 1), r->sec->alen, r->sec->blen, f->max, nr.alen, nr.blen);
			fprintf(stderr, "len(%lu)\n", strlen(nr.path));

			char const *a = kv_at(p.seq, i * 2), *b = kv_at(p.seq, i * 2 + 1);
			for(uint64_t j = 0, k = 0; j < strlen(nr.path); j++) {
				if(nr.path[j] != 'I') {
					fprintf(stderr, "%c", " AC G   T"[a[k++]]);
				} else {
					fprintf(stderr, "-");
				}
			}
			fprintf(stderr, "\n");
			for(uint64_t j = 0, k = 0; j < strlen(nr.path); j++) {
				if(nr.path[j] != 'D') {
					fprintf(stderr, "%c", " AC G   T"[b[k++]]);
				} else {
					fprintf(stderr, "-");
				}
			}
			fprintf(stderr, "\n");

			fprintf(stderr, "%s\n", nr.path);
			print_path(r);

			free(nr.path);
#endif

		}


		if(r == NULL) {
			fprintf(stderr, "%s\n%s\n", kv_at(p.seq, i * 2), kv_at(p.seq, i * 2 + 1));
			continue;
		}

		bench_start(conv);
		gaba_dp_dump_cigar_forward(c, p.alen + p.blen, r->path->array, 0, r->path->len);
		bench_end(conv);
	}
	gaba_dp_clean(dp);
	gaba_clean(ctx);
	free(c);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.score = score,
		.fail = fail
	});
}

static inline
struct bench_pair_s bench_edlib(
	struct params p)
{
	/** init context */
	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	#ifdef EDLIB_COUNT_BLOCKS
	edlib_blocks = 0;
	#endif


	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		// int64_t klim = MAX2((int64_t)(1.5 * (double)(kv_at(p.len, i * 2) + kv_at(p.len, i * 2 + 1)) * 0.2 + 0.5), 10);
		int64_t klim = 10000;
		EdlibAlignConfig cf = (EdlibAlignConfig){ .k = klim, .mode = EDLIB_MODE_SHW, .task = EDLIB_TASK_DISTANCE };
		EdlibAlignConfig ct = (EdlibAlignConfig){ .k = klim, .mode = EDLIB_MODE_SHW, .task = EDLIB_TASK_PATH };

		bench_start(fill);
		EdlibAlignResult f = edlibAlign(kv_at(p.seq, i * 2), kv_at(p.len, i * 2), kv_at(p.seq, i * 2 + 1), kv_at(p.len, i * 2 + 1), cf);
		score += f.editDistance < 0;
		bench_end(fill);

		edlibFreeAlignResult(f);

		#ifndef EDLIB_COUNT_BLOCKS
		bench_start(trace);
		EdlibAlignResult t = edlibAlign(kv_at(p.seq, i * 2), kv_at(p.len, i * 2), kv_at(p.seq, i * 2 + 1), kv_at(p.len, i * 2 + 1), ct);
		bench_end(trace);
		
		bench_start(conv);
		char *cigar = edlibAlignmentToCigar(t.alignment, t.alignmentLength, EDLIB_CIGAR_STANDARD);
		bench_end(conv);
		free(cigar);

		edlibFreeAlignResult(t);
		#endif
	}

	trace.a -= fill.a;

	#ifdef EDLIB_COUNT_BLOCKS
	fprintf(stderr, "blocks(%llu), cells(%llu)\n", edlib_blocks, edlib_blocks * 64);
	#endif

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.score = score
	});
}

static inline
struct bench_pair_s bench_bwamem(
	struct params p)
{
	/** init context */
	uint64_t cells = 0;
	bench_t fill, trace, conv;
	bench_init(fill);
	bench_init(trace);
	bench_init(conv);

	int8_t const mat[16] = { M, -X, -X, -X, -X, M, -X, -X, -X, -X, M, -X, -X, -X, -X, M };
	uint32_t *cigar = malloc(sizeof(uint32_t) * 10000);
	for(int64_t i = 0; i < p.cnt; i++) {
		int max_off, rlen, qlen;
		uint64_t s[3] = {0};

		/*
		s[i+1] = ksw_extend2(
			kv_at(p.len, i * 2), (uint8_t const *)kv_at(p.seq, i * 2),
			kv_at(p.len, i * 2 + 1), (uint8_t const *)kv_at(p.seq, i * 2 + 1),
			4, mat, GI, GE, GI, GE,
			100, 0, 100, M * kv_at(p.len, i * 2), &rlen, &qlen, NULL, NULL, &max_off);
		*/
		int w = abs(kv_at(p.len, i * 2) - kv_at(p.len, i * 2 + 1)) + 3;
		int n_cigar = 0;
		ksw_global2_bench(
			kv_at(p.len, i * 2), (uint8_t const *)kv_at(p.seq, i * 2),
			kv_at(p.len, i * 2 + 1), (uint8_t const *)kv_at(p.seq, i * 2 + 1),
			4, mat, GI, GE, GI, GE, w, &n_cigar, &cigar,
			&fill, &trace);
		cells += w * kv_at(p.len, i * 2 + 1);
	}
	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.conv = conv,
		.fail = cells
	});
}


struct wavefront_work_s {
	Work_Data *w;
	Align_Spec *s;
	Path apath, bpath;
};

struct wavefront_work_s *wavefront_init_work(double id)
{
	struct wavefront_work_s *work = (struct wavefront_work_s *)malloc(sizeof(struct wavefront_work_s));

	Work_Data *w = New_Work_Data();
	enlarge_vector(w, 10000 * (6 * sizeof(int) + sizeof(uint64_t)));
	enlarge_points(w, 4 * 4096 / 100 * sizeof(uint16_t) + sizeof(Path));


	float freq[4] = { 0.25, 0.25, 0.25, 0.25 };
	Align_Spec *s = New_Align_Spec(id, 100, freq);


	work->w = w;
	work->s = s;
	work->apath.trace = malloc(4096);
	work->bpath.trace = malloc(4096);

	return(work);
}


void wavefront_clean_work(struct wavefront_work_s *_work)
{
	struct wavefront_work_s *work = (struct wavefront_work_s *)_work;

	free(work->apath.trace);
	free(work->bpath.trace);

	Free_Align_Spec(work->s);
	Free_Work_Data(work->w);

	free(work);
	return;	
}

int 
wavefront(
	void *_work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen)
{
	struct wavefront_work_s *work = (struct wavefront_work_s *)_work;
	Work_Data *w = work->w;
	Align_Spec *s = work->s;

	Alignment aln;
	aln.path = &work->apath;
	aln.flags = 0;
	aln.aseq = (char *)a;
	aln.bseq = (char *)b;
	aln.alen = alen;
	aln.blen = blen;

	int low = 0;
	int high = 0;

	forward_wave(w, s, &aln, &work->bpath, &low, high, 0, -INT32_MAX, INT32_MAX);
	// fprintf(stderr, "(%u, %u), (%u, %u)\n", alen, blen, aln.path->aepos, aln.path->bepos);
	return(aln.path->aepos + aln.path->bepos);
}

static inline
struct bench_pair_s bench_seqan(
	struct params p)
{
	/** init context */
	bench_t fill;
	bench_init(fill);

	for(int64_t i = 0; i < p.cnt; i++) {
		bench_start(fill);
		seqan_editdist((uint8_t const *)kv_at(p.seq, i * 2), (uint8_t const *)kv_at(p.seq, i * 2 + 1));
		bench_end(fill);
	}
	return((struct bench_pair_s){
		.fill = fill
	});
}

static inline
struct bench_pair_s bench_wavefront(
	struct params p)
{
	/** init context */
	bench_t fill;
	bench_init(fill);

	struct wavefront_work_s *wwork = wavefront_init_work(0.7);
	for(int64_t i = 0; i < p.cnt; i++) {
		bench_start(fill);
		wavefront(wwork, 
			(uint8_t const *)kv_at(p.seq, i * 2), kv_at(p.len, i * 2),
			(uint8_t const *)kv_at(p.seq, i * 2 + 1), kv_at(p.len, i * 2 + 1));
		bench_end(fill);
	}
	wavefront_clean_work(wwork);
	return((struct bench_pair_s){
		.fill = fill
	});
}

static inline
void print_result(
	int table,
	struct bench_pair_s p)
{
	if(table == 0) {
		printf("%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
			bench_get(p.fill),
			bench_get(p.trace),
			bench_get(p.conv),
			bench_get(p.fill) + bench_get(p.trace) + bench_get(p.conv),
			p.score,
			p.fail);
	} else {
		printf("%ld\t%ld\t%ld\t%ld\t",
			bench_get(p.fill),
			bench_get(p.trace),
			bench_get(p.conv),
			bench_get(p.fill) + bench_get(p.trace) + bench_get(p.conv));
	}
	return;
}

/**
 * @fn main
 */
int main(int argc, char *argv[])
{
	struct params p = { .max = 33000, .frac = 1.0, .mlen = 100 };

	/** parse args */
	int i;
	while((i = getopt(argc, argv, "i:l:c:f:a:th")) != -1) {
		if(parse_args(&p, i, optarg) != 0) { exit(1); }
	}

	kv_init(p.buf);
	kv_init(p.seq);
	kv_init(p.len);

	for(i = 0; i < 64; i++) {
		kv_push(p.buf, '\0');
	}

	#define _finish(_base) { \
		uint64_t l = (uint64_t)MIN2(p.max, (p.frac * (kv_size(p.buf) - base))); \
		for(i = 0; i < p.mlen; i++) { \
			kv_push(p.buf, random_base()); \
		} \
		if(p.cnt++ % 2) { \
			p.blen = MAX2(p.blen, l + p.mlen); \
		} else { \
			p.alen = MAX2(p.alen, l + p.mlen); \
		} \
		kv_push(p.len, l + p.mlen); \
		kv_push(p.seq, (char *)base); \
		p.buf.n = base + l + p.mlen; \
		kv_push(p.buf, '\0'); \
	}

	FILE *fp = strcmp(p.file, "-") == 0 ? stdin : fopen(p.file, "r");
	int c;
	uint64_t base = kv_size(p.buf);
	while((c = getc(fp)) != EOF) {
		if(c == '\n') {
			_finish(base);
			base = kv_size(p.buf);
		} else {
			kv_push(p.buf, c);
		}
	}
	_finish(base);
	p.cnt /= 2;
	for(i = 0; i < 64; i++) {
		kv_push(p.buf, '\0');
	}
	if(fp != stdin) {
		fclose(fp);
	}

	for(i = 0; i < kv_size(p.seq); i++) {
		kv_at(p.seq, i) += (uint64_t)p.buf.a;
	}
	for(i = 0; i < kv_size(p.buf); i++) {
		if(kv_at(p.buf, i) == '\0') { continue; }
		switch(kv_at(p.buf, i)) {
			default:
			case 'A': case 'a': c = 1; break;
			case 'C': case 'c': c = 2; break;
			case 'G': case 'g': c = 4; break;
			case 'T': case 't': c = 8; break;
		}
		kv_at(p.buf, i) = c;
	}

	/* our implementations prefer 4-bit encoding for the default, but 2-bit is also allowed when compilation flag is changed, calculation speeds do not depend on the input encodings */
	print_result(p.table, bench_adaptive_editdist(p));
	// print_result(p.table, bench_ddiag_linear(p));
	print_result(p.table, bench_ddiag_affine(p));
	// print_result(p.table, bench_diff_linear(p));
	print_result(p.table, bench_diff_affine(p));
	// print_result(p.table, bench_gaba_linear(p));
	print_result(p.table, bench_gaba_affine(p));
	print_result(p.table, bench_edlib(p));		/* edlib allows any encoding since it transforms input sequences to internal representations */

	print_result(p.table, bench_seqan(p));
	print_result(p.table, bench_wavefront(p));

	/* convert to 2bit since score profile calculation overhead will be minimized with 2-bit encoding for the bwamem (ksw.c) implementation */
	for(i = 0; i < kv_size(p.buf); i++) {
		if(kv_at(p.buf, i) == '\0') { continue; }
		uint8_t const trans[16] = { [1] = 0, [2] = 1, [4] = 2, [8] = 3 };
		kv_at(p.buf, i) = trans[kv_at(p.buf, i)];
	}

	print_result(p.table, bench_bwamem(p));

	if(p.table != 0) {
		printf("\n");
	}

	kv_destroy(p.buf);
	kv_destroy(p.seq);
	kv_destroy(p.len);
	return 0;
}

/**
 * end of bench.c
 */
