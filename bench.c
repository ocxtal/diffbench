
/**
 * @file bench.c
 * @brief speed benchmark of libsea
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "bench.h"
#include "aed.h"
#include "gaba/gaba.h"
#include "edlib.h"

#define BIT_WIDTH 			8
#define BAND_WIDTH 			32

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

struct params {
	int64_t len;
	int64_t cnt;
	double x;
	double d;
	char **pa;
	char **pb;
	int table;
};

int parse_args(struct params *p, int c, char *arg)
{
	switch(c) {
		/**
		 * sequence generation parameters.
		 */
		case 'l': p->len = atoi((char *)arg); return 0;
		case 'x': p->x = atof((char *)arg); return 0;
		case 'd': p->d = atof((char *)arg); return 0;
		/**
		 * benchmarking options
		 */
		case 'c': p->cnt = atoi((char *)arg); return 0;
		case 'a': printf("%s\n", arg); return 0;
		case 't': p->table = 1; return 0;
		/**
		 * the others: print help message
		 */
		case 'h':
		default: print_usage(); return -1;
	}
}

struct bench_pair_s {
	bench_t fill, trace;
	int64_t score;
};

static inline
struct bench_pair_s bench_adaptive_editdist(
	struct params p,
	char const *a,
	int64_t alen,
	char const *b,
	int64_t blen)
{
	void *base = aligned_malloc(2 * 4 * sizeof(uint64_t) * (alen + blen + 65));
	memset(base, 0, 4 * sizeof(uint64_t));
	void *ptr = base + 4 * sizeof(uint64_t);

	char *buf = (char *)aligned_malloc(alen + blen);


	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		bench_start(fill);
		struct aed_fill_s f = aed_fill(ptr, (uint8_t const *)a, alen, (uint8_t const *)b, blen);
		score += f.score;
		bench_end(fill);

		bench_start(trace);
		int64_t len = aed_trace(buf, alen + blen, ptr, f);
		bench_end(trace);
	}
	free(base);
	free(buf);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_ddiag_linear(
	struct params p,
	char const *a,
	int64_t alen,
	char const *b,
	int64_t blen)
{
	struct sea_params param = { 0, 1, -2, -1, 0, 0, 0, 50, 32 };


	struct sea_result *aln = (struct sea_result *)aligned_malloc(sizeof(char) * (alen + blen + 32 + 1) + sizeof(struct sea_result));
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diag_linear_dynamic_banded_matsize(alen, blen, 32);
	void *base = aligned_malloc(alnsize);
	memset(base, 0, 32 * 2 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 2 * sizeof(int16_t);


	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		aln->score = 0;
		aln->a = a; aln->apos = 0; aln->alen = alen;
		aln->b = b; aln->bpos = 0; aln->blen = blen;
		aln->len = alnsize;


		bench_start(fill);
		struct mpos o = diag_linear_dynamic_banded_fill(aln, param, mat);
		o = diag_linear_dynamic_banded_search(aln, param, mat, o);
		score += aln->score;
		bench_end(fill);


		bench_start(trace);
		diag_linear_dynamic_banded_trace(aln, param, mat, o);
		bench_end(trace);
	}
	free(aln);
	free(base);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_ddiag_affine(
	struct params p,
	char const *a,
	int64_t alen,
	char const *b,
	int64_t blen)
{
	struct sea_params param = { 0, 1, -2, -2, -1, 0, 0, 50, 32 };


	struct sea_result *aln = (struct sea_result *)aligned_malloc(sizeof(char) * (alen + blen + 32 + 1) + sizeof(struct sea_result));
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diag_affine_dynamic_banded_matsize(alen, blen, 32);
	void *base = aligned_malloc(alnsize);
	memset(base, 0, 32 * 6 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 6 * sizeof(int16_t);


	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		aln->score = 0;
		aln->a = a; aln->apos = 0; aln->alen = alen;
		aln->b = b; aln->bpos = 0; aln->blen = blen;
		aln->len = alnsize;


		bench_start(fill);
		struct mpos o = diag_affine_dynamic_banded_fill(aln, param, mat);
		o = diag_affine_dynamic_banded_search(aln, param, mat, o);
		score += aln->score;
		bench_end(fill);


		bench_start(trace);
		diag_affine_dynamic_banded_trace(aln, param, mat, o);
		bench_end(trace);
	}
	free(aln);
	free(base);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_diff_linear(
	struct params p,
	char const *a,
	int64_t alen,
	char const *b,
	int64_t blen)
{
	struct sea_params param = { 0, 1, -2, -1, 0, 0, 0, 50, 32 };


	struct sea_result *aln = (struct sea_result *)aligned_malloc(sizeof(char) * (alen + blen + 32 + 1) + sizeof(struct sea_result) + 1000000);
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diff_linear_dynamic_banded_matsize(alen, blen, 32);
	void *base = aligned_malloc(alnsize);
	memset(base, 0, 32 * 2 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 2 * sizeof(int16_t);


	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		aln->score = 0;
		aln->a = a; aln->apos = 0; aln->alen = alen;
		aln->b = b; aln->bpos = 0; aln->blen = blen;
		aln->len = alnsize;


		bench_start(fill);
		struct mpos o = diff_linear_dynamic_banded_fill(aln, param, mat);
		bench_end(fill);


		bench_start(trace);
		o = diff_linear_dynamic_banded_search(aln, param, mat, o);
		score += aln->score;
		diff_linear_dynamic_banded_trace(aln, param, mat, o);
		bench_end(trace);
	}
	free(aln);
	free(base);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_diff_affine(
	struct params p,
	char const *a,
	int64_t alen,
	char const *b,
	int64_t blen)
{
	struct sea_params param = { 0, 1, -2, -2, -1, 0, 0, 50, 32 };


	struct sea_result *aln = (struct sea_result *)aligned_malloc(sizeof(char) * (alen + blen + 32 + 1) + sizeof(struct sea_result));
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diff_affine_dynamic_banded_matsize(alen, blen, 32);
	void *base = aligned_malloc(alnsize);
	memset(base, 0, 32 * 6 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 6 * sizeof(int16_t);


	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		aln->score = 0;
		aln->a = a; aln->apos = 0; aln->alen = alen;
		aln->b = b; aln->bpos = 0; aln->blen = blen;
		aln->len = alnsize;


		bench_start(fill);
		struct mpos o = diff_affine_dynamic_banded_fill(aln, param, mat);
		bench_end(fill);


		bench_start(trace);
		o = diff_affine_dynamic_banded_search(aln, param, mat, o);
		score += aln->score;
		diff_affine_dynamic_banded_trace(aln, param, mat, o);
		bench_end(trace);
	}
	free(aln);
	free(base);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_gaba_linear(
	struct params p,
	char const *a,
	int64_t alen,
	char const *b,
	int64_t blen)
{
	char *c = (char *)malloc(p.len);

	/** init context */
	gaba_t *ctx = gaba_init(GABA_PARAMS(
		.xdrop = 50,
		GABA_SCORE_SIMPLE(1, 2, 0, 1)));
	struct gaba_section_s asec = gaba_build_section(0, (uint8_t const *)a, alen);
	struct gaba_section_s bsec = gaba_build_section(2, (uint8_t const *)b, blen);

	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	void const *lim = (void const *)0x800000000000;

	gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);
	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		gaba_dp_flush(dp, lim, lim);

		bench_start(fill);
		struct gaba_fill_s *f = gaba_dp_fill_root(dp, &asec, 0, &bsec, 0);
		score += f->max;
		bench_end(fill);
		
		bench_start(trace);
		struct gaba_alignment_s *r = gaba_dp_trace(dp, f, NULL, NULL);
		gaba_dp_dump_cigar_forward(c, p.len, r->path->array, 0, r->path->len);
		bench_end(trace);
	}
	gaba_dp_clean(dp);
	gaba_clean(ctx);
	free(c);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_gaba_affine(
	struct params p,
	char const *a,
	int64_t alen,
	char const *b,
	int64_t blen)
{
	char *c = (char *)malloc(p.len);

	/** init context */
	gaba_t *ctx = gaba_init(GABA_PARAMS(
		.xdrop = 50,
		GABA_SCORE_SIMPLE(1, 2, 1, 1)));
	struct gaba_section_s asec = gaba_build_section(0, (uint8_t const *)a, alen);
	struct gaba_section_s bsec = gaba_build_section(2, (uint8_t const *)b, blen);

	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	void const *lim = (void const *)0x800000000000;

	gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);
	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		gaba_dp_flush(dp, lim, lim);

		bench_start(fill);
		struct gaba_fill_s *f = gaba_dp_fill_root(dp, &asec, 0, &bsec, 0);
		score += f->max;
		bench_end(fill);
		
		bench_start(trace);
		struct gaba_alignment_s *r = gaba_dp_trace(dp, f, NULL, NULL);
		gaba_dp_dump_cigar_forward(c, p.len, r->path->array, 0, r->path->len);
		bench_end(trace);
	}
	gaba_dp_clean(dp);
	gaba_clean(ctx);
	free(c);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_edlib(
	struct params p,
	char const *a,
	int64_t alen,
	char const *b,
	int64_t blen)
{
	/** init context */
	EdlibAlignConfig cf = (EdlibAlignConfig){ .k = -1, .mode = EDLIB_MODE_SHW, .task = EDLIB_TASK_DISTANCE };
	EdlibAlignConfig ct = (EdlibAlignConfig){ .k = -1, .mode = EDLIB_MODE_SHW, .task = EDLIB_TASK_PATH };

	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		bench_start(fill);
		EdlibAlignResult f = edlibAlign(a, alen, b, blen, cf);
		score += f.editDistance;
		bench_end(fill);

		edlibFreeAlignResult(f);

		bench_start(trace);
		EdlibAlignResult t = edlibAlign(a, alen, b, blen, ct);
		score += t.editDistance;
		bench_end(trace);
		
		edlibFreeAlignResult(t);
	}

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}


static inline
void print_result(
	int table,
	struct bench_pair_s p)
{
	if(table == 0) {
		printf("%lld\t%lld\t%lld\t%lld\n",
			bench_get(p.fill),
			bench_get(p.trace),
			bench_get(p.fill) + bench_get(p.trace),
			p.score);
	} else {
		printf("%lld\t%lld\t%lld\t%lld\t",
			bench_get(p.fill),
			bench_get(p.trace),
			bench_get(p.fill) + bench_get(p.trace),
			p.score);
	}
	return;
}

/**
 * @fn main
 */
int main(int argc, char *argv[])
{
	struct params p = (struct params){
		.len = 10000,
		.cnt = 10000,
		.x = 0.1,
		.d = 0.1,
		.pa = p.pb = NULL
	};

	/** parse args */
	int i;
	while((i = getopt(argc, argv, "l:x:d:c:a:th")) != -1) {
		if(parse_args(&p, i, optarg) != 0) { exit(1); }
	}

	if(p.table == 0) {
		printf("len\t%lld\ncnt\t%lld\nx\t%f\nd\t%f\n", p.len, p.cnt, p.x, p.d);
	} else {
		printf("%lld\t", p.len);
	}


	char *a = generate_random_sequence(p.len);
	char *b = generate_mutated_sequence(a, p.len, p.x, p.d, 1024);

	a = add_tail(a, 0, 500);
	b = add_tail(b, 0, 500);

	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	encode(a, alen);
	encode(b, blen);

	a = add_margin((uint8_t *)a, alen, 32, 32);
	b = add_margin((uint8_t *)b, blen, 32, 32);

	print_result(p.table, bench_adaptive_editdist(p, a + 32, alen, b + 32, blen));
	if(p.len < 35000) {
		print_result(p.table, bench_ddiag_linear(p, a + 32, alen, b + 32, blen));
		print_result(p.table, bench_ddiag_affine(p, a + 32, alen, b + 32, blen));
	}
	print_result(p.table, bench_diff_linear(p, a + 32, alen, b + 32, blen));
	print_result(p.table, bench_diff_affine(p, a + 32, alen, b + 32, blen));
	print_result(p.table, bench_gaba_linear(p, a + 32, alen, b + 32, blen));
	print_result(p.table, bench_gaba_affine(p, a + 32, alen, b + 32, blen));
	// print_result(bench_edlib(p, a + 32, alen, b + 32, blen));

	if(p.table != 0) {
		printf("\n");
	}

	free(a);
	free(b);
	return 0;
}

/**
 * end of bench.c
 */
