
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

#define BIT_WIDTH 			8
#define BAND_WIDTH 			32

#include "ddiag.h"

/**
 * @fn print_usage
 */
void print_usage(void)
{
	fprintf(stderr, "usage: bench -l <len> -c <cnt> -x <mismatch rate> -d <indel rate>\n");
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
uint8_t encode_base(
	char c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases { A = 0x00, C = 0x01, G = 0x02, T = 0x03 };
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
	struct params p)
{
	char *a = generate_random_sequence(p.len);
	char *b = generate_mutated_sequence(a, p.len, p.x, p.d, 1024);
	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	encode(a, alen);
	encode(b, blen);

	void *ptr = malloc(2 * 4 * sizeof(uint64_t) * (alen + blen + 65));
	char *buf = (char *)malloc(alen + blen);


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

	free(a);
	free(b);
	free(ptr);
	free(buf);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_ddiag_linear(
	struct params p)
{
	char *a = generate_random_sequence(p.len);
	char *b = generate_mutated_sequence(a, p.len, p.x, p.d, 1024);
	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	encode(a, alen);
	encode(b, blen);

	struct sea_params param = { 0, 2, -3, -5, -1, 0, 0, 0, 32 };


	struct sea_result *aln = (struct sea_result *)malloc(sizeof(char) * (alen + blen + 32 + 1) + sizeof(struct sea_result));
	aln->aln = (void *)((struct sea_result *)aln + 1);
	((char *)(aln->aln))[0] = 0;
	aln->score = 0;
	aln->ctx = NULL;


	int64_t alnsize = diag_linear_dynamic_banded_matsize(alen, blen, 32);
	void *mat = malloc(alnsize);


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
		struct mpos o = diag_linear_dynamic_banded_fill(
			aln, param,
			(char *)mat + 32 * sizeof(int16_t));
		o = diag_linear_dynamic_banded_search(
			aln, param,
			(char *)mat + 32 * sizeof(int16_t),
			o);
		score += aln->score;
		bench_end(fill);
		

		bench_start(trace);
		diag_linear_dynamic_banded_trace(
			aln, param,
			(char *)mat + 32 * sizeof(int16_t),
			o);
		bench_end(trace);
	}

	/**
	 * clean malloc'd memories
	 */
	free(a);
	free(b);
	free(aln);
	free(mat);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
struct bench_pair_s bench_gaba_linear(
	struct params p)
{
	char *a = generate_random_sequence(p.len);
	char *b = generate_mutated_sequence(a, p.len, p.x, p.d, 1024);
	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	encode(a, alen);
	encode(b, blen);

	char *c = (char *)malloc(p.len);

	/** init context */
	gaba_t *ctx = gaba_init(GABA_PARAMS(
		.xdrop = 100,
		.score_matrix = GABA_SCORE_SIMPLE(2, 3, 5, 1)));
	struct gaba_section_s asec = gaba_build_section(0, (uint8_t const *)a, alen);
	struct gaba_section_s bsec = gaba_build_section(2, (uint8_t const *)b, blen);

	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);

	void const *lim = (void const *)0x800000000000;

	int64_t score = 0;
	for(int64_t i = 0; i < p.cnt; i++) {
		// gaba_dp_flush(dp, lim, lim);
		gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);

		bench_start(fill);
		struct gaba_fill_s *f = gaba_dp_fill_root(dp, &asec, 0, &bsec, 0);
		score += f->max;
		bench_end(fill);
		
		bench_start(trace);
		struct gaba_result_s *r = gaba_dp_trace(dp, f, NULL, NULL);
		gaba_dp_dump_cigar(c, p.len, r->path->array, r->path->offset, r->path->len);
		bench_end(trace);
		gaba_dp_clean(dp);
	}
	gaba_clean(ctx);

	/**
	 * clean malloc'd memories
	 */
	free(a);
	free(b);
	free(c);

	return((struct bench_pair_s){
		.fill = fill,
		.trace = trace,
		.score = score
	});
}

static inline
void print_result(
	struct bench_pair_s p)
{
	printf("%lld\t%lld\t%lld\t%lld\n",
		bench_get(p.fill),
		bench_get(p.trace),
		bench_get(p.fill) + bench_get(p.trace),
		p.score);
	return;
}

/**
 * @fn main
 */
int main(int argc, char *argv[])
{
	int64_t i;
	char *a;
	char *b;
	char *c;
	struct params p;
	bench_t fill, trace, parse;

	/** set defaults */
	p.len = 10000;
	p.cnt = 10000;
	p.x = 0.1;
	p.d = 0.1;
	p.pa = p.pb = NULL;

	/** parse args */
	while((i = getopt(argc, argv, "q:t:o:l:x:d:c:a:seb:h")) != -1) {
		if(parse_args(&p, i, optarg) != 0) { exit(1); }
	}

	fprintf(stderr, "len\t%lld\ncnt\t%lld\nx\t%f\nd\t%f\n", p.len, p.cnt, p.x, p.d);


	print_result(bench_adaptive_editdist(p));
	print_result(bench_ddiag_linear(p));
	print_result(bench_gaba_linear(p));

	return 0;
}

/**
 * end of bench.c
 */
