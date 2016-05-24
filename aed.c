
/**
 * @file aed.c
 *
 * @brief adaptive banded edit-distance alogrithm
 */
#define UNITTEST_UNIQUE_ID		111
#include "unittest.h"
#include "bench.h"
#include "log.h"
#include "aed.h"
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

/* constants */
#define BW				( (uint64_t)64 )
#define ONES			( (uint64_t)0xffffffffffffffff )


struct aed_vec_s {
	uint64_t ph, mh, pv, mv;
};

struct aed_block_s {
	struct aed_vec_s vec[64];
	uint64_t dir;
	int64_t rem;
};

static inline
struct aed_fill_s aed_fill_search_minpos(
	int64_t arem,
	int64_t brem,
	int64_t score,
	struct aed_block_s *blk)
{

	debug("enter search_minpos arem(%lld), brem(%lld), score(%lld), blk(%p)", arem, brem, score, blk);

	uint64_t aidx = BW - arem;
	uint64_t bidx = BW - brem;

	int64_t j = blk->rem;
	uint64_t dir = blk->dir<<(64 - j);

	int64_t const score_max = INT32_MAX;
	struct aed_fill_s hmin = { score_max }, vmin = { score_max };

	int64_t hscore = score_max, vscore = score_max;
	for(int64_t i = 0; i < 128; i++) {
		debug("i(%lld), j(%lld), aidx(%lld), bidx(%lld), hscore(%lld), vscore(%lld), dir(%llx)",
			i, j, (int64_t)aidx, (int64_t)bidx, hscore, vscore, dir);

		if(j-- == 0) {
			dir = (--blk)->dir;
			j = 63;

			debug("reload blk(%p), dir(%llx)", blk, dir);
		}

		if(aidx + bidx == BW) {
			debug("aidx + bidx == 64, aidx(%llu), bidx(%llu), hscore(%lld), vscore(%lld)", aidx, bidx, hscore, vscore);
			int64_t diff = vscore - hscore;
			hscore += diff;
			hmin.score += diff;
		}

		if(dir & 0x8000000000000000) {
			/* go up */
			if(aidx < BW) {
				int64_t pd = 0x01 & (blk->vec[j].pv>>aidx);
				int64_t md = 0x01 & (blk->vec[j].mv>>aidx);
				vscore -= (pd - md);

				debug("go up pv(%llx), mv(%llx), pd(%lld), md(%lld)", blk->vec[j].pv, blk->vec[j].mv, pd, md);

				int64_t pc = 0x01 & (blk->vec[j].pv>>32);
				int64_t mc = 0x01 & (blk->vec[j].mv>>32);
				score -= (pc - mc);

				debug("pd(%lld), mc(%lld)", pc, mc);

				if(aidx + bidx < 64 && vscore < vmin.score) {
					vmin = (struct aed_fill_s){ vscore, blk - (j == 0), 0x3f & (j - 1), aidx };
				}
			}
			bidx--;

		} else {
			/* go left */
			if(bidx < BW) {
				int64_t pd = 0x01 & (blk->vec[j].ph>>(63 - bidx));
				int64_t md = 0x01 & (blk->vec[j].mh>>(63 - bidx));
				hscore -= (pd - md);

				debug("go left ph(%llx), mh(%llx), pd(%lld), md(%lld)", blk->vec[j].ph, blk->vec[j].mh, pd, md);

				int64_t pc = 0x01 & (blk->vec[j].ph>>32);
				int64_t mc = 0x01 & (blk->vec[j].mh>>32);
				score -= (pc - mc);

				debug("pd(%lld), mc(%lld)", pc, mc);

				if(aidx + bidx < 64 && hscore < hmin.score) {
					hmin = (struct aed_fill_s){ hscore, blk - (j == 0), 0x3f & (j - 1), 63 - bidx };
				}
			}
			aidx--;

			if(aidx == 32) {
				int64_t diff = score - hscore;
				hscore += diff;
				vscore += diff;
				hmin.score += diff;
				vmin.score += diff;
			}

		}

		dir <<= 1;
	}

	debug("finish, score(%lld)", score);
	debug("hmin, score(%lld), idx(%llu)", hmin.score, hmin.idx);
	debug("vmin, score(%lld), idx(%llu)", vmin.score, vmin.idx);

	return((hmin.score < vmin.score) ? hmin : vmin);
}

struct aed_fill_s aed_fill(
	void *base,
	uint8_t const *a,
	int64_t alen,
	uint8_t const *b,
	int64_t blen)
{
	debug("a(%p), alen(%lld), b(%p), blen(%lld)", a, alen, b, blen);

	/* init sequence pointers */
	uint8_t const *alim = a + alen;
	uint8_t const *blim = b + blen;
	uint8_t const *ablim = alim, *bblim = blim;
	uint8_t const *abase = a;
	uint8_t const *bbase = b;

	uint8_t const atail[BW] = { 0 };
	uint8_t const btail[BW] = {
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
	};
	uint8_t const *atlim = atail + BW, *btlim = btail + BW;

	/* filled lengths */
	int64_t aflen = 0, bflen = 0;

	/* init sequence vector */
	uint64_t register al = 0, ah = 0, bl = ONES, bh = ONES;
	for(int64_t i = 0; i < BW/2; i++) {
		al = (al<<1) | (uint64_t)(*a & 0x01);
		ah = (ah<<1) | ((uint64_t)(*a & 0x02)>>1);
		a++;
	}
	for(int64_t i = 0; i < BW/2 - 1; i++) {
		bl = (bl>>1) | ((uint64_t)(*b & 0x01)<<(BW - 1));
		bh = (bh>>1) | ((uint64_t)(*b & 0x02)<<(BW - 2));
		b++;
	}

	/* init diff vector */
	uint64_t register ph = ONES>>(BW/2);
	uint64_t register mh = ONES<<(BW/2);
	uint64_t register pv = ONES<<(BW/2 + 1);
	uint64_t register mv = ONES>>(BW/2 - 1);

	/* score accumulator */
	int64_t register score = 0;
	int64_t register acc = -2;

	/* direction vector */
	uint64_t register dir = 0;
	int64_t register cnt = 0;

	#define _update_vectors() { \
		uint64_t eq = ~((al ^ bl) | (ah ^ bh)); \
		uint64_t xh = eq | mh; \
		uint64_t xv = eq | mv; \
		uint64_t _ph = mv | ~(xh | pv); \
		uint64_t _pv = mh | ~(xv | ph); \
		uint64_t _mh = pv & xh; \
		uint64_t _mv = ph & xv; \
		ph = _ph; mh = _mh; pv = _pv; mv = _mv; \
		debug("ph(%llx), mh(%llx), pv(%llx), mv(%llx)", ph, mh, pv, mv); \
	}

	#define _update_acc(_p, _m) { \
		acc += (_p)>>(BW - 1); \
		acc -= (_p) & 0x01; \
		acc -= (_m)>>(BW - 1); \
		acc += (_m) & 0x01; \
	}

	#define _update_score(_p, _m) { \
		debug("update score, p(%llx), m(%llx), score(%lld)", _p, _m, score); \
		score += ((_p)>>(BW/2)) & 0x01; \
		score -= ((_m)>>(BW/2)) & 0x01; \
	}

	struct aed_block_s *blk = (struct aed_block_s *)base;
	while(1) {
		for(int64_t i = 0; i < 64; i++) {

			if(a >= alim || b >= blim) {
				if(a != atlim && b != btlim) {
					debug("a(%p), alim(%p), b(%p), blim(%p)", a, alim, b, blim);

					aflen += a - abase; abase = a;
					bflen += b - bbase; bbase = b;

					debug("update tail, aflen(%lld), bflen(%lld)", aflen, bflen);
					if(a == ablim) { a = abase = atail; alim = atail + BW; }
					if(b == bblim) { b = bbase = btail; blim = btail + BW; }

				} else {
					debug("dir(%llx), i(%lld), comp_dir(%llx)", dir, i, dir>>(64 - i));
					blk->dir = dir>>(64 - i);
					blk->rem = i;
					return(aed_fill_search_minpos(
						atlim - a, btlim - b, score, blk));
				}
			}

			if(acc < 0) {
				debug("go down, acc(%lld), score(%lld)", acc, score);
				dir = (dir>>1) | ((uint64_t)0x01<<63);

				debug("fetch next base on b(%p)", b);
				bl = (bl>>1) | ((uint64_t)(*b & 0x01)<<(BW - 1));
				bh = (bh>>1) | ((uint64_t)(*b & 0x02)<<(BW - 2));
				b++;

				debug("shift vectors right");
				pv = (pv>>1) | 0x8000000000000000;
				mv = mv>>1;

				debug("update vectors");
				_update_vectors();
				_update_acc(pv, mv);
				_update_score(pv, mv);

			} else {
				debug("go right, acc(%lld), score(%lld)", acc, score);
				dir >>= 1;

				debug("fetch next base on a(%p)", a);
				al = (al<<1) | (uint64_t)(*a & 0x01);
				ah = (ah<<1) | ((uint64_t)(*a & 0x02)>>1);
				a++;

				debug("shift vectors left");
				ph = (ph<<1) | 0x01;
				mh = mh<<1;

				debug("update vectors");
				_update_vectors();
				_update_acc(ph, mh);
				_update_score(ph, mh);

			}

			blk->vec[i] = (struct aed_vec_s){ ph, mh, pv, mv };
		}

		debug("save dir(%llx), blk(%p)", dir, blk);
		blk->dir = dir;
		blk->rem = 0;

		dir = 0;
		blk++;
	}

	/* never reaches here */
	return((struct aed_fill_s){ 0 });
}

int64_t aed_trace(
	char *buf,
	int64_t size,
	void *base,
	struct aed_fill_s f)
{
	struct aed_block_s *blk = f.blk;
	int64_t rem = f.rem;
	uint64_t idx = f.idx;
	uint64_t dir = blk->dir<<(63 - rem);

	char *p = buf + size - 1;
	char type = '\0';
	int64_t len = 0;

	#define _aed_trace_pop(_blk, _rem, _dir) { \
		(_dir) <<= 1; \
		if((_rem)-- == 0) { \
			(_dir) = (--(_blk))->dir; \
			(_rem) = 63; \
		} \
	}

	#define _aed_trace_push(_p, _type, _len, _c) { \
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

	while(blk >= (struct aed_block_s *)base) {
		int64_t pv = 0x01 & (blk->vec[rem].pv>>idx);
		int64_t ph = 0x01 & (blk->vec[rem].ph>>idx);
		int64_t d = (dir & 0x8000000000000000) ? 0x01 : 0;

		debug("ph(%llx), mh(%llx), pv(%llx), mv(%llx)", blk->vec[rem].ph, blk->vec[rem].mh, blk->vec[rem].pv, blk->vec[rem].mv);
		debug("blk(%p), rem(%lld), idx(%llu), dir(%llx)", blk, rem, idx, dir);
		debug("pv(%lld), ph(%lld), d(%lld)", pv, ph, d);

		if(pv == 1) {
			debug("go up, ins, type(%c)", type);
			idx -= (1 - d);
			_aed_trace_push(p, type, len, 'I');
		} else if(ph == 1) {
			debug("go left, del, type(%c)", type);
			idx += d;
			_aed_trace_push(p, type, len, 'D');
		} else {
			debug("go diag, match, type(%c)", type);
			_aed_trace_pop(blk, rem, dir);
			int64_t d2 = (dir & 0x8000000000000000) ? 0x01 : 0;
			idx += (d + d2 - 1);
			_aed_trace_push(p, type, len, 'M');
		}
		_aed_trace_pop(blk, rem, dir);
	}

	_aed_trace_push(p, type, len, ' ');

	debug("%s", p + 1);

	char *b = buf;
	int64_t l = 0;
	while(p < buf + size - 1) {
		*b++ = *++p;
		l++;
	}
	return(l);
}

#if 0
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

static inline
int print_cigar(char *cigar)
{
	char *p = cigar;
	while(*p) {
		char *q = p;
		while(isdigit(*q)) { q++; }
		switch(*q) {
			case 'M': printf("%s", "\x1b[33m"); break;
			case 'I': printf("%s", "\x1b[32m"); break;
			case 'D': printf("%s", "\x1b[35m"); break;
		}

		for(; p <= q; p++) {
			putchar(*p);
		}
		printf("%s", "\x1b[39m");
	}
	printf("\n");
	return(0);
}

unittest()
{
	uint8_t const a[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	uint8_t const b[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
 	struct aed_vec_s *ptr = (struct aed_vec_s *)
		malloc(2 * sizeof(struct aed_vec_s) * (sizeof(a) + sizeof(b) + 65));
	char *buf = (char *)malloc(sizeof(a) + sizeof(b));

	struct aed_fill_s f = aed_fill(ptr, a, sizeof(a), b, sizeof(b));
	int64_t len = aed_trace(buf, sizeof(a) + sizeof(b), ptr, f);
	printf("%s\n%lld\n", buf, len);

	free(ptr);
	free(buf);
}

unittest()
{
	int64_t const len = 200;
	int64_t const cnt = 1;

	srand(time(NULL));

	char *a = generate_random_sequence(len);
	char *b = generate_mutated_sequence(a, len, 0.1, 0.1, 1024);
	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	debug("%s\n%s", a, b);

	encode(a, alen);
	encode(b, blen);


 	struct aed_vec_s *ptr = (struct aed_vec_s *)
		malloc(2 * sizeof(struct aed_vec_s) * (alen + blen + 65));
	char *buf = (char *)malloc(alen + blen);


	bench_t fill, trace;
	bench_init(fill);
	bench_init(trace);


	for(int64_t i = 0; i < cnt; i++) {
		bench_start(fill);
		struct aed_fill_s f = aed_fill(ptr, (uint8_t const *)a, alen, (uint8_t const *)b, blen);
		bench_end(fill);

		bench_start(trace);
		int64_t len = aed_trace(buf, alen + blen, ptr, f);
		bench_end(trace);

		print_cigar(buf);
	}


	log("%lld\t%lld\t%lld",
		bench_get(fill),
		bench_get(trace),
		bench_get(fill) + bench_get(trace));

	free(a);
	free(b);
	free(ptr);
	free(buf);
}
#endif


/**
 * end of aed.c
 */
