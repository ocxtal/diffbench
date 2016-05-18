
/**
 * @file aed.c
 *
 * @brief adaptive banded edit-distance alogrithm
 */
#define UNITTEST_UNIQUE_ID		111
#include "unittest.h"
#include "bench.h"
#include "log.h"
#include <stdint.h>

/* constants */
#define BW				( (uint64_t)64 )
#define ONES			( (uint64_t)0xffffffffffffffff )


struct aed_vec_s {
	uint64_t ph, mh, pv, mv;
};

struct aed_fill_s {
	int64_t score;
	int64_t acc;
	struct aed_vec_s const *mat;
};


struct aed_fill_s aed_fill(
	uint8_t const *a,
	int64_t alen,
	uint8_t const *b,
	int64_t blen)
{
	debug("a(%p), alen(%lld), b(%p), blen(%lld)", a, alen, b, blen);

	/* malloc mem */
	struct aed_vec_s *ptr = (struct aed_vec_s *)
		malloc(2 * sizeof(struct aed_vec_s) * (alen + blen + 65));
	struct aed_vec_s *mat = ptr;

	/* init sequence pointers */
	uint8_t const *alim = a + alen;
	uint8_t const *blim = b + blen;
	uint8_t const *ablim = alim, *bblim = blim;
	uint8_t const *abase = a;
	uint8_t const *bbase = b;

	uint8_t const atail[32] = { 0 };
	uint8_t const btail[32] = { 0x03 };
	uint8_t const *atlim = atail + 32, *btlim = btail + 32;

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
	}

	#define _update_acc(_p, _m) { \
		acc += (_p)>>(BW - 1); \
		acc -= (_p) & 0x01; \
		acc -= (_m)>>(BW - 1); \
		acc += (_m) & 0x01; \
	}

	#define _update_score(_p, _m) { \
		score += ((_p)>>(BW/2)) & 0x01; \
		score -= ((_m)>>(BW/2)) & 0x01; \
	}

	while(a != atlim && b != btlim) {
		while(a < alim && b < blim) {
			if(acc < 0) {
				debug("go down, acc(%lld), score(%lld)", acc, score);
				dir = (dir>>1) | ((uint64_t)0x01<<63);

				debug("fetch next base on b(%p)", b);
				bl = (bl>>1) | ((uint64_t)(*b & 0x01)<<(BW - 1));
				bh = (bh>>1) | ((uint64_t)(*b & 0x02)<<(BW - 2));
				b++;

				debug("shift vectors right");
				pv >>= 1;
				mv >>= 1;

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
				ph <<= 1;
				mh <<= 1;

				debug("update vectors");
				_update_vectors();
				_update_acc(ph, mh);
				_update_score(ph, mh);

			}

			*mat++ = (struct aed_vec_s){ ph, mh, pv, mv };
			if(++cnt == 64) {
				*mat++ = (struct aed_vec_s){ dir, cnt, 0, 0 };
				cnt = 0;
			}
		}

		debug("a(%p), alim(%p), b(%p), blim(%p)", a, alim, b, blim);

		aflen += a - abase; abase = a;
		bflen += b - bbase; bbase = b;

		debug("update tail, aflen(%lld), bflen(%lld)", aflen, bflen);
		if(a == ablim) { a = abase = atail; alim = atail + 32; }
		if(b == bblim) { b = bbase = btail; blim = btail + 32; }
	}

	debug("finish, acc(%lld), score(%lld)", acc, score);
	*mat++ = (struct aed_vec_s){ dir, cnt, 0, 0 };
	cnt = 0;

	free(ptr);

	return((struct aed_fill_s){
		.score = score,
		.acc = acc,
		.mat = ptr
	});
}


void aed_trace(
	struct aed_fill_s f)
{
	/* search min q-coordinate */


}


char random_base(void)
{
	// char const table[4] = {'A', 'C', 'G', 'T'};
	char const table[4] = { 0x01, 0x02, 0x04, 0x08 };
	return(table[rand() % 4]);
}

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

void encode(char *ptr, int64_t len)
{
	for(int64_t i = 0; i < len; i++) {
		ptr[i] = encode_base(ptr[i]);
	}
	return;
}

unittest()
{
	uint8_t const a[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	uint8_t const b[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

	aed_fill(a, sizeof(a), b, sizeof(b));
}


unittest()
{
	char *a = generate_random_sequence(50000);
	char *b = generate_mutated_sequence(a, 50000, 0.1, 0.1, 1024);
	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	encode(a, alen);
	encode(b, blen);


	bench_t t;
	bench_init(t);

	bench_start(t);
	for(int64_t i = 0; i < 1000; i++) {
		aed_fill((uint8_t const *)a, alen, (uint8_t const *)b, blen);
	}
	bench_end(t);

	log("%lld us", bench_get(t));

	free(a);
	free(b);
}

/**
 * end of aed.c
 */
