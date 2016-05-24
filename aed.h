
/**
 * @file aed.h
 */
#ifndef _AED_H_INCLUDED
#define _AED_H_INCLUDED

struct aed_fill_s {
	int64_t score;
	struct aed_block_s *blk;
	int64_t rem;
	uint64_t idx;
};

struct aed_fill_s aed_fill(
	void *base,
	uint8_t const *a,
	int64_t alen,
	uint8_t const *b,
	int64_t blen);

int64_t aed_trace(
	char *buf,
	int64_t size,
	void *base,
	struct aed_fill_s f);

#endif
/**
 * end of aed.h
 */
