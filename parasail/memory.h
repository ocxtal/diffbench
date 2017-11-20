/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_INTERNAL_H_
#define _PARASAIL_INTERNAL_H_

#include <stdint.h>

#include "parasail.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void * parasail_memalign(size_t alignment, size_t size);
extern int * parasail_memalign_int(size_t alignment, size_t size);
extern int8_t * parasail_memalign_int8_t(size_t alignment, size_t size);
extern int16_t * parasail_memalign_int16_t(size_t alignment, size_t size);
extern int32_t * parasail_memalign_int32_t(size_t alignment, size_t size);
extern int64_t * parasail_memalign_int64_t(size_t alignment, size_t size);

extern void parasail_free(void *ptr);

extern void parasail_memset(void *b, int c, size_t len);
extern void parasail_memset_int(int *b, int c, size_t len);
extern void parasail_memset_int8_t(int8_t *b, int8_t c, size_t len);
extern void parasail_memset_int16_t(int16_t *b, int16_t c, size_t len);
extern void parasail_memset_int32_t(int32_t *b, int32_t c, size_t len);
extern void parasail_memset_int64_t(int64_t *b, int64_t c, size_t len);

extern parasail_result_t* parasail_result_new();
extern parasail_result_t* parasail_result_new_table1(const int a, const int b);
extern parasail_result_t* parasail_result_new_table3(const int a, const int b);
extern parasail_result_t* parasail_result_new_rowcol1(const int a, const int b);
extern parasail_result_t* parasail_result_new_rowcol3(const int a, const int b);

extern parasail_profile_t* parasail_profile_new(
        const char * s1, const int s1Len, const parasail_matrix_t *matrix);

extern char* parasail_reverse(const char *s, int end);

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_INTERNAL_H_ */