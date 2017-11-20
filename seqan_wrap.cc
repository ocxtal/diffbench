
#define NDEBUG		/* disable debugging macros in seqan.h */

#include <string>           // std::string
#include "seqan/align.h"    // base header of SeqAn::align
#include <stdint.h>

using namespace seqan;

extern "C" {

	#include "bench.h"
	extern bench_t seqan_fill;
	int seqan_editdist(uint8_t const *a, uint8_t const *b)
	{
		typedef String<Dna, CStyle> TSequence;
		TSequence seq1(a);
		TSequence seq2(b);
		bench_start(seqan_fill);
		int score = globalAlignmentScore(seq1, seq2, MyersBitVector());
		bench_end(seqan_fill);
		return(score);
	}
};

