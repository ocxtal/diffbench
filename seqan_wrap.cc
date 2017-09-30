
#include <string>           // std::string
#include <seqan/align.h>    // base header of SeqAn::align
#include <stdint.h>

using namespace seqan;

extern "C" {

	int seqan_editdist(uint8_t const *a, uint8_t const *b)
	{
		typedef String<Dna, CStyle> TSequence;
		TSequence seq1(a);
		TSequence seq2(b);
		int score = globalAlignmentScore(seq1, seq2, MyersBitVector());
		return(score);
	}
};

