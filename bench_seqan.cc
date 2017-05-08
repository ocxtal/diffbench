
///
// @file bemch_seqan.cc
//
// @detail
// Compile:
//   $ g++ -Wall -O3 -I<path/to/seqan> -o bench_seqan bench_seqan.cc
//
#include <string>           // std::string
#include <seqan/align.h>    // base header of SeqAn::align
#include <sys/time.h>       // for gettimeofday

using namespace seqan;

char random_base(void);
std::string generate_random_sequence(int len);
std::string generate_mutated_sequence(std::string seq, double x, double d, int bw);

int main(void)
{
    long const len[6] = {100, 300, 1000, 3000, 10000, 30000};
    long const cnt = 1000;
    struct timeval tvs, tve;

    for(int j = 0; j < 6; j++) {
        std::string a = generate_random_sequence(len[j]);
        std::string b = generate_mutated_sequence(a, 0.03, 0.1, 1024);

        // SeqAn localAlignment
        typedef String<Dna, CStyle> TSequence;
        TSequence seq1(a);
        TSequence seq2(b);

        gettimeofday(&tvs, NULL);
        for(int i = 0; i < cnt; i++) {
            volatile int score = globalAlignmentScore(seq1, seq2, MyersBitVector());
        }
        gettimeofday(&tve, NULL);

        long us = (tve.tv_sec - tvs.tv_sec) * 1000000 + (tve.tv_usec - tvs.tv_usec);
        std::cout << us << ", " << (double)(len[j] * len[j]) / (double)us * cnt << std::endl;
    }
    return 0;
}

// random sequence generation functoins. see benchmark.c for details.
char random_base(void)
{
    switch(rand() % 4) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'A';
    }
}

std::string generate_random_sequence(int len)
{
    std::string seq;
    seq.reserve(len);
    seq.resize(len);
    for(int i = 0; i < len; i++) {
        seq[i] = random_base();
    }
    return seq;
}

std::string generate_mutated_sequence(std::string seq, double x, double d, int bw)
{
    int const len = seq.length();
    std::string mutated_seq;

    mutated_seq.reserve(len);
    mutated_seq.resize(len);
    for(int i = 0, j = 0, wave = 0; i < len; i++) {
        if(((double)rand() / (double)RAND_MAX) < x) {
            mutated_seq[i] = random_base(); j++;    /** mismatch */
        } else if(((double)rand() / (double)RAND_MAX) < d) {
            if(rand() & 0x01 && wave > -bw+1) {
                mutated_seq[i] = (j < len) ? seq[j++] : random_base();
                j++; wave--;                        /** deletion */
            } else if(wave < bw-2) {
                mutated_seq[i] = random_base();
                wave++;                             /** insertion */
            } else {
                mutated_seq[i] = (j < len) ? seq[j++] : random_base();
            }
        } else {
            mutated_seq[i] = (j < len) ? seq[j++] : random_base();
        }
    }
    return mutated_seq;
}
