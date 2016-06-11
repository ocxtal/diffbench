# Benchmark on difference recurrence relation


This repository contains benchmarking programs for the difference-recurrence DP algorithm, which is adopted in [the GABA nucleotide sequence alignment library](https://github.com/ocxtal/libgaba). For an overview and assessment results of the other algorithm, adaptive banded DP, that is adopted in the GABA library, see [https://github.com/ocxtal/adaptivebandbench](https://github.com/ocxtal/adaptivebandbench).


## Derivation of the difference recurrence

The semi-global variant of the Gotoh's algorithm is shown below. It calculates an alignment path with its left end fixed (aligned) and right end free.

![swg](https://github.com/ocxtal/diffbench/blob/master/fig/swg.png)

The four difference variables, shown below, are introduced to represent difference of horizontally and vertically adjacent scores.

![vals](https://github.com/ocxtal/diffbench/blob/master/fig/diffvals.png)

The original recurrence is transformed into a difference form.

![rec](https://github.com/ocxtal/diffbench/blob/master/fig/diffrec.png)

The four difference variables are bounded by constants which are determined by the scoring parameters.

![bound](https://github.com/ocxtal/diffbench/blob/master/fig/diffbound.png)

Adding gap-penalty offsets and modifying E and F differences make the formula simple, enabling faster calculation on super-scaler processors.

![vals2](https://github.com/ocxtal/diffbench/blob/master/fig/diffvals2.png)

![rec](https://github.com/ocxtal/diffbench/blob/master/fig/diffrec2.png)

The bounding constants are also modified.

![bound2](https://github.com/ocxtal/diffbench/blob/master/fig/diffbound2.png)

This transformation is inspired by the [Loving's bit-parallel global alignment algorithm](http://bioinformatics.oxfordjournals.org/content/30/22/3166.short). In contrast to the Loving's algorithm that keeps difference values in multiple bit vectors, ours keeps difference variables in contiguous cells in a set of SIMD registers. The GABA library implements the last recurrence relation and its linear-gap penalty variant combined with the [adaptive banded heuristic](https://github.com/ocxtal/adaptivebandbench).


## Speed benchmark

The matrix calculation and traceback speed of the linear-gap and affine-gap penalty implementation, found in the GABA library, is measured and compared to the adaptive banded without difference recurrence, the adaptive banded variant of the Myers' bit-parallel edit distance algorithm (see [Myers, 1999](http://dl.acm.org/citation.cfm?id=316550) and [Kimura, 2010](http://www.ncbi.nlm.nih.gov/pubmed/22809415)).


The benchmark program is built with the [waf](https://github.com/waf-project/waf) build framework. The following command line make waf to build the program with the Intel C Compiler (icc).

```
$ ./waf configure CC=icc build
```

The result, shown in the table and figure below, is the result of the AVX2-enabled build of the GABA library, taken on 3.6GHz Intel Haswell processor. The fastest was the linear-gap penalty GABA library, being 5% faster than the adaptive-banded bit-parallel edit distance implementation. The difference algorithm accelerated the non-difference implementation by a factor 2. The affine-gap penalty implementations exhibited a similar trend. The average clock cycles of the fill-in loop were 18 in the linear-gap implementation and 25 in the affine-gap one, respectively. The traceback run in 7 cycle per loop in the both implementation.




### Performance notices

The fill-in loop of the affine-gap implementation in the GABA library uses at most 16 SIMD vectors at the same time when it is compiled with SSE4.1 option. The Intel C Compiler (icc) backend is very smart to keep all the vectors on the 16 xmm registers though the other compilers, gcc and clang, generates one or two register spills and reloads. It is confirmed that this (and some other factors) causes 10~20% performance degregation in the both linear and affine implementations.


## License

Apache v2.

Copyright (c) 2016 Hajime Suzuki
