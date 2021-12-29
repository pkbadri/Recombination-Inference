# Recombination-Inference
The C++ program fixeds.c outputs 20 different summary statistics (described in Padhukasahasram et al. 2006 and references there) under models with constant population size and uniform or nonuniform crossing-over and gene-conversion rates. The positions of segregating sites, population crossing-over and gene-conversion rates can be specified by the user at the top of the file fixeds.c. Program can be compiled using: g++ -O3 -o fixeds fixeds.c mtrand.cpp -Wno-deprecated and run by typing ./fixeds. The summaries that are output are respectively:
1. Frequency of D' < 1.0 (30% acceptance error)
2. Frequency of D' < 0.75
3. Frequency of D' < 0.50
4. Average number of distinct haplotypes for 5 kb sliding windows with 4 kb overlap. (15% acceptance error)
5. Average number of distinct haplotypes for 10 kb sliding windows with 9 kb overlap. (15% acceptance error)
6. Total number of distinct haplotypes H (15% acceptance error).
7. Frequency of D'(AB) < D'(AC) and D'(BC) < D'(AC) for 5 kb range. (30% acceptance error)
8. Frequency of D'(AB) < D'(AC) or D'(BC) < D'(AC) for 5 kb range. (30% acceptance error)
9. Frequency of D'(AB) < D'(AC) and D'(BC) < D'(AC) for 10 kb range. (30% acceptance error)
10. Frequency of D'(AB) < D'(AC) or D'(BC) < D'(AC) for 10 kb range. (30% acceptance error)
11. Frequency of D'(AB) < 1.00 and D'(BC) < 1.00 for 50 kb range.
12. Frequency of D'(AB) < 0.75 and D'(BC) < 0.75 for 50 kb range.
13. Frequency of D'(AB) < 0.50 and D'(BC) < 0.50 for 50 kb range. (30% acceptance error)
14. Frequency of D'(AB) < 0.25 and D'(BC) < 0.25 for 50 kb range.
15. Frequency of D'(AB) < 0.10 and D'(BC) < 0.10 for 50 kb range.
16. Frequency of D'(AB) < 0.50 and D'(BC) < 0.50 for 5 kb range.
17. Frequency of D'(AB) < 0.25 and D'(BC) < 0.25 for 5 kb range.
18. Frequency of D'(AB) < 0.50 and D'(BC) < 0.50 for 10 kb range.
19. Frequency of D'(AB) < 0.25 and D'(BC) < 0.25 for 10 kb range.
20. Frequency of D'(AB) < 0.10 and D'(BC) < 0.10 for 10 kb range.

Summaries that are most accurate for joint crossing-over and gene-conversion estimation are 1, 6, 13 for long-range data and 4, 5, 7, 8, 9, 10 for short-range. Suggested acceptance errors are given in brackets. These choices are slightly different from those used in Padhukasahasram et al 2006 and work better under fixed segregating sites model. Including the bounds on the minimum number of recombination events (Myers and Griffiths 2003, Song, Wu and Gusfield 2005 or Bafna and Bansal 2005) in the rejection scheme and smoothing the likelihood surfaces can result in further improvements in accuracy. For sequence length different from 50 kb, it may be necessary to change the distance cutoffs i.e. choose values other than 5 kb, 10 kb and 50 kb. Distance cutoffs can be changed at the end of the code. 
