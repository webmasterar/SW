SW: Smith-Waterman
===

Two implementations of the Smith-Waterman algorithm.

The first implementation (LS) is a linear space algorithm and just prints out
the similarity score for the longest shared factor between two sequences. The
second algorithm (TB) can be used for pairwise circular sequence alignment but
we recommend you use Circular Needleman-Wunsch (see the cNW project) instead.
TB requires O(nm) space but does a traceback to find the best position to
rotate the sequence at and also prints out the similarity score.

If you plan to use circular sequences you need to double up your sequences so
the longest factor of *x* is found in *y* regardless of where its starts.

GNU GPLv3 License; Copyright (C) 2015 Solon P. Pissis, Ahmad Retha and Fatima Vayani.

