# FeasChar
MAGMA routines for determining feasible characters of groups in exceptional algebraic groups.

This code is being released to allow verification of the tables found in [_On non-generic finite subgroups of exceptional algebraic groups_](https://doi.org/10.1090/memo/1207) [(ArXiv version)](https://arxiv.org/abs/1511.03356). This follows an error found in Table 6.298 there, where the feasible characters should each have a 29-dimensional constituent instead of one 28-dimensional constituent. This error occurred during shortening of the tables for space purposes; the present algorithm, a version of which was used for almost all cases in the _Memoir_, returns the correct table.

### Files included:
- FeasChar.M - Defines FeasChar(G,LIE_TYPE,p) which outputs a table of feasible characters for the finite group G on the adjoint module for an exceptional group of type LIE_TYPE and also on a minimal non-trivial module when LIE_TYPE is not E8.
    
- EltTraces.M - Defines EFOs_to_file(n) which writes a file named "n.M" containing the eigenvalues of elements of order n in groups of exceptional type.
  
-  2.M, 3.M, ..., 17.M - Defines pre-calculated eigenvalues for elements of these orders.

- ModsByInduction.M - Defines ModsByInduction, a routine while allows computation of absolutely irreducible modules for a given group. Importantly, by specifying an optional parameter "DimLim", this will produce only modules of this dimension or less, in such a way as to produce significant speed-ups in some situations where MAGMA's in-built routines will fail (because they attempt to compute _all_ irreducible modules, including some whose dimension is too large for us to care about here).

### Usage:
- run MAGMA
- [if not pre-computed] load "EltTraces.M" and run EFOs_to_file(n) for each order n of an element in the finite group of interest
- load "n.M" for each relevent element order n
- Run FeasChar(G,LIE_TYPE,p) where LIE_TYPE is one of "G2", "F4", "E6", "E7", "E8" and p is 0 or prime

### Customisation

Some groups have elements of very large order, for which generating or storing the relevant element orders may be prohibitive (in terms of time or space). It is therefore possible to set an optional parameter LIMITING_ORDER in FeasChar, which will cause elements of larger orders to be ignored in the calculations. This comes with the risk of finding more "feasible" characters which will not actually correspond to an embedding of the group G.

Similarly, to avoid calculating all irreducible modules for a group (because, for instance, ModsByInduction has been used to find those of dimension at most 248), setting the optional variable MODULES tells FeasChar to only consider characters coming from these particular modules.

### Example usage:

    > G := Alt(17);
    > FeasChar(G,"E8",2 : LIMITING_ORDER := 17, MODULES := ModsByInduction(G,2 : DimLim := 248) );
    
    *** WARNING: Elements of order(s) [ 21, 33, 35, 39, 45, 55, 63, 105 ] not used.
    \begin{tabular}{r|*{5}{c}}
    & \multicolumn{5}{c}{$V_{248}$} \\ \hline
            & 1 &16 &118 &128 &128 \\ \hline
    %
     1)     & 2 & 0 & 1 & 0 & 1 \\
     2)     & 2 & 0 & 1 & 1 & 0  \\ \hline
    \end{tabular}
