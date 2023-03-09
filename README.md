# FeasChar
MAGMA routines for determining feasible characters of groups in exceptional algebraic groups

### Files included:
- FeasChar.M - Defines FeasChar(G,LIE_TYPE,p) which outputs a table of feasible characters for the finite group G on the adjoint module for an exceptional group of type LIE_TYPE and also on a minimal non-trivial module when LIE_TYPE is not E8
    
- EltTraces.M - Defines EFOs_to_file(n) which writes a file named "n" containing the eigenvalues of elements of order n in groups of exceptional type.
  
-  2.M, 3.M, ... - Defines pre-calculated eigenvalues for elements of these orders.

- ModsByInduction.M - Defines ModsByInduction, a routine while allows computation of all low-dimensional irreducible modules for a given group. This routine will finish (quickly) in some situations where MAGMA's in-built routines will fail (because they attempt to compute _all_ irreducible modules, including some whose dimension is too large for us to care about here).

### Usage:
- run MAGMA
- [if not pre-computed] load "EltTraces.M" and run EFOs_to_file(n) for each order n of an element in the finite group of interest.
- load "n.M" for each relevent element order n
- Run FeasChar(G,LIE_TYPE,p) where LIE_TYPE is one of "G2", "F4", "E6", "E7", "E8" and p is 0 or prime.

### Customisation

Some groups have elements of very large order, for which generating or storing the relevant element orders may be prohibitive (in terms of time or space). It is therefore possible to set an optional parameter LIMITING_ORDER in FeasChar, which will cause elements of larger orders to be ignored in the calculations. This comes with the risk of finding more "feasible" characters which will not actually correspond to an embedding of the group G.
- Example usage: FeasChar(G,"E8",2 : LIMITING_ORDER := 25);

### Modules for groups
