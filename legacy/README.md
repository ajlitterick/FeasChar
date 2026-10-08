# Legacy Magma code

*The code in this folder, and `README-original.md`, were written by A. J. Litterick. This
README was written by an AI system (Claude), under the author's direction, in October 2026.
See "Who wrote what" in the main README.*

This folder contains the original Magma implementation of FeasChar. It is the code that was used
to generate the tables of feasible characters in Chapter 6 of

> A. J. Litterick, *On non-generic finite subgroups of exceptional algebraic groups*,
> Mem. Amer. Math. Soc. **253** (2018), no. 1207.
> [doi:10.1090/memo/1207](https://doi.org/10.1090/memo/1207),
> [arXiv:1511.03356](https://arxiv.org/abs/1511.03356).

It is kept unchanged for reference. For reproducing or checking the tables, use the GAP code in
the top-level folder (see the main `README.md`), which supersedes this.

## Files

- `FeasChar.M` defines `FeasChar(G, LIE_TYPE, p)`. It gives the feasible characters of a finite
  group `G` on the adjoint module of an exceptional group of type `LIE_TYPE`, and also on a
  minimal module when `LIE_TYPE` is not E8.
- `EltTraces.M` defines `EFOs_to_file(n)`. This writes `n.M`, holding the traces of semisimple
  elements of order `n` on the relevant modules.
- `2.M`, …, `17.M` hold the precomputed element data for orders 2–17.
- `ModsByInduction.M` computes the absolutely irreducible modules of dimension up to a bound.
- `README-original.md` is the original README of this repository, with usage notes and an example.

## Notes on this version

These points matter when comparing its output with the Memoir:

1. **Order limit.** `FeasChar` ignores elements of order above `LIMITING_ORDER`. This defaults
   to 17, while the Memoir states that orders 2–37 were used. A few Memoir rows fail only at
   elements of order 18–37, which suggests that the default was in force for those tables.
   See `../tables/DIFFERENCES.md`.
2. **The 3875-dimensional module for E8.** For E8 the released code tests only the adjoint
   module. The Memoir's E8 tables agree, with very few exceptions, with the characters that are
   also compatible with the Weyl module V(λ₁) of dimension 3875. So that test was in effect used
   for the Memoir, but it is not in this code.
3. **Element data.** The files `n.M` cover orders up to 17 only.

The GAP implementation in the top-level folder:
- tests every element order;
- includes all end-node modules;
- checks power maps between classes.

It reproduces these files' element data exactly, for every type and n = 2–17.
