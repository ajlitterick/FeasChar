# FeasChar

GAP code for computing **feasible characters** of finite groups in simple algebraic groups, of any
Dynkin type.

It reproduces and updates the tables of Chapter 6 of

> A. J. Litterick, *On non-generic finite subgroups of exceptional algebraic groups*,
> Mem. Amer. Math. Soc. **253** (2018), no. 1207.
> [doi:10.1090/memo/1207](https://doi.org/10.1090/memo/1207),
> [arXiv:1511.03356](https://arxiv.org/abs/1511.03356).

The original Magma code, which generated the tables in the Memoir, is kept unchanged in
[`legacy/`](legacy/).

## What is computed

Let G be a simply connected simple algebraic group in characteristic p ≥ 0, and H a finite group,
for example a minimal preimage in G of a finite simple subgroup of G/Z(G). A **feasible
character** of H on a G-module V is a Brauer character of H which is consistent with H being a
subgroup of G. That is, there is an assignment of semisimple classes of G to the p-regular
classes of H such that:

- every class of H of order n goes to a class of G of order n, central classes going to central
  classes;
- the assignment is compatible with the power maps of H;
- for each module V considered, the resulting trace function is a non-negative integer
  combination of the irreducible Brauer characters of H.

The search runs over these class assignments ("value-first") rather than over multiplicity
vectors. All modules are tested together, with interval bounds pruning the search.

**Modules tested by default**
- the adjoint module V(θ);
- the Weyl modules V(ωᵢ) for each end node i of the Dynkin diagram, keeping only one of each
  dual pair, so E₆ uses the 27 but not the 27*;
- in characteristic p, also the irreducible module L(λ) whenever it is smaller than the Weyl
  module V(λ). For example:
  - F₄, p = 3: L(ω₄) = 25;
  - F₄, p = 2: L(ω₁) = 26̃;
  - E₇, p = 2: L(ω₁) = 132;
  - E₈, p = 2: L(ω₁) = 3626.

  The weight multiplicities of L(λ) are F. Lübeck's (see `lubeck_data.g`).

**Classes of G**
- The semisimple classes of G of order dividing N are the points of the fundamental alcove
  (Kac coordinates), so central elements are handled correctly.
- Every element order of H is used, with no cut-off. Very large modules (dimension over 100) are
  used for element orders up to 37 by default; see `MaxOrder`.

## Quick start

You need [GAP](https://www.gap-system.org) 4.12 or later, with the CTblLib package, which is
included in standard installations. From this folder:

```gap
gap> Read("feasgen.g");
gap> tbl := CharacterTable("A6");;
gap> res := FeasibleCharacters("F", 4, tbl, 0);;      # Alt6 < F4, p = 0  (Memoir Table 6.5)
gap> FGPrintTable(res, tbl);
```

The output is the Memoir's six rows, on V₅₂ and V₂₆, listed up to automorphisms of H. For other
characteristics, pass a prime as the last argument; the Brauer table `tbl mod p` is used. For
other types:

```gap
gap> FeasibleCharacters("E", 8, CharacterTable("M11"), 2);   # adjoint, 3875, 147250, L(w1), L(w2)
gap> FeasibleCharacters("A", 2, CharacterTable("A5"), 0);    # A5 in SL3: adjoint 8, natural 3
gap> FeasibleCharacters("C", 3, CharacterTable("L2(7)"), 0); # L2(7) in Sp6
```

More are in [`examples/`](examples/).

### Main function

`FeasibleCharacters(type, rank, tbl, p : options)`
- `type`, `rank`: the Dynkin type of G, for example `"E", 8`, `"F", 4` or `"B", 3`.
- `tbl`: the ordinary character table of H.
- `p`: 0, or a prime.

| Option | Meaning |
|---|---|
| `Modules` | the list of modules to test (`FGModule(G, hw)`, `FGIrrModule(G, hw, p)`); the default is described above |
| `PowerMaps` | `"full"` (the default): one class of G per class of H, compatible with all power maps. `"divisors"`: the original FeasChar test, comparing only traces of x^(n/j) for j dividing n, which reproduces the Memoir's setting |
| `MaxOrder` | an element-order limit per module; modules that cannot be used for all elements of H are dropped and reported |
| `Strategy` | `"auto"` (the default), `"value"` or `"hybrid"`; the result is the same in each case (see below) |

### How the search works

- **One step per cyclic subgroup.** Classes of H that generate conjugate cyclic subgroups are
  Galois conjugate, so they are assigned together, in one step.
- **Candidates are roots.** Steps are taken in increasing element order. When the step for a
  class c of order n is reached, the classes of c^p are already assigned for every prime p
  dividing n, so only the classes of G that are p-th roots of those are tried.
- **Hybrid strategy.** Power maps do not determine classes of composite order: an involution
  usually has several classes of square roots. So, when the smallest module (V_min, or the
  adjoint for E₈) has few degree decompositions, the search first lists the possible characters
  of H on it. Then, for each, it searches with every class restricted to the classes of G with
  the right trace. For example, Sp₆(2) < E₇ takes about 1 s this way, against over an hour
  without it.
- **Pruning.** Interval bounds on the partial multiplicities, for every module and every
  irreducible character, prune branches early. Each solution is then checked exactly.

The result is a record:
- `res.sols` lists the compatible collections of multiplicity vectors, one per module, each with
  a witness assignment of classes of G (Kac coordinates).
- `FGProject(res, j)` gives the distinct characters on module `j`.
- `FGPrintTable(res, tbl)` prints a Memoir-style table.

## The Memoir's tables

[`tables/`](tables/) contains the **updated tables**, numbered exactly as in the Memoir
(Table 6.1 – Table 6.352), with the same captions and P/N labels:

- `tables/tables.pdf` (compiled from `tables/tables.tex`): all tables, with
  - the rows that are compatible with **all** the modules above;
  - the full power-map test;
  - every element order.
- `tables/DIFFERENCES.md`: every difference from the Memoir, with its reason.
- `tables/F4.json`, `E6.json`, `E7.json`, `E8.json`: the same tables in machine-readable form,
  plus `tables/*.g`, the same data as GAP files. Each table records:
  - its number, caption, group and characteristic;
  - the irreducible characters, as positions in GAP's character-table library;
  - every module tested, with its dominant weights and multiplicities;
  - for each row, the multiplicities on **every** module tested, not only L(G) and V_min.
  - for each row, a **witness**: the Kac coordinates of an element of G for each p-regular class
    of H.

### Checking a row independently

[`tools/verify_row.g`](tools/verify_row.g) checks any row from its witness alone. It does not
use the search code, so it is independent of `feasgen.g`.

```gap
gap> Read("tables/E8.g");; Read("tools/verify_row.g");
gap> VerifyRow(FEASCHAR_TABLES.E8, "6.243", 1);      # Table 6.243, row 1
gap> VerifyAll(FEASCHAR_TABLES.E8);                  # every row of every E8 table
```

It checks that:
- each witness is a point of the fundamental alcove of G, of the right order, and central exactly
  when the class of H is central;
- the eigenvalues on every module are compatible with all power maps of H;
- the resulting trace functions decompose over the irreducible (Brauer) characters of H with
  exactly the stated, non-negative multiplicities.

The same checks are easy to reproduce in other systems from the JSON files: Magma, Sage or
Python.

**Limitation.** This shows that every listed row **is** feasible. It does not show that a table
is **complete**, that is, that no other feasible rows exist. Completeness rests on the search
(`feasgen.g`) and on how it was checked against the legacy code and the Memoir (see
"Verification" below).

Most tables are unchanged. The differences fall into a few kinds:
- rows removed by the extra modules: V(λ₂) for E₇ and E₈, and L(λ) in bad characteristic;
- rows that fail at elements of order 18–37 (the legacy code's default order limit was 17);
- a few rows missing from the Memoir;
- typographical errors.

## Verification

`tools/` contains the scripts used to check this code against the legacy code and the Memoir:

- **Class data** (`tools/fg_dump.g`, `tools/check_classdata.py`): for types G₂, F₄, E₆, E₇ and E₈
  and element orders 2–17, the traces of all semisimple elements on the adjoint and minimal modules
  agree exactly with the legacy files `legacy/n.M`.
- **Search:** in `"divisors"` mode, the GAP search returns exactly the characters of the legacy
  E₈ search, on every Memoir table tested.
- **Tables** (`tools/memoir_numbers.py`, `tools/build_index.py`, `tools/parse_tables_gen.py`,
  `tools/check_dims.py`): these turn the Memoir's LaTeX tables into data, number them as in the
  Memoir, and check that every row has the right dimension.

## Who wrote what: human and AI contributions

This repository mixes work done by people and work done by an AI system.

**Human work**
- **Mathematics and data:**
  - the Memoir and its tables are by A. J. Litterick;
  - the code in [`legacy/`](legacy/) and its original README (`legacy/README-original.md`) are by
    A. J. Litterick;
  - the weight multiplicities of the irreducible modules L(λ) are F. Lübeck's published data,
    from his tables of weight multiplicities.
- **Direction of the new work** (by A. J. Litterick):
  - the requirements: a search over class assignments rather than multiplicity vectors, all
    end-node modules, and any Dynkin type;
  - the mathematical corrections made during development:
    - use the irreducible modules L(λ) in bad characteristic;
    - treat the two 27-dimensional E₆ modules as a single module;
    - account for the graph automorphisms of F₄ (p = 2) and G₂ (p = 3);
  - the layout of this repository.

**AI-generated work**, by Claude (Anthropic, model Claude Opus 5.5) in October 2026, working
under the author's direction:
- `feasgen.g`;
- `lubeck_data.g` (converted from Lübeck's web pages by `tools/lubeck2gap.py`);
- everything in [`tools/`](tools/) and [`examples/`](examples/);
- the updated tables in [`tables/`](tables/), including `DIFFERENCES.md`;
- this README and `legacy/README.md`.

**How the AI-generated code has been checked.** The checks below were designed and run by the AI,
and no person has independently audited the code.
- **Class data:** agree with the legacy Magma files `legacy/n.M` for every type and element order
  2–17.
- **Search:** in `"divisors"` mode it returns exactly the legacy-style E₈ results on every table
  tested.
- **Tables:** reproduce the Memoir except where `tables/DIFFERENCES.md` says otherwise.
- **Spot checks:** for classical types, for example SL(2,5) < SL₂ and A₅ < SL₃.

Treat any new claim in `tables/DIFFERENCES.md` (rows removed, rows added, errata) as a
computational result to be confirmed. Do not treat it as a verified theorem.

Commits made with AI assistance carry a `Co-Authored-By: Claude` line.

## Licence

See [LICENSE](LICENSE).
