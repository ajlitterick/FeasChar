# verify_row.g -- check a row of the updated tables from its witness, independently of the
# search code (feasgen.g is not used).
#
# Usage, from the repository folder:
#   gap> Read("tables/E8.g");               # defines FEASCHAR_TABLES.E8 (also F4, E6, E7)
#   gap> Read("tools/verify_row.g");
#   gap> VerifyRow(FEASCHAR_TABLES.E8, "6.243", 1);     # table 6.243, row 1
#   gap> VerifyTable(FEASCHAR_TABLES.E8, "6.243");       # all rows of table 6.243
#   gap> VerifyAll(FEASCHAR_TABLES.E8);                  # every recomputed table
#
# For a row, the witness assigns to each p-regular class of H an element t of the simply
# connected group G, given by Kac coordinates s (t = exp(2 pi i x), x = sum (s_i/N) w_i^vee).
# The checks are:
#   1. s is a point of the fundamental alcove (s >= 0, marks . s <= N), and t has exact order N
#      equal to the order of the class; t is central exactly when the class is central.
#   2. Power maps: for every class c of order N and every e = 1 .. N-1, the eigenvalues of t^e on
#      each module (computed from those of t) are those of the element assigned to c^e.
#   3. For each module, the trace function decomposes over the irreducible (Brauer) characters
#      of H with exactly the multiplicities stated in the row, which are non-negative integers.
#   4. The displayed row (L(G), V_min) has the same composition factors as the witness, up to
#      relabelling (rows are shown up to table automorphisms).
# Module weights come from the dominant weights and multiplicities stored in the entry; for the
# Weyl modules V(lambda) these are also compared with GAP's DominantCharacter.
#
# Provenance: written by an AI system (Claude, Anthropic), October 2026; see the README.

VR_Find := function(list, number)
  local e;
  e := First(list, t -> t.number = number);
  if e = fail then Error("no table ", number); fi;
  return e;
end;

VR_Group := function(type)   # "E8" -> data of the simply connected group
  local L, R, S, Sinv, pos, theta;
  L := SimpleLieAlgebra(type{[1]}, Int(type{[2 .. Length(type)]}), Rationals);
  R := RootSystem(L);
  S := SimpleSystem(R); Sinv := S^-1;
  pos := PositiveRoots(R) * Sinv;
  theta := First(pos, c -> Sum(c) = Maximum(List(pos, Sum)));
  return rec(L := L, W := WeylGroup(R), S := S, Sinv := Sinv, marks := theta,
             rank := Length(S));
end;

VR_Weights := function(G, M)   # all weights (root coordinates) with multiplicities
  local wts, mults, i, w, dc;
  if M.kind = "V" then
    dc := DominantCharacter(G.L, M.highest_weight);
    if Set(TransposedMat(dc)) <> Set(TransposedMat([M.dominant_weights, M.multiplicities])) then
      Error("module ", M.name, ": dominant weights differ from DominantCharacter");
    fi;
  fi;
  wts := []; mults := [];
  for i in [1 .. Length(M.dominant_weights)] do
    for w in WeylOrbitIterator(G.W, M.dominant_weights[i]) do
      Add(wts, w * G.Sinv); Add(mults, M.multiplicities[i]);
    od;
  od;
  if Sum(mults) <> M.dim then Error("module ", M.name, ": dimension ", Sum(mults), " <> ", M.dim); fi;
  return [wts, mults];
end;

VR_Counts := function(WM, s, N)   # multiplicity of zeta_N^k as an eigenvalue (index k+1)
  local cnt, P, i, e;
  P := WM[1] * s;
  cnt := ListWithIdenticalEntries(N, 0);
  for i in [1 .. Length(P)] do
    if not IsInt(P[i]) then Error("non-integral eigenvalue exponent"); fi;
    e := P[i] mod N + 1;
    cnt[e] := cnt[e] + WM[2][i];
  od;
  return cnt;
end;

VR_PowerCounts := function(cnt, N, e)   # counts of t^e (of order N/gcd(N,e)) from those of t
  local g, N2, c2, k;
  g := Gcd(N, e); N2 := N / g;
  c2 := ListWithIdenticalEntries(N2, 0);
  for k in [0 .. N - 1] do
    if cnt[k + 1] <> 0 then
      c2[(k * e / g) mod N2 + 1] := c2[(k * e / g) mod N2 + 1] + cnt[k + 1];
    fi;
  od;
  return c2;
end;

VR_Cache := rec();

VerifyRow := function(list, number, k)
  local T, row, G, tbl, modtbl, fus, inv, ok, msgs, W, cls, nC, j, c, N, s, v, ordH, cent,
        M, WM, cnts, e, je, chi, B, m, mults, name, prof, degs, disp, Prof, mods, ia, iv, tau;
  T := VR_Find(list, number);
  row := First(T.rows, r -> r.row = k);
  if row = fail then Error("table ", number, " has no row ", k); fi;
  ok := true; msgs := [];
  G := VR_Group(T.type);
  tbl := CharacterTable(T.character_table.name);
  if T.characteristic = 0 then modtbl := tbl; fus := [1 .. NrConjugacyClasses(tbl)];
  else modtbl := tbl mod T.characteristic; fus := GetFusionMap(modtbl, tbl); fi;
  inv := []; for j in [1 .. Length(fus)] do inv[fus[j]] := j; od;
  W := row.witness; cls := W.classes; nC := Length(fus);
  if List(cls, x -> x.class) <> ClassNames(tbl){fus} then Error("class names do not match the library table"); fi;
  ordH := OrdersClassRepresentatives(tbl){fus};
  cent := List(fus, c -> SizesConjugacyClasses(tbl)[c] = 1);
  # 1. alcove points, orders, centrality
  for j in [1 .. nC] do
    s := cls[j].kac; N := cls[j].order;
    if N <> ordH[j] then ok := false; Add(msgs, Concatenation("order mismatch at ", cls[j].class)); fi;
    if N > 1 then
      if not (ForAll(s, x -> x >= 0) and G.marks * s <= N) then
        ok := false; Add(msgs, Concatenation("not an alcove point at ", cls[j].class));
      fi;
      v := G.Sinv * s / N;
      if Lcm(List(v, DenominatorRat)) <> N then
        ok := false; Add(msgs, Concatenation("wrong element order at ", cls[j].class));
      fi;
      if cent[j] <> ForAll(s, x -> x = 0 or x = N) then
        ok := false; Add(msgs, Concatenation("centrality mismatch at ", cls[j].class));
      fi;
    fi;
  od;
  # 2. and 3., module by module
  B := List(Irr(modtbl), ValuesOfClassFunction);
  for M in T.modules_tested do
    name := Concatenation(T.type, "_", M.name);
    if not IsBound(VR_Cache.(name)) then VR_Cache.(name) := VR_Weights(G, M); fi;
    WM := VR_Cache.(name);
    cnts := List([1 .. nC], j -> VR_Counts(WM, cls[j].kac, cls[j].order));
    for j in [2 .. nC] do
      N := cls[j].order;
      for e in [2 .. N - 1] do
        je := inv[PowerMap(tbl, e)[fus[j]]];
        if VR_PowerCounts(cnts[j], N, e) <> cnts[je] then
          ok := false;
          Add(msgs, Concatenation("power map: ", M.name, " at ", cls[j].class, "^", String(e)));
        fi;
      od;
    od;
    chi := List([1 .. nC], j -> CycList(cnts[j]));
    m := chi * B^-1;
    if not ForAll(m, x -> IsInt(x) and x >= 0) then
      ok := false; Add(msgs, Concatenation(M.name, ": not a non-negative integer combination"));
    fi;
    if m <> W.multiplicities.(M.name) then
      ok := false; Add(msgs, Concatenation(M.name, ": multiplicities differ from the stated ones"));
    fi;
  od;
  # 4. displayed row against the witness (composition factors, up to relabelling)
  degs := List(T.irreducibles, x -> x.degree);
  mods := T.display_modules;
  Prof := function(vs)
    return SortedList(Filtered(List([1 .. Length(degs)], i -> Concatenation([degs[i]], List(vs, v -> v[i]))),
                               x -> ForAny(x{[2 .. Length(x)]}, y -> y <> 0)));
  end;
  disp := Prof(List(mods, nm -> row.display.(nm)));
  prof := Prof(List(mods, nm -> W.multiplicities.(nm)));
  tau := T.type = "F4" and T.characteristic = 2;
  if disp <> prof and not (tau and disp = Prof([W.multiplicities.(mods[1]),
                             W.multiplicities.(mods[1]) - W.multiplicities.(mods[2])])) then
    ok := false; Add(msgs, "displayed row and witness have different composition factors");
  fi;
  return rec(table := number, row := k, ok := ok, problems := msgs);
end;

VerifyTable := function(list, number)
  local T, out, k, r;
  T := VR_Find(list, number);
  out := [];
  for k in List(T.rows, r -> r.row) do
    r := VerifyRow(list, number, k);
    Print("Table ", number, " row ", k, ": ", ListBlist(["FAILED", "verified"], [not r.ok, r.ok])[1]);
    if r.problems <> [] then Print("  ", r.problems); fi;
    Print("\n");
    Add(out, r);
  od;
  return out;
end;

VerifyAll := function(list)
  local res, T;
  res := [];
  for T in list do
    if T.status <> "NOT RECOMPUTED" then Append(res, VerifyTable(list, T.number)); fi;
  od;
  Print(Number(res, r -> r.ok), " of ", Length(res), " rows verified\n");
  return res;
end;
