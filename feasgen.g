# feasgen.g -- feasible characters of a finite group on the adjoint and end-node modules of a
# simply connected simple algebraic group G of any Dynkin type, by a value-first search.
#
# Provenance: written by an AI system (Claude, Anthropic; Claude Opus 5.5), October 2026, under
# the direction of A. J. Litterick, who specified the method and supplied mathematical
# corrections.  Validated as described in the README; not independently audited.
#
# Instead of enumerating multiplicity vectors (FeasChar's odometer), we assign to each
# p-regular class of H a semisimple class of G of the same order, and read off the
# multiplicities m = chi * IBr^-1 for every module at once (adjoint, end-node Weyl modules,
# and in bad characteristic the irreducible L(lambda)).  Partial sums of m are bounded with
# floating-point intervals, so branches giving a negative (or too large) multiplicity on some
# module are cut early.
#
# Main entry point:
#   FeasibleCharacters(type, rank, tbl, p : Modules, PowerMaps, MaxOrder, Eps, NodeLimit)
#     type, rank  Dynkin type of G (simply connected), e.g. "E", 8 or "F", 4 or "A", 2
#     tbl         ordinary character table of H (a minimal preimage S~ of S < G/Z(G))
#     p           0 (ordinary characters) or a prime (Brauer characters of tbl mod p)
#   Options:
#     Modules     list of module records (FGModule, FGIrrModule); default FGDefaultModules(G, p):
#                 the adjoint V(theta) and V(omega_i) for each end node i (not both of a dual
#                 pair), each followed by L(lambda) when that is smaller in characteristic p
#                 (Luebeck's data, lubeck_data.g)
#     PowerMaps   "full" (default): one semisimple class of G per class of H, compatible with
#                 all power maps (including Galois conjugation); "divisors": FeasChar's
#                 test, only traces of x^(n/j), j | n, compared (reproduces FeasChar/memoir)
#     MaxOrder    integer or list (per module): modules with a smaller limit than the largest
#                 p-regular element order of H are dropped (reported in res.dropped); default
#                 unlimited for the adjoint and modules of dimension <= 100, else 37
#     Strategy    "auto" (default), "value" or "hybrid".  "value": branch over classes of G for
#                 each cyclic subgroup class of H (candidates indexed by their p-th powers, so
#                 only roots of already assigned classes are tried).  "hybrid": first enumerate
#                 the characters of H on the smallest module (degree decompositions whose values
#                 are traces of classes of G), then for each one search with every class
#                 restricted to the classes of G with that trace.  "auto" uses "hybrid" when the
#                 smallest module has at most 200000 degree decompositions.  Same result either way.
#     Eps         tolerance of the floating-point bounds (default 1e-3)
#     NodeLimit   abort after this many search nodes (res.aborted = true)
#   Returns rec(group, p, mode, modules, degs, sols, dropped, nodes, aborted, steps, time);
#   sols lists the distinct compatible collections rec(mults := [m_1, ..., m_k], fusion), one
#   multiplicity vector over IBr(H mod p) per module, with a witness fusion (in full mode,
#   [order, Kac coordinates] of the image of each p-regular class).
# Helpers: FGProject(res, j), FGModuleIndex(res, name), FGJointProfile(degs, mlist),
#   FGPrintTable(res, tbl : Modules := [indices]) (memoir-style table, condensed under table
#   automorphisms).
#
# Class data are validated against FeasChar's element files 2.M - 17.M (G2, F4, E6, E7, E8),
# and the E8 search against the earlier E8 code; see feasgen-report.md.
#
# Conventions.
#  * GAP gives weights in fundamental-weight coordinates; S := SimpleSystem(R) has the simple
#    roots as rows, so a weight w has simple-root coordinates w * S^-1.
#  * The semisimple classes of G_sc of order dividing N are the points x = sum (s_i/N) w_i^vee
#    of the closed fundamental alcove, i.e. tuples s >= 0 with marks . s <= N (Kac
#    coordinates without s_0).  The element t = exp(2 pi i x) has eigenvalue
#    zeta_N^(c(lambda) . s) on a weight lambda with root coordinates c(lambda), and its order is
#    the lcm of the denominators of (omega_j . s)/N.
#  * Bourbaki labelling is used only for Luebeck's data (lubeck_data.g); the permutation from
#    Bourbaki to GAP node labels is found automatically.

if not IsBound(FG_DIR) then FG_DIR := "."; fi;
if not IsBound(LUBECK) then
  if IsReadableFile(Concatenation(FG_DIR, "/lubeck_data.g")) then
    Read(Concatenation(FG_DIR, "/lubeck_data.g"));
  else LUBECK := rec(); fi;
fi;
FG_CACHE_DIR := Concatenation(FG_DIR, "/cache");
FG_TMP := fail;
FG := rec(groups := rec(), counts := rec(), classdata := rec(), cossin := []);

FG_BOURBAKI_EDGES := rec(
  E6 := [[1,3],[3,4],[4,5],[5,6],[2,4]],
  E7 := [[1,3],[3,4],[4,5],[5,6],[6,7],[2,4]],
  E8 := [[1,3],[3,4],[4,5],[5,6],[6,7],[7,8],[2,4]],
  F4 := [[1,2],[2,3],[3,4]],
  G2 := [[1,2]]);

# ---------------------------------------------------------------------------------------
# Floats

FGCosSin := function(N)
  if not IsBound(FG.cossin[N]) then
    FG.cossin[N] := [List([0 .. N - 1], k -> Cos(2 * FLOAT.PI * k / N)),
                     List([0 .. N - 1], k -> Sin(2 * FLOAT.PI * k / N))];
  fi;
  return FG.cossin[N];
end;

FGComplex := function(x)    # [Re x, Im x] as floats
  local N, co, cs, re, im, k;
  if IsRat(x) then return [Float(x), 0.]; fi;
  N := Conductor(x); co := CoeffsCyc(x, N); cs := FGCosSin(N);
  re := 0.; im := 0.;
  for k in [1 .. N] do
    if co[k] <> 0 then
      re := re + Float(co[k]) * cs[1][k]; im := im + Float(co[k]) * cs[2][k];
    fi;
  od;
  return [re, im];
end;

# ---------------------------------------------------------------------------------------
# The group

FGGroup := function(type, rank)
  local key, L, R, S, Sinv, C, pos, theta, ends, G;
  key := Concatenation(type, String(rank));
  if IsBound(FG.groups.(key)) then return FG.groups.(key); fi;
  L := SimpleLieAlgebra(type, rank, Rationals);
  R := RootSystem(L);
  S := SimpleSystem(R);
  Sinv := S^-1;
  C := CartanMatrix(R);
  pos := PositiveRoots(R) * Sinv;
  if not ForAll(Flat(pos), IsInt) then Error("roots not integral in root coordinates"); fi;
  theta := First(pos, c -> Sum(c) = Maximum(List(pos, Sum)));
  ends := Filtered([1 .. rank], i -> Number([1 .. rank], j -> j <> i and C[i][j] <> 0) <= 1);
  G := rec(type := type, rank := rank, name := key, L := L, R := R, W := WeylGroup(R),
           S := S, Sinv := Sinv, C := C, roots := Concatenation(pos, -pos),
           theta := theta, marks := theta, thetaFW := theta * S, ends := ends,
           tuples := rec());
  FG.groups.(key) := G;
  return G;
end;

FGDualWeight := function(G, hw)    # -w_0(hw), as a dominant weight
  local w, i;
  w := -ShallowCopy(hw);
  while true do
    i := PositionProperty(w, x -> x < 0);
    if i = fail then return w; fi;
    w := w - w[i] * G.S[i];
  od;
end;

# ---------------------------------------------------------------------------------------
# Modules: weights in root coordinates with multiplicities

# The W-orbit of a dominant weight, in root coordinates (cached per group).
FGOrbit := function(G, dw)
  local key, wts, w;
  if not IsBound(G.orbits) then G.orbits := rec(); fi;
  key := JoinStringsWithSeparator(List(dw, String), ".");
  if not IsBound(G.orbits.(key)) then
    wts := [];
    for w in WeylOrbitIterator(G.W, dw) do Add(wts, w * G.Sinv); od;
    G.orbits.(key) := rec(dw := ShallowCopy(dw), key := key, wts := wts,
                          integral := ForAll(wts, w -> ForAll(w, IsInt)));
  fi;
  return G.orbits.(key);
end;

# A module is stored as dominant weights with multiplicities; its exponent counts are the
# corresponding combination of orbit counts (see FGCounts).
FGModuleFromDominant := function(G, hw, doms, mults, kind, p)
  local M;
  M := rec(hw := ShallowCopy(hw), kind := kind, p := p, doms := List(doms, ShallowCopy),
           dmults := ShallowCopy(mults));
  M.dim := Sum([1 .. Length(doms)], i -> mults[i] * Length(FGOrbit(G, doms[i]).wts));
  M.name := Concatenation(kind, String(M.dim));
  M.key := Concatenation(G.name, "_", kind, String(p), "_", JoinStringsWithSeparator(List(hw, String), "."));
  M.isadjoint := (kind = "V" and hw = G.thetaFW);
  return M;
end;

FGModule := function(G, hw)      # Weyl module V(hw); its character is the Weyl character
  local dc;
  dc := DominantCharacter(G.L, hw);
  return FGModuleFromDominant(G, hw, dc[1], dc[2], "V", 0);
end;

# Bourbaki -> GAP node permutation: pi[i] is the GAP node of Bourbaki node i.
FGBourbakiPerm := function(G)
  local r, edges, adj, perms, pi, ok, ent, hwG, dc, cand;
  if IsBound(G.bperm) then return G.bperm; fi;
  if not IsBound(FG_BOURBAKI_EDGES.(G.name)) or not IsBound(LUBECK.(G.name)) then
    G.bperm := fail; return fail;
  fi;
  r := G.rank; edges := FG_BOURBAKI_EDGES.(G.name);
  perms := Filtered(List(SymmetricGroup(r), g -> ListPerm(g, r)), pi ->
     Number(Combinations([1 .. r], 2), e -> G.C[pi[e[1]]][pi[e[2]]] <> 0) = Length(edges) and
     ForAll(edges, e -> G.C[pi[e[1]]][pi[e[2]]] <> 0));
  for pi in perms do
    ok := true;
    for ent in LUBECK.(G.name) do
      if ent.dim > 3000 or Sum(ent.hw) <> 1 or not (ent.chars = "all" or StartsWith(ent.chars, "not")) then
        continue;
      fi;
      hwG := []; hwG{pi} := ent.hw;
      dc := DominantCharacter(G.L, hwG);
      cand := List(ent.dom, function(d) local v; v := []; v{pi} := d[1]; return [v, d[2]]; end);
      if Set(cand) <> Set(TransposedMat(dc)) then ok := false; break; fi;
    od;
    if ok then G.bperm := pi; return pi; fi;
  od;
  Error("no Bourbaki labelling consistent with Luebeck's data");
end;

FGCharsMatch := function(chars, p)   # does Luebeck's characteristic description include p?
  local list;
  if chars = "all" then return true; fi;
  if StartsWith(chars, "not") then
    list := List(SplitString(chars{[4 .. Length(chars)]}, ","), s -> Int(NormalizedWhitespace(s)));
    return not p in list;
  fi;
  list := List(SplitString(chars, ","), s -> Int(NormalizedWhitespace(s)));
  return p in list;
end;

# Luebeck's entry for L(hw) in characteristic p (hw and weights converted to GAP labels),
# or fail if no data.
FGLubeck := function(G, hw, p)
  local pi, hwB, ents, ent;
  pi := FGBourbakiPerm(G);
  if pi = fail then return fail; fi;
  hwB := hw{pi};
  ents := Filtered(LUBECK.(G.name), e -> e.hw = hwB and FGCharsMatch(e.chars, p));
  if Length(ents) <> 1 then return fail; fi;
  ent := ents[1];
  return rec(dim := ent.dim, chars := ent.chars,
             doms := List(ent.dom, function(d) local v; v := []; v{pi} := d[1]; return v; end),
             mults := List(ent.dom, d -> d[2]));
end;

FGIrrModule := function(G, hw, p)   # L(hw) in characteristic p, from Luebeck's data
  local lb;
  lb := FGLubeck(G, hw, p);
  if lb = fail then return fail; fi;
  return FGModuleFromDominant(G, hw, lb.doms, lb.mults, "L", p);
end;

# Default modules: V(theta) and V(omega_i) for end nodes i (dropping duals of modules already
# present), each followed by L(lambda) when that is smaller in characteristic p.
FGDefaultModules := function(G, p)
  local hws, kept, hw, mods, V, Lm, e;
  hws := [G.thetaFW];
  for e in G.ends do Add(hws, List([1 .. G.rank], j -> 0)); hws[Length(hws)][e] := 1; od;
  kept := [];
  for hw in hws do
    if not hw in kept and not FGDualWeight(G, hw) in kept then Add(kept, hw); fi;
  od;
  mods := [];
  for hw in kept do
    V := FGModule(G, hw);
    Add(mods, V);
    if p > 0 then
      Lm := FGIrrModule(G, hw, p);
      if Lm <> fail and Lm.dim < V.dim then Add(mods, Lm); fi;
    fi;
  od;
  return mods;
end;

# ---------------------------------------------------------------------------------------
# Classes of G_sc of exact order N

FGTuples := function(G, N)
  local key, out, cen, marks, r, rec_;
  key := String(N);
  if IsBound(G.tuples.(key)) then return G.tuples.(key); fi;
  marks := G.marks; r := G.rank; out := []; cen := [];
  rec_ := function(i, s, used)
    local k, v;
    if i > r then
      v := G.Sinv * s / N;
      if Lcm(List(v, DenominatorRat)) = N then
        Add(out, ShallowCopy(s)); Add(cen, ForAll(s, x -> x = 0 or x = N));
      fi;
      return;
    fi;
    for k in [0 .. QuoInt(N - used, marks[i])] do
      s[i] := k;
      rec_(i + 1, s, used + k * marks[i]);
    od;
    Unbind(s[i]);
  end;
  rec_(1, [], 0);
  G.tuples.(key) := rec(tuples := out, central := cen);
  return G.tuples.(key);
end;

# Exponent counts of one W-orbit: for each tuple of exact order N, the number of orbit weights
# on which t acts as zeta_N^k (index k+1).  Cached in memory, and on disk for large orbits.
FGOrbitCounts := function(G, dw, N)
  local O, key, file, T, out, j, P, cnt, i, e, tmp;
  O := FGOrbit(G, dw);
  key := Concatenation(G.name, "_orb_", O.key, "_", String(N));
  if IsBound(FG.counts.(key)) then return FG.counts.(key); fi;
  file := Concatenation(FG_CACHE_DIR, "/", key, ".g");
  T := FGTuples(G, N).tuples;
  if IsReadableFile(file) then
    FG_TMP := fail; Read(file);
    if IsList(FG_TMP) and Length(FG_TMP) = Length(T) then   # else unreadable: recompute
      FG.counts.(key) := FG_TMP; return FG_TMP;
    fi;
  fi;
  out := [];
  for j in [1 .. Length(T)] do
    P := O.wts * T[j];
    if not O.integral and not ForAll(P, IsInt) then Error("non-integral exponent"); fi;
    cnt := ListWithIdenticalEntries(N, 0);
    for i in [1 .. Length(P)] do
      e := P[i] mod N + 1;
      cnt[e] := cnt[e] + 1;
    od;
    Add(out, cnt);
  od;
  if Length(O.wts) * Length(T) > 10^6 and IsDirectoryPath(FG_CACHE_DIR) then
    # write to a private temporary file, then rename, so that parallel runs never see a
    # partly written cache file
    tmp := Filename(DirectoryTemporary(), "counts.g");
    PrintTo(tmp, "FG_TMP := ", out, ";\n");
    Exec(Concatenation("mv -f ", tmp, " ", file));
  fi;
  FG.counts.(key) := out;
  return out;
end;

# Exponent counts of module M: for each tuple of exact order N, the multiplicity of zeta_N^k.
FGCounts := function(G, M, N)
  local key, out, i, oc, j;
  key := Concatenation(M.key, "_", String(N));
  if IsBound(FG.counts.(key)) then return FG.counts.(key); fi;
  out := List(FGTuples(G, N).tuples, t -> ListWithIdenticalEntries(N, 0));
  for i in [1 .. Length(M.doms)] do
    oc := FGOrbitCounts(G, M.doms[i], N);
    for j in [1 .. Length(out)] do out[j] := out[j] + M.dmults[i] * oc[j]; od;
  od;
  FG.counts.(key) := out;
  return out;
end;

# Signature classes: tuples with equal exponent counts on all modules in mods.
FGClassData := function(G, mods, N)
  local key, T, C, sigs, keys, cls, j, pos, m, D;
  key := Concatenation(JoinStringsWithSeparator(List(mods, m -> m.key), "+"), "_", String(N));
  if IsBound(FG.classdata.(key)) then return FG.classdata.(key); fi;
  T := FGTuples(G, N);
  C := List(mods, m -> FGCounts(G, m, N));
  sigs := List([1 .. Length(T.tuples)], j -> Concatenation(List(C, c -> c[j])));
  keys := Set(sigs);
  cls := List(keys, k -> rec(ntuples := 0));
  for j in [1 .. Length(sigs)] do
    pos := PositionSorted(keys, sigs[j]);
    if cls[pos].ntuples = 0 then
      cls[pos].tuple := T.tuples[j]; cls[pos].central := T.central[j];
    elif cls[pos].central <> T.central[j] then
      Error("central and non-central elements with the same signature");
    fi;
    cls[pos].ntuples := cls[pos].ntuples + 1;
  od;
  for pos in [1 .. Length(keys)] do
    cls[pos].counts := List([1 .. Length(mods)], m -> keys[pos]{[(m - 1) * N + 1 .. m * N]});
  od;
  D := rec(N := N, keys := keys, cls := cls, pow := List(cls, c -> []), tr := List(mods, m -> []));
  FG.classdata.(key) := D;
  return D;
end;

# Index (in the order N/gcd(N,e) class data) of the e-th power of class idx.
FGPower := function(G, mods, D, idx, e)
  local N, g, N2, f, sig, c, c2, k, D2, pos;
  N := D.N; e := e mod N;
  if IsBound(D.pow[idx][e + 1]) then return D.pow[idx][e + 1]; fi;
  g := Gcd(N, e); N2 := N / g; f := e / g;
  sig := [];
  for c in D.cls[idx].counts do
    c2 := ListWithIdenticalEntries(N2, 0);
    for k in [0 .. N - 1] do
      if c[k + 1] <> 0 then c2[(k * f) mod N2 + 1] := c2[(k * f) mod N2 + 1] + c[k + 1]; fi;
    od;
    Append(sig, c2);
  od;
  D2 := FGClassData(G, mods, N2);
  pos := PositionSorted(D2.keys, sig);
  if pos > Length(D2.keys) or D2.keys[pos] <> sig then Error("power class not found"); fi;
  D.pow[idx][e + 1] := pos;
  return pos;
end;

FGTrace := function(D, idx, m, e)    # trace of t^e on module m, for class idx of D
  local N, c, c2, k, j;
  N := D.N; c := D.cls[idx].counts[m];
  if e mod N = 1 then return CycList(c); fi;
  c2 := ListWithIdenticalEntries(N, 0);
  for k in [1 .. N] do
    if c[k] <> 0 then j := ((k - 1) * e) mod N + 1; c2[j] := c2[j] + c[k]; fi;
  od;
  return CycList(c2);
end;

FGTrace1 := function(D, idx, m)
  if not IsBound(D.tr[m][idx]) then D.tr[m][idx] := CycList(D.cls[idx].counts[m]); fi;
  return D.tr[m][idx];
end;

FGTraceF := function(D, idx, m)     # [Re, Im] of the trace of t, as floats
  local cs, c;
  cs := FGCosSin(D.N); c := D.cls[idx].counts[m];
  return [c * cs[1], c * cs[2]];
end;

# ---------------------------------------------------------------------------------------
# The search

FeasibleCharacters := function(type, rank, tbl, p)
  local G, mode, eps, mods, maxord, nodelimit, modtbl, fus, inv, nC, B, nI, A, ordH, cent,
        pmcache, PMi, dropped, keep, lim, m, nM, degs, dims, Are, Aim, ub, order, seen, steps,
        j, n, orbit, D, st, cands, idx, img, ok, e, k, g, checks, divs, pcs, vid, vidrec,
        VID, key, vals, contrib, cplx, c, lo, hi, s, nS, sufLo, sufHi, assigned, cur, sols,
        solkeys, nodes, t0, leaf, search, aborted, zero, nvid, es, prs, pcls, index,
        s1, allsteps, RunSteps, strategy, small, ndec, npsi, allowed, psis, dfs, psi, i;
  t0 := Runtime();
  G := FGGroup(type, rank);
  mode := ValueOption("PowerMaps"); if mode = fail then mode := "full"; fi;
  eps := ValueOption("Eps"); if eps = fail then eps := 1.e-3; fi;
  mods := ValueOption("Modules"); if mods = fail then mods := FGDefaultModules(G, p); fi;
  maxord := ValueOption("MaxOrder");
  nodelimit := ValueOption("NodeLimit"); if nodelimit = fail then nodelimit := infinity; fi;
  if p = 0 then modtbl := tbl; fus := [1 .. NrConjugacyClasses(tbl)];
  else
    modtbl := tbl mod p;
    if modtbl = fail then Error("no Brauer table"); fi;
    fus := GetFusionMap(modtbl, tbl);
  fi;
  inv := []; for j in [1 .. Length(fus)] do inv[fus[j]] := j; od;
  nC := Length(fus);
  B := List(Irr(modtbl), ValuesOfClassFunction);
  nI := Length(B);
  if p = 0 then
    A := List([1 .. nC], j -> List([1 .. nI], i ->
           SizesConjugacyClasses(tbl)[j] * ComplexConjugate(B[i][j]) / Size(tbl)));
  else
    A := B^-1;
  fi;
  ordH := OrdersClassRepresentatives(tbl){fus};
  if ordH[1] <> 1 then Error("class 1 is not the identity"); fi;
  cent := List(fus, c -> SizesConjugacyClasses(tbl)[c] = 1);
  pmcache := [];
  PMi := function(j, e)
    if not IsBound(pmcache[e]) then pmcache[e] := PowerMap(tbl, e); fi;
    return inv[pmcache[e][fus[j]]];
  end;
  # module order limits
  dropped := []; keep := [];
  for m in [1 .. Length(mods)] do
    if IsList(maxord) then lim := maxord[m];
    elif IsInt(maxord) then lim := maxord;
    elif mods[m].dim <= 100 then lim := infinity;
    else lim := 37; fi;
    if mods[m].isadjoint then lim := infinity; fi;
    if Maximum(ordH) > lim then
      Add(dropped, Concatenation(mods[m].name, " (p-regular elements of order ",
          String(Maximum(ordH)), " > ", String(lim), ")"));
    else Add(keep, mods[m]); fi;
  od;
  mods := keep; nM := Length(mods);
  degs := List(B, b -> b[1]);
  dims := List(mods, m -> m.dim);
  Are := []; Aim := [];
  for j in [1 .. nC] do
    cplx := List(A[j], FGComplex);
    Are[j] := List(cplx, z -> z[1]); Aim[j] := List(cplx, z -> z[2]);
  od;
  ub := Concatenation(List(dims, d -> List(degs, x -> Float(d / x))));
  s1 := Position(dims, Minimum(dims));    # the smallest module (stage 1 of the hybrid strategy)
  zero := ListWithIdenticalEntries(nM * nI, 0.);
  # value ids (divisors mode): an integer per tuple of traces over the modules, via a sorted list
  vidrec := [[], []]; nvid := 0;
  VID := function(v)
    local pos;
    pos := PositionSorted(vidrec[1], v);
    if pos <= Length(vidrec[1]) and vidrec[1][pos] = v then return vidrec[2][pos]; fi;
    nvid := nvid + 1;
    Add(vidrec[1], v, pos); Add(vidrec[2], nvid, pos);
    return nvid;
  end;
  contrib := function(classes, zs)   # zs[c][m] = [Re, Im] of the trace at classes[c] on module m
    local out, m, c, z, part;
    out := [];
    for m in [1 .. nM] do
      part := ListWithIdenticalEntries(nI, 0.);
      for c in [1 .. Length(classes)] do
        z := zs[c][m];
        part := part + z[1] * Are[classes[c]] - z[2] * Aim[classes[c]];
      od;
      Append(out, part);
    od;
    return out;
  end;
  # steps, in increasing element order
  order := [2 .. nC]; SortBy(order, j -> ordH[j]);
  seen := []; steps := [];
  for j in order do
    if j in seen then continue; fi;
    n := ordH[j];
    D := FGClassData(G, mods, n);
    cands := [];
    if mode = "full" then
      orbit := Set(PrimeResidues(n), e -> PMi(j, e));
      prs := Filtered(Set(Factors(n)), p -> p < n);     # prime divisors p with c^p non-trivial
      pcls := List(prs, p -> PMi(j, p));
      es := PrimeResidues(n);
      SortBy(es, e -> PMi(j, e) <> j);    # stabiliser of class j first: most candidates fail there
      for idx in [1 .. Length(D.cls)] do
        if D.cls[idx].central <> cent[j] then continue; fi;
        img := []; ok := true;
        for e in es do
          k := PMi(j, e); g := FGPower(G, mods, D, idx, e);
          if IsBound(img[k]) and img[k] <> g then ok := false; break; fi;
          img[k] := g;
        od;
        if not ok then continue; fi;
        # The classes c^p (p a prime divisor of n, p < n) have lower order and are assigned
        # before this step.  Consistency with them implies consistency with every c^e,
        # gcd(e, n) > 1, since c^e = (c^p)^(e/p) and the earlier assignments are themselves
        # consistent with power maps.  So candidates are indexed by the classes of t^p: during
        # the search only the p-th roots of the classes already assigned are tried.
        Add(cands, rec(assign := img{orbit},
                       key := String(List(prs, p -> FGPower(G, mods, D, idx, p))),
                       trs := List(orbit, k -> FGTrace1(D, img[k], s1)),   # exact traces on module s1
                       contrib := contrib(orbit, List(orbit, k -> List([1 .. nM], m -> FGTraceF(D, img[k], m))))));
      od;
      index := true;     # built in RunSteps
    else
      orbit := [j]; pcls := []; index := fail;
      divs := Filtered(DivisorsInt(n), d -> d > 1);
      pcs := List(divs, d -> PMi(j, n / d));
      keep := [];       # sorted list of power-trace tuples already used
      for idx in [1 .. Length(D.cls)] do
        if D.cls[idx].central <> cent[j] then continue; fi;
        vals := List(divs, d -> List([1 .. nM], m -> FGTrace(D, idx, m, n / d)));
        if vals in keep then continue; fi;
        AddSet(keep, vals);
        checks := List([1 .. Length(divs) - 1], i -> [pcs[i], VID(vals[i])]);
        Add(cands, rec(assign := [VID(vals[Length(divs)])], checks := checks,
                       vals := [vals[Length(divs)]],
                       contrib := contrib([j], [List(vals[Length(divs)], FGComplex)])));
      od;
    fi;
    UniteSet(seen, orbit);
    Add(steps, rec(n := n, orbit := orbit, cands := cands, D := D, pcls := pcls, index := index));
  od;
  SortBy(steps, st -> [st.n, Length(st.cands)]);
  allsteps := steps;
  assigned := []; cur := []; sols := []; solkeys := []; nodes := 0; aborted := false;
  leaf := function()
    local chi, mults, m, s, i, mm;
    mults := [];
    for m in [1 .. nM] do
      chi := []; chi[1] := dims[m];
      for s in [1 .. nS] do
        for i in [1 .. Length(steps[s].orbit)] do
          if mode = "full" then chi[steps[s].orbit[i]] := FGTrace1(steps[s].D, cur[s].assign[i], m);
          else chi[steps[s].orbit[i]] := cur[s].vals[i][m]; fi;
        od;
      od;
      mm := chi * A;
      if not ForAll(mm, x -> IsInt(x) and x >= 0) then return; fi;
      Add(mults, mm);
    od;
    if not mults in solkeys then
      AddSet(solkeys, mults);
      if mode = "full" then     # witness: Kac coordinates (and order) of the image of each class
        Add(sols, rec(mults := mults, fusion := Concatenation([[1, List([1 .. G.rank], i -> 0)]],
          List([2 .. nC], k -> [ordH[k], FGClassData(G, mods, ordH[k]).cls[assigned[k]].tuple]))));
      else
        Add(sols, rec(mults := mults, fusion := ShallowCopy(assigned)));
      fi;
    fi;
  end;
  search := function(s, partial)
    local st, c, np, i, cs, key;
    if s > nS then leaf(); return; fi;
    st := steps[s];
    if st.index <> fail then      # full mode: only the roots of the classes assigned to c^p
      key := String(List(st.pcls, k -> assigned[k]));
      if not IsBound(st.index.(key)) then return; fi;
      cs := st.index.(key);
    else
      cs := st.cands;
    fi;
    for c in cs do
      nodes := nodes + 1;
      if nodes > nodelimit then aborted := true; return; fi;
      if IsBound(c.checks) and not ForAll(c.checks, x -> assigned[x[1]] = x[2]) then continue; fi;
      np := partial + c.contrib;
      if Minimum(np + sufHi[s + 1]) < -eps then continue; fi;
      if Maximum(np + sufLo[s + 1] - ub) > eps then continue; fi;
      for i in [1 .. Length(st.orbit)] do assigned[st.orbit[i]] := c.assign[i]; od;
      cur[s] := c;
      search(s + 1, np);
      if aborted then return; fi;
    od;
    for i in st.orbit do Unbind(assigned[i]); od;
  end;
  # Run the search over the given steps (each with its candidate list): pruning bounds from
  # these candidates, root index (full mode), then depth-first search from the identity class.
  RunSteps := function(stp)
    local s, lo, hi, st, c;
    steps := stp; nS := Length(steps);
    for st in steps do
      if st.index <> fail then
        st.index := rec();
        for c in st.cands do
          if not IsBound(st.index.(c.key)) then st.index.(c.key) := []; fi;
          Add(st.index.(c.key), c);
        od;
      fi;
    od;
    sufLo := []; sufHi := [];
    sufLo[nS + 1] := zero; sufHi[nS + 1] := zero;
    for s in [nS, nS - 1 .. 1] do
      if steps[s].cands = [] then
        lo := zero; hi := zero;
      else
        lo := List([1 .. nM * nI], r -> Minimum(List(steps[s].cands, c -> c.contrib[r])));
        hi := List([1 .. nM * nI], r -> Maximum(List(steps[s].cands, c -> c.contrib[r])));
      fi;
      sufLo[s] := sufLo[s + 1] + lo; sufHi[s] := sufHi[s + 1] + hi;
    od;
    if ForAny(steps, st -> st.cands = []) then return; fi;
    search(1, Concatenation(List(dims, d -> Float(d) * Are[1])));
  end;
  # Strategy.  "value": search over class assignments directly.  "hybrid": first enumerate the
  # characters psi of H on the smallest module (degree decompositions whose value at each class
  # is the trace of some candidate class of G), then search with each class restricted to the
  # classes of G whose trace on that module is psi(c).  The union over psi is the full search.
  # "auto": hybrid when the number of degree decompositions of the smallest module is small.
  strategy := ValueOption("Strategy"); if strategy = fail then strategy := "auto"; fi;
  if mode <> "full" then strategy := "value"; fi;
  small := Filtered([1 .. nI], i -> degs[i] <= dims[s1]);
  if strategy = "auto" then
    ndec := ListWithIdenticalEntries(dims[s1] + 1, 0); ndec[1] := 1;   # degree decompositions
    for i in small do
      for k in [degs[i] .. dims[s1]] do ndec[k + 1] := ndec[k + 1] + ndec[k - degs[i] + 1]; od;
    od;
    if ndec[dims[s1] + 1] <= 200000 then strategy := "hybrid"; else strategy := "value"; fi;
  fi;
  npsi := fail;
  if strategy = "value" then
    RunSteps(allsteps);
  else
    # stage 1: characters on module s1
    allowed := [];
    for st in allsteps do
      for i in [1 .. Length(st.orbit)] do allowed[st.orbit[i]] := Set(st.cands, c -> c.trs[i]); od;
    od;
    SortBy(small, i -> -degs[i]);
    psis := [];
    dfs := function(k, mult, rem)
      local i, mm, psi;
      if rem = 0 then
        psi := List([1 .. nC], c -> Sum(small, i -> mult[i] * B[i][c]));
        if ForAll([2 .. nC], c -> psi[c] in allowed[c]) then Add(psis, psi); fi;
        return;
      fi;
      if k > Length(small) then return; fi;
      i := small[k];
      for mm in [0 .. QuoInt(rem, degs[i])] do
        mult[i] := mm;
        dfs(k + 1, mult, rem - mm * degs[i]);
      od;
      mult[i] := 0;
    end;
    dfs(1, ListWithIdenticalEntries(nI, 0), dims[s1]);
    npsi := Length(psis);
    # stage 2: one restricted search per psi
    for psi in psis do
      RunSteps(List(allsteps, st -> rec(n := st.n, orbit := st.orbit, D := st.D, pcls := st.pcls,
        index := true, cands := Filtered(st.cands, c -> ForAll([1 .. Length(st.orbit)], i -> c.trs[i] = psi[st.orbit[i]])))));
      if aborted then break; fi;
    od;
  fi;
  return rec(type := type, rank := rank, group := Identifier(tbl), p := p, mode := mode,
             strategy := strategy, npsi := npsi,
             modules := List(mods, m -> rec(name := m.name, kind := m.kind, hw := m.hw, dim := m.dim)),
             degs := degs, sols := sols, dropped := dropped, nodes := nodes, aborted := aborted,
             steps := List(allsteps, st -> [st.n, Length(st.orbit), Length(st.cands)]),
             time := (Runtime() - t0) / 1000.);
end;

# ---------------------------------------------------------------------------------------
# Output helpers

FGProject := function(res, j)       # distinct multiplicity vectors on module j
  return Set(res.sols, s -> s.mults[j]);
end;

FGModuleIndex := function(res, name)
  return PositionProperty(res.modules, m -> m.name = name);
end;

# Joint profile: sorted multiset of [deg, m_1, ..., m_k] over irreducibles with a non-zero entry.
FGJointProfile := function(degs, mlist)
  return SortedList(Filtered(List([1 .. Length(degs)], i ->
           Concatenation([degs[i]], List(mlist, m -> m[i]))), x -> ForAny(x{[2 .. Length(x)]}, y -> y <> 0)));
end;

FGIrrLabels := function(degs)   # "8a", "8b", ... in GAP order; plain degree if unique
  local out, i, d, k;
  out := [];
  for i in [1 .. Length(degs)] do
    d := degs[i];
    if Number(degs, x -> x = d) = 1 then out[i] := String(d);
    else
      k := Number([1 .. i], j -> degs[j] = d);
      out[i] := Concatenation(String(d), [CHARS_LALPHA[k]]);
    fi;
  od;
  return out;
end;

# Permutations of the irreducibles induced by the table automorphisms of tbl (mod p).
FGIrrPerms := function(tbl, p)
  local modtbl, irr, aut, out, g, perm;
  if p = 0 then modtbl := tbl; else modtbl := tbl mod p; fi;
  irr := List(Irr(modtbl), ValuesOfClassFunction);
  aut := AutomorphismsOfTable(modtbl);
  out := [];
  for g in Elements(aut) do
    perm := List(irr, x -> Position(irr, Permuted(x, g)));
    Add(out, perm);
  od;
  return out;
end;

# Memoir-style table on the modules with indices mods (default: all), condensed under table
# automorphisms (orbit representatives = lexicographically largest images).
FGPrintTable := function(res, tbl)
  local mods, perms, labs, rows, canon, s, row, i, m, line, w;
  mods := ValueOption("Modules");
  if mods = fail then mods := [1 .. Length(res.modules)]; fi;
  perms := FGIrrPerms(tbl, res.p);
  labs := FGIrrLabels(res.degs);
  canon := function(ml)
    return Maximum(List(perms, pp -> Concatenation(List(ml, v -> List([1 .. Length(v)], i -> v[pp[i]])))));
  end;
  rows := Set(res.sols, s -> canon(s.mults{mods}));
  rows := Reversed(rows);
  w := Maximum(List(labs, Length)) + 1;
  Print(res.group, " < ", res.type, res.rank, ", p = ", res.p, ": ", Length(rows),
        " row(s) up to table automorphisms (", Length(res.sols), " in all)\n");
  line := "     ";
  for m in mods do
    Append(line, Concatenation("| ", res.modules[m].name, " "));
    Append(line, ListWithIdenticalEntries(Maximum(0, w * Length(labs) - Length(res.modules[m].name) - 1), ' '));
  od;
  Print(line, "\n");
  line := "     ";
  for m in mods do
    Append(line, "| ");
    for i in [1 .. Length(labs)] do Append(line, String(labs[i], w)); od;
  od;
  Print(line, "\n");
  for i in [1 .. Length(rows)] do
    line := String(Concatenation(String(i), ")"), -5);
    for m in [1 .. Length(mods)] do
      Append(line, "| ");
      for s in [1 .. Length(labs)] do
        Append(line, String(rows[i][(m - 1) * Length(labs) + s], w));
      od;
    od;
    Print(line, "\n");
  od;
  if res.dropped <> [] then Print("  modules not used: ", res.dropped, "\n"); fi;
end;
