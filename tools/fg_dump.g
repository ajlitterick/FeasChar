# Validation (a): dump the (V_min, adjoint) power tuples of all elements of G_sc of exact order
# n = 2..17 (adjoint only for E8), as Python data, for comparison with FeasChar's n.M files.
Read("feasgen.g");
TraceCnt := function(cnt, n, e)
  return Sum([1 .. n], k -> cnt[k] * E(n)^(((k - 1) * e) mod n));
end;
FloatStr := function(x)
  local z;
  z := FGComplex(x);
  return Concatenation("(", String(Round(z[1] * 10^6) / 10^6), ",", String(Round(z[2] * 10^6) / 10^6), ")");
end;
out := "classdata.py";
PrintTo(out, "data = {}\n");
for t in [["G",2],["F",4],["E",6],["E",7],["E",8]] do
  G := FGGroup(t[1], t[2]);
  adj := FGModule(G, G.thetaFW);
  vmin := List(Filtered(G.ends, i -> ForAny([1 .. G.rank], j -> j <> i)), function(i)
            local hw; hw := List([1 .. G.rank], j -> 0); hw[i] := 1; return FGModule(G, hw); end);
  vmin := Filtered(vmin, m -> not m.isadjoint);
  SortBy(vmin, m -> m.dim);
  if G.name = "E8" then mods := [adj]; else mods := [vmin[1], adj]; fi;
  Print(G.name, ": modules ", List(mods, m -> m.name), "\n");
  for n in [2 .. 17] do
    T := FGTuples(G, n).tuples;
    C := List(mods, m -> FGCounts(G, m, n));
    divs := Filtered(DivisorsInt(n), d -> d > 1);
    AppendTo(out, "data[('", G.name, "',", n, ")] = [\n");
    for j in [1 .. Length(T)] do
      AppendTo(out, " (", JoinStringsWithSeparator(List(C, c -> Concatenation("(",
         JoinStringsWithSeparator(List(divs, d -> FloatStr(TraceCnt(c[j], n, n / d))), ","), ",)")), ","), ",),\n");
    od;
    AppendTo(out, "]\n");
    Print("  n = ", n, ": ", Length(T), " elements\n");
  od;
od;
QUIT;
