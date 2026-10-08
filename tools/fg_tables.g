# Updated version of Memoir table MEMOIR[IDX] (memoir_index.g), as LaTeX, plus a one-line
# summary of the differences.  Writes out/tables/<number>.tex and out/tables/<number>.txt.
#
# Rows: feasible characters on L(G) (and V_min for F4, E6, E7) that are compatible with ALL
# default modules (adjoint and end-node Weyl modules; in bad characteristic also L(lambda)),
# with full power maps and all element orders; listed up to automorphisms of the character table
# (for (F4, 2) also up to the graph automorphism of F4).  Flags P/N are carried over from the
# Memoir row with the same composition factors; rows not in the Memoir are marked "new".
Read("feasgen.g");
Read("memoir_index.g");
T := MEMOIR[IDX];
OUTD := "out/tables";
base := Concatenation(OUTD, "/", T.number);
typ := T.type{[1]}; rank := Int(T.type{[2]});
isE8 := T.type = "E8";
iscover := IsDigitChar(T.group[1]) and Length(T.group) > 1 and T.group[2] = '.';
tau := (T.type = "F4" and T.p = 2);
WriteText := function(name, str)   # write without GAP's line wrapping
  local out;
  out := OutputTextFile(name, false);
  SetPrintFormattingStatus(out, false);
  WriteAll(out, str);
  CloseStream(out);
end;
FGToJSON := function(x)   # minimal JSON writer: records, lists, strings, integers, booleans, fail
  local s, c, nm;
  if x = true then return "true"; elif x = false then return "false"; elif x = fail then return "null";
  elif IsInt(x) then return String(x);
  elif IsStringRep(x) or (IsString(x) and x <> []) then
    s := "\"";
    for c in x do
      if c = '\\' then Append(s, "\\\\"); elif c = '"' then Append(s, "\\\"");
      elif c = '\n' then Append(s, "\\n"); else Add(s, c); fi;
    od;
    Add(s, '"'); return s;
  elif IsRecord(x) then
    return Concatenation("{", JoinStringsWithSeparator(List(RecNames(x), nm ->
             Concatenation(FGToJSON(nm), ": ", FGToJSON(x.(nm)))), ", "), "}");
  elif IsList(x) then
    return Concatenation("[", JoinStringsWithSeparator(List(x, FGToJSON), ", "), "]");
  fi;
  Error("FGToJSON: cannot serialise ", x);
end;
Summary := function(status, rest)
  WriteText(Concatenation(base, ".txt"), Concatenation(T.number, "\t", T.type, "\t", T.group, "\t",
            String(T.p), "\t", status, "\t", rest, "\n"));
end;
Stub := function(reason)   # table not recomputed: keep the caption, say why
  WriteText(Concatenation(base, ".tex"), Concatenation(
    "\\setcounter{table}{", String(Int(SplitString(T.number, ".")[2]) - 1), "}\n",
    "\\begin{longtable}{l}\n\\caption{", T.caption, "}\\\\\n",
    "Not recomputed: ", reason, ". See the Memoir.\n\\end{longtable}\n\n"));
  Summary("NOT RECOMPUTED", reason);
end;
tbl := fail;
# GAP names where the Memoir's name is ambiguous: of the two triple covers of U4(3), only 3_1
# has faithful 27-dimensional modules (as the Memoir notes under Table 6.104)
GAPNAME := rec(); GAPNAME.("3.U4(3)") := "3_1.U4(3)";
if IsBound(GAPNAME.(T.group)) then T.group := GAPNAME.(T.group); fi;
if T.group <> "An" then tbl := CharacterTable(T.group); fi;
if tbl = fail then Stub("no character table available in GAP"); QuitGap(); fi;
if T.p > 0 and (tbl mod T.p) = fail then Stub("no Brauer table available in GAP"); QuitGap(); fi;

DegOf := function(lab)
  local i;
  i := 1;
  while i <= Length(lab) and IsDigitChar(lab[i]) do i := i + 1; od;
  return Int(lab{[1 .. i - 1]});
end;
nG := Length(T.mods);
# Memoir profile of a row: sorted [deg, m_1, ..., m_nG] over the irreducibles (labels such as
# "17_a,b,c" count three times; for covers the two column groups are disjoint).
MemoirProf := function(row)
  local d, g, i, lab, k, lk, v;
  d := rec();
  for g in [1 .. nG] do
    for i in [1 .. Length(T.mods[g].labels)] do
      lab := T.mods[g].labels[i];
      if iscover then lab := Concatenation(String(g), ":", lab); fi;
      for k in [1 .. 1 + Number(lab, ch -> ch = ',')] do
        lk := Concatenation(lab, "#", String(k));
        if not IsBound(d.(lk)) then
          d.(lk) := Concatenation([DegOf(T.mods[g].labels[i])], ListWithIdenticalEntries(nG, 0));
        fi;
        d.(lk)[g + 1] := d.(lk)[g + 1] + row.mult[g][i];
      od;
    od;
  od;
  v := List(RecNames(d), l -> d.(l));
  # published erratum (FeasChar README): in Table 6.298 each row has a 29-dimensional
  # constituent instead of one of its 28-dimensional constituents
  if T.number = "6.298" then
    i := PositionProperty(v, x -> x[1] = 28 and x[2] > 0);
    v := List(v, ShallowCopy);
    v[i][2] := v[i][2] - 1;
    Add(v, [29, 1]);
  fi;
  return v;
end;
SplitP := function(prof)   # per-module profiles (insensitive to pairing labels across groups)
  return List([1 .. nG], g -> SortedList(Filtered(List(prof, x -> [x[1], x[g + 1]]), y -> y[2] <> 0)));
end;
Canon := function(prof)
  local sw;
  prof := SortedList(Filtered(prof, x -> ForAny(x{[2 .. Length(x)]}, y -> y <> 0)));
  if not tau then return prof; fi;
  sw := SortedList(Filtered(List(prof, x -> [x[1], x[2], x[2] - x[3]]), x -> x[2] <> 0 or x[3] <> 0));
  return Maximum(prof, sw);
end;

# --- compute
G := FGGroup(typ, rank);
VA := FGModule(G, G.thetaFW);
mods := FGDefaultModules(G, T.p);
iA := 1;
if isE8 then iV := fail;
else
  dimV := Int(T.mods[2].name{[2 .. Length(T.mods[2].name)]});
  iV := PositionProperty(mods, m -> m.kind = "V" and m.dim = dimV);
fi;
t0 := Runtime();
res := FeasibleCharacters(typ, rank, tbl, T.p : Modules := mods);
proj := function(r, iv) if isE8 then return Set(r.sols, s -> [s.mults[1]]);
                        else return Set(r.sols, s -> [s.mults[1], s.mults[iv]]); fi; end;
# position of V_min among the modules actually used (large modules may have been dropped)
if isE8 then iVr := fail; else iVr := Position(List(res.modules, m -> m.name), mods[iV].name); fi;
rows := proj(res, iVr);
# the Memoir's own modules, full power maps: to give reasons for removed rows
if isE8 then base2 := [VA]; else base2 := [mods[1], mods[iV]]; fi;
res0 := FeasibleCharacters(typ, rank, tbl, T.p : Modules := base2);
rows0 := proj(res0, 2);
secs := (Runtime() - t0) / 1000.;
degs := res.degs;
ProfOf := function(r)   # r = [mA] or [mA, mV]
  return List([1 .. Length(degs)], i -> Concatenation([degs[i]], List(r, m -> m[i])));
end;
mp := List(T.rows, r -> Canon(MemoirProf(r)));
msp := List(T.rows, r -> SplitP(Canon(MemoirProf(r))));

# --- condense under table automorphisms (and tau)
perms := FGIrrPerms(tbl, T.p);
CanonRow := function(r)
  local imgs;
  imgs := List(perms, pp -> List(r, m -> List([1 .. Length(m)], i -> m[pp[i]])));
  if tau then Append(imgs, List(perms, pp -> List([r[1], r[1] - r[2]], m -> List([1 .. Length(m)], i -> m[pp[i]])))); fi;
  return Maximum(imgs);
end;
reps := Reversed(Set(rows, CanonRow));
# --- match to Memoir rows
matched := [];      # Memoir row indices reproduced
info := [];
for r in reps do
  pr := Canon(ProfOf(r));
  ks := Filtered([1 .. Length(mp)], k -> mp[k] = pr);
  if ks = [] then ks := Filtered([1 .. Length(msp)], k -> msp[k] = SplitP(pr)); fi;
  UniteSet(matched, ks);
  Add(info, rec(row := r, memoir := ks,
     P := ForAny(ks, k -> T.rows[k].P), N := ForAny(ks, k -> T.rows[k].N)));
od;
# list rows in the Memoir's order; new rows last
SortBy(info, x -> [Minimum(Concatenation(x.memoir, [infinity])), Position(reps, x.row)]);
# also count GAP rows that are automorphic images of reproduced ones (profile matches)
allprofs := Set(rows, r -> Canon(ProfOf(r)));
allsplit := Set(allprofs, SplitP);
removed := Filtered([1 .. Length(mp)], k -> not mp[k] in allprofs and not msp[k] in allsplit);
profs0 := Set(rows0, r -> Canon(ProfOf(r)));
# for each extra module M: the profiles feasible on the Memoir's modules together with M alone
single := rec();
SingleProfs := function(M)
  local r;
  if not IsBound(single.(M.name)) then
    r := FeasibleCharacters(typ, rank, tbl, T.p : Modules := Concatenation(base2, [M]));
    single.(M.name) := Set(proj(r, 2), x -> Canon(ProfOf(x)));
  fi;
  return single.(M.name);
end;
reason := function(k)
  local extra, killers, inS;
  if mp[k] in profs0 or msp[k] in Set(profs0, SplitP) then
    extra := Filtered(mods, m -> not m.name in List(base2, b -> b.name) and
                                 ForAny(res.modules, x -> x.name = m.name));
    inS := function(M) local S; S := SingleProfs(M); return mp[k] in S or msp[k] in Set(S, SplitP); end;
    killers := Filtered(extra, M -> not inS(M));
    if killers = [] then   # removed only by a combination of modules
      return Concatenation("not compatible with ", JoinStringsWithSeparator(List(extra, m -> m.name), " + "), " together");
    fi;
    return Concatenation("not compatible with ", JoinStringsWithSeparator(List(killers, m -> m.name), "; nor with "));
  fi;
  return "not feasible on the Memoir's modules (all element orders, full power maps)";
end;

# --- columns: irreducibles occurring in some representative row, per module group
groups := [[iA, Concatenation("V_{", String(VA.dim), "}")]];
if not isE8 then Add(groups, [2, Concatenation("V_{", String(dimV), "}")]); fi;
labs := FGIrrLabels(degs);
TexLab := function(l)
  local i;
  i := 1;
  while i <= Length(l) and IsDigitChar(l[i]) do i := i + 1; od;
  if i > Length(l) then return l; fi;
  return Concatenation(l{[1 .. i - 1]}, "_{", l{[i .. Length(l)]}, "}");
end;
cols := List([1 .. Length(groups)], g -> Filtered([1 .. Length(degs)], i -> ForAny(reps, r -> r[g][i] <> 0)));
spec := Concatenation("rc", Concatenation(List(cols, c -> Concatenation("|", Concatenation(List(c, x -> "c"))))));
s := "";
Append(s, Concatenation("\\setcounter{table}{", String(Int(SplitString(T.number, ".")[2]) - 1), "}\n"));
Append(s, Concatenation("\\begin{longtable}{", spec, "}\n\\caption{", T.caption, "}\\\\\n"));
hdr := "& ";
for g in [1 .. Length(groups)] do
  Append(hdr, Concatenation("& \\multicolumn{", String(Length(cols[g])), "}{c", ListWithIdenticalEntries(Minimum(1, Length(groups) - g), '|'), "}{$", groups[g][2], "$} "));
od;
Append(hdr, "\\\\\n& ");
for g in [1 .. Length(groups)] do
  for i in cols[g] do Append(hdr, Concatenation("& $", TexLab(labs[i]), "$ ")); od;
od;
Append(hdr, "\\\\ \\hline\n");
Append(s, Concatenation(hdr, "\\endfirsthead\n", hdr, "\\endhead\n"));
for k in [1 .. Length(info)] do
  x := info[k];
  fl := [];
  if x.P then Add(fl, "\\possprim"); fi;
  if x.N then Add(fl, "\\nongcr"); fi;
  if x.memoir = [] then Add(fl, "\\newrow"); fi;
  Append(s, Concatenation(String(k), ") & ", JoinStringsWithSeparator(fl, ", ")));
  for g in [1 .. Length(groups)] do
    for i in cols[g] do Append(s, Concatenation(" & ", String(x.row[g][i]))); od;
  od;
  Append(s, " \\\\\n");
od;
if reps = [] then Append(s, "\\multicolumn{2}{l}{(no feasible characters)} \\\\\n"); fi;
Append(s, "\\end{longtable}\n");
notes := [];
if T.number = "6.298" then
  Add(notes, "Memoir rows read with the published correction: one constituent 28 replaced by 29 in each row.");
fi;
if removed <> [] then
  Add(notes, Concatenation("Memoir rows not reproduced: ", JoinStringsWithSeparator(List(removed, k ->
      Concatenation(String(T.rows[k].n), ") (", reason(k), ")")), "; "), "."));
fi;
if ForAny(info, x -> x.memoir = []) then Add(notes, "Rows marked \\newrow\\ do not appear in the Memoir."); fi;
if res.dropped <> [] then
  Add(notes, Concatenation("Modules not used (elements of large order): ",
      JoinStringsWithSeparator(List(res.dropped, d -> SplitString(d, " ")[1]), ", "), "."));
fi;
Add(notes, Concatenation("Modules tested: ", JoinStringsWithSeparator(List(res.modules, m -> m.name), ", "), "."));
Append(s, "\\noindent ");
for n in notes do Append(s, Concatenation(n, "\\\\\n")); od;
Append(s, "\n");
WriteText(Concatenation(base, ".tex"), s);
status := "MATCH";
if removed <> [] or ForAny(info, x -> x.memoir = []) then status := "DIFFERS"; fi;
Summary(status, Concatenation("rows ", String(Length(reps)), "; memoir rows ", String(Length(T.rows)),
  "; new ", String(Number(info, x -> x.memoir = [])), "; removed ",
  String(List(removed, k -> [T.rows[k].n, reason(k)])), "; modules ",
  String(List(res.modules, m -> m.name)), "; dropped ", String(res.dropped), "; ", String(secs), " s"));

# --- machine-readable record (JSON), with a witness class assignment for each row
modsused := Filtered(mods, m -> ForAny(res.modules, x -> x.name = m.name));
bperm := FGBourbakiPerm(G);
ModRec := function(m)
  local r;
  r := rec(name := m.name, kind := m.kind, dim := m.dim, highest_weight := m.hw,
           dominant_weights := m.doms, multiplicities := m.dmults);
  if m.kind = "L" then r.characteristic := m.p; fi;
  if bperm <> fail then r.highest_weight_bourbaki := m.hw{bperm}; fi;
  return r;
end;
if T.p = 0 then modtbl := tbl; fusp := [1 .. NrConjugacyClasses(tbl)];
else modtbl := tbl mod T.p; fusp := GetFusionMap(modtbl, tbl); fi;
witkeys := []; witsols := [];   # canonical display row -> a solution in its orbit
for x in res.sols do
  if isE8 then key := CanonRow([x.mults[1]]); else key := CanonRow([x.mults[1], x.mults[iVr]]); fi;
  if not key in witkeys then Add(witkeys, key); Add(witsols, x); fi;
od;
modnames := List(res.modules, m -> m.name);
rowsj := [];
for k in [1 .. Length(info)] do
  x := info[k];
  w := witsols[Position(witkeys, x.row)];
  r := rec(row := k, flags := Concatenation(List(Filtered([[x.P, "P"], [x.N, "N"], [x.memoir = [], "new"]],
               y -> y[1]), y -> y[2])),
           memoir_rows := List(x.memoir, j -> T.rows[j].n),
           display := rec(),
           witness := rec(multiplicities := rec(),
             classes := List([1 .. Length(fusp)], j -> rec(class := ClassNames(tbl)[fusp[j]],
                          order := w.fusion[j][1], kac := w.fusion[j][2]))));
  r.display.(modnames[1]) := x.row[1];
  if not isE8 then r.display.(modnames[iVr]) := x.row[2]; fi;
  for j in [1 .. Length(modnames)] do r.witness.multiplicities.(modnames[j]) := w.mults[j]; od;
  Add(rowsj, r);
od;
J := rec(number := T.number, caption := T.caption, type := T.type, group := T.group,
  characteristic := T.p, status := status,
  character_table := rec(library := "GAP CTblLib", version := InstalledPackageVersion("ctbllib"),
    name := T.group, brauer := T.p > 0,
    note := "irreducible i is Irr(CharacterTable(name) mod characteristic)[i] (Irr(CharacterTable(name)) if characteristic 0)"),
  irreducibles := List([1 .. Length(degs)], i -> rec(index := i, degree := degs[i], label := labs[i])),
  modules_tested := List(modsused, ModRec),
  modules_not_used := res.dropped,
  display_modules := Concatenation([modnames[1]], List(Filtered([iVr], i -> i <> fail), i -> modnames[i])),
  rows := rowsj,
  memoir_rows_removed := List(removed, k -> rec(memoir_row := T.rows[k].n, reason := reason(k))),
  conventions := "rows up to automorphisms of the character table (and the F4 graph automorphism for p = 2); display rows are images of the witness row under such an automorphism");
WriteText(Concatenation(base, ".json"), FGToJSON(J));
QUIT;
