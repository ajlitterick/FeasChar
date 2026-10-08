# Feasible characters for classical groups: the same code works for any Dynkin type.
# The default modules are the adjoint module and V(omega_i) for each end node i
# (for type A_n: the natural module and its dual are a dual pair, so only one is used).
# Run from the repository folder:   gap examples/classical.g
Read("feasgen.g");

# SL(2,5) = 2.A5 in SL_2: the centre must act as -1 on the natural module.
tbl := CharacterTable("2.A5");
res := FeasibleCharacters("A", 1, tbl, 0);
FGPrintTable(res, tbl);

# A5 in SL_3, ordinary characters and in characteristic 2.
tbl := CharacterTable("A5");
FGPrintTable(FeasibleCharacters("A", 2, tbl, 0), tbl);
FGPrintTable(FeasibleCharacters("A", 2, tbl, 2), tbl);

# L2(7) in Sp_6 = C_3 (modules: adjoint 21, natural 6, and V(omega_3) of dimension 14).
tbl := CharacterTable("L2(7)");
FGPrintTable(FeasibleCharacters("C", 3, tbl, 0), tbl);
QUIT;
