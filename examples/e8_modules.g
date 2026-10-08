# E8: the adjoint module together with all end-node modules V(omega_1) (dim 3875) and
# V(omega_2) (dim 147250); in characteristic 2 also the irreducible modules L(omega_1) (3626)
# and L(omega_2) (143376).  Run from the repository folder:   gap examples/e8_modules.g
Read("feasgen.g");

tbl := CharacterTable("M11");
res := FeasibleCharacters("E", 8, tbl, 2);
Print(List(res.modules, m -> m.name), "\n");
FGPrintTable(res, tbl : Modules := [1]);      # the adjoint module only, as in the Memoir

# Alt8, p = 0: V(omega_2) removes one character that the Memoir lists (Table 6.243).
tbl := CharacterTable("A8");
res := FeasibleCharacters("E", 8, tbl, 0);
FGPrintTable(res, tbl : Modules := [1]);

# FeasChar's original test (powers by divisors of the element order,
# adjoint module only), for comparison.
G := FGGroup("E", 8);
res := FeasibleCharacters("E", 8, tbl, 0 : Modules := [FGModule(G, G.thetaFW)], PowerMaps := "divisors");
Print("adjoint only, divisors test: ", Length(res.sols), " characters\n");
QUIT;
