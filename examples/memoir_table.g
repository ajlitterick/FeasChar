# Reproduce a table of the Memoir: Alt6 < F4, p = 0 (Memoir Table 6.5).
# Run from the repository folder:   gap examples/memoir_table.g
Read("feasgen.g");

res := FeasibleCharacters("F", 4, CharacterTable("A6"), 0);
FGPrintTable(res, CharacterTable("A6"));
# Six rows, condensed under the automorphisms of the character table, as in the Memoir.
QUIT;
