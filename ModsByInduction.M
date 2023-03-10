PruneByIsomorphism := function(MODULES : debug := false)
  // Intermediate routine. Input a sequence of absolutely irreducible modules I,
  // output one representative for each isomorphism type of module in I.
  //if debug then printf "Pruning..."; end if;
  for i in [1..#MODULES] do
    for j in [(i+1)..#MODULES] do
      if IsDefined(MODULES,i) and IsDefined(MODULES,j) then
        X, Y := ModulesOverCommonField(MODULES[i],MODULES[j]);
        if IsIsomorphic(X,Y) then MODULES := Undefine(MODULES,j); end if;
      end if;
    end for;
  end for;
  //if debug then printf "Done.\n"; end if;

  return Sort(	[MODULES[i] : i in [1..#MODULES] | IsDefined(MODULES,i)],
  				func<x,y | Dimension(x) - Dimension(y)>);

end function;



ModsByInduction := function(G,CHAR : DimLim := 0, SUBGRP_SIZE := 1, debug := false)
/*  Input:
    - a group G for which we can calculate maximal subgroups,
    - a characteristic CHAR.
    - optionally, a parameter DimLim
    - optionally, a parameter SUBGRP_SIZE which determines the point at which we
      use pre-existing routines to calculate modules for a subgroup

    Output: a list of absolutely irreducible modules for G in characteristic CHAR,
    of dimension at most DimLim (if specified)

    Example:
> G := Alt(17);
> time ModsByInduction(G,2 : DimLim := 248);
[
    GModule of dimension 1 over GF(2),
    GModule of dimension 16 over GF(2),
    GModule of dimension 118 over GF(2),
    GModule of dimension 128 over GF(2),
    GModule of dimension 128 over GF(2)
]
Time: 59.440
*/

// Initialise chain of subgroups for induction
H := [G];
if debug then ORDS := [#G]; end if;

// if no DimLim is set then make it large enough not to matter
if DimLim eq 0 then DimLim := #G; end if;

// Calculate a chain of maximal subgroups from G to a p-subgroup of G.
while (#(H[#H]) gt SUBGRP_SIZE and CHAR gt 0) or (not IsSoluble(H[#H]) and CHAR eq 0) do
  H[#H + 1] := Sort(MaximalSubgroups(H[#H]),func<x,y | y`order - x`order>)[1]`subgroup;
  if PrimeDivisors(#H[#H]) in { [], [CHAR]} then break; end if;
  if debug then ORDS[#ORDS + 1] := [#(H[#H])]; end if;
end while;

if debug then
  printf "Subgroup chain length: %o.\n", #H;
  printf "Indices: %o\n",[ ORDS[i]/ORDS[i+1] : i in [1..#ORDS-1] ];
end if;

// Base case
if CHAR gt 0 then
	if PrimeDivisors(#H[#H]) in { [], [CHAR] }  then
  		CURRENT_MODULES := [TrivialModule(H[#H],GF(CHAR))];
	else
  		CURRENT_MODULES := [x : x in AbsolutelyIrreducibleModules(H[#H],GF(CHAR)) | Dimension(x) le DimLim];
	end if;
else
	exponent := Exponent(G);
	CURRENT_MODULES := [x : x in IrreducibleModules(H[#H],CyclotomicField(Exponent(G) : Sparse := true)) | Dimension(x) le DimLim];
end if;


// Loop from H[#H] = Base subgroup,  up to H[1] = G, inducing
// modules then pruning out redundant isomorphism types.

if CHAR gt 0 then
	while #H gt 1 do
  
	H := Exclude(H,H[#H]);

	NEXT_MODULES := [];
	for i in CURRENT_MODULES do
  		NEXT_MODULES := NEXT_MODULES cat &cat[AbsolutelyIrreducibleConstituents(x) : x in Constituents(Induction(i,H[#H]))];
    	if debug then printf "."; end if;  
	end for;

// We only want I to contain one of each isomorphism type of
// abs. irreducible module.
	CURRENT_MODULES := [AbsoluteModuleOverMinimalField(i) : i in NEXT_MODULES];  
  	CURRENT_MODULES := PruneByIsomorphism([i : i in CURRENT_MODULES | Dimension(i) le DimLim]);

	if debug then printf "\nLoop done, %o remaining. Next loop: %o modules. Max dim: %o.\n",
		#H - 1, #CURRENT_MODULES, Dimension(CURRENT_MODULES[#CURRENT_MODULES]); end if;

	end while;
else
	while #H gt 1 do
  
	H := Exclude(H,H[#H]);

	NEXT_MODULES := [];
	for i in CURRENT_MODULES do
  		NEXT_MODULES := NEXT_MODULES cat &cat[CompositionFactors(Induction(i,H[#H]))];
    if debug then printf "."; end if;  
	end for;

	CURRENT_MODULES := PruneByIsomorphism([i : i in CURRENT_MODULES | Dimension(i) le DimLim]);

	if debug then printf "\nLoop done, %o remaining. Next loop: %o modules. Max dim: %o.\n",
		#H - 1, #CURRENT_MODULES, Dimension(CURRENT_MODULES[#CURRENT_MODULES]); end if;

	end while;
end if;

if debug then "Reducing to minimal field..."; end if;

return [AbsoluteModuleOverMinimalField(i) : i in CURRENT_MODULES];

end function;
