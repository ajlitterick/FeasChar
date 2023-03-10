FeasChar := procedure(G,LIE_TYPE,CHAR_K : debug := false, MODULES := [], ADJOINT_MODULES := [], MINIMAL_MODULES := [], LIMITING_ORDER := 17, ELEMENT_ORDERS := []);

ASSIGNED_LETTERS := [<"G2","G2_ELTS",7,14>,<"F4","F4_ELTS",26,52>,<"E6","E6_ELTS",27,78>,<"E7","E7_ELTS",56,133>,<"E8","E8_ELTS",0,248>];

LIE_TYPE := CartanName(RootDatum(LIE_TYPE));

if not LIE_TYPE in [i[1] : i in ASSIGNED_LETTERS] then
	"Invalid Lie type!";
	return;
end if;

if not ((IntegerRing()!CHAR_K eq IntegerRing()!0) or (IsPrime(IntegerRing()!CHAR_K))) then
	"Invalid characteristic! Must be zero or prime.";
	return;
end if;

DIM_MIN := [i[3] : i in ASSIGNED_LETTERS | LIE_TYPE eq i[1]][1];
DIM_AD := [i[4] : i in ASSIGNED_LETTERS | LIE_TYPE eq i[1]][1];

if CHAR_K eq 0 then
	CHAR_K := 2;
	repeat
		CHAR_K := NextPrime(CHAR_K);
	until not (CHAR_K in PrimeDivisors(#G));
end if;

if debug then printf "Calculating conjugacy classes and power map of G... "; end if;
Clss := Classes(G);
p := PowerMap(G);


if (IsEmpty(ADJOINT_MODULES) xor IsEmpty(MINIMAL_MODULES)) then
	"Error: Either specify MODULES, or both ADJOINT_MODULES and MINIMAL_MODULES.";
	return;
end if;

if IsEmpty(ADJOINT_MODULES) and IsEmpty(MINIMAL_MODULES) then
	if IsEmpty(MODULES) then
		if debug then printf "Done.\nCalculating absolutely irreducible G-modules... "; end if;
		MODULES := AbsolutelyIrreducibleModules(G,GF(CHAR_K));
	end if;
	ADJOINT_MODULES := MODULES;
	MINIMAL_MODULES := MODULES;
	if debug then printf " Done.\n"; end if;
end if;

DIMS_AD := [Dimension(i) : i in ADJOINT_MODULES | Dimension(i) le DIM_AD];
DIMS_MIN := [Dimension(i) : i in MINIMAL_MODULES | Dimension(i) le DIM_MIN];

FIELD_OF_DEFINITION := MinimalField([ C[i][j] : i,j in [1..#C] ] ) where C := CharacterTable(G);

if debug then FIELD_OF_DEFINITION; end if;

MODULI_AD := DIMS_AD;
for i in [1..#MODULI_AD] do
  MODULI_AD[i] := Gcd([MODULI_AD[j] : j in [i..#MODULI_AD]]);
end for;
MODULI_MIN := DIMS_MIN;
for i in [1..#MODULI_MIN] do
  MODULI_MIN[i] := Gcd([ IntegerRing() | MODULI_MIN[j] : j in [i..#MODULI_MIN]]);
end for;

B_AD := [BrauerCharacter(i) : i in ADJOINT_MODULES | Dimension(i) le DIM_AD];
B_MIN := [BrauerCharacter(i) : i in MINIMAL_MODULES | Dimension(i) le DIM_MIN];
if debug then printf "Done.\n"; end if;

if IsEmpty(ELEMENT_ORDERS) then
  ELEMENT_ORDERS := [c[1] : c in Clss | (c[1] gt 1) and (Gcd(c[1],CHAR_K) eq 1)];
end if;
EXCESS_ELEMENT_ORDERS := Sort(SetToSequence({x : x in ELEMENT_ORDERS | x gt LIMITING_ORDER}),func<x,y | x - y>);
if #EXCESS_ELEMENT_ORDERS gt 0 then
  printf "*** WARNING: Elements of order(s) %o not used.\n",EXCESS_ELEMENT_ORDERS;
  ELEMENT_ORDERS := [x : x in ELEMENT_ORDERS | x le LIMITING_ORDER];
end if;

// Store tuples giving the class numbers of order-i elements and their powers
T := [* 0 : i in [1..Max(ELEMENT_ORDERS)] *]; on := [];
for n in ELEMENT_ORDERS do
	on[n] := [[p(i,n/j) : j in Divisors(n) | j gt 1 ] : i in [1..#Clss] | Clss[i][1] eq n];
	
	T[n] := eval (LIE_TYPE cat "_ELTS" cat IntegerToString(n));
end for;

// D is the transpose of B, i.e. the traces of each class on the G-modules.
D_AD := [[FIELD_OF_DEFINITION | i[j] : i in B_AD] : j in [1..#B_AD[1]]];
if not (LIE_TYPE eq "E8") then
 	D_MIN := [[FIELD_OF_DEFINITION | i[j] : i in B_MIN] : j in [1..#B_MIN[1]]];
end if;

// Pair up classes with their traces on irred. modules.
ELT_CHAR_ARRAY_AD := [<Clss[i][1],D_AD[i]> : i in [1..#Clss]];
COEFFS_AD := [IntegerRing() | x : x in ELT_CHAR_ARRAY_AD[1][2]];
if not (LIE_TYPE eq "E8") then
  ELT_CHAR_ARRAY_MIN := [<Clss[i][1],D_MIN[i]> : i in [1..#Clss]];
  COEFFS_MIN := [IntegerRing() | x : x in ELT_CHAR_ARRAY_MIN[1][2]];
end if;


index_tick := function(integer_array,coefficients,maximum_sum,MODULI);
  // Increments the Brauer character coefficients
  // The first two arguments should be arrays of the same length,
  // the last should be the maximum allowed values of integer_array .* coeffients.
  
  // Returning 'true' means that integer_array can't be incremented without hitting a maximum_value somewhere.

  max_index := #integer_array;

  repeat
  	current_index := max_index;
 
    repeat

      DIFF := maximum_sum-&+[ IntegerRing() | integer_array[i]*coefficients[i] : i in [1..current_index-1] ];

      if
        (DIFF mod MODULI[current_index] eq 0)
      and
        (integer_array[current_index] lt DIFF div coefficients[current_index])
      then
        integer_array[current_index] +:= 1;
        break;
      else
        integer_array[current_index] := 0;
        current_index -:= 1;
      end if;

    until (current_index eq 0); // if this is reached, we are unable to tick;

    if debug and current_index eq 1 then printf "."; end if;

    if current_index eq 0 then
      return integer_array,true;
    end if;

  until IntegerRing()!(&+[coefficients[index]*integer_array[index] : index in [1..max_index]]) eq maximum_sum;

  // So here, the character degrees sum appropriately

  return integer_array,false;
	
end function;


if debug then printf "Beginning main loop.\n"; end if;

printf "\\begin{tabular}{r|*{%o}{c}", #DIMS_AD;
if not (LIE_TYPE eq "E_8") then printf "|*{%o}{c}", #DIMS_MIN; end if;
printf "}\n";

printf "& \\multicolumn{%o}{c}{$V_{%o}$}", #DIMS_AD, DIM_AD;
if not (LIE_TYPE eq "E_8") then printf "& \\multicolumn{%o}{|c}{$V_{%o}$}", #DIMS_MIN, DIM_MIN; end if;
printf " \\\\ \\hline\n\t";

for i in DIMS_AD do printf "&%2o ", i; end for;
if not (LIE_TYPE eq "E_8") then printf " "; for i in DIMS_MIN do printf "&%2o ", i; end for; end if;
printf "\\\\ \\hline\n";
print "%";

if LIE_TYPE eq "E8" then // only need to consider the adjoint module.
	thingAd := [0 : i in [1..#B_AD]];

	done := false;
	num := 0;

	while not done do

		charAd := &+[B_AD[i]*thingAd[i] : i in [1..#B_AD]]; 
		
		
		// Check if the adjoint traces match up
		if (&and[&and{[FIELD_OF_DEFINITION!charAd[i[j]] : j in [1..#i]] in T[n] : i in on[n] } : n in ELEMENT_ORDERS])
		then
		
		// Adjoint traces match up.
		num +:= 1;
		if num gt 1 then
			printf "\\\\\n";
		end if;
		printf "%2o)\t",num;
		for item in thingAd do printf "&%2o ", item; end for;
		end if;
		
		thingAd,done := index_tick(thingAd,COEFFS_AD,DIM_AD,MODULI_AD);
	end while;

else // we're not in E8
	thingAd := [0 : i in [1..#B_AD]];
	thingMin := [0 : i in [1..#B_MIN]];
		
	doneAd := false;
	num := 0;

	while not doneAd do
		charAd := &+[B_AD[i]*thingAd[i] : i in [1..#B_AD]];
		
		if (DIM_AD eq charAd[1])
		and
			(&and[&and{[FIELD_OF_DEFINITION!charAd[i[j]] : j in [1..#i]] in [z[2] : z in T[n]] : i in on[n] } : n in ELEMENT_ORDERS])
		then
		
			if debug then printf "\nAdjoint module possibility found.\n"; end if;
		
			doneMin := false;
		
			while not doneMin do
				
				charMin := &+[B_MIN[i]*thingMin[i] : i in [1..#thingMin]];
			
				// Here, the adjoint traces match up. Now check if the min traces match up
				if (charMin[1] eq DIM_MIN)
				and
					(&and[&and{[[FIELD_OF_DEFINITION!charMin[i[j]] : j in [1..#i]],[charAd[i[j]] : j in [1..#i]]] in T[n] : i in on[n] } : n in ELEMENT_ORDERS])
				then
				
				// Now, min and adjoint traces match up
				num +:= 1;
				if num gt 1 then printf "\\\\\n"; end if;
				printf "%2o)\t",num;
				for item in thingAd do printf "&%2o ", item; end for;
				printf " ";
				for item in thingMin do printf "&%2o ", item; end for;
				end if;
				
				thingMin,doneMin := index_tick(thingMin,COEFFS_MIN,DIM_MIN,MODULI_MIN);
			end while;
		end if;
		
		thingAd,doneAd := index_tick(thingAd,COEFFS_AD,DIM_AD,MODULI_AD);
		
	end while;
		
end if;

// If we have printed some possibility, need to terminate the line. If not, no need.
if num gt 0 then printf " \\\\ \\hline\n\\end{tabular}"; end if;
	
end procedure;
