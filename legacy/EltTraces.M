AJL_MoodyPatera := function(LIE_TYPE,n)
// Auxiliary function computing the number of classes of elements
// of a order n in a complex group of exceptional type LIE_TYPE

  LIE_TYPE := CartanName(RootDatum(LIE_TYPE));

  Q<x,y> := PolynomialRing(RationalField(),2);
  R<X,Y>,r := Q/ideal<Q|[1 - x^3,1 - y^2]>;
  P<t> := PowerSeriesRing(R : Precision := n+2);
  
  G2 := (P!((1-t)*(1-(t^2))*(1-(t^3))))^(-1);
  F4 := (P!((1-t)*((1-(t^2))^2)*(1-(t^3))*(1-t^4)))^(-1);
  E6 := (P!( (1-t)*(1-(X*t))*(1-(X^2)*t)
        * (1-(t^2))*(1-(X*(t^2)))*(1-((X^2)*(t^2)))
        * (1-(t^3)) ))^(-1);
  E7 := (P!( (1-t)*(1-(Y*t))
          * ((1-(t^2))^2)*(1-(Y*(t^2)))
          * (1-(t^3))*(1-(Y*(t^3)))
          * (1-(t^4)) ))^(-1);
  E8 := (P!( (1-t)
          * ((1-(t^2))^2)
          * ((1-(t^3))^2)
          * ((1-(t^4))^2)
          * (1-(t^5))
          * (1-(t^6)) ))^(-1);
  
  if LIE_TYPE in ["G2","G_2"] then G := G2;
  elif LIE_TYPE in ["F4", "F_4"] then G := F4;
  elif LIE_TYPE in ["E6", "E_6"] then G := E6;
  elif LIE_TYPE in ["E7", "E_7"] then G := E7;
  elif LIE_TYPE in ["E8", "E_8"] then G := E8;
  else
    printf "Invalid Lie type.\n";
    return 0;
  end if;

  COEFFS := [Coefficient(G,i) :  i in [1..n]];

  NUMBEROF := COEFFS;
  for ORDER in [2..n] do
    if IsPrime(ORDER) then
      NUMBEROF[ORDER] := COEFFS[ORDER] - 1;
    else
      NUMBEROF[ORDER] := NUMBEROF[ORDER] - &+[NUMBEROF[k] : k in Divisors(ORDER) | not k eq ORDER];
    end if;
  end for;

  if LIE_TYPE eq CartanName(RootDatum("E6")) then
    for ORDER in [1..n] do
      NUMBEROF[ORDER] := Coefficient(NUMBEROF[ORDER],X,0);
    end for;
  elif LIE_TYPE eq CartanName(RootDatum("E7")) then
    for ORDER in [1..n] do
      NUMBEROF[ORDER] := Coefficient(NUMBEROF[ORDER],Y,0);
    end for;
  end if;
  
  return NUMBEROF;

end function;



exelts_mp := function(LIE_TYPE, n : debug := false);
/*
Given a Lie type, a sequence of high weights and an integer n > 1, returns
the traces of all elements of G of order n on the modules of these high
weights, where G is the simple algebraic group of type LIE_TYPE over
an alg. closed field of characteristic not dividing n.
*/

// Based on an algorithm described in https://doi.org/10.1137/0605037

LIE_TYPE := CartanName(RootDatum(LIE_TYPE));

DATUM := RootDatum(LIE_TYPE: Isogeny := "SC");
RANK := Rank(DATUM);
limit := AJL_MoodyPatera(LIE_TYPE,n)[n];
FUND_INDEX := #FundamentalGroup(DATUM);
if debug then printf "Looking for %o classes.\n",limit; end if;
if LIE_TYPE eq "G2" then MARKS := [3,2,1];
elif LIE_TYPE eq "F4" then MARKS := [2,3,4,2,1];
elif LIE_TYPE eq "E6" then MARKS := [1,2,2,3,2,1,1];
elif LIE_TYPE eq "E7" then MARKS := [2,2,3,4,3,2,1,1];
elif LIE_TYPE eq "E8" then MARKS := [2,3,4,6,5,4,3,2,1];
end if;

TIPS := [i : i in [1..#MARKS] | MARKS[i] eq 1 ];
REST := [i : i in [1..#MARKS] | not i in TIPS ];

  if LIE_TYPE eq "G2" then
    W := [&join[SetToMultiset(WeightOrbit(DATUM,x)) : x in [[0,0],[1,0]]],
        &join[SetToMultiset(WeightOrbit(DATUM,x)) : x in [[0,0]^^2,[1,0],[0,1]]]];
  elif LIE_TYPE eq "F4" then
    W :=  [&join[SetToMultiset(WeightOrbit(DATUM,x)) : x in [[0,0,0,0]^^2,[0,0,0,1]] ],
        &join[SetToMultiset(WeightOrbit(DATUM,x)) : x in [[0,0,0,0]^^4,[0,0,0,1],[1,0,0,0]] ] ];
  elif LIE_TYPE eq "E6" then
    W :=  [SetToMultiset(WeightOrbit(DATUM,[1,0,0,0,0,0])),
        &join[SetToMultiset(WeightOrbit(DATUM,x)) : x in [[0,1,0,0,0,0],[0,0,0,0,0,0]^^6]] ];
  elif LIE_TYPE eq "E7" then
    W :=  [SetToMultiset(WeightOrbit(DATUM,[0,0,0,0,0,0,1])),
        &join[SetToMultiset(WeightOrbit(DATUM,x)) : x in [[1,0,0,0,0,0,0],[0,0,0,0,0,0,0]^^7]] ];
  elif LIE_TYPE eq "E8" then
    W :=  [&join[SetToMultiset(WeightOrbit(DATUM,x)) : x in [[1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1]^^7,[0,0,0,0,0,0,0,0]^^35] ],
        &join[SetToMultiset(WeightOrbit(DATUM,x)) : x in [[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0]^^8]]];
  else
    "Error: Please give an exceptional Lie type.";
    return -1;
end if;

INV_CARTAN := GL(RANK,RationalField())!CartanMatrix(LIE_TYPE)^(-1);
W := [{* w*INV_CARTAN : w in v *} : v in W];
EIGS := {};
num := 0;

if LIE_TYPE in {"G2","F4","E8"} then MULTIPLIERS := [1];
elif LIE_TYPE eq "E6" then MULTIPLIERS := [1,3];
elif LIE_TYPE eq "E7" then MULTIPLIERS := [1,2];
end if;

OrderMult := function(s)
  if LIE_TYPE in {"G2","F4","E8"} then return 1;
  elif LIE_TYPE eq "E6" then return IntegerRing()!(3/Gcd(3,s[1]+s[5]-(s[3]+s[6])));
  elif LIE_TYPE eq "E7" then return IntegerRing()!(2/Gcd(2,s[2]+s[5]+s[7]));
  end if;
end function;

if debug then "Generating first approximation to search space..."; end if;
COEFF_TUPLES := CartesianProduct([[0..(n div MARKS[j])] : j in [1..RANK] ]);
if debug then "Done. First approximation:", #COEFF_TUPLES; "Generating second approximation..."; end if;
COEFF_TUPLES := { [j : j in i] cat [r] : r in {(n div M) - &+[i[j]*MARKS[j] : j in [1..RANK]] : M in MULTIPLIERS | n mod M eq 0} meet {0..n}, i in COEFF_TUPLES }; // Sum is n, 2n or 3n; no need to call OrderMult
if debug then "Done. Unfiltered search space:", #COEFF_TUPLES; "Filtering by GCD check..."; end if;
COEFF_TUPLES := { i : i in COEFF_TUPLES | (Gcd(SequenceToSet(i)) eq 1)};
if debug then "Done. Current approximation:", #COEFF_TUPLES; "Generating final search space..."; end if;
COEFF_TUPLES := { i : i in COEFF_TUPLES | &+[i[j]*MARKS[j] : j in [1..RANK+1]]*OrderMult(i) eq n };
if debug then "Done. Final search space:", #COEFF_TUPLES; end if;

EIGS := {
    [{* [IntegerRing()!((OrderMult([i : i in C])*FUND_INDEX)*C[x]*RationalField()!j[x]) : x in [1..RANK] ] : j in W[i] *} : i in [1,2]]
    : C in COEFF_TUPLES
};

num := #EIGS;

if num lt limit then printf "***WARNING: Not all classes found: %o of %o.\n",num,limit; end if;
if num gt limit then printf "***Too many classes found! This doesn't make sense! %o of %o.\n",num,limit; end if;

return EIGS;

end function;



EFOs_to_file := procedure(n)
  FILE_NAME := IntegerToString(n) cat ".M";
  LIE_TYPES := ["G2","F4","E6","E7","E8"];

  t := Cputime();

  printf "n = %o: ",n;
  for L in [1..#LIE_TYPES] do
    if (L eq 3) and (n mod 3 eq 0) then
      VAR := RootOfUnity(3*n);
    elif (L eq 4) and (n mod 2 eq 0) then
      VAR := RootOfUnity(2*n);
    else
      VAR := RootOfUnity(n);
    end if;

    printf LIE_TYPES[L] cat " ";
    X := exelts_mp(LIE_TYPES[L],n);

      X_traces := { [ [&+{* &*[VAR^Q : Q in q]^IntegerRing()!e : q in i*} : e in [n/j : j in Divisors(n) | not (j eq 1)] ] : i in x] : x in X};

      // if E8, forget the 3875-diml module
      if L eq 5 then X_traces := {i[2] : i in X_traces}; end if;

      Write(FILE_NAME,LIE_TYPES[L] cat "_ELTS" cat IntegerToString(n) cat " := ");
      Write(FILE_NAME,X_traces, "Magma");
      Write(FILE_NAME,";");

      // If running this for many element orders, better to free up memory in between
      delete X_traces;
      delete X;

  end for;
  printf "Time taken: %o.\n",Cputime(t);

end procedure;
