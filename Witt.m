
(* Computations in the Witt ring *)

(*

See: NIRANJAN RAMACHANDRAN: ZETA FUNCTIONS, GROTHENDIECK GROUPS, AND THE WITT RING

example usage:

<<Witt`

metaZeta = ( (1 - Teichmuller[\[Alpha]]*s) (1 - Teichmuller[\[Beta]]*s) ) /
           ( (1 - Teichmuller[1]*s) (1 -Teichmuller[q]*s) )
ser = Series[metaZeta,{s,0,3}]
cf3 = Coefficient[ser,s,3]
WittToFunction[cf3]

*)

ClearAll["Witt`*"];
BeginPackage["Witt`"]
Needs["Useful`"]

Teichmuller     ::usage = "Teichmuller[a] is the Teichmuller element corresponding to a"

WittToFunction  ::usage = "WittToFunction[w,t] converts an element of the Witt to a function of t (t can be omitted)"
WittToLinearComb::usage = "WittToLinearComb[w] returns a list of {coeff,Teichmuller[a]} pairs"

ghostCoordinate ::usage = "ghostCoordinate[w,n] returns the n-th ghost coordinate"
ghostCoordinates::usage = "ghostCoordinates[w,n] returns the first n ghost coordinates"

(* ---------------------------------------- *)

(* reserve the variable "t" to be used as the formal variable of the Witt ring *)

ClearAll[Global`t]
SetAttributes[ Global`t , {Protected, ReadProtected} ]; 

(* ---------------------------------------- *)

Begin["`Private`"]

SetAttributes[ Teichmuller , {} ]; 
ClearAll[Teichmuller];

Teichmuller /: Teichmuller[a_]*Teichmuller[b_] := Teichmuller[a*b]
Teichmuller /: Teichmuller[a_]^n_ := Teichmuller[a^n]

Teichmuller /: Format[Teichmuller[a_]] := 
  DisplayForm[InterpretationBox[RowBox[{"[", a, "]"}], Teichmuller[a]]]

SetAttributes[ Teichmuller , {Protected, ReadProtected, Locked} ]; 

TeichmullerFun[a_   ] := TeichmullerFun[a,Global`t]
TeichmullerFun[a_,t_] := 1/(1-a*t)

isAllZero[es_List] := And @@ Table[e == 0, {e, es}]

(*returns the index of the only 1 exponent if there is exactly 1 and \
the rest is zero*)
expoIsLinear[es_List] := expoIsLinearWorker[1, es]
expoIsLinearWorker[_, {}] := Null;
expoIsLinearWorker[i_, es_List] := 
 Module[{this, rest}, this = First[es];
  rest = Rest[es];
  If[this == 1, If[isAllZero[rest], i, Null], 
   If[this == 0, expoIsLinearWorker[i + 1, rest], Null]]]

(* returns a list of {coeff,Teichmuller[a]} pairs *)
WittToLinearComb[X_] := Module[{coeffs, vars, fun, fun1, list},
  vars = Variables[X];
  vars = Select[vars, MatchQ[#, Teichmuller[_]] &];
  coeffs = MCoeffs[X, vars];
  fun[{cf_, expos_}] := 
   Module[{i = expoIsLinear[expos]}, 
    If[UnsameQ[i, Null], {cf, vars[[i]]}, 
     Throw["WittToFunction: error"]]];
  Map[fun, coeffs]
  ]

(* returns the Witt element as a function of t *)
WittToFunction[X_] := WittToFunction[X, Global`t]
WittToFunction[X_, t_] := Module[{list, fun1},
  fun1[{cf_, Teichmuller[a_]}] := (1/(1 - a*t))^cf;
  list = Map[fun1, WittToLinearComb[X]];
  Product[list[[i]], {i, 1, Length[list]}]
  ]

(* the n-th ghost coordinate *)
ghostCoordinate[X_, n_Integer] :=
 Module[{f},
  f[{cf_, Teichmuller[a_]}] := cf*a^n;
  SumList[Map[f, WittToLinearComb[X]]]
  ]

(* the first N ghost coordinates *)
ghostCoordinates[X_, n_Integer] :=
 Module[{f, list, g},
  list = WittToLinearComb[X];
  f[k_, {cf_, Teichmuller[a_]}] := cf*a^k;
  g[k_] := SumList[Map[f[k, #] &, list]];
  Table[g[k], {k, 1, n}]
  ]

(* ---------------------------------------- *)

End[]

EndPackage[]
