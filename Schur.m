
(* Symmetric polynomials, Schur functions and so on *)


(* example usage:

<< Schur`

poly = 3 ss[3, 2] + (8 - c) ss[2, 2, 1]
sToH[%, s, h]
hToE[%, h, e]
eToX[%, e, {x, 5}];
xToS[%, {x, 5}, s]
sToH[%, s, h]

*)


Clear["Schur`*"];
BeginPackage["Schur`"]
Needs["Useful`"]
Needs["Partitions`"]

ss::usage = "ss[3,2,1] - create standard Schur polynomials"
ee::usage = "ee[3,2,1] - create standard elementary symmetric polynomials"
hh::usage = "hh[3,2,1] - create standard homogeneous symmetric polynomials"
pp::usage = "pp[3,2,1] - create standard power symmetric polynomials"

mkSubscript::usage = "mkSubscript[var,3,2,1] - create subscripted variables. All kinds of input is accepted."
mkSchur    ::usage = "mkSchur[var,{3,2,1}] is like mkSubscript, but takes an arbitrary integer sequence, and normalizes the resulting Schur polynomial (with a sign or zero)"

lexSort     ::usage = "lexicographical sort of a list of lists"
lexCompare  ::usage = "lexicographical comparison of lists (-1 = LT, 0 = EQ , +1 = GT)"
lexMaximumBy::usage = "lexMaximumBy[list,fun] returns the lexicographically largest element after applying fun"

JacobiTrudiE::usage = "JacobiTrudiE[lam,e] returns the Jacobi-Trudi determinant for s_lambda in elementary symm polys"
JacobiTrudiH::usage = "JacobiTrudiH[lam,h] returns the Jacobi-Trudi determinant for s_lambda in complete homogeneous symm polys"

expandE::usage = "expandE[expr,e] expands terms like e[3,2,2] into e[3]*e[2]^2"
expandH::usage = "expandH[expr,h] expands terms like h[3,2,2] into h[3]*h[2]^2"
expandP::usage = "expandP[expr,p] expands terms like p[3,2,2] into p[3]*p[2]^2"

collapseE::usage = "collapseE[expr,e] collapes terms like e[3]*e[2]^2 into e[3,2,2]"
collapseH::usage = "collapseH[expr,h] collapes terms like h[3]*h[2]^2 into h[3,2,2]"
collapseP::usage = "collapseP[expr,p] collapes terms like p[3]*p[2]^2 into p[3,2,2]"

collectS::usage = "collectS[expr,s] collects the coefficients of variables like s[3,2,2] together"
collectE::usage = "collectE[expr,e] collects the coefficients of variables like e[3,2,2] together"
collectH::usage = "collectH[expr,h] collects the coefficients of variables like h[3,2,2] together"
collectP::usage = "collectP[expr,p] collects the coefficients of variables like p[3,2,2] together"

transposeSchur::usage = "transposeSchur[X,s] transposes the indices of the s[...] using partition duality"

expandGeneric  ::usage = "expandGeneric[expr,var] expands terms like var[3,2,2] into var[3]*var[2]^2"
collapseGeneric::usage = "collapseGeneric[expr,var] collapes terms like var[3]*var[2]^2 into var[3,2,2]"
collectGeneric ::usage = "collectGeneric[expr,var] collects the coefficients of variables like var[3,2,2] together"
coeffsGeneric  ::usage = "coeffsGeneric[expr,var] returns the list of {idxs,coefficient} pairs where idxs means var[idxs]"
coeffsGenericBi::usage = "coeffsGenericBi[expr,var1,var2] returns the list of {idxs1,idxs2,coefficient} triple"

sToE::usage = "sToE[expr,s,e] converts Schur polynomials to elementary symmetric polynomials"
eToS::usage = "eToS[expr,e,s] converts elementary symmetric polynomials to Schur polynomials"

sToH::usage = "sToH[expr,s,h] converts Schur polynomials to complete homogeneous symmetric polynomials"
hToS::usage = "hToS[expr,h,s] converts complete homogeneous symmetric polynomials to Schur polynomials"

eToH::usage = "eToH[expr,e,h] converts elementary symmetric polynomials to complete homogeneous symmetric polynomials"
hToE::usage = "hToE[expr,h,e] converts complete homogeneous symmetric polynomials to elementary symmetric polynomials"

sToX::usage = "sToX[expr,s,{x,n}] expand Schur polynomials in terms of x1...xn"
eToX::usage = "eToX[expr,e,{x,n}] expand elementary symmetric polynomials in terms of x1...xn"
hToX::usage = "hToX[expr,h,{x,n}] expand complete homogeneous polynomials in terms of x1...xn"

xToS::usage = "xToS[expr,{x,n},s] converts symmetric polynomials in x1..xn to Schur polynomials"
xToE::usage = "xToE[expr,{x,n},e] converts symmetric polynomials in x1..xn to elementary symmetric polynomials"
xToH::usage = "xToH[expr,{x,n},h] converts symmetric polynomials in x1..xn to complete homogoneous polynomials"

lookupEinS::usage = "lookupEinS[lam,s] looks up the expansion of e[lam] in Schur polynomials"
lookupHinS::usage = "lookupHinS[lam,s] looks up the expansion of h[lam] in Schur polynomials"

lookupEinH::usage = "lookupEinH[lam,h] looks up the expansion of e[lam] in complete homogeneous symmetric polynomials"
lookupHinE::usage = "lookupHinE[lam,e] looks up the expansion of h[lam] in elementary symmetric polynomials"

lookupSinH::usage = "lookupSinH[lam,h] looks up the expansion of s[lam] in complete homogeneous symmetric polynomials"
lookupSinE::usage = "lookupSinE[lam,e] looks up the expansion of s[lam] in elementary symmetric polynomials"

EMatrix::usage = "EMatrix[lam,mu,n] is the matrix E_{lam/mu}(n), whose determinant appear in the expansion of total Chern class of tensor products"
EDet   ::usage = "EDet[lam,mu,n] is the determinant E_{lam/mu}(n) which appears in the expansion of total Chern class of tensor products"

schurWeylCharacter::usage = "schurWeylCharacter[lam,{x,n}] implements the Weyl character formula for Schur polynomials"

extractSubscripts::usage = "debugging"
filterSubscripts::usage  = "debugging"

addGrading::usage = "addGrading[X,s,h] multiplies each term s[...] by h^deg"

eToDiffX::usage = "eToDiffX[X,e,{x,n},{y,m}] makes the subsitution c(E) -> c(x-y)"
eToDiffE::usage = "eToDiffE[X,e,{c,n},{d,m}] makes the subsitution c(E) -> c(A-B)" 
eToDiffS::usage = "eToDiffS[X,e,{s,n},{t,m}] makes the subsitution c(E) -> s(A-B)" 
sToDiffX::usage = "eToDiffE[X,s,{x,n},{y,m}] makes the subsitution s(E) -> s(x-y)"
sToDiffE::usage = "sToDiffE[X,s,{c,n},{d,m}] makes the subsitution s(E) -> c(A-B)" 
sToDiffS::usage = "sToDiffS[X,s,{s,n},{s,m}] makes the subsitution s(E) -> s(A-B)" 

fromDiffS::usage = "fromDiffS[X,{s,n},{t,m},S] tries to convert something in (s,t) to a polynomial in the difference alphabet S"
fromDiffX::usage = "fromDiffX[X,{x,n},{y,m},S,deg] tries to convert something in (x,y) to a polynomial in the difference alphabet S"

(* ========================================================================== *)

Begin[ "`Private`"]

SetAttributes[ Global`s , {Protected, ReadProtected} ]; 
SetAttributes[ Global`e , {Protected, ReadProtected} ]; 
SetAttributes[ Global`h , {Protected, ReadProtected} ]; 
SetAttributes[ Global`p , {Protected, ReadProtected} ]; 

ss[i__] := mkSubscript[ Global`s, i ];
ee[i__] := mkSubscript[ Global`e, i ];
hh[i__] := mkSubscript[ Global`h, i ];
pp[i__] := mkSubscript[ Global`p, i ];

(* ------ internal free variables used in tables ------ *)

SetAttributes[ Schur`Private`S , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Schur`Private`E , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Schur`Private`H , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Schur`Private`P , {Protected, ReadProtected, Locked} ]; 

SS = Schur`Private`S
EE = Schur`Private`E
HH = Schur`Private`H
PP = Schur`Private`P


(* accept all three forms *)
Clear[mkSubscript]
mkSubscript[var_, {}   ]  := 1   (* hackety hack *)
mkSubscript[var_, idxs_]  := If[ListQ[idxs], Subscript @@ Prepend[idxs, var], Subscript[var, idxs]]
mkSubscript[var_, idxs__] := mkSubscript[var, {idxs}]

(* --------------- Schur stufff------------------- *)

Clear[mkSchur]
mkSchur[var_, {}   ]  := 1   (* hackety hack *)
mkSchur[var_, idxs_]  := If[ListQ[idxs], mkSchurCF[var,idxs] , mkSchurCF[var,{idxs}] ]
mkSchur[var_, idxs__] := mkSchurCF[var, {idxs}]

(* like mkSubscript, but normalizes to a signed Schur polynomial (or zero) *)
mkSchurCF[ps_List] := mkSchurCF[Global`s, ps, 1]
mkSchurCF[s_, ps_] := mkSchurCF[s, ps, 1]
mkSchurCF[s_, ps0_, cf0_] := Module[
  {ps = ps0, cf = cf0, n = Length[ps0], ok, i, a, b},
  ps = ps0;
  ok = True;
  For[i = 1, i <= n, i++,
   a = partitionIndex[ps, i];
   b = partitionIndex[ps, i + 1];
   If[a == b - 1, Return[0]];
   If[a < b - 1,
    ok = False;
    cf = -cf;
    ps = ReplacePart[ps, {i -> b - 1, i + 1 -> a + 1}]];
   ];
  If[ok, 
   cf*mkSubscript[s, removeZeros[ps]],
   mkSchurCF[s, ps, cf]
   ]
  ]

(* the Weyl character formula for Schur polynomials *)
schurWeylCharacter[part_List, x_Symbol] := schurWeylCharacter[part, {x, Length[part]}]
schurWeylCharacter[part_List, {x_Symbol, n_Integer}] := Module[
  {A, V, M},
  M = Table[ 
    Subscript[x, j]^(partitionIndex[part, i] + n - i) , {i, 1, n}, {j,
      1, n}];
  A = Det[M];
  V = Product[
    Product[Subscript[x, i] - Subscript[x, j], {j, i + 1, n}], {i, 1, 
     n - 1}];
  Expand[Factor[A/V]]
  ]

(* for ex.: Map[schurNormalizationSanityCheck, Combinatorica`Compositions[7, 4]] *)
schurNormalizationSanityCheck[p_] := schurNormalizationSanityCheck[p,False]
schurNormalizationSanityCheck[p_,debug_] := Module[
  {n = Length[p] + 1, lhs, rhs, tmp, x, s},
  x = Global`x;
  s = Global`s;
  lhs = schurWeylCharacter[p, {x, n}];
  tmp = mkSchur[s, p];
  rhs = sToX[tmp, s, {x, n}];
  If[debug,
    Print[lhs];
    Print[rhs];
  ];
  Expand[lhs - rhs]
  ]

(* --------------- misc ------------------- *)

removeZeros[L_] := Select[L, # != 0 &]

constantTerm[X_, vars_] := mcoeff[X, vars, replicate[0, Length[vars]]]

EMatrix[lam_, mu_, n_] := Table[
  Binomial[
   partitionIndex[lam, i] + n - i,
   partitionIndex[mu , j] + n - j],
  {i, 1, n}, {j, 1, n}]

EDet[lam_, mu_, n_] := Det[EMatrix[lam, mu, n]]

(* --------- lexicographic ordering -------- *)

(* lexicographical sort. Taken from:
   <https://mathematica.stackexchange.com/questions/85606/lexicographic-\
   ordering-of-lists-of-lists>
*)

lexSort$helper[x_, lev_] := 
 If[Length[x] == 1, x, 
  With[{smallest = LengthWhile[x, Length[#] == lev &]}, 
   Take[x, smallest]~Join~lexSort1[Drop[x, smallest], lev + 1]]]

(* This code by a Mathematica expert does not handle the empty list... see, thats *EXACTLY* why Mathematica programming sucks... *) 
ClearAll[lexSort1];
lexSort1[{}, _ : None] := {};
lexSort1[{x_}, _ : None] := {x};
lexSort1[l_List] := lexSort1[l, 1];
lexSort1[l_List, lev_] := With[{min = Min[Length /@ l]}, Composition[
    Flatten[#, 1] &,
    Map[lexSort$helper[#, lev] &, #] &,
    SplitBy[#, Take[#, {lev, min}] &] &,
    SortBy[#, Take[#, {lev, min}] &] &
    ]@l]

ClearAll[lexSort]
lexSort[list_] := Join[
  Select[list, IsEmptyList],
  lexSort1[Select[list, IsNotEmptyList]]
  ]

(* LT -> -1, EQ -> 0, GT -> 1 *)
Clear[lexCompare];
lexCompare[{}, {}] = 0;
lexCompare[{}, xs_] := -1;
lexCompare[xs_, {}] := 1
lexCompare[xs_, ys_] := Module[{x = First[xs], y = First[ys]},
  If[x < y, -1, If[x > y, 1, lexCompare[Rest[xs], Rest[ys]]]]]

lexMaximumBy[list_, fun_] := Module[
  {i, best, result, n, this, val},
  result = list[[1]];
  best = fun[result];
  n = Length[list];
  For[i = 2, i <= n, i++,
   this = list[[i]];
   val = fun[this];
   If[lexCompare[val, best] > 0, best = val; result = this, Null];
   ];
  result]

(* -------------------------------------------------------------------------- *)


(* extract the subscript lists of the given variable *)
extractSubscripts[term_, var_, maxwidth_ , maxheight_ ] := filterSubscripts[ extractSubscripts[term,var] , maxwidth , maxheight ]
extractSubscripts[term_, var_, maxwidth_]               := filterSubscripts[ extractSubscripts[term,var] , maxwidth ]
extractSubscripts[term_, var_] := Module[{xs, list, vars},
  list = Cases[Variables[term], Subscript[var, xs__] -> {xs}];
  (* we have to handle the constant term too ... *)
  
  vars = Map[mkSubscript[var, #] &, list];
  If[SameQ[Expand[constantTerm[term, vars]], 0], Null, 
   PrependTo[list, {}]];
  Sort[DeleteDuplicates[list]]
  ]

filterSubscripts[list_,maxwidth_]            := filterSubscripts[list,maxwidth,Infinity]
filterSubscripts[list_,maxwidth_,maxheight_] := Module[
  {width,height},
  width[p_] := Length[p];
  height[p_]:= Max@@p;     (* note: can be -Infinity! *)
  Select[list , width[#]<=maxwidth && height[#]<=maxheight & ]
  ]

(* expand terms like e_{1,1,2,3} to e_ 1^2*e_ 2*e_ 3 *)
expandGeneric[term_, e_] := Module[
  {list = extractSubscripts[term, e]},
  term /. 
   Table[mkSubscript[e, idxs] -> 
     Product[Subscript[e, i], {i, idxs}], {idxs, list}]
  ]

(* collapse terms like e_ 1^2*e_ 2*e_ 3 to e_{3,2,1,1} *)
collapseGeneric[term_, e_]                          := collapseGeneric[term,e,Infinity,Infinity]
collapseGeneric[term_, e_, {maxwidth_, maxheight_}] := collapseGeneric[term,e,maxwidth,maxheight]
collapseGeneric[term_, e_, maxwidth_ ]              := collapseGeneric[term,e,maxwidth,Infinity]
collapseGeneric[term_, e_, maxwidth_, maxheight_] := Module[
  {A, js, m, vars, fun, width, height},
  A = Expand[expandGeneric[term, e]];
  js = Flatten[extractSubscripts[A, e]];
  m = Max @@ Append[js, 0];
  vars = Table[Subscript[e, i], {i, 1, m}];
  width[xs_]  := Length[xs];
  height[xs_] := Max@@xs;
  fun[{c_, idxs_}] := Module[{p = fromExpoVec[idxs]},
    If[ width[p] <= maxwidth && height[p] <= maxheight
    , c * mkSubscript[e, p]
    , 0 ]];
  Sum[fun[cf], {cf, MCoeffs[A, vars]}]
  ]


expandE[term_]     := expandGeneric[term, Global`e]
expandE[term_, e_] := expandGeneric[term, e]

expandH[term_]     := expandGeneric[term, Global`h]
expandH[term_, h_] := expandGeneric[term, h]

expandP[term_]     := expandGeneric[term, Global`p]
expandP[term_, p_] := expandGeneric[term, p]

collapseE[term_]               := collapseGeneric[term, Global`e]
collapseE[term_, e_]           := collapseGeneric[term, e]
collapseE[term_, e_, wd_]      := collapseGeneric[term, e , wd]
collapseE[term_, e_, wd_, ht_] := collapseGeneric[term, e , wd, ht]

collapseH[term_]               := collapseGeneric[term, Global`h]
collapseH[term_, h_]           := collapseGeneric[term, h]
collapseH[term_, h_, wd_]      := collapseGeneric[term, h , wd]
collapseH[term_, h_, wd_, ht_] := collapseGeneric[term, h , wd, ht]

collapseP[term_]               := collapseGeneric[term, Global`p]
collapseP[term_, p_]           := collapseGeneric[term, p]
collapseP[term_, p_, wd_]      := collapseGeneric[term, p , wd]
collapseP[term_, p_, wd_, ht_] := collapseGeneric[term, p , wd, ht]

(* ---------- coeffs ----------- *)

coeffsGeneric[term_, e_ ]                  := coeffsGeneric[term, e , Infinity, Infinity ]
coeffsGeneric[term_, e_ , {maxw_ , maxh_}] := coeffsGeneric[term, e , maxw , maxh        ] 
coeffsGeneric[term_, e_ , maxw_ ]          := coeffsGeneric[term, e , maxw , Infinity   ]
coeffsGeneric[term_, e_ , maxw_ , maxh_] := Module[
  {A = collapseGeneric[term, e, maxw, maxh], idxSet, varSet, L, kst},
  idxSet = extractSubscripts[A, e];
  idxSet1 = Select[ idxSet, IsNotEmptyList[#]& ];
  varSet = Map[mkSubscript[e, #] &, idxSet1];
  L = Table[ {idx,Coefficient[A,mkSubscript[e,idx]]} , {idx,idxSet1} ];
  kst = Expand[constantTerm[term, varSet]];
  If[SameQ[kst, 0], L, PrependTo[L, {{},kst}]]
  ]

(* like coeffsGeneric but for 2 variables *)
coeffsGenericBi[expr_, sr_, sq_] := Module[
  {f, A, B},
  A = coeffsGeneric[expr, sr];
  f[{ridx_, cf_}] := 
   Table[{ridx, pair[[1]], pair[[2]]}, {pair, 
     coeffsGeneric[cf, sq]}];
  Flatten[Map[f, A], 1]
  ]

(* ---------- transpose indices of Schur polys  ---------- *)

transposeSchur[A_, s_] := Module[{T},
  T = coeffsGeneric[A, s];
  Sum[pair[[2]]*mkSchur[s, dualPart[pair[[1]]]], {pair, T}]
  ]

(* ---------- collect ---------- *)

collectGeneric[term_, e_ ]                  := collectGeneric[term, e , Infinity, Infinity ]
collectGeneric[term_, e_ , {maxw_ , maxh_}] := collectGeneric[term, e , maxw , maxh        ] 
collectGeneric[term_, e_ , maxw_ ]          := collectGeneric[term, e , maxw , Infinity   ]
collectGeneric[term_, e_ , maxw_ , maxh_] := Module[
  {A = collapseGeneric[term, e, maxw, maxh], idxSet, varSet},
  idxSet = extractSubscripts[A, e];
  varSet = Map[mkSubscript[e, #] &, idxSet];
  varSet = Select[varSet, UnsameQ[#, 1] &];
  Collect[A, varSet]
  ]

collectS[term_]               := collectGeneric[term, Global`s]
collectS[term_, s_]           := collectGeneric[term, s]
collectS[term_, s_, wd_]      := collectGeneric[term, s , wd]
collectS[term_, s_, wd_, ht_] := collectGeneric[term, s , wd, ht]

collectE[term_]               := collectGeneric[term, Global`e]
collectE[term_, e_]           := collectGeneric[term, e]
collectE[term_, e_, wd_]      := collectGeneric[term, e , wd]
collectE[term_, e_, wd_, ht_] := collectGeneric[term, e , wd, ht]

collectH[term_]               := collectGeneric[term, Global`h]
collectH[term_, h_]           := collectGeneric[term, h]
collectH[term_, h_, wd_]      := collectGeneric[term, h , wd]
collectH[term_, h_, wd_, ht_] := collectGeneric[term, h , wd, ht]

collectP[term_]               := collectGeneric[term, Global`p]
collectP[term_, p_]           := collectGeneric[term, p]
collectp[term_, p_, wd_]      := collectGeneric[term, p , wd]
collectp[term_, p_, wd_, ht_] := collectGeneric[term, p , wd, ht]

(* ------------------ conversions ------------------- *)

Clear[JacobiTrudiCached]
JacobiTrudiCached[{}] = 1;
JacobiTrudiCached[lam_] := JacobiTrudiCached[lam] = Module[
   {n = Length[lam], hh, A},
   hh[i_] := If[i < 0, 0, If[i == 0, 1, Subscript[HH, i]]];
   A = Table[Table[hh[lam[[i]] + j - i], {j, 1, n}], {i, 1, n}];
   Det[A]
   ]

JacobiTrudiH[lam_] := JacobiTrudiH[lam, Global`h]
JacobiTrudiH[lam_, h_] := JacobiTrudiCached[lam] /. {HH -> h}

JacobiTrudiE[lam_] := JacobiTrudiE[lam, Global`e]
JacobiTrudiE[lam_, e_] := JacobiTrudiCached[DualPart[lam]] /. {HH -> e}

lookupSinH[lam_] := lookupSinH[lam, Global`h]
lookupSinH[lam_, h_] := JacobiTrudiCached[lam] /. {HH -> h}

lookupSinE[lam_] := lookupSinE[lam, Global`e]
lookupSinE[lam_, e_] := JacobiTrudiCached[DualPart[lam]] /. {HH -> e}


sToE[term_]                           := sToE[term, Global`s, Global`e]
sToE[term_, s_Symbol]                 := sToE[term, s, Global`e]
sToE[term_, s_Symbol,  e_Symbol]      := sToE[term, s, {e,Infinity,Infinity}]
sToE[term_, s_Symbol, {e_Symbol,wd_}] := sToE[term, s, {e,wd,Infinity}]
sToE[term_, s_Symbol, {e_Symbol,wd_,ht_}] := Module[
  {list = extractSubscripts[term, s], B},
  B = Expand[
    term /. Table[
      mkSubscript[s, idxs] -> JacobiTrudiE[idxs, e], {idxs, list}]];
  collectE[B, e , wd, ht]]

sToH[term_]                            := sToH[term, Global`s, Global`h]
sToH[term_, s_Symbol]                  := sToH[term, s, Global`h]
sToH[term_, s_Symbol,  h_Symbol]       := sToH[term, s, {h,Infinity,Infinity}]
sToH[term_, s_Symbol, {h_Symbol,wd_}]  := sToH[term, s, {h,wd,Infinity}]
sToH[term_, s_Symbol, {h_Symbol,wd_,ht_}] := Module[
  {list = extractSubscripts[term, s], B},
  B = Expand[
    term /. Table[
      mkSubscript[s, idxs] -> JacobiTrudiH[idxs, h], {idxs, list}]];
  collectH[B, h, wd, ht]]
  
(* -------------------------------------------- *)

calc$eToS[term_]         := calc$eToS[term, Global`e, Global`s]
calc$eToS[term_, e_]     := calc$eToS[term, e, Global`s]
calc$eToS[term_, e_, s_] := Module[
  {A, js, m, vars, cfs, Acc, clam, c, es, lam, fun, lamSet, varSet},
  A = Expand[expandE[term, e]];
  Acc = 0;
  lamSet = {};
  fun[{c_, es_}] := {c, DualPart[fromExpoVec[es]]};
  While[UnsameQ[A, 0],
   js = Flatten[extractSubscripts[A, e]];
   m = Max @@ Append[js, 0];
   vars = Table[Subscript[e, i], {i, 1, m}];
   cfs = Map[fun, MCoeffs[A, vars]];
   clam = lexMaximumBy[cfs, Snd];
   c = Fst[clam];
   lam = Snd[clam];
   Acc = Acc + c*mkSubscript[s, lam];
   A = Expand[A - c*JacobiTrudiE[lam, e]];
   AppendTo[lamSet, lam]
   (* Print[{c,lam}]; *)
   ];
  lamSet = lexSort[DeleteDuplicates[lamSet]]; 
  Acc = Expand[Acc];
  varSet = Map[mkSubscript[s, #] &, lamSet];
  varSet = Select[varSet, UnsameQ[#, 1] &];
  Collect[Acc, varSet]
  ]

calc$hToS[term_]         := calc$hToS[term, Global`h, Global`s]
calc$hToS[term_, h_]     := calc$hToS[term, h, Global`s]
calc$hToS[term_, h_, s_] := Module[
  {A, js, m, vars, cfs, Acc, clam, c, es, lam, fun, lamSet, varSet},
  A = Expand[expandH[term, h]];
  Acc = 0;
  lamSet = {};
  fun[{c_, es_}] := {c, DualPart[fromExpoVec[es]]};
  While[UnsameQ[A, 0],
   js = Flatten[extractSubscripts[A, h]];
   m = Max @@ Append[js, 0];
   vars = Table[Subscript[h, i], {i, 1, m}];
   cfs = Map[fun, MCoeffs[A, vars]];
   clam = lexMaximumBy[cfs, Snd];
   c = Fst[clam];
   lam = DualPart[Snd[clam]];
   Acc = Acc + c*mkSubscript[s, lam];
   A = Expand[A - c*JacobiTrudiH[lam, h]];
   AppendTo[lamSet, lam]
   (* Print[{c,lam}];  *)
   ];
  lamSet = lexSort[DeleteDuplicates[lamSet]]; 
  Acc = Expand[Acc];
  varSet = Map[mkSubscript[s, #] &, lamSet];
  varSet = Select[varSet, UnsameQ[#, 1] &];
  Collect[Acc, varSet]
  ]

Clear[cached$EinS];
cached$EinS[lam_] := cached$EinS[lam] = Module[ {e} ,
  calc$eToS[ mkSubscript[e, lam], e, SS] ]

Clear[cached$HinS];
cached$HinS[lam_] := cached$HinS[lam] = Module[ {h},
  calc$hToS[ mkSubscript[h, lam], h, SS] ]

lookupEinS[lam_, s_] := cached$EinS[lam] /. {SS -> s}
lookupHinS[lam_, s_] := cached$HinS[lam] /. {SS -> s}

eToS[term_]                      := eToS[term, Global`e, Global`s]
eToS[term_,  e_]                 := eToS[term, e, Global`s]
eToS[term_,  e_, s_Symbol]       := eToS[term, e, {s,Infinity,Infinity}]
eToS[term_,  e_, {s_ ,wd_ }]     := eToS[term, e, {s,wd,Infinity}]
eToS[term0_, e_, {s_, wd_, ht_}] := Module[
  {list, term, B},
  term = collapseE[term0, e];
  list = extractSubscripts[term, e];
  B = Expand[
    term /. Table[
      mkSubscript[e, idxs] -> lookupEinS[idxs, s], {idxs, list}]];
  collectS[B, s, wd, ht]
  ]

hToS[term_]                       := hToS[term, Global`h, Global`s]
hToS[term_,  h_]                  := hToS[term, h, Global`s]
hToS[term_,  h_, s_Symbol]        := hToS[term, h, {s,Infinity,Infinity}]
hToS[term_,  h_, {s_ ,wd_ }]      := hToS[term, h, {s,wd,Infinity}]
hToS[term0_, h_, {s_, wd_, ht_}]  := Module[
  {list, term, B},
  term = collapseH[term0, h];
  list = extractSubscripts[term, h];
  B = Expand[
    term /. Table[
      mkSubscript[h, idxs] -> lookupHinS[idxs, s], {idxs, list}]];
  collectS[B, s, wd, ht]
  ]

(* --------------------------------------- *)

Clear[cached$EinH];
cached$EinH[lam_] := 
 Module[ {e,s},
  sToH[eToS[mkSubscript[e, lam], e, s], s, HH]]

Clear[cached$HinE];
cached$HinE[lam_] := 
 Module[ {h,s},
  sToE[hToS[mkSubscript[h, lam], h, s], s, EE]]

lookupEinH[lam_, h_] := cached$EinH[lam] /. {HH -> h}
lookupHinE[lam_, e_] := cached$HinE[lam] /. {EE -> e}

eToH[term_]                     := eToH[term, Global`e, Global`h]
eToH[term_,  e_]                := eToH[term, e, Global`h]
eToH[term_,  e_, h_Symbol]      := eToH[term, e, {h,Infinity,Infinity} ]
eToH[term_,  e_, {h_,wd_}]      := eToH[term, e, {h,wd      ,Infinity} ]
eToH[term0_, e_, {h_,wd_,ht_}]  := Module[
  {list, term, B},
  term = collapseE[term0, e];
  list = extractSubscripts[term, e];
  B = Expand[
    term /. Table[
      mkSubscript[e, idxs] -> lookupEinH[idxs, h], {idxs, list}]];
  collectH[B, h, wd, ht]
  ]

hToE[term_]                  := hToE[term, Global`h, Global`e]
hToE[term_,  h_]             := hToE[term, h, Global`e]
hToE[term_,  h_, e_Symbol]   := hToE[term, h, {e, Infinity, Infinity}]
hToE[term_,  h_, {e_,wd_} ]  := hToE[term, h, {e, wd      , Infinity}]
hToE[term0_, h_, {e_,wd_,ht_}] := Module[
  {list, term, B},
  term = collapseH[term0, h];
  list = extractSubscripts[term, h];
  B = Expand[
    term /. 
     Table[mkSubscript[h, idxs] -> lookupHinE[idxs, e], {idxs, list}]];
  collectE[B, e , wd, ht]
  ]

(* ------------------------------------ *)

sToX[term_, s_, {x_, n_}] := sToX[term, s, x, n]
sToX[term_, s_, x_, n_] := Module[ {e},
  eToX[ sToE[term, s, e], e, {x, n}]]


eToX[term_, e_, {x_, n_}] := eToX[term, e, x, n]
eToX[term_, e_, x_, n_] := Module[
  {A, total, q, elist, js, efun, m},
  A = expandE[term, e];
  total = Product[1 + q*Subscript[x, i], {i, 1, n}];
  elist = CoefficientList[total, q];
  efun[k_] := 
   If[k < 0, 0, If[k == 0, 1, If[k > n, 0, elist[[k + 1]]]]];
  js = Flatten[extractSubscripts[A, e]];
  m = Max @@ Append[js, 0];
  Expand[A /. Table[Subscript[e, i] -> efun[i], {i, 1, m}]]
  ]

hToX[term_, h_, {x_, n_}] := hToX[term, h, x, n]
hToX[term_, h_, x_, n_] := Module[
  {A, total, q, hlist, js, hfun, m, B},
  A = expandE[term, h];
  js = Flatten[extractSubscripts[A, h]];
  m = Max @@ Append[js, 0];
  total = Product[1/(1 - q*Subscript[x, i]), {i, 1, n}];
  hlist = CoefficientList[Series[total, {q, 0, m}], q];
  hfun[k_] := If[k < 0, 0, If[k == 0, 1, hlist[[k + 1]]]];
  Expand[A /. Table[Subscript[h, i] -> hfun[i], {i, 1, m}]]
  ]      

xToE[term_, {x_, n_}, e_] := xToE[term, x, n, e]
xToE[term_, x_, n_, e_] := Module[
  {pair,A},
  pair = SymmetricReduction[term,
    Table[Subscript[x, i], {i, 1, n}],
    Table[Subscript[e, i], {i, 1, n}]];
  If[SameQ[Snd[pair], 0], Null, 
   Throw["xToE: polynomial is not symmetric"]];
  A=collectE[Expand[Fst[pair]], e , {Infinity, n}];
  A
  ]

(* this behaves really strangely ???
xToH[term_, {x_, n_}, h_] := xToH[term, x, n, h]
xToH[term_,  x_, n_ , h_] := Module[ {e} ,
 eToH[ xToE[term, {x, n}, e], e, h] ]
*)

xToH[term_, {x_, n_}, h_] := xToH[term, x, n, h]
xToH[term_,  x_, n_ , h_] := Module[ {e} ,
 sToH[ xToS[term, {x, n}, e], e, {h,n}] ]

xToS[term_, {x_, n_}, s_] := xToS[term, x, n, s]
xToS[term_,  x_, n_ , s_] := Module[ {e},
 eToS[ xToE[term, {x, n}, e], e, {s,n}] ]


(* ----------------------------------------------- *)
(* difference alphabets *)

eToDiffX[X0_, e_, {x_, n_}, {y_, m_}] := Module[
  {X, A, B, C, N, h},
  X = Expand[X0];
  X = expandE[X, e];
  N = Max[Join[{0},Flatten[extractSubscripts[X, e]]]];
  A = Product[1 + h*Subscript[y, j], {j, 1, m}];
  B = Product[1 + h*Subscript[x, i], {i, 1, n}];
  C = Normal[Series[A/B, {h, 0, N}]];
  Expand[ X /. Table[Subscript[e, i] -> Coefficient[C, h, i], {i, 1, N}] ]
  ]

eToDiffE[X0_, e_, {c_, n_}, {d_, m_}] := Module[
  {X, A, B, C, N, h},
  X = Expand[X0];
  X = expandE[X, e];
  N = Max[Join[{0},Flatten[extractSubscripts[X, e]]]];
  A = 1 + Sum[ h^j * Subscript[d, j], {j, 1, m}];
  B = 1 + Sum[ h^i * Subscript[c, i], {i, 1, n}];
  C = Normal[Series[A/B, {h, 0, N}]];
  Expand[ X /. Table[Subscript[e, i] -> Coefficient[C, h, i], {i, 1, N}] ]
  ]

sToDiffE[X_, s_, cn_, dm_] := Module[
  {A, e},
  A = sToE[X, s, e];
  eToDiffE[A, e, cn, dm]
  ]

eToDiffS[X_, e_, {s_, n_}, {t_, m_}] := Module[
  {A, c, d},
  A = eToDiffE[X, e, {c, n}, {d, m}];
  A = eToS[A, c, {s, n}];
  A = eToS[A, d, {t, m}];
  Expand[A]
  ]

sToDiffS[X_, s_, sn_, tm_] := Module[
  {A, e},
  A = sToE[X, s, e];
  eToDiffS[A, e, sn, tm]
  ]

sToDiffX[X_, s_, xn_, ym_] := Module[
  {A, e},
  A = sToE[X, s, e];
  eToDiffX[A, e, xn, ym]
  ]

(* multiplies each term by h^deg *)
addGrading[X0_, s_, h_] := Module[{X, A, mk},
  X = Expand[X0];
  A = coeffsGeneric[X, s];
  mk[{idx_, cf_}] := cf*h^SumList[idx] * mkSubscript[s, idx];
  Collect[Sum[mk[pair], {pair, A}], h]
  ]

(* tries to convert to a homogeneous expression in a difference alphabet *)
fromDiffHomogeneous[X_, 0   , {u_, n_}, {v_, m_}, s_] := X;
fromDiffHomogeneous[0,  deg_, {u_, n_}, {v_, m_}, s_] := 0;
fromDiffHomogeneous[X_, deg_, {u_, n_}, {v_, m_}, s_] := 
 Module[{parts, unkn, vars, goal, A, lhs, rhs, sols},
  parts = partitions[deg];
  unkn = Table[Subscript[A, p], {p, parts}];
  goal = Sum[Subscript[A, p]*mkSchur[s, p], {p, parts}];
  lhs = sToDiffS[goal, s, {u, n}, {v, m}];
  rhs = Expand[X];
  vars = Complement[Union[Variables[lhs], Variables[rhs]], unkn];
  sols = mySolveAlways[lhs, rhs, vars, unkn];
  If[Length[sols] != 1, 
   Throw["fromDiffHomogeneous: no solution (or too many)"],
   goal /. sols[[1]]
   ]
  ]

(* tries to convert to an expression given in two set of Schur polynomials into Schur polynomials in the difference alphabet *)
fromDiffS[X0_, {u_, n_}, {v_, m_}, s_] := Module[{X, deg, h},
  X = Expand[X0];
  X = addGrading[X, u, h];
  X = addGrading[X, v, h];
  deg = Exponent[X, h];
  Sum[fromDiffHomogeneous[Coefficient[X, h, d], d, {u, n}, {v, m}, s], {d, 0, deg}]
  ]

(* tries to convert to an expression given in two set of roots into Schur polynomials in the difference alphabet *)
fromDiffX[X0_, {a_, n_}, {b_, m_}, ss_, deg_] := Module[{X, Y, h, ser, s, t, i, j},
  X = X0;
  X = X /. Table[Subscript[a, i] -> h*Subscript[a, i], {i, 1, n}];
  X = X /. Table[Subscript[b, j] -> h*Subscript[b, j], {j, 1, m}];
  ser = Normal[Series[X, {h, 0, deg}]];
  Y = xToS[ser /. {h -> 1}, {a, n}, s];
  Y = xToS[Y, {b, m}, t];
  fromDiffS[Y /. {h -> 1}, {s, n}, {t, m}, ss]
  ]

(* ----------------------------------------------- *)

End[]

EndPackage[]
