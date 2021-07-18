
(* Miscellenaous useful functions *)


Clear["Useful`*"];
BeginPackage[ "Useful`"]
Needs["Combinatorica`"]

toExportForm::usage = "convert to export form"

Fst::usage = "Fst[{x,y}] returns the first element of a pair"
Snd::usage = "Snd[{x,y}] returns the second element of a pair"

replicate::usage = "replicate[A,n] replicates A into a list, n times"

extendListWithZeros::usage = "extendListWithZeros[L,n] add zeros to the end of L till it reaches length n"

SumList::usage = "SumList[L] returns the sum of the elements of the list L"

Zip::usage     = "Zip[as,bs] creates a list of pairs (triples) from two (three) lists"
ZipWith::usage = "Zip[f,as,bs] creates a list of f[a,b] applications"

IsEmptyList::usage = "returns True for the empty list, False otherwise"
IsNotEmptyList::usage = "returns False for the empty list, True otherwise"

unique::usage = "make a set from a list by removing duplicates and sorting"

reverseSort::usage = "like Sort but in reverse order"

fallingFactorial::usage = "fallingFactorial[x,k] computes x(x-1)(x-2)...(x-k+1)"

ogf2egf::usage = "ogf2egf[A,x,t] converts ordinary generating function to exponential"
egf2ogf::usage = "egf2ogf[A,t,x] converts exponential generating function to ordinary"

mcoeff::usage = "mcoeff[A,{x1,x2,x3},{e1,e2,e3}] extracts the coefficient of x1^e1*x2^e2*x3^e3" 

MCoeffs::usage = "MCoeffs[A,{x1,x2,x3}] returns a list of {coeff,exponent vector} pairs"

mySolveAlways::usage = "mySolveAlways[lhs,rhs,{p1,p2,..},{a1,a2,..}] solves lhs=rhs where p1,p2.. are paramteres"

EmptyPartQ::usage = "returns true for the empty partition"

DualPart::usage = "returns the dual partition"

toExpoVec::usage = "converts a partition to an exponential vector"
fromExpoVec::usage = "converts an exponential vector to a partition"

toExpoForm::usage = "convert a partition to {index,exponent} pairs"

partitionMonom::usage = "partitionMonom[lam] or partitionMonom[lam,x] the monom x1^e1*x2^e2*..."

partitionIndex::usage = "partitionMonom[lam,i] is the same as lam[[i]] but allows overindexing"
partitionReplace1::usage = "partitionReplace1[lam,i,y] replaces the i-th element with y; but allows overindexing"
partitionReplaceParts::usage = "partitionReplaceParts[lam,{{i1,y1},{i2,y2}}] replaces the i1-th element with y1 etc; but allows overindexing"

FitsInto::usage = "FitsInto[w, h, lam] returns true if the partition lambda fits in to a W x H box"
GrassmannBetti::usage = "GrassmannBetti[n, k, d] return 2d-th Betti number of the grassmannian of Gr_k(C^n)"

SortingPermutation::usage = "returns a permutation sorting a list" 
ReverseSortingPermutation::usage = "returns a permutation sorting a list in reverse order"

(* ========================================================================== *)

Begin[ "`Private`"]

toExportForm[X_] := 
 ToString[X, FormatType -> InputForm, PageWidth -> Infinity, 
  TotalWidth -> Infinity]

(* --------------- LISTS ----------------- *)

Fst[pair_] := pair[[1]]
Snd[pair_] := pair[[2]]

replicate[A_, n_] := Module[ {i}, Table[A, {i, 1, n}] ]

extendListWithZeros[L_, n_] := Join[ L , replicate[0,n - Length[L]] ]

reverseSort[p_] := Sort[p, Greater]

SumList[L_] := Module[{x},Sum[x, {x, L}]]

unique[L_] := Sort[DeleteDuplicates[%]]

Zip[as_, bs_, cs_] := Zip3[as,bs,cs]
Zip[as_, bs_] := Module[{ la=Length[as] , lb=Length[bs] , n } , 
  If[la==lb , MapThread[{#1, #2} &, {as, bs}] ,
    n = Min[la,lb];
    MapThread[{#1, #2} &, {Take[as,n], Take[bs,n]}] ] ]

ZipWith[f_, as_, bs_     ] := Map[ f[#[[1]], #[[2]]] &         , Zip[as, bs]     ]
ZipWith[f_, as_, bs_, cs_] := Map[ f[#[[1]], #[[2]], #[[3]]] & , Zip[as, bs, cs] ]

Zip3[as_, bs_, cs_] := Module[{ la=Length[as] , lb=Length[bs] , lc=Length[cs], n } , 
  If[la==lb && lb==lc , MapThread[{#1, #2, #3} &, {as, bs, cs}] ,
    n = Min[la,lb,lc];
    MapThread[{#1, #2, #3} &, {Take[as,n], Take[bs,n], Take[cs,n]}] ] ]

IsEmptyList[L_] := SameQ[L, {}]
IsNotEmptyList[L_] := UnsameQ[L, {}]

(* --------------- MISC ----------------- *)

fallingFactorial[X_, k_] := Module[ {i} , Product[X - i, {i, 0, k - 1}] ]

(* ----------------- EGF vs OGF ------------------------ *)

egf2ogf[G_, t_, x_] := Module[{xxx},Expand[(LaplaceTransform[G, t, xxx]*xxx) /. {xxx -> 1/x}]]
ogf2egf[F_, x_, t_] := Module[{xxx},Expand[InverseLaplaceTransform[(F*x) /. {x -> 1/xxx}, xxx, t]]]

(*
Sum[Subscript[a, i]*y^i, {i, 0, 6}]
ogf2egf[%, y, s]
egf2ogf[%, s, y]
*)

(* ------------- COEFFICIENT EXTRACTION ------------ *)

mcoeff[expr_, vars_, expos_] := 
 If[Length[vars] == 0, expr, 
  mcoeff[Coefficient[expr, First[vars], First[expos]], Rest[vars], 
   Rest[expos]]]

(* helper function *)
InsertExpo[e_, pair_] := {pair[[1]], Join[{e}, pair[[2]]]}

(* returns a list of {coeff,exponent vector} pairs *)
MCoeffs[expr_, vars_] := 
 If[Length[vars] == 0, {{expr, {}}}, 
  Module[{var, list, ii, jj, kk, prev, new, rest, pair, cf, expo, 
    big}, var = First[vars]; rest = Rest[vars];
   list = CoefficientList[expr, var]; kk = Length[list]; big = {};
   For[ii = 1, ii <= kk, ii++, cf = list[[ii]];(*Print[{ii,cf}];*)
    If[cf =!= 0, prev = MCoeffs[cf, rest];
     new = Map[InsertExpo[ii - 1, #] &, prev];(*Print["new",new];
     Print["big old",big];*)big = Join[big, new];(*Print["big new",
     big]*), 0]]; big]]

(* -------------------- SOLVE ALWAYS ---------------- *)

(* find given exponent in MCoeffs; basically lookup of (value,key) pairs *)
FindSnd[list_, snd_, deflt_] := 
 Module[{sel = Select[list, Snd[#] == snd &]}, 
  If[Length[sel] == 0, deflt, Fst[sel[[1]]]]]

(* actually working SolveAlways,for multivariate polynomials at least *)
mySolveAlways[P_, Q_, vars_] := 
 Module[{list1, list2, expos, eqs}, list1 = MCoeffs[P, vars];
  list2 = MCoeffs[Q, vars];
  expos = DeleteDuplicates[Join[Map[Snd, list1], Map[Snd, list2]]];
  eqs = Table[
    FindSnd[list1, expo, 0] == FindSnd[list2, expo, 0], {expo, expos}];
  (*Print[eqs];*)
  Solve[eqs]]

mySolveAlways[P_, Q_, vars_, solveVars_] := 
 Module[{list1, list2, expos, eqs}, list1 = MCoeffs[P, vars];
  list2 = MCoeffs[Q, vars];
  expos = DeleteDuplicates[Join[Map[Snd, list1], Map[Snd, list2]]];
  eqs = Table[
    FindSnd[list1, expo, 0] == FindSnd[list2, expo, 0], {expo, expos}];
  (*Print[eqs];*)
  Solve[eqs, solveVars]]

(* ------------ PARTITIONS -------------- *)

EmptyPartQ[part_] := Length[part] == 0

Clear[DualPart];
DualPart[{}  ] = {};
DualPart[lam_] := 
 Module[{i,m = lam[[1]]}, Table[Length[Select[lam, # >= i &]], {i, 1, m}]]

Clear[toExpoVec];
toExpoVec[{}   ] = {};
toExpoVec[part_] := 
 Module[{j,k = Max[part]}, 
  Table[Length[Select[part, # == j &]], {j, 1, k}]]

Clear[fromExpoVec];
fromExpoVec[{} ] = {};
fromExpoVec[es_] := Module[ {i},
  Reverse[Join @@ Table[replicate[i, es[[i]]], {i, 1, Length[es]}]] ]

toExpoForm[lam_] := 
  Module[{es = toExpoVec[lam], n}, n = Length[es]; Zip[Range[n], es]]

partitionMonom[p_   ] := partitionMonom[p,x]
partitionMonom[p_,x_] := Module[{fun,ie},
  fun[{i_, e_}] := Subscript[x,i]^e;
  Product[fun[ie], {ie, toExpoForm[p]}]
  ]

(* allows over-indexing *)
partitionIndex[list_, i_] := If[i <= Length[list], list[[i]], 0]

(* allows over-indexing *)
partitionReplace1[list_, i_, y_] := Module[{n = Length[list]},
  If[i <= n,
   ReplacePart[list, i -> y],
   ReplacePart[Join[list, replicate[0, i - n]], i -> y]
   ]]

partitionReplace$pair[list_, {i_, y_}] := partitionReplace1[list, i, y]

partitionReplaceParts[list_, {}] := list
partitionReplaceParts[list_, iys_] := If[ ListQ[iys] , 
  partitionReplaceParts[partitionReplace$pair[list, First[iys]], Rest[iys]] , Throw["partitionReplaceParts: expecting a list"]]

(* ---------- BETTI NUMBERS OF GRASSMANNIANS ------------ *)

FitsInto[w_, h_, part_] :=  Length[part] <= w && (EmptyPartQ[part] || part[[1]] <= h)

GrassmannBetti[n_, k_, j_] := Length[Select[Combinatorica`Partitions[j], FitsInto[n - k, k, #] &]]

(* ----------- SORTING PERMUTATION ------------------------------ *)

SortingPermutation[L_] := 
 Module[{n = Length[L], A}, A = Zip[Range[n], L]; A = SortBy[A, Snd]; 
  Map[Fst, A]]

ReverseSortingPermutation[L_] :=
  Module[{n = Length[L], A}, A = Zip[Range[n], L]; 
  A = Reverse[SortBy[A, Snd]]; Map[Fst, A]]

(*
testList = {7, 5, 99, 3, 4, 2, 100, 10}
testPerm = ReverseSortingPermutation[testList]
invPerm = InversePermutation[testPerm]
testSorted = Table[testList[[i]], {i, testPerm}]
Permute[testList, testPerm]
testCheck = Table[testSorted[[i]], {i, invPerm}]
Permute[testSorted, invPerm]
*)

End[]

EndPackage[]
