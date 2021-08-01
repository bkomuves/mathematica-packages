
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
MaxList::usage = "SumList[L] returns the maximum of the elements of the list L (for empty list, minus infinity)"

Zip::usage     = "Zip[as,bs] creates a list of pairs (triples) from two (three) lists"
ZipWith::usage = "Zip[f,as,bs] creates a list of f[a,b] applications"

IsEmptyList::usage = "returns True for the empty list, False otherwise"
IsNotEmptyList::usage = "returns False for the empty list, True otherwise"

unique::usage = "make a set from a list by removing duplicates and sorting"

reverseSort::usage = "like Sort but in reverse order"

fallingFactorial::usage = "fallingFactorial[x,k] computes x(x-1)(x-2)...(x-k+1)"

ogf2egf::usage = "ogf2egf[A,x,t] converts ordinary generating functions to exponential"
egf2ogf::usage = "egf2ogf[A,t,x] converts exponential generating functions to ordinary"

mcoeff ::usage = "mcoeff[A,{x1,x2,x3},{e1,e2,e3}] extracts the coefficient of x1^e1*x2^e2*x3^e3" 
MCoeffs::usage = "MCoeffs[A,{x1,x2,x3}] returns a list of {coeff,exponent vector} pairs"

mySolveAlways::usage = "mySolveAlways[lhs,rhs,{p1,p2,..},{a1,a2,..}] solves lhs=rhs where p1,p2.. are parameters"

(* -------- random stuff ---------- *)

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
MaxList[L_] := Max@@L

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

(* ---------- BETTI NUMBERS OF GRASSMANNIANS ------------ *)

(* TODO: more efficient impl *)
GrassmannBetti[n_, k_, j_]   := Length[Select[Combinatorica`Partitions[j], fitsIntoBlock[k , n - k, #] &]]

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
