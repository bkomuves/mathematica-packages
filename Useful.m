
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

mcoeff::usage = "mcoeff[A,{x1,x2,x3},{e1,e2,e3}] extracts the coefficient of x1^e1*x2^e2*x3^e3" 

MCoeffs::usage = "MCoeffs[A,{x1,x2,x3}] returns a list of {coeff,exponent vector} pairs"

mySolveAlways::usage = "mySolveAlways[lhs,rhs,{p1,p2,..},{a1,a2,..}] solves lhs=rhs where p1,p2.. are parameters"

(* ---------- partitions --------- *)

emptyPartQ::usage = "returns true for the empty partition"
EmptyPartQ::usage = "returns true for the empty partition"

dualPart::usage = "returns the dual partition"
DualPart::usage = "returns the dual partition"

samePartitionQ::usage = "samePartitionQ[lambda,mu] returns True if lambda = mu (discarding zeros)"

subPartitionQ::usage = "subPartitionQ[lambda,mu] returns True if mu is a sub-partition of lambda"

toExpoVec  ::usage = "converts a partition to an exponential vector"
fromExpoVec::usage = "converts an exponential vector to a partition"

toExpoForm  ::usage = "convert a partition to a list of {index,exponent} pairs"
fromExpoForm::usage = "crate a partition from a list of {index,exponent} pairs"

partitionWidth ::usage = "width (length) of a partition"
partitionHeight::usage = "height (maximum) of a partition"
partitionWeight::usage = "weight (sum) of a partition"

partitionMonom::usage = "partitionMonom[lam] or partitionMonom[lam,x] the monom x1^e1*x2^e2*... where lam=(1^e1,2^e2...)"

partitionIndex::usage        = "partitionIndex[lam,i] is the same as lam[[i]] but allows overindexing"
partitionReplace1::usage     = "partitionReplace1[lam,i,y] replaces the i-th element with y; but allows overindexing"
partitionReplaceParts::usage = "partitionReplaceParts[lam,{{i1,y1},{i2,y2}}] replaces the i1-th element with y1 etc; but allows overindexing"

blockPartition     ::usage = "blockPartition[height,width] or blockPartition[{height,width}] is the partition (height^width)"
partitionComplement::usage = "partitionComplement[{n,k},lam] is the complement of lambda in the block partition (n^k)"

fitsIntoBlock::usage = "fitsIntoBlock[{H, W}, lam] returns true if the partition lambda fits in to a H x W block"

partitions       ::usage = "partitions[n] returns the list of all partitions of weight n"
allPartitions    ::usage = "allPartitions[n] returns the list of all partitions of weight at most n"
partitionsInBlock::usage = "partitionsInBlock[{n,m},k] returns all partitions (optionally: of weight k) fitting into the block (n^m)"

partitionsInBlock$slow::usage = "another implementation of partitionsInBlock"

subPartitionsOf    ::usage  = "subPartitionsOf[lam] (resp. subParitionsOf[lam,k]) returns the sub-partitions of lam (of weight k)"
subPartitionsOfSize::usage  = "subPartitionsOfSize[lam,k] returns the sub-partitions of lam with weigth k"

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

(* ------------ PARTITIONS -------------- *)

emptyPartQ[part_] := Length[part] == 0
EmptyPartQ[part_] := Length[part] == 0

Clear[dualPart];
dualPart[{}  ] = {};
dualPart[lam_] := 
 Module[{i,m = lam[[1]]}, Table[Length[Select[lam, # >= i &]], {i, 1, m}]]

DualPart[lam_] := dualPart[lam];

samePartitionQ[lam_,mu_] := Module[
  {n = Max[ Length[lam], Length[mu]] },
    And @@ Table[ partitionIndex[lam,i] == partitionIndex[mu,i] , {i,1,n} ]
  ]

blockPartition[{height_Integer,width_Integer}] := replicate[height,width]
blockPartition[ height_Integer,width_Integer ] := replicate[height,width]

partitionComplement[{ht_Integer,wd_Integer},lam_List] := partitionComplement[ ht,wd, lam ] 
partitionComplement[ ht_Integer,wd_Integer ,lam_List] :=
  If[ subPartitionQ[ blockPartition[{ht,wd}] , lam],
    Select[Table[ ht - partitionIndex[lam,wd+1-i] , {i,1,wd}] , #>0& ],
    Throw["partitionComplement: lambda does not fit into the block"]
    ]


Clear[subPartitionsOfSize];
subPartitionsOfSize[lam_, 0] := {{}}
subPartitionsOfSize[{}, n_Integer] := If[n == 0, {{}}, {}];
subPartitionsOfSize[{m_}, n_Integer] := If[n <= m , {{n}}, {}]
subPartitionsOfSize[list_, n_Integer] := 
 Module[{hd = First[list], tl = Rest[list], k, ps, head},
  head[{}] := 0;
  head[xs_List] := First[xs];
  Flatten[
   Table[Table[
     Join[{k}, p], {p, 
      Select[subPartitionsOfSize[tl, n - k], head[#] <= k &]}], {k, 1,
      Min[hd, n]}], 1]
  ]

Clear[subPartitionsOf];
subPartitionsOf[p_List,n_Integer] := subPartitionsOfSize[p,n] 
subPartitionsOf[{}] := {{}};
subPartitionsOf[{n_}] := Join[{{}}, Table[{i}, {i, 1, n}]]
subPartitionsOf[list_] := 
 Module[{hd = First[list], tl = Rest[list], ps}, 
  ps = Flatten[
    Table[Table[
      Join[{i}, xs], {i, Max[1, partitionIndex[xs, 1]], hd}], {xs, 
      subPartitionsOf[tl]}], 1];
  Join[{{}}, ps]]


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
  Module[{es = toExpoVec[lam], n}, n = Length[es]; 
    Select[ Zip[Range[n], es] , Snd[#]>0& ] ]

fromExpoForm[pairs_] := 
  Module[{f},
    f[{i_,e_}] := replicate[i,e];
    reverseSort[Flatten[Map[f,pairs]]]
    ]

partitionMonom[p_   ] := partitionMonom[p,x]
partitionMonom[p_,x_] := Module[{fun,ie},
  fun[{i_, e_}] := Subscript[x,i]^e;
  Product[fun[ie], {ie, toExpoForm[p]}]
  ]

(* width (length) of a partition *)
partitionWidth[list_] := Length[list];

(* height (maximum) of a partition *)
partitionHeight[list_] := Max[0,Max@@list];

(* weight (sum) of a partition *)
partitionWeight[list_] := Plus@@list;

(* is mu a subpartition of lambda? *)
subPartitionQ[lam_,mu_] := Module[
  {n},
  n = Max[Length[lam],Length[mu]];
  And @@ Table[ partitionIndex[lam,i] >= partitionIndex[mu,i] , {i,1,n} ]
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

partitions[n_]    := Combinatorica`Partitions[n]
allPartitions[n_] := Flatten[Table[Combinatorica`Partitions[k],{k,0,n}],1];

(* partitions fitting into a block - TODO: more efficient implementation *)
Clear[partitionsInBlock]
partitionsInBlock[{n_,m_}   ] := subPartitionsOf[blockPartition[{n,m}]  ]
partitionsInBlock[{n_,m_},k_] := subPartitionsOf[blockPartition[{n,m}],k]
partitionsInBlock[ n_,m_    ] := partitionsInBlock[ {n,m}  ]
partitionsInBlock[ n_,m_, k_] := partitionsInBlock[ {n,m},k]

Clear[partitionsInBlock$slow]
partitionsInBlock$slow[{n_, m_}    ] := partitionsInBlock$slow[n, m   ]
partitionsInBlock$slow[{n_, m_}, k_] := partitionsInBlock$slow[n, m, k]
partitionsInBlock$slow[n_, m_, k_] := 
  Select[Combinatorica`Partitions[k, n], Length[#] <= m &] 
partitionsInBlock$slow[n_, m_] := Flatten[
  Table[partitionsInBlock$slow[n, m, k], {k, 0, n*m}], 1]

fitsIntoBlock[ {h_, w_}, part_]             :=  fitsIntoBlock[ h, w, part ]
fitsIntoBlock[ h_Integer, w_Integer, part_] :=  Length[part] <= w && (EmptyPartQ[part] || part[[1]] <= h)

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
