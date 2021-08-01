
(* Functions related to partitions *)

Clear["Partitions`*"];
BeginPackage[ "Partitions`"]
Needs["Combinatorica`"]
Needs["Useful`"]

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

(* ========================================================================== *)

Begin[ "`Private`"]

(* ------------ PARTITIONS -------------- *)

emptyPartQ[part_] := Length[Select[part,#>0&]] == 0
EmptyPartQ[part_] := emptyPartQ[part] 

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

(* ------------------------------------ *)
End[]

EndPackage[]
