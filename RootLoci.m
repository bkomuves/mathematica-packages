
Clear["RootLoci`*"];
BeginPackage[ "RootLoci`"]
Needs["Useful`"]

PartitionClosure::usage       = "PartitionClosure[p] returns the set of coarsenings of a partition"
PartitionClosure$v1::usage    = "PartitionClosure$v1[p] is an alternative (slower?) implementation"

(* ========================================================================== *)

Begin["Private`"]

(* ---------------- closure set of partitions ------------- *)

disjPairs[n_] :=  Flatten[Table[Table[{i, j}, {j, i + 1, n}], {i, 1, n - 1}], 1]

Clear[PartitionClosureStep1]
PartitionClosureStep1[p_] := PartitionClosureStep1[p] = Module[
   {es, r, ijs, is, fun1, fun2, tmp},
   es = toExpoVec[p];
   r = Length[es];
   is = Range[r];
   ijs = disjPairs[r];
   fun1[i_] := If[es[[i]] >= 2,
     partitionReplaceParts[
      es, {{i, es[[i]] - 2}, {2 i, partitionIndex[es, 2 i] + 1}}],
     Null];
   fun2[{i_, j_}] := If[es[[i]] >= 1 && es[[j]] >= 1,
     partitionReplaceParts[
      es, {{i, es[[i]] - 1}, {j, es[[j]] - 1}, {i + j, 
        partitionIndex[es, i + j] + 1}}],
     Null];
   tmp = Select[Join[Map[fun1, is], Map[fun2, ijs]],  Not[SameQ[ #, Null]] &];
   (* Print[tmp]; *)
   Sort[DeleteDuplicates[Map[fromExpoVec, tmp]]]
   ]

Clear[PartitionClosure]
PartitionClosure[p_] := PartitionClosure[p] = Module[
   {layer, stuff},
   layer = PartitionClosureStep1[p];
   stuff = Map[PartitionClosure, layer];
   (*Print[layer];
   Print[stuff];*)
   Union @@ Prepend[stuff, {p}]
   ]

(* alternative implementation (slower?) *)
coarsen[part_, i_, j_] := Module[{part2, idxs, n = Length[part]},
  idxs = Select[Range[n], (# != i && # != j) &];
  part2 = Table[part[[k]], {k, idxs}];
  Reverse[Sort[Join[part2, {part[[i]] + part[[j]]}]]]]

closureLevel1[part_] := 
 Module[{n = Length[part]}, 
  Union[Table[
    coarsen[part, ij[[1]], ij[[2]]], {ij, Subsets[Range[n], {2}]}]]]

closureOfSet[parts_] := 
 Module[{A = Union[parts], B, C}, 
  B = Flatten[Table[closureLevel1[p], {p, A}], 1];
  C = Union[A, B]; If[A == C, A, closureOfSet[C]]]

Clear[PartitionClosure$v1]
PartitionClosure$v1[part_] := PartitionClosure$v1[part] = closureOfSet[{part}]

(* test it:
  Table[SameQ[PartitionClosure[p], PartitionClosure$v1[p]], {p, Partitions[13]}]
*)

(* ---------------------------------------------------- *)

End[]

EndPackage[]
