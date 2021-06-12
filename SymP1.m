
(* Computations in the cohomology and K-theory of the symmetric products of P^1 *)

(* example usage:

nn = 3
dd = 4
vv[i_] := Subscript[v, i]
vs[n_] := Table[vv[i], {i, 1, n}]
A = c + 5 u - (7 + a) u^2 + 6 u^3
cohomDeltaPF[%, ProjSpace[nn,u], vs[dd]]
B1 = cohomPsiPF[%, Table[ProjSpace[nn,vv[i]], {i, 1, dd}],  w]
B2 = cohomOmegaPF[A, ProjSpace[nn, u], dd, w]
Expand[B1 - B2]

*)

Clear["SymP1`*"];
BeginPackage[ "SymP1`"]
Needs["Useful`"]

Cohomology::usage      = "Cohomology (non-equivariant)"
KTheory::usage         = "K-theory (non-equivariant)"
EquivCohomology::usage = "T^2-equivariant cohomology"
EquivKTheory::usage    = "T^2-equivariant K-theory"

ProjSpace::usage = "ProjSpace[n,u] represents the projective space \[DoubleStruckCapitalP]^n with a generator u"
\[DoubleStruckCapitalP]::usage = "\[DoubleStruckCapitalP][n,u] or \[DoubleStruckCapitalP]^n[u] is the projective space of dim n with a generator u"

Dimension::usage = "the dimension of a projective space"
Variable::usage  = "the variable (generator) of a projective space"

cohomNormalize::usage "cohomNormalize[A,space(s)] normalizes the class in the (product of) the given spaces"

cohomOmegaPF::usage = "cohomOmegaPF[A,space,d,var] is the pushforward of A along the replicating map Omega^d in cohomology"
cohomDeltaPF::usage = "cohomDeltaPF[A,space,vars] is the pushforward of A along the diagonal map Delta^d in cohomology"
cohomPsiPF::usage   = "cohomPsiPF[A,spaces,var] is the pushforward of A along the merging map Psi in cohomology"
cohomCollapsePF::usage = "cohomCollapsePF[A,space] is the pushforward of A along the collapsing map pt:P^n->pt"

cohomOmegaPB::usage = "cohomOmegaPB[A,space,d,var] is the pullback of A along the replicating map Omega^d in cohomology"
cohomDeltaPB::usage = "cohomDeltaPB[A,space,vars] the pullback of A along the diagonal map Delta^d in cohomology"
cohomPsiPB::usage   = "cohomPsiPB[A,spaces,var] is the pullback along the merging map Psi in cohomology"

(* ========================================================================== *)

Begin["Private`"]

(* ---------------- SUPPORTED THEORIES ---------------------- *)

KTheory         = Symbol["KTheory"];
Cohomology      = Symbol["Cohomology"];
EquivKTheory    = Symbol["EquivariantKTheory"];
EquivCohomology = Symbol["EquivariantCohomology"];

SetAttributes[ KTheory         , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Cohomology      , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ EquivKTheory    , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ EquivCohomology , {Protected, ReadProtected, Locked} ]; 

SetAttributes[ \[Alpha]      , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ \[Beta]       , {Protected, ReadProtected, Locked} ]; 

SetAttributes[ X , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Y , {Protected, ReadProtected, Locked} ]; 

(* ---------------- Projective spaces ------------------------- *)

Clear[\[DoubleStruckCapitalP]]
\[DoubleStruckCapitalP][n_Integer, u_] := ProjSpace[n, u]
\[DoubleStruckCapitalP] /: \[DoubleStruckCapitalP]^ n_Integer[u_] := ProjSpace[n, u]

ProjSpaceForm[ProjSpace[n_, v_]] := 
 Superscript[\[DoubleStruckCapitalP], n][v]
ProjSpace /: Format[ProjSpace[n_, v_]] := 
 ProjSpaceForm[ProjSpace[n, v]]

Dimension[ProjSpace[n_Integer, _]] := n
Variable[ProjSpace[_Integer, v_]] := v

SetAttributes[ ProjSpace , {Protected, ReadProtected, Locked} ]; 

(* ---------------- COHOMOLOGY ------------------------- *)

(* normal form of cohomology classes *)
Clear[cohomNormalize]
cohomNormalize[A_, ProjSpace[n_, u_]] := Normal[Series[A, {u, 0, n}]]
cohomNormalize[A_, {}] := A
cohomNormalize[A_, {space}] := cohomNormalize[A, space]
cohomNormalize[A_, ps_List] := 
 cohomNormalize[cohomNormalize[A, First[ps]], Rest[ps]]


(* pushforward of A in P^n along Omega^d *)
Clear[cohomOmegaPF]
cohomOmegaPF[A_, ProjSpace[n_, u_], d_] := cohomOmegaPF[A, {u, n}, d, u]
cohomOmegaPF[A_, ProjSpace[n_, u_], 0, v_] := 1
cohomOmegaPF[A_, ProjSpace[n_, u_], 1, v_] := A /. {u -> v}
cohomOmegaPF[A0_, ProjSpace[n_, u_], d_, v_] := 
 Module[{A = Expand[A0]}, 
  Sum[Coefficient[A, u, n - k]*d^k*v^(d*n - k), {k, 0, n}]]

(* pushforward of A in P^n1 x P^n2 x ...along Psi_{n1,n2,...}  *)
Clear[cohomPsiPF]
cohomPsiPF[A_, {}, w_] := 1
cohomPsiPF[A_, {ProjSpace[n_, u_]}, w_] := A /. {u -> w}
cohomPsiPF[A0_, {ProjSpace[n_, u_], ProjSpace[m_, v_]}, w_] := 
 Module[{A = Expand[A0]}, 
  Sum[Coefficient[Coefficient[A, u, n - i], v, m - j]*
    Binomial[i + j, i]*w^(n + m - i - j), {i, 0, n}, {j, 0, m}]]
cohomPsiPF[A_, spaces_List, w_] := Module[
  {p1, p2, ps, y, m, tmp},
  p1 = First[spaces];
  ps = Rest[spaces];
  m = Sum[Dimension[p], {p, ps}];
  p2 = ProjSpace[m, y];
  tmp = cohomPsiPF[A, Rest[spaces], y];
  Collect[Expand[cohomPsiPF[tmp, {p1, p2}, w]], w]
  ]

(* pushforward of A in P^n along Delta^d *)
Clear[cohomDeltaPF]
cohomDeltaPF[A_, space_, d_Integer] := 
 cohomDeltaPF[A, space, 
  Table[Subscript[Variable[space], i], {i, 1, d}]]
cohomDeltaPF[A_, space_, {}] := 1
cohomDeltaPF[A_, ProjSpace[n_, u_], {v_}] := A /. {u -> v}
cohomDeltaPF[A_, ProjSpace[n_, u_], {v1_, v2_}] :=
  Sum[Coefficient[A, u, n - k]*
   Sum[v1^(n - i)*v2^(n - (k - i)), {i, 0, k}], {k, 0, n}]
cohomDeltaPF[A_, ProjSpace[n_, u_], vs_List] := Module[
  {v1, vs2, y, m, tmp},
  v1 = First[vs];
  vs2 = Rest[vs];
  tmp = cohomDeltaPF[A, ProjSpace[n, u], {v1, y}];
  Expand[cohomDeltaPF[tmp, ProjSpace[n, y], vs2]]
  ]

(* pushforward of A in P^n along pt:P^n->pt *)
Clear[cohomCollapsePF]
cohomCollapsePF[A_, ProjSpace[n_, u_]] := Coefficient[Expand[A], u, n]

(* pullback along Delta^d *)
Clear[cohomDeltaPB]
cohomDeltaPB[A_, space_, d_Integer] := 
 cohomDeltaPB[A, space, 
  Table[Subscript[Variable[space], i], {i, 1, d}]]
cohomDeltaPB[A_, ProjSpace[n_, u_], {}] := A
cohomDeltaPB[A_, ProjSpace[n_, u_], {v_}] := A /. {v -> u}
cohomDeltaPB[A_, ProjSpace[n_, u_], vs_List] := Module[{tmp, i},
  tmp = A /. Table[vs[[i]] -> u, {i, 1, Length[vs]}];
  cohomNormalize[tmp, ProjSpace[n, u]]
  ]

(* pullback along Psi *)
Clear[cohomPsiPB]
cohomPsiPB[A_, {}, w_] := A
cohomPsiPB[A_, {ProjSpace[n_, u_]}, w_] := A /. {w -> u}
cohomPsiPB[A_, spaces_List, w_] := 
 A /. {w -> Sum[Variable[p], {p, spaces}]}

(* pullback along Omega *)
Clear[cohomOmegaPB]
cohomOmegaPB[A_, ProjSpace[n_, u_], d_Integer] := 
 cohomOmegaPB[A, ProjSpace[n, u], d, u]
cohomOmegaPB[A_, ProjSpace[n_, u_], 0, v_] := A
cohomOmegaPB[A_, ProjSpace[n_, u_], 1, v_] := A /. {v -> u}
cohomOmegaPB[A0_, ProjSpace[n_, u_], d_Integer, v_] := 
 Module[{A = Expand[A0]},
  Sum[Coefficient[A, v, k]*(d*u)^k, {k, 0, n}]]

(* -------------------------------------- *)

End[]

EndPackage[]
