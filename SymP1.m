
(* Computations in the cohomology and K-theory of the symmetric products of P^1 *)


(* conventions:

- in cohomology, the generator is -c_1(L), -1 times the first Chern class of the tautological line bundle on P^n
- in K-theory, the generator is the class of the tautological line bundle on P^n
- in equivariant cohomology, \[Alpha] and \[Beta] are the universal Chern roots
- in equivariant K-theory, the corresponding universal bundles are X and Y, with the convention that \[Alpha] = -c_1(X) and \[Beta] = -c_1(Y)

*)


(* example usage:

<<SymP1`

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
AbstractTheory::usage  = "Relative Grothendieck group of varietes"

ProjSpace::usage = "ProjSpace[n,u] represents the projective space \[DoubleStruckCapitalP]^n with a generator u"
\[DoubleStruckCapitalP]::usage = "\[DoubleStruckCapitalP][n,u] or \[DoubleStruckCapitalP]^n[u] is the projective space of dim n with a generator u"

Dimension::usage = "the dimension of a projective space"
Variable::usage  = "the variable (generator) of a projective space"

(* ----- theory-polymorphic API ---- *)

relation  ::usage = "For example relation[EquivKTheory,ProjSpace[n,L] returns the relation in the equivariant K-theory of P^n"

normalize ::usage = "normalize[theory,A,space(s)] normalizes the class in the (product of) the given spaces"

OmegaPF   ::usage = "OmegaPF[theory,A,space,d,var] is the pushforward of A along the replicating map Omega^d in cohomology"
DeltaPF   ::usage = "DeltaPF[theory,A,space,vars] is the pushforward of A along the diagonal map Delta^d in cohomology"
PsiPF     ::usage = "PsiPF[theory,A,spaces,var] is the pushforward of A along the merging map Psi in cohomology"
CollapsePF::usage = "CollapsePF[theory,A,space] is the pushforward of A along the collapsing map pt:P^n->pt"

OmegaPB   ::usage = "cohomOmegaPB[theory,A,space,d,var] is the pullback of A along the replicating map Omega^d in cohomology"
DeltaPB   ::usage = "cohomDeltaPB[theory,A,space,vars] the pullback of A along the diagonal map Delta^d in cohomology"
PsiPB     ::usage = "cohomPsiPB[theory,A,spaces,var] is the pullback along the merging map Psi in cohomology"

(* ----- cohomology (non-equivariant) ----- *)

cohomNormalize ::usage = "cohomNormalize[A,space(s)] normalizes the class in the (product of) the given spaces"

cohomOmegaPF   ::usage = "cohomOmegaPF[A,space,d,var] is the pushforward of A along the replicating map Omega^d in cohomology"
cohomDeltaPF   ::usage = "cohomDeltaPF[A,space,vars] is the pushforward of A along the diagonal map Delta^d in cohomology"
cohomPsiPF     ::usage = "cohomPsiPF[A,spaces,var] is the pushforward of A along the merging map Psi in cohomology"
cohomCollapsePF::usage = "cohomCollapsePF[A,space] is the pushforward of A along the collapsing map pt:P^n->pt"

cohomOmegaPB   ::usage = "cohomOmegaPB[A,space,d,var] is the pullback of A along the replicating map Omega^d in cohomology"
cohomDeltaPB   ::usage = "cohomDeltaPB[A,space,vars] the pullback of A along the diagonal map Delta^d in cohomology"
cohomPsiPB     ::usage = "cohomPsiPB[A,spaces,var] is the pullback along the merging map Psi in cohomology"

(* ----- K-theory (non-equivariant) ----- *)

ktheoryNormalize ::usage = "ktheoryNormalize[A,space(s)] normalizes the class in the (product of) the given spaces"

ktheoryOmegaPF   ::usage = "ktheoryOmegaPF[A,space,d,var] is the pushforward of A along the replicating map Omega^d in K-theory"
ktheoryDeltaPF   ::usage = "ktheoryDeltaPF[A,space,vars] is the pushforward of A along the diagonal map Delta^d in K-theory"
ktheoryPsiPF     ::usage = "ktheoryPsiPF[A,spaces,var] is the pushforward of A along the merging map Psi in K-theory"
ktheoryCollapsePF::usage = "ktheoryCollapsePF[A,space] is the pushforward of A along the collapsing map pt:P^n->pt"

ktheorySingleDeltaBangL ::usage = "ktheorySingleDeltaBangL[n,d,k]"
ktheorySingleOmegaBangL ::usage = "ktheorySingleOmegaBangL[n,d,k]"
ktheorySinglePsiBangL   ::usage = "ktheorySinglePsiBangL[n,m,i,j]"

(* ----- equivariant cohomology ----- *)

equivCohomNormalize ::usage = "equivCohomNormalize[A,space(s)] normalizes the class in the (product of) the given spaces"

equivCohomOmegaPF   ::usage = "equivCohomOmegaPF[A,space,d,var] is the pushforward of A along the replicating map Omega^d in equivariant cohomology"
equivCohomDeltaPF   ::usage = "equivCohomDeltaPF[A,space,vars] is the pushforward of A along the diagonal map Delta^d in equivariant cohomology"
equivCohomPsiPF     ::usage = "equivCohomPsiPF[A,spaces,var] is the pushforward of A along the merging map Psi in equivariant cohomology"
equivCohomCollapsePF::usage = "equivCohomCollapsePF[A,space] is the pushforward of A along the collapsing map pt:P^n->pt"

(* ========================================================================== *)

Begin["`Private`"]
(* Print[Context[]]; *)

(* ---------------- SUPPORTED THEORIES ---------------------- *)

KTheory             = Symbol["KTheory"];
Cohomology          = Symbol["Cohomology"];
EquivKTheory        = Symbol["EquivariantKTheory"];
EquivCohomology     = Symbol["EquivariantCohomology"];
AbstractTheory      = Symbol["AbstractTheory"];

(*
EquivCohomologyT2   = Symbol["EquivariantCohomologyT2"];
EquivCohomologyGL2  = Symbol["EquivariantCohomologyGL2"];
*)

SetAttributes[ KTheory         , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Cohomology      , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ EquivKTheory    , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ EquivCohomology , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ AbstractTheory  , {Protected, ReadProtected, Locked} ]; 

ClearAll[ Global`\[Alpha] , Global`\[Alpha] , Global`X , Global`Y  ]

(* Chern roots *)
SetAttributes[ Global`\[Alpha] , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Global`\[Beta]  , {Protected, ReadProtected, Locked} ]; 

(* X = Exp[-\[Alpha]], Y = Exp[-\[Beta]]  *)
SetAttributes[ Global`X , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Global`Y , {Protected, ReadProtected, Locked} ]; 

\[Alpha] = Global`\[Alpha] ;
\[Beta]  = Global`\[Beta]  ;

X = Global`X;
Y = Global`Y;

ClearAll[ 
  Global`\[CapitalOmega]   , 
  Global`\[CapitalDelta]   , 
  Global`\[CapitalPsi]     , 
  Global`\[ScriptCapitalA] , 
  Global`\[ScriptCapitalS]   ];

\[CapitalDelta]   = Global`\[CapitalDelta]   
\[CapitalOmega]   = Global`\[CapitalOmega]   
\[CapitalPsi]     = Global`\[CapitalPsi]   
\[ScriptCapitalA] = Global`\[ScriptCapitalA]
\[ScriptCapitalS] = Global`\[ScriptCapitalS]

(* --------- internal free variables --------- *)

(* these are used in precomputed tables *)
SetAttributes[ H , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ L , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ U , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ C , {Protected, ReadProtected, Locked} ]; 

HH = SymP1`Private`H;
LL = SymP1`Private`L;
UU = SymP1`Private`U;
CC = SymP1`Private`C;

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

(* -------------- theory-polymorphic API ----------------------- *)

normalize[ Cohomology      , A_ , space_ ] := cohomNormalize        [ A , space ]
normalize[ KTheory         , A_ , space_ ] := ktheoryNormalize      [ A , space ]
normalize[ EquivCohomology , A_ , space_ ] := equivCohomNormalize   [ A , space ]
normalize[ EquivKTheory    , A_ , space_ ] := equivKTheoryNormalize [ A , space ]

OmegaPF[ Cohomology      , A_ , space_ , d_ , var_ ] := cohomOmegaPF        [ A , space , d , var ]
OmegaPF[ KTheory         , A_ , space_ , d_ , var_ ] := ktheoryOmegaPF      [ A , space , d , var ]
OmegaPF[ EquivCohomology , A_ , space_ , d_ , var_ ] := equivCohomOmegaPF   [ A , space , d , var ]
OmegaPF[ EquivKTheory    , A_ , space_ , d_ , var_ ] := equivKTheoryOmegaPF [ A , space , d , var ]

DeltaPF[ Cohomology      , A_ , space_ , vars_ ] := cohomDeltaPF        [ A , space , vars ]
DeltaPF[ KTheory         , A_ , space_ , vars_ ] := ktheoryDeltaPF      [ A , space , vars ]
DeltaPF[ EquivCohomology , A_ , space_ , vars_ ] := equivCohomDeltaPF   [ A , space , vars ]
DeltaPF[ EquivKTheory    , A_ , space_ , vars_ ] := equivKTheoryDeltaPF [ A , space , vars ]

PsiPF[ Cohomology      , A_ , spaces_ , var_ ] := cohomPsiPF        [ A , spaces , var ]
PsiPF[ KTheory         , A_ , spaces_ , var_ ] := ktheoryPsiPF      [ A , spaces , var ]
PsiPF[ EquivCohomology , A_ , spaces_ , var_ ] := equivCohomPsiPF   [ A , spaces , var ]
PsiPF[ EquivKTheory    , A_ , spaces_ , var_ ] := equivKTheoryPsiPF [ A , spaces , var ]

CollpasePF[ Cohomology      , A_ , space_ ] := cohomCollapsePF        [ A , space ]
CollpasePF[ KTheory         , A_ , space_ ] := ktheoryCollapsePF      [ A , space ]
CollpasePF[ EquivCohomology , A_ , space_ ] := equivCohomCollapsePF   [ A , space ]
CollpasePF[ EquivKTheory    , A_ , space_ ] := equivKTheoryCollapsePF [ A , space ]

OmegaPB[ Cohomology      , A_ , space_ , d_ , var_ ] := cohomOmegaPB        [ A , space , d , var ]
OmegaPB[ KTheory         , A_ , space_ , d_ , var_ ] := ktheoryOmegaPB      [ A , space , d , var ]
OmegaPB[ EquivCohomology , A_ , space_ , d_ , var_ ] := equivCohomOmegaPB   [ A , space , d , var ]
OmegaPB[ EquivKTheory    , A_ , space_ , d_ , var_ ] := equivKTheoryOmegaPB [ A , space , d , var ]

DeltaPB[ Cohomology      , A_ , space_ , vars_ ] := cohomDeltaPB        [ A , space , vars ]
DeltaPB[ KTheory         , A_ , space_ , vars_ ] := ktheoryDeltaPB      [ A , space , vars ]
DeltaPB[ EquivCohomology , A_ , space_ , vars_ ] := equivCohomDeltaPB   [ A , space , vars ]
DeltaPB[ EquivKTheory    , A_ , space_ , vars_ ] := equivKTheoryDeltaPB [ A , space , vars ]

PsiPB[ Cohomology      , A_ , spaces_ , var_ ] := cohomPsiPB        [ A , spaces , var ]
PsiPB[ KTheory         , A_ , spaces_ , var_ ] := ktheoryPsiPB      [ A , spaces , var ]
PsiPB[ EquivCohomology , A_ , spaces_ , var_ ] := equivCohomPsiPB   [ A , spaces , var ]
PsiPB[ EquivKTheory    , A_ , spaces_ , var_ ] := equivKTheoryPsiPB [ A , spaces , var ]

(* ---------------- Relations ------------------------ *)

weight[n_,i_]  := (n-i)*Global`\[Alpha] + i*Global`\[Beta]
kWeight[n_,i_] := Global`X^(n-i) * Global`Y^i

Clear[relation]

relation[ Cohomology, ProjSpace[n_,u_] ] := u^(n+1) 

relation[ KTheory , ProjSpace[n_,L_] ] :=
   relation[ KTheory , ProjSpace[n,L] ] = Expand [ (1-L)^(n+1) ]

relation[ EquivCohomology, ProjSpace[n_,u_] ] := 
  relation[ EquivCohomology, ProjSpace[n,u] ] =
    Expand[ Product[u+weight[n,i], {i,0,n} ] ]

relation[ EquivKTheory   , ProjSpace[n_,L_] ] :=
  relation[ EquivKTheory   , ProjSpace[n,L] ] =
    Expand[ Product[1-L*kWeight[n,i], {i,0,n} ] ]

(* ---------------- COHOMOLOGY ------------------------- *)

(* normal form of cohomology classes *)
Clear[cohomNormalize]
cohomNormalize[A_, ProjSpace[n_, u_]] := Normal[Series[A, {u, 0, n}]]
cohomNormalize[A_, {}] := A
cohomNormalize[A_, {space_}] := cohomNormalize[A, space]
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

(* ---------------- K-theory ---------------------- *)

(* normal form of K-theory classes *)
Clear[ktheoryNormalize]
ktheoryNormalize[A_, ProjSpace[n_, L_]] := Module[ {H} , 
  Expand[ Normal[Series[A/.{L->1-H}, {H, 0, n}]] /. {H->1-L} ] ]
ktheoryNormalize[A_, {}] := A
ktheoryNormalize[A_, {space_}] := ktheoryNormalize[A, space]
ktheoryNormalize[A_, ps_List] := 
 ktheoryNormalize[ktheoryNormalize[A, First[ps]], Rest[ps]]

(* schur polynomials of partitions of length = 2 *) 
schur2[a_, b_, x1_, x2_] := 
 If[b <= a, Sum[x1^(a - i)*x2^(b + i), {i, 0, a - b}],
  If[b == a + 1, 0, -schur2[b - 1, a + 1, x1, x2]]]

(* pushforward of H^k in P^n along Delta^d *)
ktheorySingleDeltaBangH[n_, 2, k_] := 
 Sum[(-1)^i*
   schur2[n, k + i, Subscript[HH, 1],  Subscript[HH, 2]], {i, 0, Min[1, n - k]}]

(* pushforward of L^k in P^n along Delta^d *)
Clear[ktheorySingleDeltaBangL];
ktheorySingleDeltaBangL[n_, 2, k_] := Module[{A, B, H, d},
  d = 2;
  A = Expand[(1 - H)^k];
  B = Sum[ Coefficient[A, H, i]*ktheorySingleDeltaBangH[n, d, i], {i, 0,  k}];
  Expand[B /. Table[Subscript[HH, i] -> (1 - Subscript[LL, i]), {i, 1, d}]]
  ]

(* pushforward of A in P^n along Delta^d *)
Clear[ktheoryDeltaPF]
ktheoryDeltaPF[A_, space_, d_Integer] := 
 ktheoryDeltaPF[A, space, 
  Table[Subscript[Variable[space], i], {i, 1, d}]]
ktheoryDeltaPF[A_, space_, {}] := 1
ktheoryDeltaPF[A_, ProjSpace[n_, L_], {L1_}] := A /. {L -> L1}
ktheoryDeltaPF[A_, ProjSpace[n_, L_], {L1_, L2_}] := 
 Sum[Coefficient[A, L, k] * (ktheorySingleDeltaBangL[n, 2, k] /.
   {Subscript[LL, 1] -> L1, 
    Subscript[LL, 2] -> L2}), {k, 0, n}]
ktheoryDeltaPF[A_, ProjSpace[n_, u_], vs_List] := 
 Module[{v1, vs2, y, m, tmp}, v1 = First[vs];
  vs2 = Rest[vs];
  tmp = ktheoryDeltaPF[A, ProjSpace[n, u], {v1, y}];
  Expand[ktheoryDeltaPF[tmp, ProjSpace[n, y], vs2]]]


(* pushforward of H1^i*H2^j *)
ktheorySinglePsiBangH[n_, m_, i_, j_] := 
 Module[{p = Min[n - i, m - j]}, 
  Sum[(-1)^k*Binomial[p, k]*Binomial[n + m - i - j - k, p]*
    HH^(i + j + k), {k, 0, p}]]

(* pushforward of L1^i*L2^j *)
ClearAll[ktheorySinglePsiBangL]
ktheorySinglePsiBangL[n_, m_, ii_, jj_] :=
 ktheorySinglePsiBangL[n, m, ii, jj] =
  Module[{H1, H2, A, B},
   A = Expand[(1 - H1)^ii*(1 - H2)^jj];
   B = Sum[
     mcoeff[A, {H1, H2}, {i, j}]*
      ktheorySinglePsiBangH[n, m, i, j], {i, 0, ii}, {j, 0, jj}];
   Expand[B /. {HH -> (1 - LL)}]
   ]

(* pushforward of A in P^n1 x P^n2 x ...along Psi_{n1,n2,...}  *)
Clear[ktheoryPsiPF]
ktheoryPsiPF[A_ , {}, L_] := 1
ktheoryPsiPF[A_ , {ProjSpace[n_, L1_]}, L_] := A /. {L1 -> L}
ktheoryPsiPF[A0_, {ProjSpace[n_, L1_], ProjSpace[m_, L2_]}, L_] := 
 Module[{A = Expand[A0]}, 
  Sum[Coefficient[Coefficient[A, L1, i], L2, j] * (ktheorySinglePsiBangL[n, m, i, j] /. {LL -> L}) , {i, 0, n}, {j, 0, m}]]
ktheoryPsiPF[A_, spaces_List, w_] := Module[
  {p1, p2, ps, y, m, tmp},
  p1 = First[spaces];
  ps = Rest[spaces];
  m = Sum[Dimension[p], {p, ps}];
  p2 = ProjSpace[m, y];
  tmp = ktheoryPsiPF[A, Rest[spaces], y];
  Collect[Expand[ktheoryPsiPF[tmp, {p1, p2}, w]], w]
  ]

ktheoryCollapsePF[A_,ProjSpace[n_,L_]] := Coefficient[Expand[A],L,0]

ClearAll[ktheorySingleOmegaBangL]
ktheorySingleOmegaBangL[n_, d_, k_] := 
 ktheorySingleOmegaBangL[n, d, k] = Module[
   {A, B, spaces, vars, L},
   vars = Table[Subscript[LL, i], {i, 1, d}];
   spaces = Table[ProjSpace[n, Subscript[LL, i]], {i, 1, d}];
   A = ktheoryDeltaPF[L^k, ProjSpace[n, L], vars];
   Expand[ktheoryPsiPF[A, spaces, LL]]
   ]

(* pushforward of A in P^n along Omega^d *)
Clear[ktheoryOmegaPF]
ktheoryOmegaPF[A_, space_           , d_Integer] :=  ktheoryOmegaPF[A, space, d , Variable[space] ]
ktheoryOmegaPF[A_, space_           , 0  , Lout_] := 1
ktheoryOmegaPF[A_, ProjSpace[n_, L_], 1  , Lout_] := A /. {L -> Lout}
ktheoryOmegaPF[A_, ProjSpace[n_, L_], d_ , Lout_] := 
  Sum[ Coefficient[A, L, k] * (ktheorySingleOmegaBangL[n, d, k] /. {LL -> Lout}) , {k, 0, n} ]

(* --------------------- EQUIVARIANT COHOMOLOGY ---------------------- *)

Begin["`EquivCohom`"]
(* Print[Context[]]; *) 

SetAttributes[ U , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ c , {Protected, ReadProtected, Locked} ]; 

UU = SymP1`Private`EquivCohom`U;
CC = SymP1`Private`EquivCohom`c;

CC1  = Subscript[CC, 1];
CC2  = Subscript[CC, 2];
CC12 = { CC1, CC2 };

UU1  = Subscript[UU, 1];
UU2  = Subscript[UU, 2];

\[Alpha] = Global`\[Alpha]; 
\[Beta]  = Global`\[Beta] ;

weight[n_,i_]  := (n-i)*Global`\[Alpha] + i*Global`\[Beta]

rel[n_] := Product[UU + weight[n, i], {i, 0, n}]

Clear[upow$nplus1];
upow$nplus1[n_] := upow$nplus1[n] = Expand[ UU^(n + 1) - rel[n]]

(* (unefficiently) normalize a/polynomial/in u *)
Clear[normalizeSlow];
normalizeSlow[n_, X0_] := 
 Module[{X = Expand[X0], m}, m = Exponent[X, UU];
  If[m <= n, X, 
   normalizeSlow[n, X /. {UU^m -> UU^(m - n - 1)*upow$nplus1[n]}]]]

(*table of normalized powers of u*)
Clear[UPowAB];
UPowAB[n_, k_] := UPowAB[n, k] = Expand[normalizeSlow[n, UU^k]]

normalizeVarAB[n_, u_, X_] := 
 Module[{m = Exponent[X, u]}, 
  Expand[Sum[
    Coefficient[X, u, k]*(UPowAB[n, k] /. {UU -> u}), {k, 0, m}]]]

(* Psi_*(u^k*v^l) *)
Clear[psiStarAB];
psiStarAB[n_, m_, 0 , 0 ] := psiStarAB[n, m, 0, 0] = Binomial[n + m, n]
psiStarAB[n_, m_, k_, l_] := psiStarAB[n, m, k, l] = 
  Module[{A, B, AB, F, R, IJ, sel, fun, cft, u, v}, 
   A = Product[u + weight[n, i], {i, 0, k - 1}];
   B = Product[v + weight[m, j], {j, 0, l - 1}];
   AB = Expand[A*B];
   sel[{i_, j_}] := (i < k) || (j < l);
   cft[{i_, j_}] := Coefficient[Coefficient[AB, u, i], v, j];
   fun[{i_, j_}] := psiStarAB[n, m, i, j];
   IJ = Flatten[Table[{i, j}, {i, 0, k}, {j, 0, l}], 1];
   IJ = Select[IJ, sel];
   F = Binomial[n + m - k - l, n - k] * Product[UU + weight[n + m, i], {i, 0, k + l - 1}];
   R = Sum[cft[ij]*fun[ij], {ij, IJ}];
   Expand[F - R]
  ]

Clear[psiStarChern];
psiStarChern[n_, m_, k_, l_] := 
 psiStarChern[n, m, k, l] = 
  SymmetricReduction[  psiStarAB[n, m, k, l], {\[Alpha], \[Beta]}, CC12 ][[1]] 
     
notk[n_, k_] := Select[Table[i, {i, 0, n}], # != k &]

(* class of the diagonal in P^n x P^n *)
Clear[deltaClassAB];
deltaClassAB[n_] := deltaClassAB[n] = 
  Expand[Factor[  
    Sum[Product[(UU1 +  weight[n, i])*(UU2 + weight[n, i]) / (weight[n, i] - weight[n, k]), {i, notk[n, k]}], {k, 0, n}]]]

Clear[deltaClassChern];
deltaClassChern[n_] := deltaClassChern[n] = 
  Expand[SymmetricReduction[
     deltaClassAB[n], {\[Alpha], \[Beta]}, CC12 ][[1]]]

(* pushforward of u^k along Delta:P^n->P^n x P^n *)
Clear[deltaStarAB];
deltaStarAB[n_, 0 ] := deltaStarAB[n, 0] = deltaClassAB[n]
deltaStarAB[n_, k_] := 
 deltaStarAB[n, k ] = 
  Module[{prod, Y, preY, Delta, rest , uxx  }, 
   Y = Expand[Product[ UU1 + weight[n, i], {i, 0, k - 1}]];
   preY = Y /. { UU1 -> uxx };
   Delta = deltaClassAB[n];
   prod = normalizeVarAB[n, UU1, Y*Delta];
   rest = 
    Sum[Coefficient[preY, uxx, i]*deltaStarAB[n, i], {i, 0, k - 1}];
   Expand[prod - rest]]

Clear[deltaStarChern];
deltaStarChern[n_,k_] :=  deltaStarChern[n,k] = 
  Expand[SymmetricReduction[  deltaStarAB[n,k], {\[Alpha], \[Beta]}, CC12 ][[1]]]

(* (Omega^d)_* u^k *)
Clear[omegaStarAB];
omegaStarAB[n_, d_, k_] := 
 Module[ {idxs = Select[ Table[i, {i, 0, n*d} ], Not[Divisible[#, d]] &]}, 
  UU^k*d^(n - k)*  Product[ UU + (n*d - i)*\[Alpha] + i*\[Beta], {i, idxs}]]

Clear[omegaStarChern];
omegaStarChern[n_,d_,k_] := omegaStarChern[n,d,k] = 
  Expand[SymmetricReduction[
     omegaStarAB[n,d,k], {\[Alpha], \[Beta]}, CC12 ][[1]]]

End[]   (* --------  `EquivCohom` ------- *)

(* normal form of cohomology classes *)
Clear[equivCohomNormalize];
equivCohomNormalize[A_, ProjSpace[n_, u_]] := SymP1`Private`EquivCohom`normalizeVarAB[n,u,A]
equivCohomNormalize[A_, {}      ] := A
equivCohomNormalize[A_, {space_}] := equivCohomNormalize[A, space]
equivCohomNormalize[A_, ps_List ] := 
  equivCohomNormalize[equivCohomNormalize[A, First[ps]], Rest[ps]]


(* pushforward of A in P^n along Omega^d *)
Clear[equivCohomOmegaPF];
equivCohomOmegaPF[A_ , ProjSpace[n_, u_], d_] := equivCohomOmegaPF[A, {u, n}, d, u]
equivCohomOmegaPF[A_ , ProjSpace[n_, u_], 0, v_] := 1;
equivCohomOmegaPF[A_ , ProjSpace[n_, u_], 1, v_] := A /. {u -> v};
equivCohomOmegaPF[A0_, ProjSpace[n_, u_], d_, v_] := 
 Module[{A = Expand[A0]}, 
   Expand[ Sum[Coefficient[A, u, k] * ( SymP1`Private`EquivCohom`omegaStarAB[n,d,k] /. {SymP1`Private`EquivCohom`UU->v} ), {k, 0, n}] ] ] 

(* pushforward of A in P^n1 x P^n2 x ...along Psi_{n1,n2,...}  *)
Clear[equivCohomPsiPF];
equivCohomPsiPF[A_ , {}, w_] := 1;
equivCohomPsiPF[A_ , {ProjSpace[n_, u_]}, w_] := A /. {u -> w};
equivCohomPsiPF[A0_, {ProjSpace[n_, u_], ProjSpace[m_, v_]}, w_] := 
 Module[{A = Expand[A0]}, 
   Expand[ Sum[Coefficient[Coefficient[A, u, i], v, j] * ( SymP1`Private`EquivCohom`psiStarAB[n,m,i,j] /. {SymP1`Private`EquivCohom`UU->w} ) , {i, 0, n}, {j, 0, m}]] ]
equivCohomPsiPF[A_, spaces_List, w_] := Module[
  {p1, p2, ps, y, m, tmp},
  p1 = First[spaces];
  ps = Rest[spaces];
  m = Sum[Dimension[p], {p, ps}];
  p2 = ProjSpace[m, y];
  tmp = equivCohomPsiPF[A, Rest[spaces], y];
  Collect[Expand[equivCohomPsiPF[tmp, {p1, p2}, w]], w]
  ]

(* pushforward of A in P^n along Delta^d *)
Clear[equivCohomDeltaPF];
equivCohomDeltaPF[A_, space_, d_Integer] := equivCohomDeltaPF[A, space, Table[Subscript[Variable[space], i], {i, 1, d}]]
equivCohomDeltaPF[A_, space_, {}] := 1;
equivCohomDeltaPF[A_, ProjSpace[n_, u_], {v_}] := A /. {u -> v};
equivCohomDeltaPF[A_, ProjSpace[n_, u_], {v1_, v2_}] :=
  Sum[ Coefficient[A, u, k] * ( SymP1`Private`EquivCohom`deltaStarAB[n,k] /. {SymP1`Private`EquivCohom`UU1->v1 , SymP1`Private`EquivCohom`UU2->v2} ) , {k, 0, n} ]
equivCohomDeltaPF[A_, ProjSpace[n_, u_], vs_List] := Module[
  {v1, vs2, y, m, tmp},
  v1 = First[vs];
  vs2 = Rest[vs];
  tmp =  equivCohomDeltaPF[A,   ProjSpace[n, u], {v1, y}];
  Expand[equivCohomDeltaPF[tmp, ProjSpace[n, y], vs2]]
  ]


(* ---------------------------------------- *)

End[]

EndPackage[]
