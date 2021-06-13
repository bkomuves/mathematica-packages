
(* computing classes of coincident root loci *)

(* example usage:

<<SymP1`
<<RootLoci`

classOfPn[ EquivMotivicChern, ProjSpace[5,L] ]

*)

Clear["RootLoci`*"];
BeginPackage[ "RootLoci`"]
Needs["Combinatorica`"]
Needs["Useful`"]
Needs["SymP1`"]

FundClass              ::usage = "The Fundamental class in cohomology"
CSM                    ::usage = "The Chern-Schwartz-MacPherson class"
ToddClass              ::usage = "The (motivic) Todd class"           
Hirzebruch             ::usage = "The (motivic) Hirzebruch T_y class" 
UnnormalizedHirzebruch ::usage = "The unnormalized Hirzebruch class"  
MotivicChern           ::usage = "The motivic Chern class K-theory"    
EquivFundClass         ::usage = "The equivariant fundamental class in equivariant cohomology"  
EquivCSM               ::usage = "The equivariant CSM class in equivariant cohomology"        
EquivMotivicChern      ::usage = "The equivariant motivic Chern class in equivariant K-theory"  
AbstractClass          ::usage = "The abstract class"       

Theory ::usage = "Theory[class] gives back the theory the given class lives in. For example Theory[MotivicChern] == KTheory"

EulerChar      ::usage = "The Euler characterestic"     
ToddGenus      ::usage = "The Todd genus"     
ChiY           ::usage = "The Hirzebruch Chi_y genus"          
PoincarePoly   ::usage = "The virtual Poincare polynomial"  
HodgeDeligneE  ::usage = "The Hodge-Deligne E-polynomial" 
NumberOfPoints ::usage = "The number of points over a finite field F_q"
HasseZeta      ::usage = "The Hasse-Weil zeta function of a variety over F_q"     

classOfPn         ::usage = "For example classOfPn[CSM,ProjSpace[n,u]] returns the CSM class of Sym^nP1 = P^n"
abstractGhostClass::usage = "abstractGhostClass[m,S] returns the m-th abstract ghost class in terms of S_n = Sym^n"
genericGhostClass ::usage = "genericGhostClass[class,space] returns the m-th ghost class computed from the classes of P^n"
ghostClass        ::usage = "ghostClass[class,space] returns the m-th ghost class, computed via hardcoded formulas"
umbralGhostClass  ::usage = "umbralGhostClass[class,m,z] returns the m-th ghost clsas in the `umbral z-coordinates'"

PartitionClosure   ::usage  = "PartitionClosure[p] returns the set of coarsenings of a partition"
PartitionClosure$v1::usage  = "PartitionClosure$v1[p] is an alternative (slower?) implementation"

(* ========================================================================== *)

Begin["`Private`"]

(* reserve y as the "y" parameter in Hirzebruch and Motivic Chern classes *)
(* and q as the order of a finite field F_q *)

ClearAll[ Global`y , Global`q ];

SetAttributes[ Global`y , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Global`q , {Protected, ReadProtected, Locked} ]; 

y = Global`y;
q = Global`q;

(* --- internal free variables used in tables --- *)

SetAttributes[ RootLoci`Private`L , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ RootLoci`Private`H , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ RootLoci`Private`U , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ RootLoci`Private`N , {Protected, ReadProtected, Locked} ]; 

LL = RootLoci`Private`L
HH = RootLoci`Private`H
UU = RootLoci`Private`U
NN = RootLoci`Private`N

(* ------------- supported characteristic classes ------------- *)

(* characteristic classes *)

FundClass              = Symbol["FundClass"        ]; 
CSM                    = Symbol["CSM"              ]; 
ToddClass              = Symbol["ToddClass"        ]; 
Hirzebruch             = Symbol["Hirzebruch"       ]; 
UnnormalizedHirzebruch = Symbol["UnnormHirzebruch" ]; 
MotivicChern           = Symbol["MotivicChern"     ]; 
EquivFundClass         = Symbol["EquivFundClass"   ]; 
EquivCSM               = Symbol["EquivCSM"         ]; 
EquivMotivicChern      = Symbol["EquivMC"          ]; 

SetAttributes[ FundClass         , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ CSM               , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ ToddClass         , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ HirzebruchClass   , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ UnnormHirzebruch  , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ MotivicChern      , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ EquivFundClass    , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ EquivCSM          , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ EquivMotivicChern , {Protected, ReadProtected, Locked} ]; 

Theory[FundClass]          = Cohomology;
Theory[CSM]                = Cohomology;
Theory[ToddClass]          = Cohomology;
Theory[HirzebruchClass]    = Cohomology;
Theory[UnnormHirzebruch]   = Cohomology;
Theory[MotivicChern]       = KTheory;
Theory[EquivFundClass]     = EquivCohomology;
Theory[EquivCSM]           = EquivCohomology;
Theory[EquivMotivicChern]  = EquivKTheory;
Theory[AbstractClass]      = AbstractTheory;

(* genera *)

EulerChar      = Symbol["EulerChar"       ];
ToddGenus      = Symbol["ToddGenus"       ];
ChiY           = Symbol["ChiY"            ];
PoincarePoly   = Symbol["PoincarePoly"    ];
HodgeDeligneE  = Symbol["HodgeDeligneE"   ];
NumberOfPoints = Symbol["NumberOfPoints"  ];
HasseZeta      = Symbol["HasseZeta"       ];

SetAttributes[ EulerChar      , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ ToddGenus      , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ ChiY           , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ PoincarePoly   , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ HodgeDeligneE  , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ NumberOfPoints , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ HasseZeta      , {Protected, ReadProtected, Locked} ]; 

(* --------------------------------------------------------- *)

weight[n_,i_]  := (n-i)*Global`\[Alpha] + i*Global`\[Beta]
kWeight[n_,i_] := Global`X^(n-i) * Global`Y^i

(* ------------ classes of Sym^nP1 = P^n ------------------- *)

Clear[classOfPn]

classOfPn[ FundClass      , ProjSpace[n_,u_] ] := 1
classOfPn[ EquivFundClass , ProjSpace[n_,u_] ] := 1

classOfPn[ CSM      , ProjSpace[n_,u_] ] := Expand[ (1+u)^(n+1) - u^(n+1) ]
classOfPn[ EquivCSM , ProjSpace[n_,u_] ] :=
  classOfPn[ EquivCSM , ProjSpace[n,u] ] = 
    Expand[ Product[1+u+weight[n,i],{i,0,n}] - Product[u+weight[n,i],{i,0,n}] ]

classOfPn[ ToddClass , ProjSpace[n_,u_] ]  :=
  classOfPn[ ToddClass , ProjSpace[n,u] ] =
    Normal[Series[ ( u / (1-Exp[-u]) )^(n+1) , {u,0,n} ]]

classOfPn[ HirzebruchClass , ProjSpace[n_,u_] ] :=
  classOfPn[ HirzebruchClass , ProjSpace[n,u] ] =
    Normal[Series[ ( u*(1+y) / (1-Exp[-u*(1+y)]) - u*y )^(n+1) , {u,0,n} ]]

classOfPn[ UnnormHirzebruch , ProjSpace[n_,u_] ] := 
  classOfPn[ UnnormHirzebruch , ProjSpace[n,u] ] =
    Normal[Series[ u*(1+y*Exp[-u]) / (1-Exp[-u]) , {u,0,n} ]]

classOfPn[ MotivicChern , ProjSpace[n_,L_] ] := 
  classOfPn[ MotivicChern , ProjSpace[n,L] ] = 
    Expand[Factor[ ( (1+y*L)^(n+1) - (-y)^(n+1)*(1-L)^(n+1) ) / (1 + y) ]]

classOfPn[ EquivMotivicChern , ProjSpace[n_,L_] ] :=
  classOfPn[ EquivMotivicChern , ProjSpace[n,L] ] =
    Expand[ Factor [ 
        ( Product[1+y*L*kWeight[n,i],{i,0,n}] - 
            (-y)^(n+1)*Product[1-L*kWeight[n,i],{i,0,n}] ) / (1 + y) ] ]

(* ----------------------- umbral bases ---------------------- *)

Clear[ umbralBasis ];

umbralBasis[ Cohomology      , ProjSpace[n_,u_] , k_ ] := k!*u^(n-k)
umbralBasis[ KTheory         , ProjSpace[n_,L_] , k_ ] := umbralBasis[ KTheory , ProjSpace[n,L] , k ] = Expand[ ktheoryUmbralBasisH [ n,  k , (1-L) ] ]
umbralBasis[ EquivCohomology , ProjSpace[n_,u_] , k_ ] := equivCohomUmbralBasis    [ n , k , u ]
umbralBasis[ EquivKTheory    , ProjSpace[n_,L_] , k_ ] := equivKTheoryUmbralBasisL [ n , k ] /. { LL -> L }

ktheoryUmbralBasisH   [ n_, k_ , H_ ] := Sum[ (-1)^(k-i) * i! * StirlingS2[k+1,i+1] * H^(n-i) , {i,0,k} ]
equivCohomUmbralBasis [ n_, k_ , u_ ] := PhatVar[ n , n - k , u ] * k!

Clear[Phat]
Phat[-1] = 0;
Phat[0]  = 1;
Phat[k1_] := Phat[k1] =  Expand[If[k1 < 0, 0, 
    Module[{k = k1 - 1}, (UU + k (\[Alpha] + \[Beta])) Phat[k] +  k (NN - k + 1) (\[Alpha]*\[Beta])*Phat[k - 1]]]]

PhatVar[n_, k_, u_] := Phat[k] /. { NN -> n , UU -> u }

(* The Pontrjagin product of k copies of L and (n-k) copies of (1-L) *)
equivKTheoryUmbralBasisL [ n_ , k_ ] := equivKTheoryUmbralBasisL [ n , k ] = Module[ {L2},
  If[ k == 0 , 
    PsiPF[ EquivKTheory , equivKTheoryUmbralBasisL [ n-1 , 0   ] * (1-L2) , { ProjSpace[n-1,LL] , ProjSpace[1,L2] } , LL ]
    PsiPF[ EquivKTheory , equivKTheoryUmbralBasisL [ n-1 , k-1 ] *   (L2) , { ProjSpace[n-1,LL] , ProjSpace[1,L2] } , LL ]
    ]]

(* to test ktheoryUmbralBasisH against *)
computeKTheoryUmbralBasisL [ n_ , k_ ] := computeKTheoryUmbralBasisL [ n , k ] = Module[ {L2},
  If[ k == 0 , 
    PsiPF[ KTheory , computeKTheoryUmbralBasisL [ n-1 , 0   ] * (1-L2) , { ProjSpace[n-1,LL] , ProjSpace[1,L2] } , LL ]
    PsiPF[ KTheory , computeKTheoryUmbralBasisL [ n-1 , k-1 ] *   (L2) , { ProjSpace[n-1,LL] , ProjSpace[1,L2] } , LL ]
    ]]

(* --------------------- ghost classes ---------------------- *)

(* compute the ghost class as a formal polynomial *)
Clear[ abstractGhostClass ]
abstractGhostClass[ m_ ] := abstractGhostClass[ m , Global`\[ScriptCapitalS] ] 
abstractGhostClass[ m_ , S_ ] := abstractGhostClass[ m ]  =
  Coefficient[ Series[ D[ Log[ 1 + Sum[ h^n*Subscript[S,n] , {n,1,m} ] ] , h ] , {h,0,m} ] , h , m-1 ]

(* returns a list of {coeff,exponents} pairs *)
Clear[ abstractGhostCoeffs ]
abstractGhostCoeffs[ m_ ] := abstractGhostCoeffs[ m ]  =
  Module[{S} , MCoeffs [ abstractGhostClass[m,S] , Table[ Subscript[S,i] , {i,1,m} ] ] ]

(* compute the ghost class from the classes of Pn *)
Clear[genericGhostClass]
genericGhostClass[class_, ProjSpace[m_, var_]] := 
 genericGhostClass[class, space] = 
  Module[{S, syms, coeffs, fun, theory, product},
   theory = Theory[class];
   syms = Table[classOfPn[class, ProjSpace[i, var]], {i, 1, m}];
   coeffs = abstractGhostCoeffs[m];
   product[dims_] := Module[
     {n = Length[dims], A, spaces, uu, xxx},
     uu[i_] := Subscript[xxx, i];
     spaces = Table[ProjSpace[dims[[i]], uu[i]], {i, 1, n}];
     A = Product[syms[[dims[[i]]]] /. {var -> uu[i]}, {i, 1, n}];
     PsiPF[theory, A, spaces, var]
     ];
   fun[{c_, expos_}] := c*product[fromExpoVec[expos]];
   SumList[Map[fun, coeffs]]
   ]
  
(* the "ghost classes" are the coefficients of the logarithmic differential *)
ghostClass      [ CSM , ProjSpace[m_, u_] ] := 2*u^m + m*u^(m-1)
umbralGhostClass[ CSM , m_ , z_ ] := 2 + m*z

ghostClass      [ MotivicChern , ProjSpace[m_, L_] ] := 
  (1 + (-y)^m)*(1-L)^m + (1 - (-y)^m) * L*(1-L)^(m-1)
umbralGhostClass[ MotivicChern , m_ , z_ ] := (1+z) + (-y)^m*(1-z)

ghostClass      [ EquivCSM , ProjSpace[m_, u_] ] := 
  (2*u + m(1+\[Alpha]+\[Beta])) * Product[ u + weight[m,j] {j,1,m-1} ]
umbralGhostClass[ EquivCSM , m_ , z_ ] := 
  (1+z*\[Alpha])^m + (1+z*\[Beta])^m + ( (1+z*\[Alpha])^m - (1+z*\[Beta])^m ) / (\[Alpha]+\[Beta])

ghostClass      [ EquivMotivicChern , ProjSpace[m_, L_] ] := 
  (1 + (-y)^m*(1-L*(X^m+Y^m))) * Product[ 1-L*kWeight[m,j] , {j,1,m-1} ]
umbralGhostClass[ EquivMotivicChern , m_ , z_ ] := 
           ( X^m*(1+z*(1-Y))^m - Y^m*(1+z*(1-X))^m ) / (X^m -Y^m) + 
  (-y)^m * ( X^m*(1+z*(1-X))^m - Y^m*(1+z*(1-Y))^m ) / (X^m -Y^m)
  

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
