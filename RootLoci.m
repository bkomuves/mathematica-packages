
(* computing classes of coincident root loci

example usage:

<<SymP1`
<<RootLoci`

classOfPn[ EquivMotivicChern, ProjSpace[5,L] ]

*)

Clear["RootLoci`*"];
BeginPackage[ "RootLoci`"]
Needs["Combinatorica`"]
Needs["Useful`"]
Needs["Partitions`"]
Needs["SymP1`"]
Needs["Witt`"]

(* ------------- supported characteristic classes ------------- *)

(* characteristic classes *)

FundClass              = Symbol["FundClass"        ]; 
CSM                    = Symbol["CSM"              ]; 
ToddClass              = Symbol["ToddClass"        ]; 
HirzebruchClass        = Symbol["Hirzebruch"       ]; 
UnnormHirzebruch       = Symbol["UnnormHirzebruch" ]; 
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

(* ... and genera ... *)

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

FundClass              ::usage = "The Fundamental class in cohomology"
CSM                    ::usage = "The Chern-Schwartz-MacPherson class"
ToddClass              ::usage = "The (motivic) Todd class"           
HirzebruchClass        ::usage = "The (motivic) Hirzebruch \*SubscriptBox[T,y] class" 
UnnormHirzebruch       ::usage = "The unnormalized Hirzebruch class"  
MotivicChern           ::usage = "The motivic Chern class K-theory"    
EquivFundClass         ::usage = "The equivariant fundamental class in equivariant cohomology"  
EquivCSM               ::usage = "The equivariant CSM class in equivariant cohomology"        
EquivMotivicChern      ::usage = "The equivariant motivic Chern class in equivariant K-theory"  
AbstractClass          ::usage = "The 'abstract' class (in the relative Grothendieck group of varieties)"       

EulerChar      ::usage = "The Euler characterestic"     
ToddGenus      ::usage = "The Todd genus"     
ChiY           ::usage = "The Hirzebruch \*SubscriptBox[\[Chi],y] genus"          
PoincarePoly   ::usage = "The virtual Poincare polynomial"  
HodgeDeligneE  ::usage = "The Hodge-Deligne E-polynomial" 
NumberOfPoints ::usage = "The number of points over a finite field \*SubscriptBox[\[DoubleStruckCapitalF],q]"
HasseZeta      ::usage = "The Hasse-Weil zeta function of a variety over \*SubscriptBox[\[DoubleStruckCapitalF],q]"     

Theory    ::usage = "Theory[class] gives back the theory the given class lives in. For example Theory[MotivicChern] == KTheory"
shortName ::usage = "shortName[class] gives back a short string, which can be used as name for example when generating tables"

umbralBasis    ::usage = "umbralBasis[theory,ProjSpace[n,u],k] returns the k-th umbral basis element in \*SuperscriptBox[\[DoubleStruckCapitalP],n]"
fromUmbralBasis::usage = "fromUmbralBasis[theory,A,z,ProjSpace[n,u]] converts from the umbral basis to the usual representation"
toUmbralBasis  ::usage = "toUmbralBasis[theory,A,ProjSpace[n,u],z] converts from the usual representation to the umbral basis" 

(* ---------- various classes --------- *)

classOfPn         ::usage = "For example classOfPn[CSM,ProjSpace[n,u]] returns the CSM class of Sym^nP1 = P^n"
abstractGhostClass::usage = "abstractGhostClass[m,S] returns the m-th abstract ghost class in terms of S_n = Sym^n"
genericGhostClass ::usage = "genericGhostClass[class,space] returns the m-th ghost class computed from the classes of P^n"
ghostClass        ::usage = "ghostClass[class,space] returns the m-th ghost class, computed via hardcoded formulas"
umbralGhostClass  ::usage = "umbralGhostClass[class,m,z] returns the m-th ghost class in the 'umbral z-coordinates'"

LSeries::usage = "LSeries[n,d] is the \[ScriptCapitalL](\*SuperscriptBox[x,d];\*SuperscriptBox[h,d]) power series up to degree n (d can be omitted, and defaults to 1)"
PSeries::usage = "PSeries[n,d] is the \[ScriptCapitalP](\*SuperscriptBox[x,d];\*SuperscriptBox[h,d]) power series up to degree n (d can be omitted, and defaults to 1))"
LPoly  ::usage = "LPoly[k] is the degree k term of \[ScriptCapitalL](x;h); LPoly[k,d] is the after the substitution by x -> \*SuperscriptBox[x,d] (d defaults to 1)"
PPoly  ::usage = "PPoly[k] is the degree k term of \[ScriptCapitalP](x;h); PPoly[k,d] is the after the substitution by x -> \*SuperscriptBox[x,d] (d defaults to 1)"

recClassOfRootLoci     ::usage = "recClassOfRootLoci[class, lambda, gen] returns the class of \*SubscriptBox[X,\[Lambda]]"
recClassOfDistinctLoci ::usage = "recClassOfDistinctLoci[class, dvec, gen] returns the class of D(n1,n2...nk)"
ExportRootLoci         ::usage = "ExportRootLoci[class, fname, n (,var)] writes a text file with a table of the classes of root loci up to n"

PartitionClosure   ::usage  = "PartitionClosure[p] returns the set of coarsenings of a partition"
PartitionClosure$v1::usage  = "PartitionClosure$v1[p] is an alternative (slower?) implementation"

(* ========================================================================== *)

Begin["`Private`"]

(* reserve y as the "y" parameter in the Chi_y genus, and Hirzebruch and Motivic Chern classes *)
(* and q as the order of a finite field F_q *)
(* and \[DoubleStruckCapitalL] as the Lefschetz motive *)
(* and x as the variable of the Poincare polynomial, and the L/P series *)
(* and u,v as the variables of the E-polynomial *)

ClearAll[ Global`x , Global`y , Global`q , Global`\[DoubleStruckCapitalL] ];

\[Alpha] = Global`\[Alpha];
\[Beta]  = Global`\[Beta];
X = Global`X;
Y = Global`Y;

SetAttributes[ Global`y                 , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Global`q                 , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Global`\[ScriptCapitalL] , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Global`x                 , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ Global`u                 , {Protected, ReadProtected        } ]; 
SetAttributes[ Global`v                 , {Protected, ReadProtected        } ]; 

x = Global`x;
y = Global`y;
q = Global`q;
\[DoubleStruckCapitalL] = Global`\[DoubleStruckCapitalL];

(* --- internal free variables used in tables --- *)

SetAttributes[ RootLoci`Private`L , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ RootLoci`Private`H , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ RootLoci`Private`U , {Protected, ReadProtected, Locked} ]; 
SetAttributes[ RootLoci`Private`N , {Protected, ReadProtected, Locked} ]; 

LL = RootLoci`Private`L
HH = RootLoci`Private`H
UU = RootLoci`Private`U
NN = RootLoci`Private`N

(* --------------------------------------------------------- *)

Clear[Theory]

Theory[ FundClass         ]   = Cohomology;
Theory[ CSM               ]   = Cohomology;
Theory[ ToddClass         ]   = Cohomology;
Theory[ HirzebruchClass   ]   = Cohomology;
Theory[ UnnormHirzebruch  ]   = Cohomology;
Theory[ MotivicChern      ]   = KTheory;
Theory[ EquivFundClass    ]   = EquivCohomology;
Theory[ EquivCSM          ]   = EquivCohomology;
Theory[ EquivMotivicChern ]   = EquivKTheory;
Theory[ AbstractClass     ]   = AbstractTheory;

Theory[ EulerChar      ]   = TrivialTheory;
Theory[ ToddGenus      ]   = TrivialTheory;
Theory[ ChiY           ]   = TrivialTheory;
Theory[ PoincarePoly   ]   = TrivialTheory;
Theory[ HodgeDeligneE  ]   = TrivialTheory;
Theory[ NumberOfPoints ]   = TrivialTheory;
Theory[ HasseZeta      ]   = TrivialTheory;

(* this can be used for example when exporting tables of classes *)
Clear[shortName]

shortName[FundClass]          = "Fund";
shortName[CSM]                = "CSM";
shortName[ToddClass]          = "Todd";
shortName[HirzebruchClass]    = "Hirz";
shortName[UnnormHirzebruch]   = "UHirz";
shortName[MotivicChern]       = "MC";
shortName[EquivFundClass]     = "eFund";
shortName[EquivCSM]           = "eCSM";
shortName[EquivMotivicChern]  = "eMC";
shortName[AbstractClass]      = "Mot";

shortName[EulerChar     ] = "Chi" 
shortName[ToddGenus     ] = "ToddG"
shortName[ChiY          ] = "ChiY"
shortName[PoincarePoly  ] = "Poin"
shortName[HodgeDeligneE ] = "EPoly"
shortName[NumberOfPoints] = "Cnt"
shortName[HasseZeta     ] = "Zeta"

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
    Normal[Series[ ( u*(1+y*Exp[-u]) / (1-Exp[-u]) )^(n+1) / (1+y) , {u,0,n} ]]

classOfPn[ MotivicChern , ProjSpace[n_,L_] ] := 
  classOfPn[ MotivicChern , ProjSpace[n,L] ] = 
    Expand[Factor[ ( (1+y*L)^(n+1) - (-y)^(n+1)*(1-L)^(n+1) ) / (1 + y) ]]

classOfPn[ EquivMotivicChern , ProjSpace[n_,L_] ] :=
  classOfPn[ EquivMotivicChern , ProjSpace[n,L] ] =
    Expand[ Factor [ 
        ( Product[1+y*L*kWeight[n,i],{i,0,n}] - 
            (-y)^(n+1)*Product[1-L*kWeight[n,i],{i,0,n}] ) / (1 + y) ] ]

classOfPn[ EulerChar      , ProjSpace[n_ , _] ] := n+1
classOfPn[ ToddGenus      , ProjSpace[n_ , _] ] := 1
classOfPn[ ChiY           , ProjSpace[n_ , _] ] := Sum[ (-Global`y)^i           , {i,0,n} ]
classOfPn[ PoincarePoly   , ProjSpace[n_ , _] ] := Sum[   Global`x^(2*i)        , {i,0,n} ]
classOfPn[ HodgeDeligneE  , ProjSpace[n_ , _] ] := Sum[  (Global`u*Global`v)^i  , {i,0,n} ]
classOfPn[ NumberOfPoints , ProjSpace[n_ , _] ] := Sum[   Global`q^i            , {i,0,n} ]
classOfPn[ HasseZeta      , ProjSpace[n_ , _] ] := Sum[ Teichmuller[ \[DoubleStruckCapitalL]^i ] , {i,0,n} ]

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
Clear[equivKTheoryUmbralBasisL];
equivKTheoryUmbralBasisL [ 0  , 0  ]  = 1;
equivKTheoryUmbralBasisL [ n_ , k_ ] := equivKTheoryUmbralBasisL [ n , k ] = Module[ {L2},
  If[ k == 0 , 
    PsiPF[ EquivKTheory , equivKTheoryUmbralBasisL [ n-1 , 0   ] * (1-L2) , { ProjSpace[n-1,LL] , ProjSpace[1,L2] } , LL ] ,
    PsiPF[ EquivKTheory , equivKTheoryUmbralBasisL [ n-1 , k-1 ] *   (L2) , { ProjSpace[n-1,LL] , ProjSpace[1,L2] } , LL ]
    ]]

(* to test ktheoryUmbralBasisH against *)
Clear[computeKTheoryUmbralBasisL]
computeKTheoryUmbralBasisL [ 0  , 0  ]  = 1;
computeKTheoryUmbralBasisL [ n_ , k_ ] := computeKTheoryUmbralBasisL [ n , k ] = Module[ {L2},
  If[ k == 0 , 
    PsiPF[ KTheory , computeKTheoryUmbralBasisL [ n-1 , 0   ] * (1-L2) , { ProjSpace[n-1,LL] , ProjSpace[1,L2] } , LL ] ,
    PsiPF[ KTheory , computeKTheoryUmbralBasisL [ n-1 , k-1 ] *   (L2) , { ProjSpace[n-1,LL] , ProjSpace[1,L2] } , LL ]
    ]]

(* convert from the umbral basis *)
fromUmbralBasis[theory_, A0_, z_, ProjSpace[n_, u_]] := 
 Module[{A = Expand[A0]},
  Sum[Coefficient[A, z, k]*
    umbralBasis[theory, ProjSpace[n, u], k], {k, 0, n}]]

(* convert to the umbral basis *)
toUmbralBasis[theory_, A0_, ProjSpace[n_, u_], z_] := Module[
  {A = Expand[A0],
   unkn, lhs, rhs, a, vars1, vars2, sols, sol, theoryvars},
  unkn = Sum[Subscript[a, i]*z^i, {i, 0, n}];
  lhs = Sum[
    Subscript[a, i]*umbralBasis[theory, ProjSpace[n, u], i], {i, 0, 
     n}];
  vars1 = Join[{u}, Select[Variables[A0], # != z &]];
  vars2 = Table[Subscript[a, i], {i, 0, n}];
  sols = mySolveAlways[lhs, A, vars1, vars2];
  sol = sols[[1]];
  (*
  Print["============"];
  Print[ProjSpace[n,u]];
  Print[vars1];
  Print[vars2];
  Print[lhs];
  Print[A0];
  Print[sols];
  *)
  unkn /. sol
  ]

(* --------------- ghost classes (or "Adams operations" on [P^1]) -------------- *)

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

ghostClass      [ HirzebruchClass , ProjSpace[m_, u_] ] := (1+(-y)^m)*u^m + u^(m-1)*(1-(-y)^m)/(1+y)
umbralGhostClass[ HirzebruchClass , m_ , z_ ] := (1+(-y)^m) + z*(1-(-y)^m)/(1+y)

ghostClass      [ UnnormHirzebruch , ProjSpace[m_, u_] ] := (1+(-y)^m)*u^m + u^(m-1)*(1-(-y)^m)
umbralGhostClass[ UnnormHirzebruch , m_ , z_ ] := (1+(-y)^m) + z*(1-(-y)^m)

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
  
(* ---------------- The L and P power series ------------- *)

xx[i_] := Subscript[Global`x, i]
xxs[n_] := Table[xx[i], {i, 1, n}]

Clear[LSeries, PSeries, LSeries$check, LPoly, PPoly];

Clear[Global`h]
h = Global`h

LSeries[n_] := LSeries[n] = 
  Normal[Series[Log[1 + Sum[xx[i]*Global`h^i, {i, 1, n}]], {h, 0, n}]]
PSeries[n_] := PSeries[n] =
  Normal[Series[
    Sum[MoebiusMu[d]/d*LSeries[n, d], {d, 1, n}], {h, 0, n}]]
LSeries$check[n_] := LSeries$check[n] =
  Normal[Series[Sum[1/d*PSeries[n, d], {d, 1, n}], {h, 0, n}]]

LSeries[n_, d_] := Module[{m = Quotient[n, d]},
  LSeries[m] /. {h -> h^d} /. Table[xx[i] -> xx[i]^d, {i, 1, m}]]
PSeries[n_, d_] := Module[{m = Quotient[n, d]},
  PSeries[m] /. {h -> h^d} /. Table[xx[i] -> xx[i]^d, {i, 1, m}]]

LPoly[k_] := Coefficient[LSeries[k], h, k]
PPoly[k_] := Coefficient[PSeries[k], h, k]
LPoly[k_, d_] := LPoly[k] /. Table[xx[i] -> xx[i]^d, {i, 1, k}]
PPoly[k_, d_] := PPoly[k] /. Table[xx[i] -> xx[i]^d, {i, 1, k}]

(*
lhs = Sum[PSeries[10, m]*Subscript[A, m]/m, {m, 1, 10}];
rhs = Sum[Subscript[B, k]/k*LSeries[10, k], {k, 1, 10}];
solB = mySolveAlways$v2[lhs, rhs, Prepend[xxs[10], h], Table[Subscript[B, i], {i, 1, 10}]]
solA = mySolveAlways$v2[lhs, rhs, Prepend[xxs[10], h], Table[Subscript[A, i], {i, 1, 10}]]
*)

(* ----------- computing classes of root loci via the recursive algo ----------- *)


(* g is the cohomology theory generator; z is a power series variable *)
ClearAll[g, z, G, Z]
SetAttributes[g, {Protected, ReadProtected}];
SetAttributes[z, {Protected, ReadProtected}];

Z = z;
G = g;
GG[i_] := Subscript[g, i];
GGs[n_] := Table[g[i], {i, 1, n}]

(* class of the root loci X_lambda, computed using the recursive algorithm *)
recClassOfRootLoci[class_, lambda_, var_] := classOfXLam[class, lambda] /. {G -> var}

(* class of the distinct loci D(n1,n2...,nr), computed using the recursive algorithm *)
recClassOfDistinctLoci[class_, n_Integer, var_] := classOfDisj1[class, n] /. {G -> var}
recClassOfDistinctLoci[class_, ns_List  , var_] := classOfDisj[class, ns] /. 
  Table[GG[i] -> Subscript[var, i], {i, 1, Length[ns]}]

posVectorQ[as_] := Map[# >= 0 &, as] /. {List -> And};

kdeTriples[p_, ns_] := 
 Module[{m = Length[ns], posQ, oneK, A}, 
  oneK[k_] := 
   Table[{k, ns - es, es}, {es, Combinatorica`Compositions[p - k, m]}];
  A = Table[oneK[k], {k, 0, p - 1}];
  A = Select[Flatten[A, 1], posVectorQ[Snd[#]] &];
  A]

Clear[classOfXLam, classOfDisj1, classOfDisj, classOfDisjSorted]

(* class of D(n) = X(1^n) *)
classOfDisj1[class_, 0 ] := classOfDisj1[class, 0] = 1
classOfDisj1[class_, 1 ] := classOfDisj1[class, 1] = classOfPn[class, ProjSpace[1, G]]
classOfDisj1[class_, n_] := classOfDisj1[class, n] = 
  Module[{parts = 
     Select[Combinatorica`Partitions[n], Length[#] < n &]}, 
   Expand[classOfPn[class, ProjSpace[n, G]] - 
     Sum[classOfXLam[class, p], {p, parts}]]]

(* class of the coincident root loci X(lambda) *)
classOfXLam[class_, {} ] := classOfXLam[class, {}] = 1
classOfXLam[class_, {1}] := classOfXLam[class, {1}] = classOfPn[class, ProjSpace[1, G]]
classOfXLam[class_, lambda_] := classOfXLam[class, lambda] = 
  Module[{es, theory, m, m1, ns, pairs, spaces, A},
   theory = Theory[class];
   es = toExpoVec[lambda];
   m = Length[es];
   ns = Range[m];
   pairs = Zip[ns, es];(*i^e_*)
   
   pairs = Select[pairs, Snd[#] > 0 &];(*!!!*)
   m1 = Length[pairs];
   ns = Map[Fst, pairs];
   es = Map[Snd, pairs];
   m = Length[ns];
   A = classOfDisj[class, es];
   spaces = Table[ ProjSpace[es[[i]], GG[i] ], {i, 1, m1}];
   (* Print[{lambda,es,ns,spaces}]; *)
   
   Expand[OmegaVecPF[theory, A, spaces, ns, G]]
   ]

(* class of D(n1,n2,...) *)
classOfDisj[class_, {}  ] := classOfDisj[class, {} ] = 1
classOfDisj[class_, {n_}] := classOfDisj[class, {n}] = classOfDisj1[class, n] /. {G -> GG[1]}
classOfDisj[class_, ns0_] := classOfDisj[class, ns0] = 
  Module[{m = Length[ns0], nis0, nis1, ns1, idxs, X, ttt}, 
   nis0 = Zip[ns0, Range[m]];
   nis1 = SortBy[nis0, -Fst[#] &];
   idxs = Map[Snd, nis1];
   ns1 = Select[Map[Fst, nis1], # > 0 &];
   X = classOfDisjSorted[class, ns1];
   X = X /. Table[GG[i] -> Subscript[ttt, i], {i, 1, m}];
   X /. Table[Subscript[ttt, i] -> GG[idxs[[i]]], {i, 1, m}]
   ]

(* a single term corresponding to a triple (k,ds,es) *)
Clear[singleKDE]
singleKDE[class_, {k_, ds_, es_}] := singleKDE[class, {k, ds, es}] = 
  Module[{theory, A, B, m, i, vars, dims, p, q, r, s, pp, qq, rr, ss, pps, qqs, rrs, sss, Zz, spaces},
   theory = Theory[class];
   m = Length[ds];
   pp[i_] := Subscript[p, i];
   qq[i_] := Subscript[q, i];
   rr[i_] := Subscript[r, i];
   ss[i_] := Subscript[s, i];
   pps = Table[pp[i], {i, 1, m}];
   qqs = Table[qq[i], {i, 1, m}];
   rrs = Table[rr[i], {i, 1, m}];
   sss = Table[ss[i], {i, 1, m}];
   vars = Join[{Zz}, pps, qqs];
   dims = Join[{k}, ds, es];
   A = classOfDisj[class, dims] /. 
     Table[GG[i] -> vars[[i]], {i, 1, 2 m + 1}];
   B = A;
   For[i = 1, i <= m, i++, 
    B = Expand[
       DeltaPF[theory, B, 
        ProjSpace[es[[i]], qq[i]], {rr[i], ss[i]}]];
    ];
   For[i = 1, i <= m, i++,
    B = Expand[PsiPF[theory, B,
        {ProjSpace[ds[[i]], pp[i]], ProjSpace[es[[i]], ss[i]]},
        GG[i]]];
    ];
   dims = Join[{k}, es];
   vars = Join[{Zz}, rrs];
   spaces = ZipWith[ProjSpace, dims, vars];
   B = PsiPF[theory, B, spaces, Z]; (* Z ???? *)
   Expand[B]
   ]

(* equiv mc of D(d1,d2,...), but we require d1>=d2>=d3>=...>=dn>0 *)
classOfDisjSorted[class_, {}  ] := classOfDisjSorted[class, {}] = 1
classOfDisjSorted[class_, {n_}] := classOfDisjSorted[{class, n}] = classOfDisj1[class, n] /. {G -> GG[1]}
classOfDisjSorted[class_, pns_] := classOfDisjSorted[class, pns ] = 
  Module[{p = pns[[1]], ns = Drop[pns, 1], A, B, rest, KDE},
   KDE = kdeTriples[p, ns];
   A = (classOfDisj1[class, p] /. {G -> Z})*classOfDisj[class, ns];
   rest = Sum[singleKDE[class, kde], {kde, KDE}];
   B = Expand[A - rest];
   B = B /. Table[GG[i] -> GG[i + 1], {i, 1, Length[ns]}];
   B = B /. {Z -> GG[1]}
   ]


(* export the classes of X(lambda) for |lambda|<=n *)
ExportRootLoci[class_, fname_, n_] := ExportRootLoci[class, fname, n, standardVar[Theory[class]]]
ExportRootLoci[class_, fname_, n_, var_] := 
 Module[{h, i, p, parts, k, m, s, j, A, name, theory},
  theory = Theory[class];
  name = shortName[class];
  h = OpenWrite[fname];
  WriteString[h, 
   "\n(* --- " <> ToString[class] <> 
    " classes of coincident root loci in " <> ToString[theory] <> 
    " up to n=" <> ToString[n] <> " --- *)"];
  For[i = 1, i <= n, i++, Print["\nn = ", i];
    WriteString[h, "\n(* =================== *)"];
    WriteString[h, 
     "\n(* ----   n = " <> ToString[i] <> "   ---- *)\n\n"];
    parts = Combinatorica`Partitions[i];
    m = Length[parts];
    For[j = 1, j <= m, j++, p = parts[[j]];
     Print["part = ", p];
     A = Expand[recClassOfRootLoci[class, p, var]];
     s = ToString[A, FormatType -> InputForm, PageWidth -> Infinity, 
       TotalWidth -> Infinity];
     s = StringJoin[name, "[", ToString[p], "] = ", s, " ;\n\n"];
     WriteString[h, s];]] 
   Close[h]
  ;]

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
