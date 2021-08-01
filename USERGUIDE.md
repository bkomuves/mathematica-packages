
Documentation for my Mathematica packages
-----------------------------------------

WARNING: this codebase is very much in flux and probably is full of bugs.
If you find a bug please send me an email to _username_ (at) gmail.

Table of content:

* `Useful` - all kind of basic functions, which are missing from Mathematica
* `Partitions` - functions related to (integer) partitions
* `Schur` - working with symmetric polynomials
* `SymP1` - computations in the cohomology and K-theory of Sym<sup>n</sup>&Popf;<sup>1</sup> = &Popf;<sup>n</sup>
* `RootLoci` - computing characteristic classes of coincident root loci 


Package `Useful`
================

Assorted useful functions.

### Lists and tuples

    SumList[L]
    MaxList[L]

Returns the sum (resp. maximum) of the elements of the list `L`. Note: the 
maximum will return minus infinity for the empty list

    Fst[{x,y}] 
    Snd[{x,y}] 
    
Returns the first and second element of a pair (or list), respectively.

    replicate[A,n] 

Replicates `A`, resulting a list `{A,A,...,A}` of length `n`.

    Zip[as,bs]
    Zip[as,bs,cs]

Create a list of pairs (or triples) from two (or three) lists.

    ZipWith[f,as,bs]
    ZipWith[f,as,bs,cs]

Creates a new list by applying the function `f` to the pairs or triples made from
the lists `as` and `bs` (and `cs`). For example:

    ZipWith[Plus, {100, 200, 300}, {4, 5, 6}] == {104, 205, 306}
    
The functions

    IsEmptyList[L]
    IsNotEmptyList[L]

returns `True` (resp. `False`) if `L` is empty, and `False` (reps. `True`) otherwise.
TODO: maybe rename these so that they conform to the Mathematica convention of 
ending with `Q`?
 
    reverseSort[L]

is equivalent to `Reverse[Sort[L]]`

    unique[L]

Sorts and removes the duplicates from a list (making it usable as a representation of a set).

    extendListWithZeros[L,n] 

Add zeros to the end of the list `L` till it reaches length `n`.

### Coefficient extraction

    mcoeff[A,{x1,x2,x3},{e1,e2,e3}] 

extracts the coefficient of `x1^e1*x2^e2*x3^e3` in `A`.

    MCoeffs[A,{x1,x2,x3}] 

Returns a list of `{coeff,expovec}` pairs, representing all the monomials in the
variables `{x1,x2,x3}`

### Equation solving

    mySolveAlways[lhs,rhs,{p1,p2,..}] 
    mySolveAlways[lhs,rhs,{p1,p2,..},{a1,a2,..}] 

tries to solve `lhs==rhs` where `p1,p2,...` are parameters and `a1,a2,...` are unknowns.
This is similar to the built-in function `SolveAlways`, but works with more than
one parameter.

### Generating functions    
    
    fallingFactorial[x,k] 

Computes the falling factorial `x*(x-1)*(x-2)*...*(x-k+1)`.
    
    ogf2egf[F,x,t]
    egf2ogf[G,t,x] 

Converts between an ordinary generating function F(x) and an exponential generating
function G(t), using Laplace transformation. For example:

    ogf2egf[ 1 / (1-a*x) , x , t ] == Exp[ a*t ]
    egf2ogf[ Exp[a*t]    , t , x ] == 1 / (1-a*x)

Of course this does not always work.

### Permutations

    SortingPermutation[L]
    ReverseSortingPermutation[L]

Given a list `L` of length `n`, these functions returns a permutation of `{1,2,...,n}`, which,
applied to the list `L` result in that list sorted in ascending (resp. descending) order.

Note: this permutation is not unique if the list contains repeated elements.

### Misc

    toExportForm[expr]

Converts an expression to "export form". This is useful when you want to export
the results of calculations into text files.
    

Package `Partitions`
====================

An (integer) partition is encoded as a non-increasing list of positive integers (same
as the built-in Mathematica convention). Sometimes we allow extra zeros at
the end.

Available functions
-------------------

    emptyPartQ[mu]

Returns `True` for the empty partition, `False` otherwise

    samePartitionQ[lam,mu]

Returns True if `lam` and `mu` are the same partitions (after discarding zeros).

    subPartitionQ[lam,mu]

Returns True if `mu` is a sub-partition of `lam` (note the order of arguments!).

    toExpoVec[mu]  
    fromExpoVec[expos]

Converts a partition to an exponent vector, and vica versa. For example

    toExpoVec[{4,4,2,1,1,1}] = {3,1,0,2}

because (4,4,2,1,1,1) = (1<sup>3</sup>, 2<sup>1</sup>, 3<sup>0</sup>, 4<sup>2</sup>).

    toExpoForm[mu]  
    fromExpoForm[pairs]

Converts a partition to a list of (part,exponent) pairs, and vica versa. This is
very similar to the above, but it lacks the zero exponents:

    toExpoForm[{4,4,2,1,1,1}] = {{1,3},{2,1},{4,2}}

Basic size queries:

    partitionWidth[mu] 
    partitionHeight[mu]
    partitionWeight[mu]

These return the width (that is, the number of parts), height (that is, the maximum of parts) and 
weight (that is the sum of parts) of a partition.

    dualPart[mu]

Returns the dual partition of `mu`

    partitionMonom[mu]
    partitionMonom[mu,x]

Returns the associated monomial 
(x<sub>1</sub><sup>e<sub>1</sub></sup> x<sub>2</sub><sup>e<sub>2</sub></sup> x<sub>3</sub><sup>e<sub>3</sub></sup>, ... )
where `mu` = (1<sup>e<sub>1</sub></sup>, 2<sup>e<sub>2</sub></sup>, 3<sup>e<sub>3</sub></sup>, ... )

    partitionIndex[lam,i]

This is the same as `lam[[i]]` but allows overindexing (in which case it returns 0)

    partitionReplace1[lam,i,y] 
    partitionReplaceParts[lam,{{i1,y1},{i2,y2}}] 

The first one replaces the `i`-th element with `y`, the second one does the same
but with multiple replacements. Again these handle overindexing (in which case
the missing parts will be filled by zeros).
    
    blockPartition[{n,k}]
    blockPartition[ n,k ]

The partition (n<sup>k</sup>)=(n,n,...,n) (k times). Note the order of arguments: the first is the
_height_ (maximum) and the second is the _width_ (length).

    fitsIntoBlock[{n,k}, lam] 

Returns `True` if the partition `lam` fits into (is a subpartition of) the block 
partition (n<sup>k</sup>).

    partitionComplement[{n,k},lam]

The the complement of `lam` in the block partition (n<sup>k</sup>). This will throw an error
if `lam` does not fits in the block.
        
    partitions[n] 
    allPartitions[n] 
    partitionsInBlock[{n,m}]
    partitionsInBlock[{n,m},k]

These lists :

* the partitions of weight n
* all partitions of weights at most n
* the partitions fitting into a block (n<sup>m</sup>)
* the partitions of weight k fitting into a block (n<sup>m</sup>)

---

    subPartitionsOf[lam]
    subPartitionsOf[lam,k]

These list the sub-partitions of `lam` (resp. those of weight k)

    
Package `Schur`
===============
    
This package implements conversion between symmetric polynomial in the standard bases

* elementary symmetric polynomials: e<sub>i</sub>
* complete homogeneous symmetric polynomials: h<sub>i</sub>
* power symmetric polynomials: p<sub>i</sub>
* Schur polynomials: s<sub>&lambda;</sub>

You can use other names than e,h,p,s (for example when you work for with 
multiple variable sets); these are just the default ones.

The notation we use is subscripts. In case of elementary, complete and power symmetric
polynomials, multiple subscripts means simply the product of the corresponding
symmetric polynomials.

You can create subscript using `mkSubscript`

    mkSubscript[var, i]
    mkSubscript[var, i1,i2,...]
    mkSubscript[var,{i1,i2,...}]

For the standard names we have convenient shorthands:

    ss[3,2,1] 
    ee[3,2,1] 
    hh[3,2,1] 
    pp[3,2,1] 

For Schur polynomials, there is a function `mkSchur` which works as `mkSubscript`
but normalizes the input to a partition (possible adding a sign or returning 0).

    mkSchur[var,i]
    mkSchur[var,i1,i2,...]
    mkSchur[var,{i1,i2,...}]

This works so that (hopefully...) the following identity is true:

    toX(mkSchur[s,idxs],s,{x,n}) == schurWeylCharacter[idxs,{x,n}]


### Conversion function

TODO: implement conversions to/from power symmetric functions!

In most of the functions below you can omit the variables name, in which case 
they use the default names (eg. `e` for elementary symmetric and `s` for Schur)

    sToE[expr,s,e]
    eToS[expr,e,s] 

Converts Schur polynomials to elementary symmetric polynomials and vica versa.
Our convention with the Schur polynomial indexing is that 
s<sub>n</sub>=h<sub>n</sub> and s<sub>1,1,...,1</sub>=e<sub>n</sub>.

    sToH[expr,s,h]
    hToS[expr,h,s]

Converts Schur polynomials to complete homogeneous symmetric polynomials and vica versa.
    
    eToH[expr,e,h] 
    hToE[expr,h,e]

Converts elementary symmetric polynomials to complete homogeneous symmetric polynomials and vica versa.
    
    sToX[expr,s,{x,n}]
    eToX[expr,e,{x,n}]
    hToX[expr,h,{x,n}]

Expand Schur, elementary symmetric and complete homogenous polynomials in terms of 
the (symmetric) variables x<sub>1</sub>, x<sub>2</sub>, ... x<sub>n</sub>.
    
    xToS[expr,{x,n},s]
    xToE[expr,{x,n},e]
    xToH[expr,{x,n},h]

Converts symmetric polynomials in the 
variables x<sub>1</sub>, x<sub>2</sub>, ... x<sub>n</sub>
to Schur, elementary symmetric and complete homogenous polynomials. 


### Low-level manipulation of formulas containing symmetric polynomials
    
As usual, the variable names can be omitted below if you use the standard names:

    expandE[expr,e]
    expandH[expr,h]
    expandP[expr,p]
    expandGeneric[expr,var]
    
Expands terms like e<sub>3,2,2</sub> into e<sub>3</sub> e<sub>2</sub><sup>2</sup>.

    collapseE[expr,e]
    collapseH[expr,h]
    collapseP[expr,p]
    collapseGeneric[expr,var]

Collapses terms like e<sub>3</sub> e<sub>2</sub><sup>2</sup> into e<sub>3,2,2</sub>.

    collectS[expr,s]
    collectE[expr,e]
    collectH[expr,h]
    collectP[expr,p]
    collectGeneric[expr,var]

Collects coefficients of variables like e<sub>3,2,2</sub> together.

    coeffsGeneric[expr,var] 

Returns the list of `{idxs,coefficient}` pairs where `idxs` means `mkSubscript[var,idxs]`

    coeffsGenericBi[expr,var1,var2] 

Similar to the previous, but works with two variable sets at the same times, returning a list of
`{idxs1,idxs2,coefficient}` triples


### Individual symmetric functions

    schurWeylCharacter[lam,{x,n}] 

This implements the Weyl character formula for Schur polynomials.
    
    JacobiTrudiE[lam,e] 
    JacobiTrudiH[lam,h] 

Returns the Jacobi-Trudi determinant for the Schur polynomial s<sub>&lambda;</sub> in terms of elementary 
(resp. complete homogeneous) symmetric polynomials.
    
    lookupEinS[lam,s] 
    lookupHinS[lam,s] 

Looks up the expansion of e<sub>&lambda;</sub> resp. h<sub>&lambda;</sub> in terms of Schur polynomials 
(these are cached, so it's only slow when called at the first time).

    lookupSinH[lam,h] 
    lookupSinE[lam,e]

Looks up the expansion of s<sub>&lambda;</sub> terms of in complete homogeneous resp. elementary symmetric polynomials.

    lookupEinH[lam,h] 
    lookupHinE[lam,e] 

Looks up the expansion of e<sub>&lambda;</sub> in terms of complete homogeneous symmetric polynomials
and vica versa.


### Branching rule coefficients

TODO: implement Littlewood-Richardson and other coefficients
    
    EMatrix[lam,mu,n] 

This is the matrix E<sub>&lambda;/&mu;</sub>(n), whose determinant appear in the expansion of total Chern class of tensor products.

    EDet[lam,mu,n] 

The the determinant det[E<sub>&lambda;/&mu;</sub>(n)] which appears in the expansion of total Chern class of tensor products.


### Miscelleneous

    lexSort[Llist] 

Lexicographical sort of a list of lists.

    lexCompare[list1,list1] 

Lexicographical comparison of two lists (-1 = LT, 0 = EQ, +1 = GT).

    lexMaximumBy[list,fun] 

Returns the lexicographically largest element after applying fun.
    
    
Package `SymP1`
===============

This package implements computations in the cohomology and K-theory of symmetric products 
Sym<sup>n</sup>&Popf;<sup>1</sup> = &Popf;<sup>n</sup> of
the complex projective line &Popf;<sup>1</sup>.

The projective space &Popf;<sup>n</sup> with the generator of the cohomology theory is denoted by

    ProjSpace[n,g]

where `n` is the dimension and `g` is the generator. In case of the K-theory the
generator is the class [L] of the tautological line bundle L&rarr;&Popf;<sup>n</sup>; and in case of cohomology
it's _minus 1 times_ the first Chern class of L: u=-c<sub>1</sub>(L). You can also use the syntax
<tt>&Popf;[n,g]</tt> in Mathematica.

The functions

    Dimension[space]
    Variable[space]

return the dimension and the generator of a projective space `ProjSpace[n,u]`.

### Available cohomology theories

The following theories are supported:

    TrivialTheory
    Cohomology
    KTheory
    EquivCohomology
    EquivKTheory   

The "trivial theory" associates the same ring to any space. The equivariant version 
mean T<sup>2</sup> or GL<sub>2</sub> -equivariant cohomology or K-theory, and the
convention is that weigths of the T<sup>2</sup>-representation are denoted by 
&alpha;,&beta; with the corresponding representations (line bundles) are
X<sup>-1</sup> being Y<sup>-1</sup> (note the signs!). So c<sub>1</sub>(X) = -&alpha;.
In the GL<sub>2</sub> case the Chern classes are then c<sub>1</sub>=&alpha;+&beta;
and c<sub>2</sub>=&alpha;&beta;.

Because Mathematica is very dynamically typed, we represent elements of these
rings simply by polynomials; and variables which are not the generator (nor
the equivariant generators) are considered part of the coefficient ring.

    theoryStandardVar[theory] 
    theoryGlobalVars[theory] 

`theoryStandardVar` returns the standard generator name for that theory (which is `u` for cohomology and `L` for K-theory;
while `theoryGlobalVars` returns the set of global names associated with that theory
(this is relevant for the equivariant theories, where it returns {&alpha;,&beta;} resp. {X,Y}).
Note: These symbols are reserved for this use in the global namespace.

The function

    relation[theory,space]  

returns the relation in the cohomology / K-theory ring of the given projective space, while

    normalize[theory,A,space] 

normalizes a class A using the relation, such that only powers 0..n of the generator appear in the results.

### The three maps

There are four important maps between these spaces:

* the diagonal map &Delta;<sup>d</sup> : &Popf;<sup>n</sup> &rarr; &Popf;<sup>n</sup> &times; &Popf;<sup>n</sup> &times; ...  &times; &Popf;<sup>n</sup>;
* the merging map &Psi;<sub>n1,n2,...nk</sub> : &Popf;<sup>n1</sup> &times; &Popf;<sup>n2</sup> &times; ...  &times; &Popf;<sup>nk</sup> &rarr; &Popf;<sup>n1+n2+...+nk</sup> 
* the replicating map &Omega;<sup>d</sup> = &Psi; &compfn; &Delta; : &Popf;<sup>n</sup> &rarr; &Popf;<sup>dn</sup>;
* and the collapsing map &pi; : &Popf;<sup>n</sup> &rarr; pt.

The essence of this package is a set of functions computing pushforward and pullbacks along these maps.

#### Pushforwards

WARNING: These function assume that the class you pass them is already normalized,
that is, it only contains powers 0..n of the generator. If this is not the case,
you have to use the `normalize` function above! 

TODO: maybe do this automatically.

    CollapsePF[theory,A,space] 

computes pushforward of `A` along &pi; from `space` to the point. So `A` lives in `space` which is a `ProjSpace[n,g]`.

    PsiPF[theory,A,spaces,var] 

computes the pushforward of `A` along &Psi;. Here `spaces` should be a list of projective spaces, 
corresponding to their cartesian product. Finally `var` is the name of the generator of the target space,
whose dimension is already determined to be the sum of the dimensions of `spaces`.

    DeltaPF[theory,A,space,vars]
    DeltaPF[theory,A,space,d   ]

computes the pushforward of `A` along the diagonal map &Delta;<sup>d</sup>. Here `vars` is the list
of generators of the different copies of &Popf;<sup>n</sup> in the target; the number `d` is
simply the length of this list. In the second version, the number `d` is specified, and 
the output generators will be the standard generator name with subscripts 1..d.

    OmegaPF[theory,A,space,d    ] 
    OmegaPF[theory,A,space,d,var] 

computes the pushforward of `A` along the replicating map &Omega;<sup>d</sup>. 
If the output generator name `var` is not specified, it will recycle the generator of `space`.

    OmegaVecPF[theory,A,spaces,ds,var]

This computes the pushforward along the composition 
&Psi; &compfn; &Omega;<sup>d1</sup> &times; &Omega;<sup>dk</sup> &times; ... &times; &Omega;<sup>dk</sup>.
It is assumed that `Length[spaces] == Length[ds]`.

#### Pullbacks

WARNING: these are not really tested (and maybe not even implemented in some cases)!
The conventions are the same as with the pushforwards above.

    OmegaPB[theory,A,space,d,var] 

computes the pullback of `A` along the replicating (or power) map
&Omega;<sup>d</sup> = &Psi; &compfn; &Delta; : &Popf;<sup>n</sup> &rarr; &Popf;<sup>dn</sup>;
Note: the conventions are the same as with the pushforwards, so `space` should have dimension `n`
and *not* `d*n` !!

    DeltaPB[theory,A,space,vars] 

computes the pullback of `A` along the diagonal map
&Delta;<sup>d</sup> : &Popf;<sup>n</sup> &rarr; &Popf;<sup>n</sup> &times; &Popf;<sup>n</sup> &times; ...  &times; &Popf;<sup>n</sup>;

    PsiPB[theory,A,spaces,var] 

computes the pullback of `A` along the merging map 
&Psi;<sub>n1,n2,...nk</sub> : &Popf;<sup>n1</sup> &times; &Popf;<sup>n2</sup> &times; ...  &times; &Popf;<sup>nk</sup> &rarr; &Popf;<sup>n1+n2+...+nk</sup> 

Remark: there are also versions of these functions for the individual theories.
But it's recommended that you use the generic versions.


Package `RootLoci`
==================

The goal of this package is to compute characteristic classes of _coincident root loci_.
These are locally closed subvarieties X<sub>&lambda;</sub> &subset; &Popf;<sup>n</sup>,
indexed by partitions &lambda; of n.

### Available classes and invariants

WARNING: not everything is implemented yet!

Invariants:

    EulerChar     
    ToddGenus     
    ChiY          
    PoincarePoly  
    HodgeDeligneE 
    NumberOfPoints
    HasseZeta     

Characteristic classes:

    FundClass             
    CSM                   
    ToddClass             
    HirzebruchClass       
    UnnormHirzebruch      
    MotivicChern          
    EquivFundClass        
    EquivCSM              
    EquivMotivicChern     

Some functions operating on classes:

    Theory[class] 

gives back the theory the given class lives in. For example `Theory[MotivicChern] == KTheory`.

    shortName[class] 

gives back a short string, which can be used as name for example when generating 
tables of these invariants.

### Umbral bases

The "umbral basis" associated to a cohomology theory is a basis in the (generalized)
cohomology ring of &Popf;<sup>n</sup> in which the pushforwards &Psi;<sub>!</sub> are
trivial to compute. This is very useful because it maps the big ring where the multiplication
is induced by &Psi;<sub>!</sub> on the sum (for all) of the cohomology rings of &Popf;<sup>n</sup>
into a normal power series ring in which we can compute as usual.

    umbralBasis[theory,ProjSpace[n,g],k] 

returns the k-th umbral basis element in &Popf;<sup>n</sup> (k should be between 0 and n).

    toUmbralBasis[theory,A,ProjSpace[n,g],z] 

converts from the usual (generalized) cohomology class representation to the umbral basis,
where the umbral variable z<sup>k</sup> corresponds to the k-th umbral basis.

    fromUmbralBasis[theory,A,z,ProjSpace[n,g]] 

converts back from the umbral basis to the usual representation.

### Various classes 

    classOfPn[class,ProjSpace[n,g]] 

The class of Sym<sup>n</sup>&Popf;<sup>1</sup> = &Popf;<sup>n</sup>

    ghostClass[class,ProjSpace[m,g]] 

This returns the m-th ghost class (computed via hardcoded formulas). These are
the coefficients of the logarithmic differential of the "relative zeta" (or exponential),
that is, the power series whose coefficients are the symmetric products.

    abstractGhostClass[m,S] 

Returns the m-th "abstract" ghost class in terms of the variables S<sub>n</sub> which
corresponds to the (class of) Sym<sup>n</sup>X.

    genericGhostClass[class,space] 

Returns the m-th ghost class computed using pushforwards from the classes 
of &Popf;<sup>n</sup> using the "abstract ghost class" above. This should give
the same result as `ghostClass` above.

    umbralGhostClass[class,m,z] 

Returns the m-th ghost class in the "umbral z-coordinates" (via hardcoded formulas)

### The &Lscr; and &Pscr; power series 

    LSeries[n]
    PSeries[n]

The &Lscr; and &Pscr; power series up to degree n.

    LSeries[n,d]
    PSeries[n,d]

The same but after the substitution x<sub>i</sub> x&map; (x<sub>i</sub>)<sup>d</sup>.

    LPoly[k]
    PPoly[k]

The degree k part (weighted homogeneous polynomial) of the above. Note: deg(x<sub>i</sub>)=i !

    LPoly[k,d]
    PPoly[k,d]

The same but after the substitution x<sub>i</sub> x&map; (x<sub>i</sub>)<sup>d</sup>.

### Classes and invariants of coincident root loci

These functions compute the classes using the recursive algorithm.

    recClassOfRootLoci[class, lambda, gen] 

computes the class of X<sub>&lambda;</sub> &subset; &Popf;<sup>n</sup>.

    recClassOfDistinctLoci[class, nvec, gen] 

computes the class of 
D(n1,n2,...,nk) &subset; &Popf;<sup>n1</sup> &times; &Popf;<sup>n2</sup> &times; ...  &times; &Popf;<sup>nk</sup>
which is the set of n1+n2+...+nk points which are all distinct from each other.

    ExportRootLoci[class, fname, n      ] 
    ExportRootLoci[class, fname, n ,var ]

exports the classes of coincident root loci (up to partitions of weight at most n) 
into a text file.

### Closure of the strata

The coincident root loci X<sub>&lambda;</sub> stratifies &Popf;<sup>n</sup>, and
the closure of a stratum is the union of it and other strata. On the level of
partitions, this corresponds to merging some parts into one.

    PartitionClosure[lambda] 
    PartitionClosure$v1[lambda]

both returns the set of coarsenings of the partition `lambda`, which corresponds to
the stratification of the closure of the strata X<sub>&lambda;</sub>. These are 
two different implementations, which should give the same but the `$v1` is slower.

