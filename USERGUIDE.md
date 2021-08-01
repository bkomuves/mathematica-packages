
Documentation for my Mathematica packages
-----------------------------------------

WARNING: this codebase is very much in flux and probably is full of bugs.
If you find a bug please send me an email to _username_ (at) gmail.

Table of content:

* `Useful` - all kind of basic functions, which are missing from Mathematica
* `Partitions` - functions related to (integer) partitions
* `Schur` - symmetric polynomials
* `SymP1` - computations in the cohomology and K-theory of Sym^n(P^1) = P^n
* `RootLoci` - computing characteristic classes of coincident root loci 


Package `Useful`
================


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


Package `RootLoci`
==================
