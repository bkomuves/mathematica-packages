
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
    
    
Package `SymP1`
===============


Package `RootLoci`
==================
