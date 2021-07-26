
My personal Mathematica packages
--------------------------------

I thought it may be a good idea to finally convert my random pieces of code
into a reusable format...


### Quickstart

Add these to `~/Library/Mathematica/Kernel/init.m` so that we are in the path:

    AppendTo[$Path, ToFileName["/Users/bkomuves/math/Mathematica"]];
    AppendTo[$Path, ToFileName["/Users/bkomuves/math/Mathematica/MyPackages"]];

Also you can add anything else you want to auto-run at start; for example loading
certain packages:

    << Useful`
    << Schur`


### Available packages

* `Useful` - all kind of basic functions, which are missing from Mathematica
* `Schur` - symmetric polynomials
* `SymP1` - computations in the cohomology and K-theory of Sym^n(P^1) = P^n
* `RootLoci` - computing characteristic classes of coincident root loci 
* `Witt` - computations in the (big) Witt ring

