
My personal Mathematica packages
--------------------------------

I thought it may be a good idea to finally convert my random pieces of code
into a reusable set of libraries...

WARNING: this codebase is very much in flux and probably is full of bugs.
If you find a bug please send me an email to _username_ (at) gmail.


### Quickstart

Copy the `*.m` files into a subdirectory, for example 

    /Users/<username>/math/Mathematica/MyPackages/

Add the following line to `~/Library/Mathematica/Kernel/init.m`  (on MacOS;
on other operating system this file should be at a similar location) so 
that we are in the path:

    AppendTo[$Path, ToFileName["/Users/<username>/math/Mathematica/MyPackages"]];

Also you can add anything else you want to auto-run at start; for example loading
certain packages:

    << Useful`
    << Partitions`
    << Schur`


### Available packages

* `Useful` - all kind of basic functions, which are missing from Mathematica
* `Partitions` - functions related to (integer) partitions
* `Schur` - symmetric polynomials
* `SymP1` - computations in the cohomology and K-theory of Sym^n(P^1) = P^n
* `RootLoci` - computing characteristic classes of coincident root loci 
* `Witt` - computations in the (big) Witt ring


### Documentation

See the `USERGUIDE.md` file.
