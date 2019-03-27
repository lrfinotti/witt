# Witt Vectors and Canonical Liftings

This contains various routines to run on
[MAGMA](http://magma.maths.usyd.edu.au/magma/) to perform computations
with Witt vectors and canonical liftings.

## Files

* `etas.m`: Contains the auxiliary "eta functions" for computations
  with Witt vectors.

* `witt.m`: Contains functions to perform operations with Witt
  vectors.  (Depends on `etas.m`.)

* `gt.m`: Contains routines to compute the Greenberg transform of
  polynomials in two variables, allowing in particular, to evaluate
  these polynomials at a pair of Witt vectors.  It also contains
  functions to compute powers of Witt vectors.  (Depends on `witt.m`.)

* `lift.m`: Contains routines to compute canonical liftings and the
  elliptic Teichmüller lift (up to 3 coordinates) and minimal degree
  liftings of ordinary elliptic curves.  (Depends on `gt.m`.)

* `lift_j.m`: contains routines to find the coordinates of the
  j-invariant of the canonical lifting.  (Depends on `gt.m`.)


## Theory

The routines are based on the methods developed in [Computations with
Witt Vectors and the Greenberg Transform][comp].  It provides more
efficient ways to perform computations with Witt vectors.  In
particular, one does not need the (*huge*!) polynomials that define
the sum and product of Witt vectors.

The algorithm to compute the canonical lifting is described in
[Degrees of the Elliptic Teichmuller Lift][canlift], and the one to
compute minimal degree liftings is described in [Minimal Degree
Liftings of Hyperelliptic Curves][minlift].

For the liftings of the j-invariant, see [Coordinates of the
j-Invariant of the Canonical Lifting][jinv].



[comp]: http://www.math.utk.edu/~finotti/papers/witt.pdf

[canlift]: http://www.math.utk.edu/~finotti/papers/degs.pdf

[minlift]: http://www.math.utk.edu/~finotti/papers/mindeg.pdf

[jinv]: http://www.math.utk.edu/~finotti/papers/jn.pdf



## Usage

The "eta functions" (described in Section 5 of [Computations with Witt
Vectors and the Greenberg Transform][comp]) that are used in all
computations can be computed in three different ways:

1. We can store some simpler version (in two variables) in memory and
   used them in the computations.

2. Store some more basic data that take much less memory and compute
   the eta functions "on the fly".  (This is describe in Section 8 of
   [Computations with Witt Vectors and the Greenberg
   Transform][comp]).  This is usually slower, but requires less
   memory.  In some cases, where the eta functions are too large, and
   one does not need many computations, it might be actually faster
   than method 1.
   
3. Use method 2, but store some values of resulting computations of
   the eta function, so that we don't need to recompute them, which is
   likely to happen with method 2.  Of course, it uses more memory
   than method 2, but likely less memory than method 1.


### Basic Computations with Witt Vectors

The following functions are provided in `witt.m`:

* `WittSum`: adds two Witt vectors.
* `WittProd`: multiplies two Witt vectors.
* `WittNeg`: gives the additive inverse of a Witt vector.
* `WittInv`: gives the multiplicative inverse of a Witt vector (if
  exists).
  
The file `gt.m` also provides `WittPower` to compute powers of Witt
vectors.
  
Each of these functions have an optional argument `choice` that allows
you to choose which method (1, 2 or 3, as above) to use in the
computations.  It defaults to method 1, which seems to be usually the
fastest one.

**Note:** If any of the above functions are called to compute with
Witt vectors over finite fields, then the operations are performed by
converting to unramified extensions of p-adic integers, performing the
operation in that ring, and then converting back to Witt vectors.
This method is *a lot* more efficient than any of the methods above!
  
Here is one simple example:
```
> load "gt.m";
Loading "gt.m"
Loading "witt.m"
Loading "etas.m"

> p:=5; d:=3; F<a>:=GF(p^d); P<x,y>:=PolynomialRing(F,2);
> v := [ x, x*y, P!2 ];
> w := [ x^2, x-y, 1+y^2];
> WittSum(v,w);
[
    x^2 + x,
    4*x^9 + 3*x^8 + 3*x^7 + 4*x^6 + x*y + x + 4*y,
    4*x^49 + 3*x^48 + 3*x^47 + 4*x^46 + 3*x^44 + 2*x^43 + x^41 + 4*x^38 + 4*x^37*y + 3*x^37 + 3*x^36*y + 2*x^36 + x^35*y + 3*x^35 + 
        3*x^34*y + 2*x^34 + 4*x^33*y + 3*x^32*y + x^31 + 3*x^30 + 2*x^29*y^2 + x^29*y + x^29 + 3*x^28*y^2 + x^28*y + x^28 + 
        4*x^27*y^2 + 2*x^27 + 2*x^26*y + x^26 + 2*x^25*y^2 + 4*x^25*y + 4*x^24*y^2 + x^24*y + x^24 + 4*x^23*y^2 + 3*x^23*y + 
        2*x^22*y^2 + 2*x^22*y + x^22 + 3*x^21*y^3 + 4*x^21*y^2 + x^21*y + 3*x^20*y^3 + 2*x^20*y^2 + 2*x^20*y + 4*x^20 + 2*x^19*y^3 + 
        2*x^19*y^2 + 2*x^19*y + 4*x^19 + x^18*y^3 + 4*x^18*y^2 + 3*x^18*y + 4*x^17*y^3 + 4*x^17*y^2 + 2*x^17*y + 4*x^17 + x^16*y^3 + 
        2*x^16*y^2 + 4*x^16*y + 2*x^16 + 4*x^15*y^3 + 4*x^15*y^2 + 3*x^15*y + 3*x^15 + 3*x^14*y^3 + 3*x^14*y^2 + x^14*y + x^13*y^4 + 
        x^13*y^3 + 4*x^13*y + x^13 + 3*x^12*y^4 + 3*x^12*y^3 + 4*x^12*y + 2*x^12 + x^11*y^3 + 4*x^11*y^2 + 2*x^11 + x^10*y^4 + 
        4*x^10*y^2 + x^10*y + x^10 + x^9*y^4 + 4*x^9*y^3 + x^9*y + 4*x^8*y^3 + x^8*y^2 + 3*x^7*y^4 + x^7*y^3 + x^6*y^4 + 4*x^5*y^4 + 
        3*x^5*y^3 + 3*x^5*y^2 + 4*x^5*y + x^4*y^5 + 4*x^4*y^4 + x^4*y^3 + 4*x^4*y^2 + 3*x^3*y^5 + 4*x^3*y^4 + 4*x^3*y^3 + 2*x^2*y^5 +
        4*x^2*y^4 + 4*x*y^5 + y^2 + 3
]

> WittProd(v,w);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 
        4*x^35*y + 4*x^34*y^2 + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 4*x^5*y^10
]

> WittNeg(v); 
[
    4*x,
    4*x*y,
    3
]

> WittPower(v,5);
[
    x^5,
    0,
    x^105*y^5
]

```
Note that we cannot invert `v` or `w` since they are not units.
```
> WittInv(v);

WittInv(
    v: [ x, x*y, 2 ]
)
WittInv1(
    v: [ x, x*y, 2 ]
)
In file "witt.m", line 325, column 19:
>>         w2 := [ PR!res[j] : j in [1..i] ] cat [x];
                     ^
Runtime error in '!': Illegal coercion
```
But:
```
> WittInv([P!2,x,y]);
[
    3,
    x,
    2*x^10 + y
]
```

Every time these functions are called, the eta polynomials are
computed, wasting time.  If one needs to make many computations for
longer length or large characteristic, it is better to save these and
reuse them.  For example, continuing the example above, we could save
the eta polynomials (for prime `p`, defined to be 5 in the example, and
length 3) with: 
```
epols:=etapols(p,2);
```
(Note that we use `2` as the second argument to give length of 3.
This is because we often see vectors of length `n+1` as
`[a0,a1,...,an]`, so `n` corresponds to the last coordinate (when
counting from 0.)  We use this convention through out.)
Then, when we need to compute with Witt vectors of length 3 with
entries in characteristic `p`, we can pass these polynomials as an
optional argument `pols`.  The routines then use those, instead of
computing them each time:
```
> WittSum(v,w : pols:=epols);
[
    x^2 + x,
    4*x^9 + 3*x^8 + 3*x^7 + 4*x^6 + x*y + x + 4*y,
    4*x^49 + 3*x^48 + 3*x^47 + 4*x^46 + 3*x^44 + 2*x^43 + x^41 + 4*x^38 + 
        4*x^37*y + 3*x^37 + 3*x^36*y + 2*x^36 + x^35*y + 3*x^35 + 3*x^34*y + 
        2*x^34 + 4*x^33*y + 3*x^32*y + x^31 + 3*x^30 + 2*x^29*y^2 + x^29*y + 
        x^29 + 3*x^28*y^2 + x^28*y + x^28 + 4*x^27*y^2 + 2*x^27 + 2*x^26*y + 
        x^26 + 2*x^25*y^2 + 4*x^25*y + 4*x^24*y^2 + x^24*y + x^24 + 4*x^23*y^2 +
        3*x^23*y + 2*x^22*y^2 + 2*x^22*y + x^22 + 3*x^21*y^3 + 4*x^21*y^2 + 
        x^21*y + 3*x^20*y^3 + 2*x^20*y^2 + 2*x^20*y + 4*x^20 + 2*x^19*y^3 + 
        2*x^19*y^2 + 2*x^19*y + 4*x^19 + x^18*y^3 + 4*x^18*y^2 + 3*x^18*y + 
        4*x^17*y^3 + 4*x^17*y^2 + 2*x^17*y + 4*x^17 + x^16*y^3 + 2*x^16*y^2 + 
        4*x^16*y + 2*x^16 + 4*x^15*y^3 + 4*x^15*y^2 + 3*x^15*y + 3*x^15 + 
        3*x^14*y^3 + 3*x^14*y^2 + x^14*y + x^13*y^4 + x^13*y^3 + 4*x^13*y + x^13
        + 3*x^12*y^4 + 3*x^12*y^3 + 4*x^12*y + 2*x^12 + x^11*y^3 + 4*x^11*y^2 + 
        2*x^11 + x^10*y^4 + 4*x^10*y^2 + x^10*y + x^10 + x^9*y^4 + 4*x^9*y^3 + 
        x^9*y + 4*x^8*y^3 + x^8*y^2 + 3*x^7*y^4 + x^7*y^3 + x^6*y^4 + 4*x^5*y^4 
        + 3*x^5*y^3 + 3*x^5*y^2 + 4*x^5*y + x^4*y^5 + 4*x^4*y^4 + x^4*y^3 + 
        4*x^4*y^2 + 3*x^3*y^5 + 4*x^3*y^4 + 4*x^3*y^3 + 2*x^2*y^5 + 4*x^2*y^4 + 
        4*x*y^5 + y^2 + 3
]

> WittProd(v,w : pols:=epols);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]
```

Remember that by default, all these functions use method 1.  You can
specify a different method with the optional argument `choice`.  To
use method 2:
```
> WittProd(v,w : choice:=2);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]
```
For method 3:
```
> WittProd(v,w : choice:=3);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]
```

As observed above, methods 2 and 3 do not use the eta polynomials, but some
smaller data, which are in fact just some binomials.  To save time
then, we can precompute these data and reuse it.  For length 3 and
characteristic `p`, we can save the date with:
```
bt:=BinTab(p,2);
```
And then we can use it by using the optional argument `bintab`:
```
> WittProd(v,w : choice:=2, bintab:=bt);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]

> WittProd(v,w : choice:=3, bintab:=bt);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]
```

**Important:** If you are performing many computations with Witt
vectors is highly recommended that you store and reuse either these
eta polynomials (for method 1) or the binomials (for methods 2 and 3)
as above.  Otherwise, they will be computed every time you call
perform a computation, and for longer lengths or higher
characteristics, these alone can take some time.


### Conversions

The file `etas.m` also provide a few basic functions.

* `IntToWitt` converts an integer to a Witt vector of a given length
  an entries in a specified prime field.

* `WittVToSeries` converts a Witt vector over a finite field to a
  p-adic (truncated) power series in an unramified extension of the
  p-adic integers.

* `SeriesToWittV` converts a a
  p-adic (truncated) power series in an unramified extension of the
  p-adic integers to a Witt vector.


Here are some examples.

Let's convert 1728 into a Witt vector of length 3 over the field
with 5 elements.
```
> IntToWitt(1728,5,2);
[ 3, 2, 0 ]
```
(Again, note we use second argument as `2` to give length 3.)

Converting from Witt vector to power series:
```
> F<a> := GF(5^3);
> WittVToSeries([a,F!2,a^7,a^100]);
-190*aa^2 - 219*aa - 305 + O(5^4)

> s := WittVToSeries([a,F!2,a^7,a^100]); s;
-190*aa^2 - 219*aa - 305 + O(5^4)

> Parent(s);
Unramified extension defined by the polynomial (1 + O(5^4))*x^3 + O(5^4)*x^2 + 
    (3 + O(5^4))*x + 3 + O(5^4)
 over 5-adic ring

> F!s;
a
```


Now, converting from series to Witt vector:
```
> Z5 := pAdicRing(5 : Precision:=4);     
> Z125<aa> := ext<Z5 | 3>;
> Z125;
Unramified extension defined by the polynomial (1 + O(5^4))*x^3 + O(5^4)*x^2 + 
    (3 + O(5^4))*x + 3 + O(5^4)
 over 5-adic ring

> F!aa;
a

> x := Random(Z125); x;
-120*aa^2 + 59*aa - 241 + O(5^4)

> SeriesToWittV(x);
[ a^34, 4, a^40, a^65 ]
```

### Greenberg Transform


The file `gt.m` provides the following functions:

* `GT`: computes the Greenberg transform of polynomials (over Witt
  vectors) in two variables.  In particular, it provides an efficient
  way to evaluate such polynomials.

* `Pol_GT_Form`: Converts a polynomial over the integers in two
  variables to the form that is needed for the input of `GT`.

The input for `GT` (a polynomial in two variables with coefficients in
the ring of Witt vectors) must be given in the form ` [ m1, m2, ... ]`
where each `mi` corresponds to a monomial and has the form `[* vi, ki,
li *]`, where `vi` is the coefficient of the monomial, and `ki` and `li` are the
powers of the first and second variable, respectively, for the
corresponding monomial.

Here are some examples:
```
> p:=3; d:=4; F<a>:=GF(p^d);
> pol := [ [* [a,F!2,a^5], 2,0 *], [* [a,F!0,F!1], 1, 1 *], [* [F!1,a^7,F!0], 0, 0 *] ];
> GT(pol);
[
    a*x0^2 + a*x0*y0 + 1,
    2*x0^6 + a^43*x0^5*y0 + a^43*x0^4*y0^2 + a^42*x0^4 + a^43*x0^3*x1 + a^2*x0^3*y0 + a^3*x0^3*y1 + a^42*x0^2*y0^2 + a^41*x0^2 + 
        a^41*x0*y0 + a^3*x1*y0^3 + a^7,
    a^5*x0^18 + a^38*x0^17*y0 + a^72*x0^16*y0^2 + a^37*x0^16 + a^3*x0^15*x1 + a^46*x0^15*y0^3 + a^31*x0^15*y0 + a^43*x0^15*y1 + 
        a^46*x0^14*x1*y0 + a^37*x0^14*y0^4 + a^2*x0^14*y0^2 + a^6*x0^14*y0*y1 + a^70*x0^14 + a^77*x0^13*x1*y0^2 + a^45*x0^13*x1 + 
        a^49*x0^13*y0^5 + a^48*x0^13*y0^3 + a^37*x0^13*y0^2*y1 + a*x0^13*y0 + a^5*x0^13*y1 + a^6*x0^12*x1^2 + a^57*x0^12*x1*y0^3 + 
        a^36*x0^12*x1*y0 + a^6*x0^12*x1*y1 + a^49*x0^12*y0^4 + a^9*x0^12*y0^3*y1 + a^4*x0^12*y0^2 + a^76*x0^12*y0*y1 + a^6*x0^12*y1^2
        + a^34*x0^12 + a^9*x0^11*x1^2*y0 + a^10*x0^11*x1*y0^4 + a^45*x0^11*x1*y0^2 + a^9*x0^11*x1*y0*y1 + a^75*x0^11*x1 + 
        a^49*x0^11*y0^7 + a^48*x0^11*y0^5 + a^49*x0^11*y0^4*y1 + a^47*x0^11*y0^3 + a^5*x0^11*y0^2*y1 + a^9*x0^11*y0*y1^2 + 
        a^9*x0^11*y0 + a^35*x0^11*y1 + a^9*x0^10*x1^2*y0^2 + a^8*x0^10*x1^2 + a^37*x0^10*x1*y0^5 + a^5*x0^10*x1*y0^3 + 
        a^9*x0^10*x1*y0^2*y1 + a^44*x0^10*x1*y0 + a^8*x0^10*x1*y1 + a^49*x0^10*y0^8 + a^35*x0^10*y0^4 + a^9*x0^10*y0^2*y1^2 + 
        a^52*x0^10*y0^2 + a^4*x0^10*y0*y1 + a^8*x0^10*y1^2 + a^16*x0^10 + a^12*x0^9*x1^3 + a^6*x0^9*x1^2*y0^3 + a^48*x0^9*x1^2*y0 + 
        a^49*x0^9*x1^2*y1 + a^9*x0^9*x1*y0^6 + a^45*x0^9*x1*y0^4 + a^46*x0^9*x1*y0^3*y1 + a^7*x0^9*x1*y0^2 + a^48*x0^9*x1*y0*y1 + 
        a^9*x0^9*x1*y1^2 + a^9*x0^9*x1 + a^49*x0^9*x2 + x0^9*y0^9 + a^8*x0^9*y0^7 + a^8*x0^9*y0^4*y1 + a^73*x0^9*y0^3 + 
        a^47*x0^9*y0^2*y1 + a^48*x0^9*y0*y1^2 + a^35*x0^9*y0 + a^49*x0^9*y1 + a^9*x0^9*y2 + a^9*x0^8*x1^2*y0^4 + a^8*x0^8*x1^2*y0^2 +
        a^7*x0^8*x1^2 + a^49*x0^8*x1*y0^7 + a^5*x0^8*x1*y0^5 + a^49*x0^8*x1*y0^4*y1 + a^35*x0^8*x1*y0^3 + a^8*x0^8*x1*y0^2*y1 + 
        a^13*x0^8*x1*y0 + a^7*x0^8*x1*y1 + a^48*x0^8*y0^8 + a^28*x0^8*y0^4 + a^8*x0^8*y0^2*y1^2 + a^16*x0^8*y0^2 + a^53*x0^8*y0*y1 + 
        a^7*x0^8*y1^2 + a^50*x0^8 + a^9*x0^7*x1^2*y0^5 + a^8*x0^7*x1^2*y0^3 + a^7*x0^7*x1^2*y0 + a^49*x0^7*x1*y0^5*y1 + 
        a^8*x0^7*x1*y0^4 + a^48*x0^7*x1*y0^3*y1 + a^13*x0^7*x1*y0^2 + a^7*x0^7*x1*y0*y1 + a^67*x0^7*x1 + a^47*x0^7*y0^7 + 
        a^46*x0^7*y0^5 + a^47*x0^7*y0^4*y1 + a^45*x0^7*y0^3 + a^53*x0^7*y0^2*y1 + a^7*x0^7*y0*y1^2 + a^7*x0^7*y0 + a^27*x0^7*y1 + 
        a^49*x0^6*x1^3*y0^3 + a^6*x0^6*x1^2*y0^6 + a^48*x0^6*x1^2*y0^4 + a^49*x0^6*x1^2*y0^3*y1 + a^53*x0^6*x1^2 + a^8*x0^6*x1*y0^7 +
        a^47*x0^6*x1*y0^5 + a^8*x0^6*x1*y0^4*y1 + a^49*x0^6*x1*y0^3*y1^2 + a^50*x0^6*x1*y0^3 + a^27*x0^6*x1*y0 + a^53*x0^6*x1*y1 + 
        a^29*x0^6*y0^4 + a^6*x0^6*y0^3*y1 + a^51*x0^6*y0^2 + a^67*x0^6*y0*y1 + a^53*x0^6*y1^2 + a*x0^6 + a^9*x0^5*x1^2*y0^7 + 
        a^8*x0^5*x1^2*y0^5 + a^7*x0^5*x1^2*y0^3 + a^48*x0^5*x1*y0^5*y1 + a^53*x0^5*x1*y0^4 + a^47*x0^5*x1*y0^3*y1 + a^67*x0^5*x1*y0^2
        + a^11*x0^5*x1 + a^45*x0^5*y0^5 + a^44*x0^5*y0^3 + a^27*x0^5*y0^2*y1 + a^17*x0^5*y0 + a^51*x0^5*y1 + a^9*x0^4*x1^2*y0^8 + 
        a^8*x0^4*x1^2*y0^6 + a^7*x0^4*x1^2*y0^4 + a^47*x0^4*x1*y0^7 + a^53*x0^4*x1*y0^5 + a^47*x0^4*x1*y0^4*y1 + a^27*x0^4*x1*y0^3 + 
        a^11*x0^4*x1*y0 + a^26*x0^4*y0^4 + a^17*x0^4*y0^2 + a^51*x0^4*y0*y1 + a^45*x0^4 + a^9*x0^3*x1^3*y0^6 + a^48*x0^3*x1^2*y0^7 + 
        a^49*x0^3*x1^2*y0^6*y1 + a^53*x0^3*x1^2*y0^3 + a^6*x0^3*x1*y0^6 + a^67*x0^3*x1*y0^4 + a^13*x0^3*x1*y0^3*y1 + a^17*x0^3*x1 + 
        a^10*x0^3*y0^3 + a^5*x0^3*y0 + a^57*x0^3*y1 + a^8*x0^2*x1^2*y0^8 + a^7*x0^2*x1^2*y0^6 + a^27*x0^2*x1*y0^5 + a^51*x0^2*x1*y0^3
        + a^45*x0^2*y0^2 + a^7*x0^2 + a^7*x0*x1^2*y0^7 + a^51*x0*x1*y0^4 + a^7*x0*y0 + a^9*x1^6 + a^9*x1^3*y1^3 + a^53*x1^2*y0^6 + 
        a^57*x1*y0^3 + a^9*x2*y0^9
]
```

As with the other functions, you can specify the method for the
computations with the optional argument `choice`.  For the method 1
(the default), again as above, you can use the optional argument
`pols` to give the function precomputed eta polynomials.  For methods
2 or 3, you can use the optional argument `bintab` to pass along the
table of the used table of binomials used with those methods.

For example:
```
> epols := etapols(p,2);
> GT(pol : pols:=epols);
[
    a*x0^2 + a*x0*y0 + 1,
    2*x0^6 + a^43*x0^5*y0 + a^43*x0^4*y0^2 + a^42*x0^4 + a^43*x0^3*x1 + a^2*x0^3*y0 + a^3*x0^3*y1 + a^42*x0^2*y0^2 + a^41*x0^2 + 
        a^41*x0*y0 + a^3*x1*y0^3 + a^7,
    a^5*x0^18 + a^38*x0^17*y0 + a^72*x0^16*y0^2 + a^37*x0^16 + a^3*x0^15*x1 + a^46*x0^15*y0^3 + a^31*x0^15*y0 + a^43*x0^15*y1 + 
        a^46*x0^14*x1*y0 + a^37*x0^14*y0^4 + a^2*x0^14*y0^2 + a^6*x0^14*y0*y1 + a^70*x0^14 + a^77*x0^13*x1*y0^2 + a^45*x0^13*x1 + 
        a^49*x0^13*y0^5 + a^48*x0^13*y0^3 + a^37*x0^13*y0^2*y1 + a*x0^13*y0 + a^5*x0^13*y1 + a^6*x0^12*x1^2 + a^57*x0^12*x1*y0^3 + 
        a^36*x0^12*x1*y0 + a^6*x0^12*x1*y1 + a^49*x0^12*y0^4 + a^9*x0^12*y0^3*y1 + a^4*x0^12*y0^2 + a^76*x0^12*y0*y1 + a^6*x0^12*y1^2
        + a^34*x0^12 + a^9*x0^11*x1^2*y0 + a^10*x0^11*x1*y0^4 + a^45*x0^11*x1*y0^2 + a^9*x0^11*x1*y0*y1 + a^75*x0^11*x1 + 
        a^49*x0^11*y0^7 + a^48*x0^11*y0^5 + a^49*x0^11*y0^4*y1 + a^47*x0^11*y0^3 + a^5*x0^11*y0^2*y1 + a^9*x0^11*y0*y1^2 + 
        a^9*x0^11*y0 + a^35*x0^11*y1 + a^9*x0^10*x1^2*y0^2 + a^8*x0^10*x1^2 + a^37*x0^10*x1*y0^5 + a^5*x0^10*x1*y0^3 + 
        a^9*x0^10*x1*y0^2*y1 + a^44*x0^10*x1*y0 + a^8*x0^10*x1*y1 + a^49*x0^10*y0^8 + a^35*x0^10*y0^4 + a^9*x0^10*y0^2*y1^2 + 
        a^52*x0^10*y0^2 + a^4*x0^10*y0*y1 + a^8*x0^10*y1^2 + a^16*x0^10 + a^12*x0^9*x1^3 + a^6*x0^9*x1^2*y0^3 + a^48*x0^9*x1^2*y0 + 
        a^49*x0^9*x1^2*y1 + a^9*x0^9*x1*y0^6 + a^45*x0^9*x1*y0^4 + a^46*x0^9*x1*y0^3*y1 + a^7*x0^9*x1*y0^2 + a^48*x0^9*x1*y0*y1 + 
        a^9*x0^9*x1*y1^2 + a^9*x0^9*x1 + a^49*x0^9*x2 + x0^9*y0^9 + a^8*x0^9*y0^7 + a^8*x0^9*y0^4*y1 + a^73*x0^9*y0^3 + 
        a^47*x0^9*y0^2*y1 + a^48*x0^9*y0*y1^2 + a^35*x0^9*y0 + a^49*x0^9*y1 + a^9*x0^9*y2 + a^9*x0^8*x1^2*y0^4 + a^8*x0^8*x1^2*y0^2 +
        a^7*x0^8*x1^2 + a^49*x0^8*x1*y0^7 + a^5*x0^8*x1*y0^5 + a^49*x0^8*x1*y0^4*y1 + a^35*x0^8*x1*y0^3 + a^8*x0^8*x1*y0^2*y1 + 
        a^13*x0^8*x1*y0 + a^7*x0^8*x1*y1 + a^48*x0^8*y0^8 + a^28*x0^8*y0^4 + a^8*x0^8*y0^2*y1^2 + a^16*x0^8*y0^2 + a^53*x0^8*y0*y1 + 
        a^7*x0^8*y1^2 + a^50*x0^8 + a^9*x0^7*x1^2*y0^5 + a^8*x0^7*x1^2*y0^3 + a^7*x0^7*x1^2*y0 + a^49*x0^7*x1*y0^5*y1 + 
        a^8*x0^7*x1*y0^4 + a^48*x0^7*x1*y0^3*y1 + a^13*x0^7*x1*y0^2 + a^7*x0^7*x1*y0*y1 + a^67*x0^7*x1 + a^47*x0^7*y0^7 + 
        a^46*x0^7*y0^5 + a^47*x0^7*y0^4*y1 + a^45*x0^7*y0^3 + a^53*x0^7*y0^2*y1 + a^7*x0^7*y0*y1^2 + a^7*x0^7*y0 + a^27*x0^7*y1 + 
        a^49*x0^6*x1^3*y0^3 + a^6*x0^6*x1^2*y0^6 + a^48*x0^6*x1^2*y0^4 + a^49*x0^6*x1^2*y0^3*y1 + a^53*x0^6*x1^2 + a^8*x0^6*x1*y0^7 +
        a^47*x0^6*x1*y0^5 + a^8*x0^6*x1*y0^4*y1 + a^49*x0^6*x1*y0^3*y1^2 + a^50*x0^6*x1*y0^3 + a^27*x0^6*x1*y0 + a^53*x0^6*x1*y1 + 
        a^29*x0^6*y0^4 + a^6*x0^6*y0^3*y1 + a^51*x0^6*y0^2 + a^67*x0^6*y0*y1 + a^53*x0^6*y1^2 + a*x0^6 + a^9*x0^5*x1^2*y0^7 + 
        a^8*x0^5*x1^2*y0^5 + a^7*x0^5*x1^2*y0^3 + a^48*x0^5*x1*y0^5*y1 + a^53*x0^5*x1*y0^4 + a^47*x0^5*x1*y0^3*y1 + a^67*x0^5*x1*y0^2
        + a^11*x0^5*x1 + a^45*x0^5*y0^5 + a^44*x0^5*y0^3 + a^27*x0^5*y0^2*y1 + a^17*x0^5*y0 + a^51*x0^5*y1 + a^9*x0^4*x1^2*y0^8 + 
        a^8*x0^4*x1^2*y0^6 + a^7*x0^4*x1^2*y0^4 + a^47*x0^4*x1*y0^7 + a^53*x0^4*x1*y0^5 + a^47*x0^4*x1*y0^4*y1 + a^27*x0^4*x1*y0^3 + 
        a^11*x0^4*x1*y0 + a^26*x0^4*y0^4 + a^17*x0^4*y0^2 + a^51*x0^4*y0*y1 + a^45*x0^4 + a^9*x0^3*x1^3*y0^6 + a^48*x0^3*x1^2*y0^7 + 
        a^49*x0^3*x1^2*y0^6*y1 + a^53*x0^3*x1^2*y0^3 + a^6*x0^3*x1*y0^6 + a^67*x0^3*x1*y0^4 + a^13*x0^3*x1*y0^3*y1 + a^17*x0^3*x1 + 
        a^10*x0^3*y0^3 + a^5*x0^3*y0 + a^57*x0^3*y1 + a^8*x0^2*x1^2*y0^8 + a^7*x0^2*x1^2*y0^6 + a^27*x0^2*x1*y0^5 + a^51*x0^2*x1*y0^3
        + a^45*x0^2*y0^2 + a^7*x0^2 + a^7*x0*x1^2*y0^7 + a^51*x0*x1*y0^4 + a^7*x0*y0 + a^9*x1^6 + a^9*x1^3*y1^3 + a^53*x1^2*y0^6 + 
        a^57*x1*y0^3 + a^9*x2*y0^9
]

> bt := BinTab(p,2);
> GT(pol : choice :=2, bintab:=bt);
[
    a*x0^2 + a*x0*y0 + 1,
    2*x0^6 + a^43*x0^5*y0 + a^43*x0^4*y0^2 + a^42*x0^4 + a^43*x0^3*x1 + a^2*x0^3*y0 + a^3*x0^3*y1 + a^42*x0^2*y0^2 + a^41*x0^2 + 
        a^41*x0*y0 + a^3*x1*y0^3 + a^7,
    a^5*x0^18 + a^38*x0^17*y0 + a^72*x0^16*y0^2 + a^37*x0^16 + a^3*x0^15*x1 + a^46*x0^15*y0^3 + a^31*x0^15*y0 + a^43*x0^15*y1 + 
        a^46*x0^14*x1*y0 + a^37*x0^14*y0^4 + a^2*x0^14*y0^2 + a^6*x0^14*y0*y1 + a^70*x0^14 + a^77*x0^13*x1*y0^2 + a^45*x0^13*x1 + 
        a^49*x0^13*y0^5 + a^48*x0^13*y0^3 + a^37*x0^13*y0^2*y1 + a*x0^13*y0 + a^5*x0^13*y1 + a^6*x0^12*x1^2 + a^57*x0^12*x1*y0^3 + 
        a^36*x0^12*x1*y0 + a^6*x0^12*x1*y1 + a^49*x0^12*y0^4 + a^9*x0^12*y0^3*y1 + a^4*x0^12*y0^2 + a^76*x0^12*y0*y1 + a^6*x0^12*y1^2
        + a^34*x0^12 + a^9*x0^11*x1^2*y0 + a^10*x0^11*x1*y0^4 + a^45*x0^11*x1*y0^2 + a^9*x0^11*x1*y0*y1 + a^75*x0^11*x1 + 
        a^49*x0^11*y0^7 + a^48*x0^11*y0^5 + a^49*x0^11*y0^4*y1 + a^47*x0^11*y0^3 + a^5*x0^11*y0^2*y1 + a^9*x0^11*y0*y1^2 + 
        a^9*x0^11*y0 + a^35*x0^11*y1 + a^9*x0^10*x1^2*y0^2 + a^8*x0^10*x1^2 + a^37*x0^10*x1*y0^5 + a^5*x0^10*x1*y0^3 + 
        a^9*x0^10*x1*y0^2*y1 + a^44*x0^10*x1*y0 + a^8*x0^10*x1*y1 + a^49*x0^10*y0^8 + a^35*x0^10*y0^4 + a^9*x0^10*y0^2*y1^2 + 
        a^52*x0^10*y0^2 + a^4*x0^10*y0*y1 + a^8*x0^10*y1^2 + a^16*x0^10 + a^12*x0^9*x1^3 + a^6*x0^9*x1^2*y0^3 + a^48*x0^9*x1^2*y0 + 
        a^49*x0^9*x1^2*y1 + a^9*x0^9*x1*y0^6 + a^45*x0^9*x1*y0^4 + a^46*x0^9*x1*y0^3*y1 + a^7*x0^9*x1*y0^2 + a^48*x0^9*x1*y0*y1 + 
        a^9*x0^9*x1*y1^2 + a^9*x0^9*x1 + a^49*x0^9*x2 + x0^9*y0^9 + a^8*x0^9*y0^7 + a^8*x0^9*y0^4*y1 + a^73*x0^9*y0^3 + 
        a^47*x0^9*y0^2*y1 + a^48*x0^9*y0*y1^2 + a^35*x0^9*y0 + a^49*x0^9*y1 + a^9*x0^9*y2 + a^9*x0^8*x1^2*y0^4 + a^8*x0^8*x1^2*y0^2 +
        a^7*x0^8*x1^2 + a^49*x0^8*x1*y0^7 + a^5*x0^8*x1*y0^5 + a^49*x0^8*x1*y0^4*y1 + a^35*x0^8*x1*y0^3 + a^8*x0^8*x1*y0^2*y1 + 
        a^13*x0^8*x1*y0 + a^7*x0^8*x1*y1 + a^48*x0^8*y0^8 + a^28*x0^8*y0^4 + a^8*x0^8*y0^2*y1^2 + a^16*x0^8*y0^2 + a^53*x0^8*y0*y1 + 
        a^7*x0^8*y1^2 + a^50*x0^8 + a^9*x0^7*x1^2*y0^5 + a^8*x0^7*x1^2*y0^3 + a^7*x0^7*x1^2*y0 + a^49*x0^7*x1*y0^5*y1 + 
        a^8*x0^7*x1*y0^4 + a^48*x0^7*x1*y0^3*y1 + a^13*x0^7*x1*y0^2 + a^7*x0^7*x1*y0*y1 + a^67*x0^7*x1 + a^47*x0^7*y0^7 + 
        a^46*x0^7*y0^5 + a^47*x0^7*y0^4*y1 + a^45*x0^7*y0^3 + a^53*x0^7*y0^2*y1 + a^7*x0^7*y0*y1^2 + a^7*x0^7*y0 + a^27*x0^7*y1 + 
        a^49*x0^6*x1^3*y0^3 + a^6*x0^6*x1^2*y0^6 + a^48*x0^6*x1^2*y0^4 + a^49*x0^6*x1^2*y0^3*y1 + a^53*x0^6*x1^2 + a^8*x0^6*x1*y0^7 +
        a^47*x0^6*x1*y0^5 + a^8*x0^6*x1*y0^4*y1 + a^49*x0^6*x1*y0^3*y1^2 + a^50*x0^6*x1*y0^3 + a^27*x0^6*x1*y0 + a^53*x0^6*x1*y1 + 
        a^29*x0^6*y0^4 + a^6*x0^6*y0^3*y1 + a^51*x0^6*y0^2 + a^67*x0^6*y0*y1 + a^53*x0^6*y1^2 + a*x0^6 + a^9*x0^5*x1^2*y0^7 + 
        a^8*x0^5*x1^2*y0^5 + a^7*x0^5*x1^2*y0^3 + a^48*x0^5*x1*y0^5*y1 + a^53*x0^5*x1*y0^4 + a^47*x0^5*x1*y0^3*y1 + a^67*x0^5*x1*y0^2
        + a^11*x0^5*x1 + a^45*x0^5*y0^5 + a^44*x0^5*y0^3 + a^27*x0^5*y0^2*y1 + a^17*x0^5*y0 + a^51*x0^5*y1 + a^9*x0^4*x1^2*y0^8 + 
        a^8*x0^4*x1^2*y0^6 + a^7*x0^4*x1^2*y0^4 + a^47*x0^4*x1*y0^7 + a^53*x0^4*x1*y0^5 + a^47*x0^4*x1*y0^4*y1 + a^27*x0^4*x1*y0^3 + 
        a^11*x0^4*x1*y0 + a^26*x0^4*y0^4 + a^17*x0^4*y0^2 + a^51*x0^4*y0*y1 + a^45*x0^4 + a^9*x0^3*x1^3*y0^6 + a^48*x0^3*x1^2*y0^7 + 
        a^49*x0^3*x1^2*y0^6*y1 + a^53*x0^3*x1^2*y0^3 + a^6*x0^3*x1*y0^6 + a^67*x0^3*x1*y0^4 + a^13*x0^3*x1*y0^3*y1 + a^17*x0^3*x1 + 
        a^10*x0^3*y0^3 + a^5*x0^3*y0 + a^57*x0^3*y1 + a^8*x0^2*x1^2*y0^8 + a^7*x0^2*x1^2*y0^6 + a^27*x0^2*x1*y0^5 + a^51*x0^2*x1*y0^3
        + a^45*x0^2*y0^2 + a^7*x0^2 + a^7*x0*x1^2*y0^7 + a^51*x0*x1*y0^4 + a^7*x0*y0 + a^9*x1^6 + a^9*x1^3*y1^3 + a^53*x1^2*y0^6 + 
        a^57*x1*y0^3 + a^9*x2*y0^9
]
```

One can also use the function `GT` to evaluate polynomials in two
variables using the optional argument `vvars`.  For instance, if given
Witt vectors `v` and `w` we want to compute, say, `v^2*w + 2*v*w+
w^3`, one can do:
```
> P<x>:=PolynomialRing(F);
> v:=[ a^29, x + a^11, a^68 ];
> w:=[ a^25, a^48, a^9*x^2 ]; 


> vone := [ F!1, F!0, F!0 ];
> vtwo := IntToWitt(2,p,2); vtwo;
[ 2, 1, 0 ]

> pol := [ [* vone, 2, 1 *], [* vtwo, 1, 1 *], [* vone, 0, 3 *] ];
> GT(pol : pols:=epols, vvars:=v cat w);
[
    a^34,
    a^17*x + a^67,
    a^65*x^6 + a^30*x^3 + a^61*x^2 + a^16*x + 1
]
```

Note that we could also have used `Pol_GT_Form` to produce the `pol`
above:
```
> PP<X,Y>:=PolynomialRing(Integers(),2);
> zpol := X^2*Y + 2*X*Y + Y^3;
> Pol_GT_Form(zpol,p,2);
[ [*
    [ 1, 0, 0 ],
    2, 1
*], [*
    [ 2, 1, 0 ],
    1, 1
*], [*
    [ 1, 0, 0 ],
    0, 3
*] ]
```
