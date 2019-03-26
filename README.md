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



[comp]: https://mathscinet.ams.org/mathscinet/search/publdoc.html?arg3=&co4=AND&co5=AND&co6=AND&co7=AND&dr=all&pg4=AUCN&pg5=TI&pg6=PC&pg7=TI&pg8=ET&review_format=html&s4=finotti&s5=&s6=&s7=&s8=All&sort=Newest&vfpref=html&yearRangeFirst=&yearRangeSecond=&yrop=eq&r=2&mx-pid=3248165

[canlift]: https://mathscinet.ams.org/mathscinet/search/publdoc.html?arg3=&co4=AND&co5=AND&co6=AND&co7=AND&dr=all&pg4=AUCN&pg5=TI&pg6=PC&pg7=TI&pg8=ET&review_format=html&s4=finotti&s5=&s6=&s7=&s8=All&sort=Newest&vfpref=html&yearRangeFirst=&yearRangeSecond=&yrop=eq&r=12&mx-pid=1924093

[minlift]: https://mathscinet.ams.org/mathscinet/search/publdoc.html?arg3=&co4=AND&co5=AND&co6=AND&co7=AND&dr=all&pg4=AUCN&pg5=TI&pg6=PC&pg7=TI&pg8=ET&review_format=html&s4=finotti&s5=&s6=&s7=&s8=All&sort=Newest&vfpref=html&yearRangeFirst=&yearRangeSecond=&yrop=eq&r=11&mx-pid=2044910

[jinv]: https://mathscinet.ams.org/mathscinet/search/publdoc.html?arg3=&co4=AND&co5=AND&co6=AND&co7=AND&dr=all&pg4=AUCN&pg5=TI&pg6=PC&pg7=TI&pg8=ET&review_format=html&s4=finotti&s5=&s6=&s7=&s8=All&sort=Newest&vfpref=html&yearRangeFirst=&yearRangeSecond=&yrop=eq&r=3&mx-pid=3127898



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

* `WittSum`: sums two Witt vectors.
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
