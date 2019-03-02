// Computes canonical lift.


// computes quotient and remainder of a div. of polynomials at the same time

// It seems I need the function below since I am dealing with polynomials
// in several variables, but I only want the division wrt one.
MyQuotrem := function(pol1,pol2,x) // assumes that pol2 is monic!!!
     quo:=0;
     rem:=pol1;
     n:=Degree(pol1,x);
     m:=Degree(pol2,x);
     for j in [0..(n-m)] do
       l1:=Coefficient(rem,x,n-j);
       ftmp:=l1*x^(n-m-j);
       quo+:=ftmp;
       rem-:=ftmp*pol2;
     end for;
     return quo, rem;
 end function;


 // Does the long division by converting to one variable
 // IT IS SLOWER THAN THE ABOVE!!!!!!
     MyQuotrem2 := function(pol1,pol2) // assumes that pol2 is monic!!!
     P:=Parent(pol1);  // must be the same as pol2
     n:=Rank(P);
     R:=BaseRing(P);
     RR:=PolynomialRing(R,n-1);
     PP:=PolynomialRing(RR);
     evalvec:=[ PP.1 ] cat [ RR.i : i in [1..(n-1)] ];
     ppol1:=Evaluate(pol1,evalvec);
     ppol2:=Evaluate(pol2,evalvec);
     quo, rem := Quotrem(ppol1,ppol2);
     rquo:=0;
     evalvec:=[ P.i : i in [2..n] ];
     for i in Support(quo) do
         rquo+:=Evaluate(Coefficient(quo,i),evalvec)*P.1^i;
     end for;
     delete quo;
     rrem:=0;
     for i in Support(rem) do
         rrem+:=Evaluate(Coefficient(rem,i),evalvec)*P.1^i;
     end for;
     delete rem;
     return rquo, rrem;
 end function;



//  ------- fct ---------

lift := function(a0,b0 : tm:=false, test:=false, ncoord:=2, minimal:=false)
// tm = time it?
// ncoord = no. of coordinates to compute (1 -> x1,P1; 2 -> x2,P2)
// test = test result
// minimal = compute minimal instead of canonical lift?
    
if tm then
    ttime:=Cputime();
    ptime:=Cputime();
end if;

F1:=Parent(a0);

error if not F1 cmpeq Parent(b0),
     "ERROR: ",a0, " and ", b0, " are not in the same field";


p:=Characteristic(F1); // characteristic

error if p eq 0,
     "ERROR: Char. of the field is zero!";

error if 4*a0^3+27*b0^2 eq 0,
     "ERROR: singular curve";





if tm then print "Comp. 1st coord.";
end if;

Z:=RingOfIntegers();


SS1<XX0,AA0,BB0>:=PolynomialRing(ResidueClassRing(p^2),3);
S1<XX0,AA0,BB0>:=PolynomialRing(Z,3);

// compute the ``special function'' for the RHS of the equations
// IMPROVE?
AuxPol1:= SS1! (XX0^(3*p)+AA0^p*XX0^p+BB0^p) - (SS1!(XX0^3+AA0*XX0+BB0))^p;
AuxPol1:=S1!(AuxPol1);
AuxPol1 div:= p;


// we will assume c1=0;
PR1<x0,c0,a1,b1>:=PolynomialRing(F1,4);

// RHS of the equation
f:=PR1!(x0^3+a0*x0+b0);

// y0^(p-1)
y0pm1:=f^(Z!((p-1)/2));

// y0^(p+1)
y0pp1:=y0pm1*f;

// Hasse Invariant
HI:=F1!(Coefficient(y0pm1,x0,p-1));

error if HI eq 0,
  "ERROR: hasse inv. = 0!";

// dx1/dx0
tmpx1:=PR1!((HI^(-1))*y0pm1-x0^(p-1));

x1:=Integral(tmpx1,x0)+c0;

delete tmpx1;

tmp:=Evaluate(AuxPol1,<x0,a0,b0>);
        
// RHS - 2nd coord
rhs1:=PR1!(3*x0^(2*p)*x1+a0^p*x1+a1*x0^p+b1) + tmp;

delete tmp;

// Find canditate to P1 (y1=y0*P1)
// rem1 must be zero!
P1, rem1 := MyQuotrem(((F1!1/2)*rhs1),y0pp1,x0);

// separate the coefficients of rem1 -- they must be zero
vcoef:=Coefficients(rem1,x0);

len:=Degree(rem1,x0)+1;

// 3 by len matrix -- 3 variables: (c0,a1,b1)
M1:=Matrix(F1,3,len,
     [[F1!(Coefficient(vcoef[j],i+1,1)) : j in [1..len]] : i in [1..3]]);

// We now solve the linear system.
VS1:=VectorSpace(F1,len);

// term independent of the variables
// (note: the 4th variable is x0, which does not show up in the coef.)
v1:=VS1![-Evaluate(vcoef[i],<0,0,0,0>) : i in [1..len]];

// find ONE solution!
vsol:=Solution(M1,v1);

// IMPROVE?
// We could find all...  Check Nullspace
// file:///usr/local/magma/doc/html/text616.htm#5579


// Now, comput the solutions

// x1 only involves c1
x1:=Evaluate(x1,<x0,vsol[1],0,0>);


P1:=Evaluate(P1,<x0,vsol[1],vsol[2],vsol[3]>);


// ********** test ***************

if test then 
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
    if Evaluate(rem1,<x0,vsol[1],vsol[2],vsol[3]>) eq 0
	then print "test1 is OK";
    else print "test 1 in bad!"; print "rem1 = ", Evaluate(rem1,<x0,vsol[1],vsol[2],vsol[3]>);
    end if;
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
end if;


delete M1, rem1, VS1, vcoef, len, v1, S1, AuxPol1;


//  Now to red. mod p^3

if tm then print "Done with 1st coord.";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
end if;

// return a1, b1, x1, P1
// if ncoord=1, stops here!
if ncoord eq 1 then
    return vsol[2], vsol[3], x1, P1;
end if;

//////////////////////////////////////////////////
////////////// SECOND COORDINATE /////////////////
//////////////////////////////////////////////////

if tm then
    ptime:=Cputime();
    print "Doing convertions and y0psm1, y0psp1";
end if;

// introduce the new variables: a2, b2 and the coeff. of x0^(pk) in x2
PR2<x0,a2,b2>:=PolynomialRing(F1,Z!((3*p+7)/2));

// We need to convert the previous results to this new ring.
Convert := procedure(~P,M) // convert P to new magma M
    P:= &+[F1!(Coefficient(P,1,i))*(M.1)^i : i in [0..Degree(P,1)]];
end procedure;

// make the conversions
Convert(~P1,PR2);
Convert(~x1,PR2);
Convert(~y0pm1,PR2);
Convert(~y0pp1,PR2);
Convert(~f,PR2);
va1:=vsol[2]; // a1
vb1:=vsol[3]; // b1

delete vsol;

// compute new powers of y0:
// y0^(p^2-1) and y0^(p^2+1)
y0psm1:=y0pm1^p*y0pm1;
y0psp1:=f*y0psm1;

if tm then
    print "Done conv. and y0ps...";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
    print "Comp. x2";
end if;

// dx2/dx0
tmpx2:=PR2!((HI^(-1-p))*y0psm1-x0^(p^2-1)-x1^(p-1)*(HI^(-1)*y0pm1-x0^(p-1)));

// we are NOT adding the term in x0^(p^2)!!
x2:=Integral(tmpx2,x0)+ &+[PR2.(i+4)*x0^(i*p) : i in [0..Z!((3*p-1)/2)] | i ne p];
//                                                                      ^^^^^^^^^

delete tmpx2;

// add extra terms (known)
if not minimal then
    tmp:= F1!(3/4)*x1^2;
    x2+:= &+[(Coefficient(tmp,x0,Z!(i/p+p)))^p*x0^i : i in [Z!((3*p^2+p)/2) .. 2*p^2-p by p]];
    delete tmp;
end if;

if tm then
    print "Done comp. x2";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
    print "comp witt prod and sum";
end if;

SS2<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(GF(p),6);
S2p2<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(ResidueClassRing(p^2),6);
S2p3<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(ResidueClassRing(p^3),6);
S2<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(Z,6);

// computation of the RHS
// IMPROVE?
// We already have the 2nd coord.  Can't we use it?


// find the formulas for Witt Sum
WS0:=(XX0+YY0); // final

WS1:=S2!(-(S2p2!(XX0+YY0))^p);

WS1+:=XX0^p+YY0^p;

WS1 div:= p;

WS1 := S2!(SS2!(WS1));

WS1+:=XX1+YY1; // final


WS2:=S2p3!(XX0^(p^2)+YY0^(p^2))-(S2p3!(WS0))^(p^2);

WS2+:=p*(S2p3!(XX1^p+YY1^p) -S2p3!((S2p2!(WS1))^p));

WS2:=S2!(WS2);

WS2 div:=p^2;

WS2+:=XX2+YY2;  // final


WS0:=SS2!(WS0); // reduce mod p
WS1:=SS2!(WS1); // reduce mod p
WS2:=SS2!(WS2); // reduce mod p



// find the formulas for Witt Product
WP0:=XX0*YY0; // final

WP1:=XX0^p*YY1+YY0^p*XX1; // final

WP2:=S2!(-(S2p2!(WP1))^p);

WP2 +:= XX0^(p^2)*YY1^p+YY0^(p^2)*XX1^p;

WP2 div:=p;

WP2 +:= XX0^(p^2)*YY2+XX1^p*YY1^p+YY0^(p^2)*XX2; // final

WP0:=SS2!(WP0); // reduce mod p
WP1:=SS2!(WP1); // reduce mod p
WP2:=SS2!(WP2); // reduce mod p



WittSum := function(v,w)
     res:=[v[1]+w[1]];
     res:=Append(res,Evaluate(WS1,<v[1],v[2],0,w[1],w[2],0>));
     res:=Append(res,Evaluate(WS2,<v[1],v[2],v[3],w[1],w[2],w[3]>));
     return res;
     end function;


WittProd := function(v,w)
     res:=[v[1]*w[1]];
     res:=Append(res,Evaluate(WP1,<v[1],v[2],0,w[1],w[2],0>));
     res:=Append(res,Evaluate(WP2,<v[1],v[2],v[3],w[1],w[2],w[3]>));
     return res;
     end function;


if tm then
    print "Done comp. witt sum/prod";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
    print "computing rhs2";
end if;

vecX:=[x0,x1,x2];
vecA:=[a0,va1,a2];
vecB:=[b0,vb1,b2];

// compute RHS
rhs2:=WittProd(vecX,vecX); // x^2
rhs2:=WittProd(vecX,rhs2); // x^3

vt:=WittProd(vecA,vecX); // a*x

rhs2:=WittSum(rhs2,vt); // x^3+a*x

delete vt;

rhs2:=WittSum(rhs2,vecB); // x^3+a*x+b

delete vecX, vecA, vecB, WittSum, WittProd, WS0, WS1, WS2,
    WP0, WP1, WP2, SS2, S2p2, S2p3, S2;

// we only care abou the 3rd coord.!
rhs2:=rhs2[3];

// now subtract parts of the LHS that are known!
// what is left is 2*y0^(p^2+1)*P2!
rhs2-:=(f^p*P1^(2*p)+F1!((2-2^p)/p)*y0pp1^p*P1^p); // NOT rhs2...


if tm then
    print "Done comp. rhs";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
    print "computing P2 and rem2";
end if;

// rem2 has to be zero
P2, rem2 := MyQuotrem(F1!(1/2)*rhs2,y0psp1,x0);
//P2, rem2 := MyQuotrem2(F1!(1/2)*rhs2,y0psp1);

delete rhs2;

vcoef:=Coefficients(rem2,x0);

len:=Degree(rem2,x0)+1;


if tm then 
    print "Done comp. P2 and rem2";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
    print "computing M2 and v2";
end if;

r:=Z!(2+(3*p-1)/2);  // number of variables = rank!
ncol:=1;

// will work with transpose!!
j:=len;

gen:=[F1!(Coefficient(vcoef[j],i+1,1)) : i in [1..(r+1)] | i ne p+3 ];
list:=[j];
M2:=Matrix(F1,ncol,r,gen);
j-:=1;


// IIRC, I add the rows just as long as they increase the rank,
// and until the rank is maximal, to try to not deal with a huge matrix...
// Maybe MAGMA would be faster to just deal with the matrix itself????
while (j ge 1) and (ncol lt r) do
    tgen:=gen cat [F1!(Coefficient(vcoef[j],i+1,1)) : i in [1..(r+1)] | i ne p+3 ];
    tM2:=Matrix(F1,ncol+1,r,tgen);
    if Rank(tM2) eq ncol+1
	then M2:=tM2; gen:=tgen; ncol+:=1; list:=Append(list,j); j-:=1; 
    else j-:=1;
    end if;
end while;

// transpose back
M2:=Transpose(M2);

if tm then
    print "%%%%%%%%%%%%%%";
    print "M2 is ", Nrows(M2) , "x" , Ncols(M2);
    print "Number of var is ", r;
    print "%%%%%%%%%%%%%%";
end if;

// Now solve the system
VS2:=VectorSpace(F1,r);

tzero:=<0 : i in [1..(3*p+7) div 2]>;

// independent terms
v2:=VS2![-Evaluate(vcoef[i],tzero) : i in list];

if tm then
    print "done comp. M2 and v2";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
    print "Solving the system";
end if;

// find ONE solution
vsol:=Solution(M2,v2);

// IMPROVE?
// Find all solutions??

if tm then
    print "done solving the system";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
end if;

delete M2, v2;

// no. of varibales - 2 (for x0 and the coeff of x0^(p^2))
r:=Z!(2+(3*p-1)/2);

// we now create a vector to evaluate the solutions
// it must have <x0,a2,b2, ..coefs of x2.. >
tp:=<x0>;

for i in [1..p+2] do
    tp:=Append(tp,vsol[i]);
end for;

// here we are making the coefficient of x0^(p^2)
// equal to zero, as we've assumed earlier
tp:=Append(tp,0);

for i in [p+3..r] do
    tp:=Append(tp,vsol[i]);
end for;


// ********** test ***************
if test then
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
    if Evaluate(rem2,tp) eq 0
	then print "test2 is OK";
    else print "test 2 in bad!"; print "rem2 = ",  Evaluate(rem2,tp);
    end if;
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
end if;

delete rem2;

x2:=Evaluate(x2,tp);

P2:=Evaluate(P2,tp);

delete tp, r;

if tm then
    print "total time = ", Cputime(ttime);
end if;

return va1, vb1, x1, P1, vsol[1], vsol[2], x2, P2;
end function;
