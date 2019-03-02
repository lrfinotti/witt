// NEW VERSION!!!!
// Uses the new functions

load 'sec_coord.m';

// Computes canonical lift.


// computes quotient and remainder of a div. of polynomials at the same time


/*

Removed it!!!

I do it "by hand" below, since it allows me to change the RHS (the dividend)
as I do it, freeing some memory.


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

*/

    

//  ------- fct ---------

lift := function(a0,b0 : tm:=false, test:=false, ncoord:=2, minimal:=false, choice:=1)
// tm = time it?
// ncoord = no. of coordinates to compute (1 -> x1,P1; 2 -> x2,P2)
// test = test result
// minimal = compute minimal instead of canonical lift?
// choice = make a_1=a_2=0 (choice:=2) or b_1=b_2=0 (choice:=1)
    
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





if tm then
     print "Comp. 1st coord.";
end if;

Z:=RingOfIntegers();


if choice eq 1 then
    // we will assume b1=0;
    PR1<x0,c0,c1,a1>:=PolynomialRing(F1,4);
    b1:=0;
else
    // we will assume a1=0;
    PR1<x0,c0,c1,b1>:=PolynomialRing(F1,4);
    a1:=0;
end if;
    
f:=PR1!(x0^3+a0*x0+b0);
sf:=SpcFct(f);

// y0^(p-1)
y0pm1:=f^(((p-1) div 2));

// y0^(p+1)
y0pp1:=y0pm1*f;

// Hasse Invariant
HI:=F1!(Coefficient(y0pm1,x0,p-1));

error if HI eq 0,
  "ERROR: hasse inv. = 0!";

// dx1/dx0
tmpx1:=PR1!((HI^(-1))*y0pm1-x0^(p-1));

x1:=Integral(tmpx1,x0)+c0+c1*x0^p;

delete tmpx1;

// RHS - 2nd coord
rhs1:=PR1!((3*x0^(2*p)+a0^p)*x1+a1*x0^p+b1) + sf;


// Find canditate to P1 (y1=y0*P1)
// rem1 must be zero!
rhs1*:=F1!(1/2);
P1:=0;
deg1:=Degree(rhs1,x0);
deg2:=(3*((p+1) div 2));

while deg1 ge deg2 do
    lterm:=Coefficient(rhs1,x0,deg1)*x0^(deg1-deg2);
    rhs1-:=lterm*y0pp1;
    P1+:=lterm;
    deg1:=Degree(rhs1,x0);
end while;

rem1:=rhs1;

delete rhs1;


// separate the coefficients of rem1 -- they must be zero
vcoef:=Coefficients(rem1,x0);

len:=Degree(rem1,x0)+1;

// 3 by len matrix -- 3 variables: (c0,c1,a1) or (c0,c1,b1)
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

// x1 only involves c0 and c1, but not b1
x1:=Evaluate(x1,<x0,vsol[1],vsol[2],0>);


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


delete M1, rem1, VS1, vcoef, len, v1;


if tm then print "Done with 1st coord.";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
end if;

// return a1, b1, x1, P1
// if ncoord=1, stops here!
if ncoord eq 1 then
    if choice eq 1 then
        return vsol[3], 0, x1, P1;
    else
        return 0, vsol[3], x1, P1;
    end if;
end if;

//////////////////////////////////////////////////
////////////// SECOND COORDINATE /////////////////
//////////////////////////////////////////////////

//  Now to red. mod p^3

if tm then
    ptime:=Cputime();
    print "Doing convertions and y0psm1, y0psp1";
end if;

// introduce the new variables: a2 (and b2=0) or b2 (and a2=0)
//  and the coeff. of x0^(pk) in x2
// NOTE: need one less var. than if a2!=0 and b2!=0!!!
if choice eq 1 then
    PR2<x0,a2>:=PolynomialRing(F1,Z!((3*p+5)/2));
    b1:=0;
    b2:=0;
else
    PR2<x0,b2>:=PolynomialRing(F1,Z!((3*p+5)/2));
    a1:=0;
    a2:=0;
end if;
  
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
Convert(~sf,PR2);

if choice eq 1 then
    va1:=vsol[3]; // a1
    vb1:=0; // b1
else
    vb1:=vsol[3]; // b1
    va1:=0; // a1
end if;

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

// we ARE adding the term in x0^(p^2)!!
x2:=Integral(tmpx2,x0)+ &+[PR2.(i+3)*x0^(i*p) : i in [0..Z!((3*p-1)/2)] ];


delete tmpx2, HI;

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
    print "comp RHS";
end if;


// find the necessary vector for the call of SpcFctv
if choice eq 1 then // b1=0
    if a0 eq 0 then
        ev:=[ 3*x0^(2*p)*x1 , va1*x0^p ];
    else
        ev:=[ 3*x0^(2*p)*x1 , a0^p*x1, va1*x0^p ];
    end if;
else // a1=0
    if a0 eq 0 then
        ev:=[ 3*x0^(2*p)*x1 , vb1 ];
    else
        ev:=[ 3*x0^(2*p)*x1 , a0^p*x1, vb1 ];
    end if;
end if;


rhs2:= (3*x0^(2*p^2)+a0^(p^2))*x2 + 3*x0^(p^2)*x1^(2*p) +
    (-NextCoord(3,p)*x0^(2*p^2)+va1^p)*x1^p + a2*x0^(p^2) + b2;

rhs2+:= SpcFctv(ev);

rhs2+:= fSpcFct( (3*x0^(2*p)+a0^p)*x1 + va1*x0^p + vb1 , sf );

rhs2+:= Spc2Fct(f);

rhs2+:= Spc3Fct(f);
    
delete sf, ev;

// now subtract parts of the LHS that are known!
// what is left is 2*y0^(p^2+1)*P2!
rhs2-:=(f^p*P1^(2*p)+F1!((2-2^p)/p)*y0pp1^p*P1^p); // NOT rhs2...


if tm then
    // -------------------------
    print "test: ", Terms(rhs2)[3];
    // -------------------------
    print "Done comp. rhs";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
    print "computing P2 and rem2";
end if;

// rem2 has to be zero
rhs2*:=F1!(1/2);
P2:=0;
deg1:=Degree(rhs2,x0);
deg2:=(3*((p^2+1) div 2));

//print "RHS := ", rhs2 , ";";

while deg1 ge deg2 do
    lterm:=Coefficient(rhs2,x0,deg1)*x0^(deg1-deg2);
    rhs2-:=lterm*y0psp1;
    P2+:=lterm;
    deg1:=Degree(rhs2,x0);
end while;

rem2:=rhs2;

delete rhs2, y0pm1, y0pp1, y0psm1, y0psp1;


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

r:=Z!((3*p+3)/2);  // number of variables (except x0) = rank!
ncol:=1;

// Use echelon form to find solution
tzero:=<0 : i in [1..(3*p+5) div 2]>;

j:=len;


M2:=Matrix(F1,ncol,r,[F1!(Coefficient(vcoef[j],i+1,1)) : i in [1..r]]);
v2:=Matrix(F1,1,1,[F1!(-Evaluate(vcoef[j],tzero))]);
j-:=1;


// IIRC, I add the rows just as long as they increase the rank,
// and until the rank is maximal, to try to not deal with a huge matrix...
// Maybe MAGMA would be faster to just deal with the matrix itself????

// Find a way to put already in echelon form so that the comp. of rank is quicker and
// so is solving the system
while (j ge 1) and (ncol lt r) do
    // -------------------------
    // print "Test: j=", j;
    // -------------------------
    tN:=Matrix(F1,1,r,[F1!(Coefficient(vcoef[j],i+1,1)) : i in [1..r]]);
    tM2,T:=EchelonForm(VerticalJoin(M2,tN));
    if Rank(tM2) eq ncol+1 then
        M2:=tM2; ncol+:=1;
        v2:=T*(VerticalJoin(v2,Matrix(F1,1,1,[F1!(-Evaluate(vcoef[j],tzero))])));
        j-:=1;
        delete T;
    else j-:=1;
    end if;
    // -------------------------
    // print "M2=", M2;
    // -------------------------
end while;

if test then
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
    if M2 eq ScalarMatrix(r,F1!1)
	then print "Echelon gives Identity";
    else print "Echelon DOESN'T give Identity!"; print "M2 = ", M2;
    end if;
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
end if;

if tm then
    print "%%%%%%%%%%%%%%";
    print "M2 is ", Nrows(M2) , "x" , Ncols(M2);
    print "Number of var is ", r;
    print "%%%%%%%%%%%%%%";
end if;

/*
if tm then
    print "done comp. M2 and v2";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
    print "Solving the system";
end if;
*/

// After it is in echelon form
vsol:=[ v2[i][1] : i in [1..(Nrows(v2))] ];

// IMPROVE?
// Find all solutions??
if #vsol ne r then
    vsol:=vsol cat [ F1!0 : i in [(#vsol + 1)..r] ];
end if;


if tm then
    print "done solving the system";
    print "Partial time = ", Cputime(ptime);
    print "Total time = ", Cputime(ttime);
    print "**************************************************";
    ptime:=Cputime();
end if;

delete M2, v2;

// no. of varibales - 1 (for x0)
r:=Z!((3*p+3)/2);

// we now create a vector to evaluate the solutions
// it must have <x0,a2 or b2, ..coefs of x2.. >
tp:=<x0>;

for i in [1..r] do
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


if choice eq 1 then
    return va1, vb1, x1, P1, vsol[1], 0, x2, P2;
else
    return va1, vb1, x1, P1, 0, vsol[1], x2, P2;
end if;


end function;


/*

Test it!!!

p:=7;F<a0,b0>:=RationalFunctionField(GF(p),2);


load 'fctWitt4.m';
a1, b1, x1, P1, a2, b2, x2, P2 := lift(a0,b0 : tm:=true, test:=true, ncoord:=2, minimal:=false, choice:=1);


load 'fctWitt3.m';
aa1, bb1, xx1, PP1, aa2, bb2, xx2, PP2 := lift(a0,b0 : tm:=true, test:=true, ncoord:=2, minimal:=false);

Px1:=Parent(x1);
Px2:=Parent(x2);
a1 eq aa1;
b1 eq bb1;
x1 eq Px1!xx1;
P1 eq Px1!PP1;
a2 eq aa2;
b2 eq bb2;
x2 eq Px2!xx2;
P2 eq Px2!PP2;

load 'fctWitt4.m';
a1, b1, x1, P1, a2, b2, x2, P2 := lift(a0,b0 : tm:=true, test:=true, ncoord:=2, minimal:=false, choice:=2);


load 'fctWitt2.m';
aa1, bb1, xx1, PP1, aa2, bb2, xx2, PP2 := lift(a0,b0 : tm:=true, test:=true, ncoord:=2, minimal:=false);

Px1:=Parent(x1);
Px2:=Parent(x2);
a1 eq aa1;
b1 eq bb1;
x1 eq Px1!xx1;
P1 eq Px1!PP1;
a2 eq aa2;
b2 eq bb2;
x2 eq Px2!xx2;
P2 eq Px2!PP2;



*/







/*

// Let's use a separate function to compute the RHS of the 3rd coord
F3rdCoord := function(v,p)
    // v is a vector which can give us some zero coeff.
    // v=[a0,a1,a2,b0,b1,b2] : a value equal to zero makes the fct
    // assume the the corresponding term is zero.
    F<A0,B0>:=RationalFunctionField(GF(p),2);
    PR<A1,A2,B1,B2,X0,X1,X2>:=PolynomialRing(F,7);

    if v[1] eq 0 then
        a0:=0;
    else
        a0:=A0;
    end if;

    if v[2] eq 0 then
        a1:=0;
    else
        a1:=A1;
    end if;

    if v[3] eq 0 then
        a2:=0;
    else
        a2:=A2;
    end if;

    if v[4] eq 0 then
        b0:=0;
    else
        b0:=B0;
    end if;

    if v[5] eq 0 then
        b1:=0;
    else
        b1:=B1;
    end if;

    if v[6] eq 0 then
        b2:=0;
    else
        b2:=B2;
    end if;

    term1:= (3*X0^(p^2)+a0^(p^2))*X2 + 3*X0^(p^2)*X1^(2*p) +
        (-NextCoord(3,p)*X0^(2*p^2)+a1^p)*X1^p + a2*X0^(p^2) + b2;

    term2:= SpcFct( (3*X0^p+a0^p)*X1 + a1*X0^p + b1 );

    term3:= fSpcFct( (3*X0^p+a0^p)*X1 + a1*X0^p + b1 , SpcFct(X0^3+a0*X0+b0) );

    term4:= Spc2Fct(X0^3+a0*X0+b0);

    term5:= Spc3Fct(X0^3+a0*X0+b0);

    return term1 + term2 + term3 + term4 + term5;

end function;

*/

