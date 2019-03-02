

load 'sec_coord.m';


p:=13; // choose p

// Will make c1=0 here!
c1:=0;

F1<a0,b0>:=RationalFunctionField(GF(p),2);
PR1<x0,c0,a1,b1>:=PolynomialRing(F1,4);



Z:=RingOfIntegers();

    
f:=PR1!(x0^3+a0*x0+b0);
sf:=SpcFct(f);

// y0^(p-1)
y0pm1:=f^(((p-1) div 2));

// y0^(p+1)
y0pp1:=y0pm1*f;

// Hasse Invariant
HI:=F1!(Coefficient(y0pm1,x0,p-1));


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

// 3 by len matrix -- 3 variables: (c0,a1,b1)
M1:=Matrix(F1,3,len,
     [[F1!(Coefficient(vcoef[j],i+1,1)) : j in [1..len]] : i in [1..3]]);

// We now solve the linear system.
VS1:=VectorSpace(F1,len);

// term independent of the variables
// (note: the 4th variable is x0, which does not show up in the coef.)
v1:=VS1![-Evaluate(vcoef[i],<0,0,0,0>) : i in [1..len]];

// find ONE solution!
vsol, ker:=Solution(M1,v1);

kerB:=Basis(ker);

print kerB;

x1:=Evaluate(x1,<x0,vsol[1],0,0>);


P1:=Evaluate(P1,<x0,vsol[1],vsol[2],vsol[3]>);



delete M1, rem1, VS1, vcoef, len, v1;




// //////////////////////////////////////////////////////////

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
Convert(~sf,PR2);

va1:=vsol[2];
vb1:=vsol[3];


delete vsol;

// compute new powers of y0:
// y0^(p^2-1) and y0^(p^2+1)
y0psm1:=y0pm1^p*y0pm1;
y0psp1:=f*y0psm1;



// dx2/dx0
tmpx2:=PR2!((HI^(-1-p))*y0psm1-x0^(p^2-1)-x1^(p-1)*(HI^(-1)*y0pm1-x0^(p-1)));

// we ARE adding the term in x0^(p^2)!!
x2:=Integral(tmpx2,x0)+ &+[PR2.(i+4)*x0^(i*p) : i in [0..Z!((3*p-1)/2)] ];

tmp:= F1!(3/4)*x1^2;
x2+:= &+[(Coefficient(tmp,x0,Z!(i/p+p)))^p*x0^i : i in [Z!((3*p^2+p)/2) .. 2*p^2-p by p]];

delete tmpx2, tmp;


ev:=[ 3*x0^(2*p)*x1 , a0^p*x1, va1*x0^p, vb1 ];


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


// rem2 has to be zero
rhs2*:=F1!(1/2);
P2:=0;
deg1:=Degree(rhs2,x0);
deg2:=(3*((p^2+1) div 2));

RHS2:=rhs2;

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


r:=Z!((3*p+5)/2);  // number of variables (except x0) = rank!


// Use echelon form to find solution
tzero:=<0 : i in [1..(3*p+7) div 2]>;


// r by len matrix -- 3 variables: (c0,a1,b1)
M2:=Matrix(F1,r,len,
     [[F1!(Coefficient(vcoef[j],i+1,1)) : j in [1..len]] : i in [1..r]]);

// We now solve the linear system.
VS2:=VectorSpace(F1,len);

// term independent of the variables
// (note: the 4th variable is x0, which does not show up in the coef.)
v2:=VS2![-Evaluate(vcoef[i],tzero) : i in [1..len]];


// After it is in echelon form
vsol, ker2:=Solution(M2,v2);

print Basis(ker2);



delete M2, v2;

// no. of varibales - 1 (for x0)
r:=Z!((3*p+5)/2);

// we now create a vector to evaluate the solutions
// it must have <x0,a2 or b2, ..coefs of x2.. >
tp:=<x0>;

for i in [1..r] do
    tp:=Append(tp,vsol[i]);
end for;


// ********** test ***************

    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
    if Evaluate(rem2,tp) eq 0
	then print "test2 is OK";
    else print "test 2 in bad!"; print "rem2 = ",  Evaluate(rem2,tp);
    end if;
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";


delete rem2;

x2:=Evaluate(x2,tp);

P2:=Evaluate(P2,tp);

delete tp, r;

/* print x2; */

/* print ""; */

/* print "a2 = "; vsol[1]; */


/* print ""; */

/* print "b2 = "; vsol[2]; */


