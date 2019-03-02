

load 'sec_coord.m';


p:=5; // choose p

F1<a0,b0>:=RationalFunctionField(GF(p),2);
PR1<x0,c0,c1,a1,b1>:=PolynomialRing(F1,5);


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

// 3 by len matrix -- 3 variables: (c0,c1,a1) or (c0,c1,b1)
M1:=Matrix(F1,4,len,
     [[F1!(Coefficient(vcoef[j],i+1,1)) : j in [1..len]] : i in [1..4]]);

// We now solve the linear system.
VS1:=VectorSpace(F1,len);

// term independent of the variables
// (note: the 4th variable is x0, which does not show up in the coef.)
v1:=VS1![-Evaluate(vcoef[i],<0,0,0,0,0>) : i in [1..len]];

// find ONE solution!
vsol, ker:=Solution(M1,v1);

kerB:=Basis(ker);

print kerB;

