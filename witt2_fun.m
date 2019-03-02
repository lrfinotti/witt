// Compute sum of witt vectors of lenth 3 using
// the fcts from sec_coord.m

load 'sec_coord.m';

// sum two vectors using recursion
WittSum2 := function(v1,v2)

    p:=Characteristic(Parent(v1[1]));
    
    c1:= v1[1] + v2[1];
    c2:= v1[2] + v2[2] + fSpcFct(v1[1],v2[1]);
    c3:= v1[3] + v2[3] + fSpcFct(v1[2],v2[2]) +
        fSpcFct(v1[2]+v2[2],fSpcFct(v1[1],v2[1])) +
        fSpc4Fct(v1[1],v2[1]);

    return [c1,c2,c3];
    
end function;


// multiply two vectors using recursion
WittProd2 := function(v1,v2)

    p:=Characteristic(Parent(v1[1]));
    
    c1:= v1[1]*v2[1];
    c2:= v1[2]*v2[1]^p + v2[2]*v1[1]^p;
    c3:= v1[3]*v2[1]^(p^2) + v1[2]^p*v2[2]^p + v2[3]*v1[1]^(p^2)
        + fSpcFct(v2[1]^p*v1[2],v1[1]^p*v2[2]);

    return [c1,c2,c3];
    
end function;


/* Testing:

load 'witt2_fun.m';
load 'witt_gen.m';

p:=13;PR<x0,x1,x2,y0,y1,y2>:=PolynomialRing(GF(p),6);
v:=[x0,x1,x2];w:=[y0,y1,y2];

time vecS:=WittSumPol(2,p);
time a:=WittSum(v,w : vecS:=vecS) ;

time b:=WittSum2(v,w); 

a eq b;


time vecP:=WittProdPol(2,p);
time a:=WittProd(v,w : vecP:=vecP);

time b:=WittProd2(v,w); 

a eq b;


*/


// compute the polynomials for the sum using the recursion
WittSum2Pol := function(p)

    PR<X0,X1,X2,Y0,Y1,Y2>:=PolynomialRing(GF(p),6);
 
    c1:= X0 + Y0;
    c2:= X1 + Y1 + fSpcFct(X0,Y0);
    c3:= X2 + Y2 + fSpcFct(X1,Y1) +
        fSpcFct(X1+Y1,fSpcFct(X0,Y0)) +
        fSpc4Fct(X0,Y0);

    return [c1,c2,c3];
  
end function;


// compute the polynomials for the product using the recursion
WittProd2Pol := function(p)
    
    PR<X0,X1,X2,Y0,Y1,Y2>:=PolynomialRing(GF(p),6);
 
    
    c1:= X0*Y0;
    c2:= X1*Y0^p + Y1*X0^p;
    c3:= X2*Y0^(p^2) + X1^p*Y1^p + Y2*X0^(p^2)
        + fSpcFct(X0^p*Y1,Y0^p*X1);

    return [c1,c2,c3];
    
end function;



// these will compute in Z/p^n

Divp := function(f)

    P:=Parent(f);
    p:=Factorization(Characteristic(P))[1][1];

    R:=PolynomialRing(Integers(),6);

    return P!(R!f div p);

end function;


Divp2 := function(f)

    P:=Parent(f);
    p:=Factorization(Characteristic(P))[1][1];

    R:=PolynomialRing(Integers(),6);

    return P!(R!f div p^2);

end function;



WittSum2Pol2 := function(p)

    PR<X0,X1,X2,Y0,Y1,Y2>:=PolynomialRing(GF(p),6);

    c1:= X0 + Y0;

    Rp2:=ResidueClassRing(p^2);
    RR<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(Rp2,6);

    c2:=X1+Y1 + PR!(Divp(XX0^p+YY0^p - (XX0+YY0)^p));
    
    Rp3:=ResidueClassRing(p^3);
    RR<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(Rp3,6);

    c3:=X2 + Y2 +
        PR!(Divp2(p*(XX1^p+YY1^p-(XX1+YY1+Divp(XX0^p+YY0^p-(XX0+YY0)^p))^p)
        + (XX0^(p^2)+YY0^(p^2)-(XX0+YY0)^(p^2)))); 
        

    return [c1,c2,c3];
  
end function;


WittProd2Pol2 := function(p)

    PR<X0,X1,X2,Y0,Y1,Y2>:=PolynomialRing(GF(p),6);

    c1:= X0 * Y0;

    c2:=X0^p*Y1+Y0^p*X1;


    Rp2:=ResidueClassRing(p^2);
    RR<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(Rp2,6);

    c3:=X0^(p^2)*Y2 + X1^p*Y1^p + Y0^(p^2)*X2
        + PR!(Divp(XX0^(p^2)*YY1^p+YY0^(p^2)*XX1^p - (XX0^p*YY1+YY0^p*XX1)^p));

    return [c1,c2,c3];

end function;


// now compute in Z

WittSum2Pol3 := function(p)

    PR<X0,X1,X2,Y0,Y1,Y2>:=PolynomialRing(GF(p),6);

    c1:= X0 + Y0;

    RR<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(Integers(),6);

    c2:=X1+Y1 + PR!((XX0^p+YY0^p - (XX0+YY0)^p) div p);
    
    c3:=X2 + Y2 +
        PR!((p*(XX1^p+YY1^p-(XX1+YY1+(XX0^p+YY0^p-(XX0+YY0)^p) div p)^p)
        + (XX0^(p^2)+YY0^(p^2)-(XX0+YY0)^(p^2))) div p^2); 
        

    return [c1,c2,c3];
  
end function;


WittProd2Pol3 := function(p)

    PR<X0,X1,X2,Y0,Y1,Y2>:=PolynomialRing(GF(p),6);

    c1:= X0 * Y0;

    c2:=X0^p*Y1+Y0^p*X1;


    RR<XX0,XX1,XX2,YY0,YY1,YY2>:=PolynomialRing(Integers(),6);

    c3:=X0^(p^2)*Y2 + X1^p*Y1^p + Y0^(p^2)*X2
        + PR!((XX0^(p^2)*YY1^p+YY0^(p^2)*XX1^p - (XX0^p*YY1+YY0^p*XX1)^p) div p);

    return [c1,c2,c3];

end function;



// evaluate the vector of polynomials in aux
EvalVec := function(vec,aux)
    return [ Evaluate(term,aux) : term in vec ];
end function;


// compute sum and product by precomputing the formulas
// for sum and product and then evaluating
WittSum22 := function(v1,v2 : vecS:=[])

    P:=Parent(v1[1]);
    p:=Characteristic(P);
    if vecS eq [] then
	vecS:=WittSum2Pol(p);
    end if;

    aux:=< (v1 cat v2)[i] : i in [1..6] >;
    return EvalVec(vecS,aux);

end function;



WittProd22 := function(v1,v2 : vecP:=[])

    P:=Parent(v1[1]);
    p:=Characteristic(P);
    if vecP eq [] then
	vecP:=WittProd2Pol(p);
    end if;
    
    aux:=< (v1 cat v2)[i] : i in [1..6] >;
    return EvalVec(vecP,aux);

end function;
