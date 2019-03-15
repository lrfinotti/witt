// PR<y0,y1,y2,x0,x1,x2,a0,a1,a2,P1,Q1,P2,Q2>:=PolynomialRing(GF(2),13);
// QR:=quo< PR | y0^2+x0*y0+x0^3+a0 >;

// auxiliar function to evaluate vectors
EvalVec := function(vec,aux)
//     res:=[];
//     for j:=1 to #vec do
//       res:=Append(res,Evaluate(vec[j],aux));
//     end for;
//     return res;
    return [ Evaluate(term,aux) : term in vec ];
end function;


// evaluates rational function -- not builtin???
EvalRatFct := function(f,v)
     F:=Parent(f);
     num:=Evaluate(Numerator(f),v);
     den:=Evaluate(Denominator(f),v);
     return num/den;
end function;



EvalRFVec := function(vec,aux)
//     res:=[];
//     for j:=1 to #vec do
//       res:=Append(res,EvalRatFct(vec[j],aux));
//     end for;
//     return res;
    return [ EvalRatFct(term,aux) : term in vec ];
end function;


// Generates the polynomial for Witt sums
// the n is (lenght - 1)
WittSumPol := function(n,p)
     PR:=PolynomialRing(Integers(),2*(n+1));
     PRp:=PolynomialRing(GF(p),2*(n+1));
     AssignNames(~PR,
      [ "X" cat IntegerToString(i) : i in [0..n] ] cat 
      [ "Y" cat IntegerToString(i) : i in [0..n] ] );
     res:=[PR.1 + PR.(n+2)];
     resp:=[PRp.1 + PRp.(n+2)];
     for j:=1 to n do
       tmp := (&+[ p^i*(PR.(i+1))^(p^(j-i)) : i in [0..j] ])
             +(&+[ p^i*(PR.(n+i+2))^(p^(j-i)) : i in [0..j] ]);
       tmp -:=(&+[ p^i*(res[i+1])^(p^(j-i)) : i in [0..j-1] ]);
       tmp div:= p^j;
       tmp:= PRp!tmp;
       resp:=Append(resp,tmp);
       res:=Append(res,PR!tmp);
     end for;
   return resp;
end function;


function newWittSumPol(n,p)
    PP:=PolynomialRing(ResidueClassRing(p),2*(n+1));
    AssignNames(~PP,
      [ "X" cat IntegerToString(i) : i in [0..n] ] cat 
      [ "Y" cat IntegerToString(i) : i in [0..n] ] );
    res:=[* *];
    for i in [1..(n+1)] do
        P1:=PolynomialRing(ResidueClassRing(p^i),2*(n+1));
        res[i]:=p^(i-1)*(P1.i + P1.(n+1+i));
    end for;
    // now add the power of the terms which we finish
    for i in [2..(n+1)] do
        print "iteration ", i, " of ", n+1;
        t:=res[i-1];
        // add the powers of the precomputed to all
        for j in [i..(n+1)] do
            print "    iteration", j-i+1, " of ", n-i+2;
            P1:=Parent(res[j]);
            t:=(P1!t)^p;
            res[j]+:=p^(i-2)*(P1.(i-1)^(p^(j-i+1))+P1.(n+i)^(p^(j-i+1))-t);
        end for;
        res[i]:=PP!(res[i] div p^(i-1));
       // res[i] is done!
    end for;
    res[1]:=PP!res[1];
    return res;
end function;


/*

// test

p:=11;
n:=3;
time v1:=WittSumPol(n,p);
time v2:=newWittSumPol(n,p);
P:=Parent(v2[1]);
[ P!x : x in v1 ] eq [ x : x in v2 ];

*/



// Generates the polynomial for Witt products
// the n is (lenght - 1)
WittProdPol := function(n,p)
     PR:=PolynomialRing(Integers(),2*(n+1));
     PRp:=PolynomialRing(GF(p),2*(n+1));
     AssignNames(~PR,
      [ "X" cat IntegerToString(i) : i in [0..n] ] cat 
      [ "Y" cat IntegerToString(i) : i in [0..n] ] );
     res:=[PR.1 * PR.(n+2)];
     resp:=[PRp.1 * PRp.(n+2)];
     for j:=1 to n do
       tmp := (&+[ p^i*(PR.(i+1))^(p^(j-i)) : i in [0..j] ])
             *(&+[ p^i*(PR.(n+i+2))^(p^(j-i)) : i in [0..j] ]);
       tmp -:=(&+[ p^i*(res[i+1])^(p^(j-i)) : i in [0..j-1] ]);
       tmp div:= p^j;
       tmp:= PRp!tmp;
       resp:=Append(resp,tmp);
       res:=Append(res,PR!tmp);
     end for;
   return res;
end function;


function newWittProdPol(n,p)
    PP:=PolynomialRing(ResidueClassRing(p),2*(n+1));
    AssignNames(~PP,
      [ "X" cat IntegerToString(i) : i in [0..n] ] cat 
      [ "Y" cat IntegerToString(i) : i in [0..n] ] );
    P1:=PolynomialRing(ResidueClassRing(p),2*(n+1));
    res:=[* P1.1* P1.(n+2) *];
    for i in [2..(n+1)] do
        P1:=PolynomialRing(ResidueClassRing(p^(i-1)),2*(n+1));
        // AssignNames(~P1,
        //             [ "x" cat IntegerToString(i) : i in [0..n] ] cat 
        //             [ "y" cat IntegerToString(i) : i in [0..n] ] ); 
         res[i]:=p^(i-2)*( &+[ (P1.j)^(p^(i-j))*(P1.(n+2+i-j))^(p^(j-1)) : j in [1..i] ] );
    end for;
    // print res;
    res[2]:=PP!res[2];
    // now add the power of the terms which we finish
    for i in [3..(n+1)] do
        print "iteration ", i, " of ", n+1;
        t:=res[i-1];
        // add the powers of the precomputed to all
        for j in [i..(n+1)] do
            print "    iteration", j-i+1, " of ", n-i+2;
            P1:=Parent(res[j]);
            t:=(P1!t)^p;
            res[j]+:=p^(i-3)*( &+[ (P1.r)^(p^(j-r))*(P1.(n+1+i-r))^(p^(j-i+r)) : r in [1..(i-1)] ]-t);
        end for;
        res[i]:=PP!(res[i] div p^(i-2));
       // res[i] is done!
    end for;
    res[1]:=PP!res[1];
    return res;
end function;



/*

// test

p:=3;
n:=5;
time v1:=WittProdPol(n,p);
time v2:=newWittProdPol(n,p);
P:=Parent(v2[1]);
[ P!x : x in v1 ] eq [ x : x in v2 ];

*/




// Generates the polynomial for Witt squares (necessary???)
// the n is (lenght - 1)
WittSquarePol := function(n,p)
     PR:=PolynomialRing(Integers(),(n+1));
     PRp:=PolynomialRing(GF(p),(n+1));
     AssignNames(~PR,
		 [ "X" cat IntegerToString(i) : i in [0..n] ]);
     res:=[(PR.1)^2];
     resp:=[(PRp.1)^2];
     for j:=1 to n do
       tmp := (&+[ p^i*(PR.(i+1))^(p^(j-i)) : i in [0..j] ]);
       tmp := tmp^2;
       tmp -:=(&+[ p^i*(res[i+1])^(p^(j-i)) : i in [0..j-1] ]);
       tmp div:= p^j;
       aux:=< PRp.i : i in [1..(n+1)] >;
       tmp:= Evaluate(tmp,aux);
       resp:=Append(resp,tmp);
       res:=Append(res,PR!tmp);
     end for;
     return res;
end function;


// Sum two Witt vectors -- vectors over a ring in charact. p
// If the pol that gives Witt sum are already computed, pass it
// as vecS
WittSum := function(v1,v2 : vecS:=[])
    n:=#v1;
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    if vecS eq [] then
	vecS:=WittSumPol(n-1,p);
    end if;
    aux:= < (v1 cat v2)[i] : i in [1..2*n] >;
    return EvalVec(vecS,aux);
end function;



// Mult. two Witt vectors -- vectors over a ring in charact. p
// If the pol that gives Witt prod. are already computed, pass it
// as vecP
WittProd := function(v1,v2 : vecP:=[])
    n:=#v1;
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    if vecP eq [] then
	vecP:=WittProdPol(n-1,p);
    end if;
    aux:=< (v1 cat v2)[i] : i in [1..2*n] >;
    return EvalVec(vecP,aux);
end function;


// Generates the polynomial for Witt negatives
// the n is (lenght - 1)
WittNegPol := function(n,p)
     PR:=PolynomialRing(Integers(),2*(n+1));
     PRp:=PolynomialRing(GF(p),2*(n+1));
     AssignNames(~PR,
      [ "X" cat IntegerToString(i) : i in [0..n] ] cat 
      [ "Y" cat IntegerToString(i) : i in [0..n] ] );
     resp:=[-PRp.1];
     vX:=[ PRp.i : i in [1..(n+1)]];
     vY:=[ PRp.i : i in [(n+2)..2*(n+1)]];
     vSum:=WittSum(vX,vY);

     evalvec:=< PRp.i : i in [1..2*(n+1)] >;

     for j:=1 to n do
        evalvec[n+1+j]:=resp[j];
        vSum[j+1]:=Evaluate(vSum[j+1],evalvec);
        tmp := -(vSum[j+1] - PRp.(n+2+j));
        resp:=Append(resp,tmp);

     end for;

     res:=[ PR!resp[i] : i in [1..(n+1)]];

     return res;
end function;



// computes -v
WittNeg := function(v : vecNeg:=[])
     n:=#v;
     P:=Parent(v[1]);
     p:=Characteristic(P);
     if vecNeg eq [] then
       vecNeg:=WittNegPol(n-1,p);
     end if;
     w:= v cat [ 0 : i in [1..n]];
     evalvec:=< w[i] : i in [1..2*n] >;
     return EvalVec(vecNeg,evalvec);
end function;



// Generates the polynomial for Witt inverses
// the n is (lenght - 1)
WittInvPol := function(n,p)
     PR:=PolynomialRing(Integers(),2*(n+1));
     PRp:=PolynomialRing(GF(p),2*(n+1));
     F:=FieldOfFractions(PR);
     Fp:=FieldOfFractions(PRp);
     AssignNames(~F,
      [ "X" cat IntegerToString(i) : i in [0..n] ] cat 
      [ "Y" cat IntegerToString(i) : i in [0..n] ] );
     resp:=[1/Fp.1];
     vX:=[ Fp.i : i in [1..(n+1)]];
     vY:=[ Fp.i : i in [(n+2)..2*(n+1)]];
     vProd:=WittProd(vX,vY);

     evalvec:=< Fp.i : i in [1..2*(n+1)] >;

     for j:=1 to n do
        evalvec[n+1+j]:=resp[j];
        vProd[j+1]:=EvalRatFct(vProd[j+1],evalvec);
        tmp := -(vProd[j+1] - Fp.1^(p^j)*Fp.(n+2+j))/Fp.1^(p^j);
        resp:=Append(resp,tmp);

     end for;

     res:=[ F!resp[i] : i in [1..(n+1)]];

     return res;
end function;

// computes v^(-1)
WittInv := function(v : vecInv:=[])
     n:=#v;
     P:=Parent(v[1]);
     p:=Characteristic(P);
     if vecInv eq [] then
       vecInv:=WittInvPol(n-1,p);
     end if;
     w:= v cat [ 0 : i in [1..n]];
     evalvec:=< w[i] : i in [1..2*n] >;
     return EvalRFVec(vecInv,evalvec);
end function;


    
// computer v^n
// use successive squares???
WittPwr := function(v,m : vecP:=[]);
    n:=#v;
    P:=Parent(v[1]);
    p:=Characteristic(P);

    if vecP eq [] then
	vecP:=WittProdPol(n-1,p);
    end if;

    if m eq 0
	then return [ P!1 ] cat [ P!0 : i in [1..(n-1)] ];
    end if;

    res:=v;
    
    for i in [2..Abs(m)] do
	res:=WittProd(res,v : vecP:=vecP);
    end for;

    if m gt 0 then
	return res;
    else return WittInv(res);
    end if;
    
end function;


// Frobenius of Witt vectors
WittFrob := function(v)
    p:=Characteristic(Parent(v[1]));
    return [ x^p : x in v ];
end function;

// Computes v/w if w is non-zero
WittQuot := function(v,w : vecInv:=[], vecP:=[])
    n:=#v;
    P:=Parent(v[1]);
    p:=Characteristic(P);

    if vecP eq [] then
	vecP:=WittProdPol(n-1,p);
    end if;

    if vecInv eq [] then
	vecInv:=WittInvPol(n-1,p);
    end if;

    return WittProd(v,WittInv(w : vecInv:=vecInv) : vecP:=vecP);
    
end function;   



// Take p-th roots -- Inverse of Frobenius
// x has to be in FF (or in a quotient field)
// **************************
// NOT VERY WELL TESTED!!!!!!!!!!!
// **************************
// USED TO BE InvFrob -- bad name!

pthroot := function(x : n:=1)
    F:=Parent(x);
    p:=Characteristic(F);

    // let's check if we are in a function field
    if IsField(F) then
	num:=Numerator(x);
	den:=Denominator(x);
	R:=Parent(num);
    else
	num:=x;
	den:=F!1;
	R:=F;
    end if;

    r:=Rank(R);


    newnum:=0;
    for term in Terms(num) do
        tterm:=term;
	for j in [1..r] do
            if r eq 1 then
                deg:=Degree(tterm);
            else
		deg:=Degree(tterm,j);
            end if;
	    
	    error if (deg mod p^n) ne 0,
		  "Exponent of ", tterm, " not divisible by ", p^n;
	    
	    if r eq 1 then
                tterm:=Evaluate(tterm,1)*(R.1)^(deg div p^n);
            else
                tterm:=Evaluate(tterm,j,1)*(R.j)^(deg div p^n);
            end if;

	end for;
	newnum+:=tterm;
    end for;

    if Degree(den) eq 0 then
        newden:=den;
    else
        newden:=0;
        for term in Terms(den) do
            tterm:=term;
	    for j in [1..r] do
                if r eq 1 then
                    deg:=Degree(tterm);
                else
		    deg:=Degree(tterm,j);
                end if;
		
		error if (deg mod p^n) ne 0,
		      "Exponent of ", tterm, " not divisible by ", p^n;
		
		if r eq 1 then
                    tterm:=Evaluate(tterm,1)*(R.1)^(deg div p^n);
                else
                    tterm:=Evaluate(tterm,j,1)*(R.j)^(deg div p^n);
                end if;
            
	    end for;
	    newden+:=tterm;

        end for;
    end if;
    
    return F!newnum/F!newden;

end function;





// compute prim. roots of 1 to help
// change integers to Witt vectors
GenPrimRoots := function(p,n)
    if p eq 2
      then return [0,1];
    end if;

    res:=[0,1];

    ZZ:=Integers();
    PR<x>:=PolynomialRing(ZZ);
    

    R:=ResidueClassRing(p^n);
    
    f:=x^(p-1)-1;
    fp:=(p-1)*x^(p-2);

    for i in [2..(p-2)] do
	root:=R!i;
	while (Evaluate(f,root) ne R!0) do 
	    root-:=Evaluate(f,root)/(Evaluate(fp,root));
	end while;
	Append(~res,ZZ!root);
    end for;

    // Append(~res,-1);
    Append(~res,p^n-1);

    return(res);

end function;
	    


// Convert Integer to Witt Vector
IntToWitt := function(n,p,l : primr:=[])

// n = number to be converted
// p = prime
// l = length
// primr = primitive roots 

    Zp:=GF(p);
    ZZ:=Integers();

    if primr eq []
	then primr:=GenPrimRoots(p,l);
    end if;

    // print primr;
    
    res:=[];
    numb:=n mod p^l;
    
    for i in [1..l] do
	entry:=numb mod p;
	Append(~res,Zp!entry);
	// print "i = ", i;
	// print "numb = ", numb;
	// print "entry = ", entry;
	// print "primr[entry+1] = ", primr[entry+1];
	// print "(numb-primr[entry+1]) div p = ", (numb-primr[entry+1]) div p;
	// print "";
	numb:=((numb-primr[entry+1]) div p);
    end for;

    return res;

end function;


// given an integer m and a vector v, computes
// m*v.  (Previous way was BADDD.)
WittMult := function(v,m : vecP:=[], primr:=[])
    n:=#v;
    P:=Parent(v[1]);
    p:=Characteristic(P);

    if vecP eq [] then
	vecP:=WittProdPol(n-1,p);
    end if;

    if primr eq []
	then primr:=GenPrimRoots(p,n);
    end if;

    vecm:=IntToWitt(m,p,n : primr:=primr);

    return WittProd(vecm,v : vecP:=vecP);
    
end function;

/*
Previous way...  BADDD!!!

WittMult := function(v,m : vecS:=[])
    n:=#v;
    P:=Parent(v[1]);
    p:=Characteristic(P);

    if vecS eq [] then
	vecS:=WittSumPol(n-1,p);
    end if;

    if m eq 0
	then return [ P!0 : i in [1..n] ];
    end if;

    res:=v;
    
    for i in [2..Abs(m)] do
	res:=WittSum(res,v : vecS:=vecS);
    end for;

    if m gt 0 then
	return res;
    else return WittNeg(res);
    end if;
    
end function;
*/





// pol a polynomial with integer coefficients;
// v is a vector of Witt vectors;
// the function evaluates the polynomial in v
EvalPolWitt := function(pol,v : vecS:= [], vecP:= [], primr:=[])
    PR:=Parent(pol);
    r:=Rank(PR);
    Fp:=Parent(v[1][1]);
    p:=Characteristic(Fp);
    l:=#v[1];

    ZZ:=Integers();
    
    error if r gt #v,
	"Rank of the polynomial is greater than the lenth of the vector";

    // pre-compute auxilar data
    if vecS eq [] then
	vecS:=WittSumPol(l-1,p);
    end if;

    if vecP eq [] then
	vecP:=WittProdPol(l-1,p);
    end if;

    if primr eq [] then
	primr:=GenPrimRoots(p,l);
    end if;

    // vector with the degress in all variables
    if r eq 1 then
	vdeg:= [ Degree(pol) ];
    else
	vdeg:=[ Degree(pol,i) : i in [1..r] ];
    end if;

    // compute the necessary powers
    // vecpwr[i][j] has v[i]^(j-1)
    vecpwr:=[* *];

    vid:=[ Fp!1 ] cat [ Fp!0 : i in [2..l] ];
    
    for i in [1..r] do
	vtmp:=[ vid, v[i] ];
	term:=v[i];
	for j in [2..vdeg[i]] do
	    Append(~vtmp,WittProd(term,vtmp[j] : vecP:=vecP));
	end for;
	Append(~vecpwr,vtmp);
    end for;

    mon:=Monomials(pol);
    coef:=Coefficients(pol);

    nmon:=[];
    for i in [1..#mon] do
	term:=vid;
	for j in [1..r] do
	    if r eq 1 then
		deg:=Degree(mon[i]);
	    else deg:=Degree(mon[i],j);
	    end if;
	    if deg ne 0 then // try to save a few products here
		term:=WittProd(term,vecpwr[j][deg+1] : vecP:=vecP);
	    end if;
	end for;
	Append(~nmon,term);
    end for;

    res:=[ Fp!0 : i in [1..l] ];
    
    for i in [1..#coef] do
	tmp:=IntToWitt(ZZ!coef[i],p,l : primr:=primr);
	tmp:=WittProd(tmp,nmon[i] : vecP:=vecP);
	res:=WittSum(res,tmp : vecS:=vecS);
    end for;

    return res;
end function;





    
/*
F<a3,d2,d3,e0,e1,e2,e3,f1,f2,f3,g1,g2,g3,h0>:=
FieldOfFractions(PolynomialRing(GF(2),14));

va:=[0, 0, e0 + e1^2 + h0, a3];
vd:=[0, e0 + h0, d2, d3];
ve:=[e0, e1, e2, e3];
vf:=[0, f1, f2, f3];
vg:=[0, g1, g2, g3];
vh:=[h0, h0^2, g1 + h0^4, g1^4*h0^4 + h0^8 + f1^4*h0^8];




j1num:=  4*AA1^2*EE1^8 + 4*AA1*EE1^9 + 9*EE1^10;
j1den:=    8*FF1*HH1 + 4*GG1^2 + 11*HH1;


j2num:=    1;
j2den:=    8*FF1*HH1^3 + 12*GG1^2*HH1^2 + 3*HH1^3;


j3num:=    8*AA1^4*DD1^3 + 6*AA1^2*DD1^4 + 8*AA1*DD1^6 + 12*DD1^6*EE1 
              + 8*DD1^5*FF1 + 7*DD1^5 + 2*DD1^4*GG1 + 8*DD1^3*GG1^2;
j3den:=    8*FF1*HH1 + 4*GG1^2 + 11*HH1;

j4num:=    DD1^8 + 4*DD1^6*HH1 + 6*DD1^4*HH1^2 + 4*DD1^2*HH1^3 + HH1^4;
j4den:=    8*FF1*HH1 + 4*GG1^2 + 11*HH1;

*/
