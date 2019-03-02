load 'etas.m';

function bin_tab(p,n,m)
    // binomial table up to choose n and modulo m
    F:=ResidueClassRing(m);
    res:=[[F!1],[F!1,F!1]];  // res[i][j]=bin(i-1,j-1)
    for i in [2..n] do
        v:=[F!1];
        for j in [1..(i-1)] do
            Append(~v,res[i][j] + res[i][j+1]);
        end for;
        Append(~v,F!1);
        Append(~res,v);
    end for;
    return res;
end function;


load 'witt_gen.m';


function Not_Zero_Vec(v)
    // checks if a vector is the zero vector
    for i in [1..#v] do
        if v[i] ne 0 then
            return true;
        end if;
    end for;

    return false;

end function;


function GT_der(f,p,i,r,n,tab,primr,pols)
    // to use to compute the GT
    // compute f^{(i,r-i)} (assuming 0<=i<=r)
    // needs to supply: table with binomials (tab)
    //                  primitive roots (primr)
    //                  vector to perform sums of Witt vectors (vecP)

    if r eq 0 then
        derf:=f;  // no derivative
    else
        derf:=[]; // derivative
        for t in f do
            // loop over the terms of f

            coef:=t[1];
            degx:=t[2];
            degy:=t[3];

            // print "coef = ", coef;

            if (degx-i ge 0) and (degy-(r-i) ge 0) and Not_Zero_Vec(coef) then
                // make sure derivative of the term is not zero

                // print "degx, degy = ", degx, degy;

                bin1:=tab[degx+1][i+1];
                bin2:=tab[degy+1][r-i+1];

                b12:=Integers()!(bin1*bin2);

                if (b12 mod p^(n-r+1)) eq 1 then
                    // no need to multiply
                    Append(~derf, [* [ coef[i] : i in [1..(n+1-r)] ] , degx-i, degy-(r-i) *]);
                elif (b12 mod p^(n-r+1)) ne 0 then
                    // if the term is not zero...
                    v1:=IntToWitt(b12,p,n+1-r : primr:=primr);

                    // print "v1 = ", v1;
                    // print "[ coef[",i,"]: i in [1..(n+1-r)] ] = ", [ coef[i] : i in [1..(n+1-r)] ];
                    // print "prod = ", newWittProd(v1,[ coef[i] : i in [1..(n+1-r)] ] : pols:=pols);

                    Append(~derf,[* newWittProd(v1,[ coef[i] : i in [1..(n+1-r)] ] : pols:=pols) , degx-i, degy-(r-i) *]);
                end if;
            end if;
        end for;

        // derivative is zero?
        if #derf eq 0 then
            derf:=[ [* [0 : s in [0..n-r]], 0, 0 *] ];
        end if;
    end if;

    return derf;

end function;




function split_Ds(f,p,max,PR,primr)
    // given a polynomial f (likely over Z/p^n) returns
    // [ \xi_0(f) , \xi_1(f), ... , \xi_max(f) ]
    
    res:=[ PR!0 : i in [0..max] ];
    for mon in Monomials(f) do
        coef:=Integers()!(MonomialCoefficient(f,mon));
        vcoef:=IntToWitt(coef,p,max+1 : primr:=primr);
        for i in [1..(max+1)] do
            res[i]+:=PR!(vcoef[i])*PR!mon;
        end for;
    end for;

    return res;

end function;



function GT_Ds(i,r,n,p,primr,PR)
    // assumes r>=1!!!

    PR1:=PolynomialRing(ResidueClassRing(p^(n+1)),2*n+2);

    PR2<t>:=PolynomialRing(PR1);

    tmp:=Min(n,n-r+1);
    
    t1:=&+[ PR2.1^s*PR1.(s+1)^(p^(n-s)) : s in [1..tmp] ];
    t2:=&+[ PR2.1^s*PR1.(n+s+2)^(p^(n-s)) : s in [1..tmp] ];

    t:=t1^i*t2^(r-i);
    delete t1, t2;

    
    t:=[ Coefficient(t,s) : s in [r..n] ];
    // t = [ D_{r,n}^{(i,r-i)}, ... , D_{n,n}^{(i,r-i)} ]

    // now, split each D_{s,n}^{(i,r-i)} into
    // [ D_{s,n,0}^{(i,r-i)}, ..., D_{r,n,n-s}^{(i,r-i)} ]
    return [ split_Ds(f,p,n-r,PR,primr) : f in t ];

end function;



function Pol_Root(f,r)
    // takes "r-th root" of a polyn.
    // in fact, leaves coefficient the same and divide
    // the powers of variables by n

    P:=Parent(f);
    n:=Rank(P);
    res:=0;
    for mon in Monomials(f) do
        nmon:=MonomialCoefficient(f,mon);
        for i in [1..n] do
            nmon*:=(P.i)^(Degree(mon,i) div r);
        end for;
        res+:=nmon;
    end for;
    
    return res;

end function;



function GT(f : pols:=[], primr:=[], tab:=[], vvars:=[])
    // f must be given by f:=[ ... , [[a0,a1,...,an],i,j] , ... ]
    n:=#(f[1][1])-1; // length of the coefficients - 1
    P:=Parent(f[1][1][1]);
    p:=Characteristic(P);

    if #pols eq 0 then
        pols:=etapols(p,n);
    end if;

    PR:=PolynomialRing(P,2*n+2);
    AssignNames(~PR,
                [ "x" cat IntegerToString(i) : i in [0..n] ] cat 
                [ "y" cat IntegerToString(i) : i in [0..n] ] );

    maxdeg:=Max([ Max(t[2],t[3]) : t in f ]); // to find the binomial table size

    if #tab eq 0 then
        tab:=bin_tab(p,maxdeg,p^(n+1)); // table of binomials to use
    end if;

    if #primr eq 0 then
        primr:=GenPrimRoots(p,n+1); // (seems it needs length)
    end if;

    res:=[ [] : s in [0..n] ];

    for r in [0..n] do
        for i in [0..r] do
            
            // first, compute the derivative
            derf:=GT_der(f,p,i,r,n,tab,primr,pols);

            // print "derf = ", derf;
            
            // now we need the Ds
            Ds:=GT_Ds(i,r,n,p,primr,PR);
            // Ds=[ [D_{r,n,0}^{(i,r-i)} , ... , D_{r,n,n-r}^{(i,r-i)}],
            //      [D_{r+1,n,0}^{(i,r-i)} , ... , D_{r,n,n-r-1}^{(i,r-i)}],
            //      ...,
            //      [D_{n,n,0}^{(i,r-i)}] ]


            // print "Ds = ", Ds;

            for m in [r..n] do
                for j in [r..m] do
                    for k in [r..j] do

                        // print "r, i, m, j, k = ", r, i, m, j, k;

                        // split the derivative (i.e., the \xi_i's)
                        // also takes care of the Frob and evalutation
                        vt:=&+[ (t[1][m-j+1])^(p^j)*(PR.1)^(t[2]*p^m)*(PR.(n+2))^(t[3]*p^m) : t in derf ];

                        // create a vector to evaluate D_{r,n}
                        // to get D_{r,m}
                        v:=< PR!0 : s in [1..(2*n+2)] >;
                        for s in [2..(m+1)] do
                            v[s]:=PR.s;
                            v[s+n+1]:=PR.(s+n+1);
                        end for;
                        
                        // get D_{r,m} from D_{r,n}
                        tmp:=Evaluate(Ds[k-r+1][j-k+1],v);
                        tmp:=Pol_Root(tmp,p^(n-m));

                        if #vvars ne 0 then 
                            res[m+1]:=res[m+1] cat [ Evaluate(term,vvars) : term in Terms(vt*tmp)];
                        else
                            res[m+1]:=res[m+1] cat Terms(vt*tmp);
                        end if;
                        
                    end for; // k
                end for; // j
            end for; // m
        end for; // i

        delete tmp, v, vt;

    end for; // r

    // now add the etas...

    // print "pols = ", pols;

    for i in [1..n] do

        // print "res[", i, "] = ", res[i];

        ve:=vetav(p,n-i+1,res[i] : pols:=pols);
        
        // print "ve = ", ve;

        for j in [(i+1)..(n+1)] do
            if ve[j-i] ne 0 then
                Append(~res[j],ve[j-i]);
            end if;
        end for;

        // don't need all pols anymore
        pols:=pols[1..(n-i)];

    end for;

    return [ &+t : t in res ];

end function; 


/*

// Test
> F:=GF(3);
> f:=[ [* [ F!1, F!0, F!2 ], 1, 2 *], [* [ F!1, F!1, F!1 ], 3, 2 *] ];
> load 'gt.m'; GT(f);
Loading "gt.m"
Loading "witt_gen.m"
Loading "etas.m"
[
    x0^3*y0^2 + x0*y0^2,
    x0^9*y0^6 + 2*x0^9*y0^3*y1 + 2*x0^7*y0^6 + 2*x0^5*y0^6 + 2*x0^3*y0^3*y1 + x1*y0^6,
    x0^27*y0^18 + x0^27*y0^15*y1 + 2*x0^27*y0^12*y1^2 + 2*x0^27*y0^9*y2 + x0^27*y1^6 +
        x0^25*y0^15*y1 + x0^25*y0^12*y1^2 + 2*x0^23*y0^18 + 2*x0^23*y0^15*y1 +
        x0^23*y0^12*y1^2 + x0^21*y0^18 + x0^21*y0^12*y1^2 + x0^21*y0^9*y1^3 +
        x0^19*y0^18 + 2*x0^19*y0^15*y1 + 2*x0^19*y0^12*y1^2 + x0^18*x1^3*y0^18 +
        2*x0^18*x1*y0^18 + 2*x0^18*x1*y0^15*y1 + 2*x0^18*x1*y0^12*y1^2 + 2*x0^17*y0^18
        + 2*x0^17*y0^15*y1 + 2*x0^17*y0^12*y1^2 + 2*x0^16*x1*y0^18 + x0^16*x1*y0^15*y1
        + 2*x0^15*y0^15*y1 + 2*x0^15*y0^12*y1^2 + x0^15*y0^9*y1^3 + x0^14*x1*y0^18 +
        x0^14*x1*y0^15*y1 + 2*x0^13*y0^18 + x0^13*y0^15*y1 + x0^13*y0^12*y1^2 +
        x0^12*x1*y0^18 + 2*x0^12*x1*y0^15*y1 + x0^12*x1*y0^12*y1^2 + 2*x0^11*y0^18 +
        x0^11*y0^12*y1^2 + 2*x0^10*x1*y0^18 + x0^10*x1*y0^15*y1 + 2*x0^9*x1^2*y0^18 +
        x0^9*x1^2*y0^15*y1 + 2*x0^9*y0^18 + x0^9*y0^9*y1^3 + 2*x0^9*y0^9*y2 + x0^9*y1^6
        + x0^8*x1*y0^15*y1 + x0^7*x1^2*y0^18 + 2*x0^6*x1*y0^12*y1^2 + x0^5*x1^2*y0^18 +
        x0^3*x1^2*y0^15*y1 + 2*x1^3*y0^9*y1^3 + x2*y0^18
]
> pols:=etapols(3,2);
> PR<x0,x1,x2,y0,y1,y2>:=PolynomialRing(F,6);
> x:=[ x0, x1, x2 ];
> y:=[ y0, y1, y2 ];
> ff:=newWittProd([ F!1, F!0, F!2 ], x : pols:=pols);                
> y2:=newWittProd(y,y : pols:=pols);                                 
> ff:=newWittProd( ff, y2 : pols:=pols);                             
> x3:=newWittProd(x,x : pols:=pols);                                 
> x3:=newWittProd(x3,x : pols:=pols);
> tmp:=newWittProd( [ F!1, F!1, F!1 ], x3 : pols:=pols);
> tmp:=newWittProd( tmp, y2 : pols:=pols);
> ff:=newWittSum(ff, tmp : pols:=pols);
> gg:=GT(f);
> ff eq [ PR!term : term in gg ];
true


> f:=[ [* [ F!1, F!0, F!2, F!1 ], 1, 0 *], [* [ F!0, F!2, F!1, F!0 ], 2 , 2 *], 
[* [ F!1, F!1, F!1, F!2 ], 3, 2 *] ];
> gg:=GT(f);
> pols:=etapols(3,3);
> PR<x0,x1,x2,x3,y0,y1,y2,y3>:=PolynomialRing(F,8);
> x:=[ x0, x1, x2, x3 ];
> y:=[ y0, y1, y2, y3 ];
> ff:=newWittProd([ F!1, F!0, F!2, F!1 ], x : pols:=pols);
> x2:=newWittProd(x,x : pols:=pols);
> y2:=newWittProd(y,y : pols:=pols);
> tmp:=newWittProd( x2, y2 : pols:=pols);
> tmp:=newWittProd([ F!0, F!2, F!1, F!0 ], tmp : pols:=pols);
> ff:=newWittSum(ff, tmp : pols:=pols);
> tmp:=newWittProd(x, x2 : pols:=pols);
> tmp:=newWittProd(tmp, y2 : pols:=pols);
> tmp:=newWittProd(tmp, [ F!1, F!1, F!1, F!2 ] : pols:=pols);
> ff:=newWittSum(ff, tmp : pols:=pols);
> ff eq [ PR!term : term in gg ];
true


////////////////////////////

p:=5;
n:=3;

load 'witt_gen.m';

vecS:=WittSumPol(n,p);
vecP:=WittProdPol(n,p);
primr:=GenPrimRoots(p,n+1);

F:=GF(p);
P:=PolynomialRing(F,2*n+2);
AssignNames(~P,
            [ "x" cat IntegerToString(i) : i in [0..n] ] cat 
            [ "y" cat IntegerToString(i) : i in [0..n] ] );


// load 'load_js.m';

PP<X,Y>:=PolynomialRing(Integers(),2);
PPP:=PolynomialRing(ResidueClassRing(p^(n+2)),2);

mp:=PP!(PPP!ClassicalModularPolynomial(p));


vx:=[ P.i : i in [1..(n+1)] ];
vy:=[ P.i : i in [(n+2)..(2*n+2)] ];

vxs:=<[P!1] cat [P!0 : i in [1..n] ], vx>;


maxdeg:=Max(Degree(mp,1), Degree(mp,2));
// compute powers of vx
for i in [2..maxdeg] do
    Append(~vxs,WittProd(vxs[i],vx : vecP:=vecP));
end for;

evalv:=vy cat [ F!0 : i in [1..(n+1)] ];
// compute powers of vy
vys:=[ [ Evaluate(t,evalv) : t in vxs[i] ] : i in [1..(#vxs)] ];

res:=[ P!0 : i in [0..n] ];

// compute GT
for t in Terms(mp) do
    coef:=LeadingCoefficient(t);
    degX:=Degree(t,1);
    degY:=Degree(t,2);
    tmp:=IntToWitt(coef,p,n+1 : primr:=primr);
    tmp:=WittProd(tmp,vxs[degX+1] : vecP:=vecP);
    tmp:=WittProd(tmp,vys[degY+1] : vecP:=vecP);
    
    res:=WittSum(res,tmp : vecS:=vecS);
    
end for;

// NOW, NEW WAY

load 'etas.m';
load 'mp_gt.m';


pols:=etapols(p,n);


mmp:=Pol_New_Form(mp,p,n : primr:=primr );

res2:=GT(mmp : pols:=pols, primr:=primr);

res eq [ P!x : x in res2 ];

// Now, using the new form mmp

res3:=[ P!0 : i in [0..n] ];

// compute GT
for t in mmp do
    coef:=t[1];
    degX:=t[2];
    degY:=t[3];
    tmp:=newWittProd(coef,vxs[degX+1] : pols:=pols);
    tmp:=newWittProd(tmp,vys[degY+1] : pols:=pols);
    
    res3:=newWittSum(res3,tmp : pols:=pols);
    
end for;


*/


function Pol_New_Form(f,p,n : primr:=[] )
    // f a polynomial in two variables over Z
    // converts to format to use with GT
    if #primr eq 0 then
        primr:=GenPrimRoots(p,n+1);
    end if;

    res:=[];
    for mon in Monomials(f) do
        coef:=IntToWitt(MonomialCoefficient(f,mon),p,n+1 : primr:=primr);
        Append(~res,[* coef, Degree(mon,1), Degree(mon,2) *]);
    end for;
    
    return res;

end function;






function GT_der2(f,p,i,r,n,tab,primr,pols)
    // to use to compute the GT
    // compute f^{(i,r-i)} (assuming 0<=i<=r)
    // needs to supply: table with binomials (tab)
    //                  primitive roots (primr)
    //                  vector to perform sums of Witt vectors (vecP)

    if r eq 0 then
        derf:=f;  // no derivative
    else
        derf:=[]; // derivative
        for t in f do
            // loop over the terms of f

            coef:=t[1];
            degx:=t[2];
            degy:=t[3];

            // print "coef = ", coef;

            if (degx-i ge 0) and (degy-(r-i) ge 0) and Not_Zero_Vec(coef) then
                // make sure derivative of the term is not zero

                // print "degx, degy = ", degx, degy;

                bin1:=tab[degx+1][i+1];
                bin2:=tab[degy+1][r-i+1];

                b12:=Integers()!(bin1*bin2);

                if (b12 mod p^(n-r+1)) eq 1 then
                    // no need to multiply
                    Append(~derf, [* [ coef[i] : i in [1..(n+1-r)] ] , degx-i, degy-(r-i) *]);
                elif (b12 mod p^(n-r+1)) ne 0 then
                    // if the term is not zero...
                    v1:=IntToWitt(b12,p,n+1-r : primr:=primr);

                    // print "v1 = ", v1;
                    // print "[ coef[",i,"]: i in [1..(n+1-r)] ] = ", [ coef[i] : i in [1..(n+1-r)] ];
                    // print "prod = ", newWittProd(v1,[ coef[i] : i in [1..(n+1-r)] ] : pols:=pols);

                    Append(~derf,[* newWittProd(v1,[ coef[i] : i in [1..(n+1-r)] ] : pols:=pols) , degx-i, degy-(r-i) *]);
                end if;
            end if;
        end for;

        // derivative is zero?
        if #derf eq 0 then
            derf:=[ [* [0 : s in [0..n-r]], 0, 0 *] ];
        end if;
    end if;

    return derf;

end function;




function split_Ds2(f,p,max,primr)
    // given a polynomial f (likely over Z/p^n) returns
    // [ \xi_0(f) , \xi_1(f), ... , \xi_max(f) ]
    
    PR:=Parent(f);
    m:=Rank(PR);
    PR3:=PolynomialRing(GF(p),m); // will return over GF(p)

    res:=[ PR3!0 : i in [0..max] ];

    for mon in Monomials(f) do
        coef:=Integers()!(MonomialCoefficient(f,mon));
        vcoef:=IntToWitt(coef,p,max+1 : primr:=primr);
        for i in [1..(max+1)] do
            res[i]+:=PR3!(vcoef[i])*PR3!mon;
        end for;
    end for;

    return res;

end function;



function GT_Ds2(n,r,p,primr,PR)
    // assumes r>=1!!!

    vecDs:=AssociativeArray(); // will have all the D[n,r,i,k,t]

    PR1:=PolynomialRing(ResidueClassRing(p^(n+1-r)),2*n+2); // x's and y's

    PR2<t>:=PolynomialRing(PR1); // t


    tmp:=Min(n,n-r+1); // in case r=0
    t1:=&+[ PR2.1^s*PR1.(s+1)^(p^(n-s)) : s in [1..tmp] ];    

    vx:=[PR2!1,t1]; // contains tthe powers of t1

    // compute other powers
    for i in [2..r] do
        vx[i+1]:=t1*vx[i];
    end for;

    tmpeval:=[ PR1.(n+2+i) : i in [0..n] ] cat [ PR1!0 : i in [0..n] ];
    for i in [0..r] do
        tmp:=&+[ Evaluate(Coefficient(vx[r-i+1],dd),tmpeval)*t^dd : dd in [0..Degree(vx[r-i+1])] ];
        tmp1:=vx[i+1]*tmp;
        for k in [r..n] do
            // print "splid Ds";
            tmp2:=split_Ds2(Coefficient(tmp1,k),p,n-k,primr);
            for t in [0..(n-k)] do
                vecDs[[n,r,i,k,t]]:=PR!tmp2[t+1];
            end for;
        end for;
    end for;

    return vecDs;

end function;






function GT2(f : pols:=[], primr:=[], tab:=[], vvars:=[])
    // f must be given by f:=[ ... , [[a0,a1,...,an],i,j] , ... ]
    n:=#(f[1][1])-1; // length of the coefficients - 1
    P:=Parent(f[1][1][1]);
    p:=Characteristic(P);

    if #pols eq 0 then
        pols:=etapols(p,n);
    end if;

    PR:=PolynomialRing(P,2*n+2);
    AssignNames(~PR,
                [ "x" cat IntegerToString(i) : i in [0..n] ] cat 
                [ "y" cat IntegerToString(i) : i in [0..n] ] );

    maxdeg:=Max([ Max(t[2],t[3]) : t in f ]); // to find the binomial table size

    if #tab eq 0 then
        tab:=bin_tab(p,maxdeg,p^(n+1)); // table of binomials to use
    end if;

    if #primr eq 0 then
        primr:=GenPrimRoots(p,n+1); // (seems it needs length)
    end if;

    res:=[ [] : s in [0..n] ];

    for r in [0..n] do
        
        // print "compute Ds: r, n = ", r, n;
        vecDs:=GT_Ds2(n,r,p,primr,PR);

        for i in [0..r] do
            
            // first, compute the derivative
            // print "compute der.";
            derf:=GT_der2(f,p,i,r,n,tab,primr,pols);

            // print "derf = ", derf;
            
            for m in [r..n] do
                for j in [r..m] do
                    for k in [r..j] do

                        // print "r, i, m, j, k = ", r, i, m, j, k;

                        // split the derivative (i.e., the \xi_i's)
                        // also takes care of the Frob and evalutation
                        vt:=&+[ (t[1][m-j+1])^(p^j)*(PR.1)^(t[2]*p^m)*(PR.(n+2))^(t[3]*p^m) : t in derf ];

                        // create a vector to evaluate D_{r,n}
                        // to get D_{r,m}
                        v:=< PR!0 : s in [1..(2*n+2)] >;
                        for s in [2..(m+1)] do
                            v[s]:=PR.s;
                            v[s+n+1]:=PR.(s+n+1);
                        end for;
                        
                        // get D_{r,m} from D_{r,n}
                        // print "evaluate Ds";
                        tmp:=Evaluate(vecDs[[n,r,i,k,j-k]],v);
                        // print "pol root";
                        tmp:=Pol_Root(tmp,p^(n-m));

                        if #vvars ne 0 then
                            // print "evaluate";
                            res[m+1]:=res[m+1] cat [ Evaluate(term,vvars) : term in Terms(vt*tmp)];
                        else
                            res[m+1]:=res[m+1] cat Terms(vt*tmp);
                        end if;
                        
                    end for; // k
                end for; // j
            end for; // m
        end for; // i

        delete tmp, v, vt;

    end for; // r

    // now add the etas...

    // print "pols = ", pols;

    for i in [1..n] do

        // print "res[", i, "] = ", res[i];

        print "compute etas";
        time ve:=vetav(p,n-i+1,res[i] : pols:=pols);
        
        // print "ve = ", ve;

        for j in [(i+1)..(n+1)] do
            if ve[j-i] ne 0 then 
                Append(~res[j],ve[j-i]);
            end if;
        end for;

        // don't need all pols anymore
        pols:=pols[1..(n-i)];

    end for;

    return [ &+t : t in res ];

end function; 

