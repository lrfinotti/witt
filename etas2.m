// given pol(X,Y) and v=[v1,...,vn], computes pol(v)
// using the "breaking in half" idea
// VERY WASTEFUL: recomputes some parts in the recurion
//                but I see no way around it...
function multiplev(v,pols)
    lgt:=#v;
    r:=#pols;
    if lgt le 1 then
        return [ 0 : i in [1..r] ];
    end if;
    if lgt eq 2 then
        return [ Evaluate(pols[i],v) : i in [1..r] ];
    end if;
    // split v in two
    v1:=[ v[i] : i in [1..(lgt div 2)] ];
    // print "v1=",v1;
    v2:=[ v[i] : i in [((lgt div 2)+1)..lgt] ];
    // print "v2=",v2;
    x1:=$$(v1,pols);
    // print "x1=",x1;
    // read tt;
    x2:=$$(v2,pols);
    // print "x2=",x2;
    // read tt;
    x3:=$$([ &+v1, &+v2 ],pols);
    // print "x3=",x3;
    res:=[ ];
    for i in [1..r] do
        // print "i = ", i;
        // print res;
        // print Parent(x1[i]), Parent(x2[i]), Parent(x3[i]);
        // print x1[i], x2[i], x3[i];
        res[i]:=[ x1[i],x2[i],x3[i] ];
    end for;
    delete x1,x2,x3,v1,v2;
    // read tt;
    for i in [2..r] do
        pols:=pols[1..(r-i+1)];
        tmp:=$$(res[i-1],pols);
        for t in [i..r] do
            Append(~res[t],tmp[t-i+1]);
        end for;
        delete tmp;
    end for;
    return [ &+x : x in res ];
end function;


/*

load 'etas2.m';
p:=7;P<x,y>:=PolynomialRing(GF(p),2);
v:=[P!2, x*y, 2*x+y, x^2+y^2, x+y+5];
// v:=[P!2,x*y, 2*x+y];
pols:=[x,x^2,y,y^2];
// pols:=[x,y];

multiplev(v,pols);

*/



function sigma(x,n,p)
    if n eq 0 then
        return x;
    end if;
    /* R:=ResidueClassRing(p^(n+1)); */
    /* y:=R!x; */
    /* y^:=p^(n-1); */
    /* z:=y^p; */
    /* z-:=y; */
    /* z*:=-1; */
    /* return (Integers()!z) div p^n; */
    return (x^(p^(n-1)) - x^(p^n)) div p^n;
end function;




function tau(p,k,l,i)
    if (i mod p) eq 0 then
        print "ERROR: i=", i, " is div. by p=", p;
        return 0;
    end if;
    if k eq 0 then
        return (-Binomial(p^l,i) div p^l);
    end if;
    return &+[ sigma($$(p,k-j,l,i),j,p) : j in [1..k] ];
end function;


// DOESN'T WORK: for n>2...  There are overlaps between terms of the summation an
// and terms coming from the etas (NOT DISJOINT) and hence my deduction doesn't work.
function etapols2(p,n)
    P<x,y>:=PolynomialRing(GF(p),2);
    res:=[ &+[ &+[ tau(p,k-j,j,i)*x^(i*p^(k-j))*y^(p^k-i*p^(k-j)) 
                   : i in [1..(p^j-1)] | ((i mod p) ne 0) ] : j in [1..k] ] 
           : k in [1..n] ];
    print res;
    // aux = [ [n1], [n1,1, n2], [n2,1, n1,1,1 n1,2, n3], ....]
    aux:=[ [[* res[1],-1 *]] ];
    for i in [2..n] do
        aux[i]:=[];
        for j in [1..(i-1)] do
            for k in [1..#aux[j]] do
                // wasteful!  recompute some parts...
                sign:=aux[j][k][2]*(-1);
                Append(~aux[i],[* multiplev(Terms(aux[j][k][1]),res[1..(i-j)])[i-j],sign *]);
            end for;
        end for;
        res[i]+:=&+[ x[1]*x[2] : x in aux[i]]; // finished eta_i
        Append(~aux[i],[* res[i],-1 *]);
    end for;
    print aux;
    return res;
end function;


// DOESN'T WORK: it was just a guess
function etapols22(p,n)
    P<x,y>:=PolynomialRing(GF(p),2);
    res:=[ &+[ &+[ tau(p,k-j,j,i)*x^(i*p^(k-j))*y^(p^k-i*p^(k-j)) 
                   : i in [1..(p^j-1)] | ((i mod p) ne 0) ] : j in [1..k] ] 
           : k in [1..n] ];
    // print res;
    for i in [2..n] do
        for j in [1..(i-1)] do
            res[i]+:=multiplev(Terms(res[i-j]),res[1..j])[j];
        end for;
        // res[i] done!
    end for;
    return res;
end function;






function etapols3(p,n)
    P<x,y>:=PolynomialRing(GF(p),2);
    res:=[];
    for k in [1..n] do
        res[k]:=[];
        for j in [1..k] do
            for i in [x : x in [1..(p^j-1)] | (x mod p) ne 0 ] do
                Append(~res[k], tau(p,k-j,j,i)*x^(i*p^(k-j))*y^(p^k-i*p^(k-j)));
            end for;
        end for;
    end for;
    // print res;
    for i in [2..n] do
        for j in [1..(i-1)] do
            Append(~res[i],multiplev(res[j],[ &+[ x : x in res[l] ] : l in [1..(i-j)] ])[i-j]);
        end for;
        // res[i] done!
    end for;
    return [ &+[ x : x in res[i] ] : i in [1..n] ];
end function;



// WRONG!!!!
function etapols4(p,n)
    P<x,y>:=PolynomialRing(GF(p),2);
    res:=[];
    for k in [1..n] do
        res[k]:=[];
        for i in [x : x in [1..(p^k-1)] | (x mod p) ne 0 ] do
            Append(~res[k], tau(p,0,k,i)*x^(i)*y^(p^k-i));
        end for;
    end for;
    // print res;
    for i in [2..n] do
        for j in [1..(i-1)] do
            Append(~res[i],multiplev(res[j],[ &+[ x : x in res[l] ] : l in [1..(i-j)] ])[i-j]);
        end for;
        // res[i] done!
    end for;
    return [ &+[ x : x in res[i] ] : i in [1..n] ];
end function;



load 'witt_gen.m';

// WRONG!!!!
function etapols5(p,n)
    P<x,y>:=PolynomialRing(GF(p),2);
    res:=[];
    for k in [1..n] do
        res[k]:=[];
        for i in [x : x in [1..(p^k-1)] | (x mod p) ne 0 ] do
            Append(~res[k], tau(p,0,k,i)*x^(i)*y^(p^k-i));
        end for;
    end for;
    // print res;
    for i in [2..n] do
        for j in [1..(i-1)] do
            Append(~res[i],multiplev(res[j],[ &+[ x : x in res[l] ] : l in [1..(i-j)] ])[i-j]);
        end for;
        // res[i] done!
    end for;
    return [ &+[ x : x in res[i] ] : i in [1..n] ];
end function;








/*

load 'etas.m'; load 'etas2.m';
p:=3;
n:=4;

time v1:=etapols(p,n);
time v2:=etapols3(p,n);
// time v3:=etapols4(p,n);
P:=Parent(v1[1]);
v2:=[ P!x : x in v2 ];
// v3:=[ P!x : x in v3 ];
v1 eq v2;
// v1 eq v3;




P<x,y>:=Parent(v1[1]);

v1[1] eq v2[1];
v1[2] eq v2[2];
v1[3] - v2[3];


PP<X,Y>:=PolynomialRing(Rationals(),2);                                  

eta1:=(X^p+Y^p - (X+Y)^p)/p; P!eta1;                                      
eta2:=(X^(p^2)+Y^(p^2) - (X+Y)^(p^2))/p^2 - eta1^p/p; P!eta2;             
eta3:=(X^(p^3)+Y^(p^3) - (X+Y)^(p^3))/p^3 - eta1^(p^2)/p^2 - eta2^p/p; P!eta3;                                           

v1 eq [ P!eta1, P!eta2, P!eta3 ];
v2 eq [ P!eta1, P!eta2, P!eta3 ];

pols:=v1;
tmp:=vetav(p,2,Terms(v1[1]) : pols:=pols[1..2]);                                                                      
n11:=tmp[1];
n21:=tmp[2];
n111:=vetav(p,1,Terms(n11) : pols:=pols[1..1])[1];                                                                     
n12:=vetav(p,1,Terms(v1[2]) : pols:=pols[1..1])[1];                                                                       

tmp:=n12+n21-n111;

i:=13; Integers()!(1/p*(-1/(p^2)*Binomial(p^3,p*i) - (-1/p^2*Binomial(p^2,i))^p)) mod p;








*/
