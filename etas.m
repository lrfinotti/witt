function etapols(p,n)
    PP<x,y>:=PolynomialRing(ResidueClassRing(p),2);
    res:=[* *];
    // res[i], at first has p^i*eta_i
    for i in [1..n] do
        P1<x,y>:=PolynomialRing(ResidueClassRing(p^(i+1)),2);
        res[i]:=x^(p^i)+y^(p^i) - (x+y)^(p^i);
    end for;
    // fix eta_1, which is done
    res[1]:=PP!(res[1] div p);
    // now add the power of the terms which we finish
    for i in [2..n] do
        t:=res[i-1];
        // add the powers of the precomputed to all
        for j in [i..n] do
            P1:=Parent(res[j]);
            t:=(P1!t)^p;
            res[j]-:=p^(i-1)*t;
        end for;
        res[i]:=PP!(res[i] div p^i);
       // res[i] is done!
    end for;
    return [x : x in res];
end function;


/*
// done with integers to test!
function nu2(p,n)
    P<X,Y>:=PolynomialRing(Integers(),2);
    PP<x,y>:=PolynomialRing(GF(p),2);
    res:=[* *];
    for i in [1..n] do
        res[i]:=X^(p^i)+Y^(p^i) - (X+Y)^(p^i);
    end for;
    res[1]:=res[1] div p;
    for i in [2..n] do
        t:=res[i-1];
        // add the powers of the precomputed to all
        for j in [i..n] do
            t^:=p;
            res[j]-:=p^(i-1)*t;
        end for;
        res[i]:=res[i] div p^i;
       // res[i] is done!
    end for;
    return [* PP!x : x in res *];
end function;
*/

/*

// can test with
for v in [ [3,7], [5,6], [7,5], [11,4] ] do
    p:=v[1]; n:=v[2];
    P<x,y>:=PolynomialRing(GF(p),2);
    [ P!x : x in nu2(p,n) ] eq [ P!x : x in nu(p,n) ];
end for;

// Passed thr first 3...  Didn't finish the 4th.

*/


// remove zeros of a vector to simplify comp. of etas
// if v is all zeros, return [0] to avoid problems with &+[]
function vRemoveZeros(v)

    tmp:=[ x : x in v | x ne 0 ];
    if #tmp ne 0 then
        return tmp;
    else
        //return [ Parent(v[1])!0 ];
        return [ 0 ];
    end if;
end function;


vetav := function(p,k,v : pols:=[ ])
    // the input v is a vector!
    // the output is a vector of resulting
    //    polynomials eta_i(v) for i=1, ... , k

    // firts, remove zeros!
    v:=vRemoveZeros(v);

 
    lgt:=#v;
    // print "length = ", lgt;
    if lgt eq 1 then
        P:=Parent(v[1]);
        return [ P!0 : i in [1..k] ];
    end if;
    if pols eq [ ] then
        pols:=etapols(p,k);
    end if;
    if lgt eq 2 then
        // //////////////////////////////////////////////////////////////////////
        // now, the special case of #v=2
        // //////////////////////////////////////////////////////////////////////
        return [ Evaluate(pol,v) : pol in pols ]; // needs to have pols precomputed
    else
        // //////////////////////////////////////////////////////////////////////
        // now if #v>2
        // split v in two
        v1:=[ v[i] : i in [1..(lgt div 2)] ];
        v2:=[ v[i] : i in [((lgt div 2)+1)..lgt] ];
        x1:=$$(p,k,v1 : pols:=pols);
        x2:=$$(p,k,v2 : pols:=pols);
        x3:=$$(p,k,[ &+v1, &+v2 ] : pols:=pols);
        res:=[ ];
        for i in [1..k] do
            // print "i = ", i;
            // print res;
            // print Parent(x1[i]), Parent(x2[i]), Parent(x3[i]);
            // print x1[i], x2[i], x3[i];
            res[i]:=vRemoveZeros([ x1[i],x2[i],x3[i] ]);
        end for;
        delete x1,x2,x3,v1,v2;
        for i in [2..k] do
            // lim1:=#res[i-1]; // should be 2^(i-2)*3-(2^(i-2)-1);
            // // print lim1 eq 2^(i-2)*3-(2^(i-2)-1);
            // for j in [1..(lim1-1)] do
            //     temp:=$$(p,k-i+1,[res[i-1][j],&+( res[i-1][(j+1)..lim1] )] : pols:=pols[1..(k-i+1)]);
            //     for t in [1..(k-i+1)] do
            //         Append(~res[i-1+t],temp[t]);
            //     end for;
            // end for;
            // delete temp;
            pols:=pols[1..(k-i+1)];
            tmp:=$$(p,k-i+1,res[i-1] : pols:=pols);
            for t in [i..k] do
                if tmp[t-i+1] ne 0 then 
                    Append(~res[t],tmp[t-i+1]);
                end if;
            end for;
            delete tmp;
        end for;
        return [ &+x : x in res ];
    end if;
end function;


load 'witt_gen.m';




// GOOD ONE (old BinTab3)
function BinTab(p,k : primr:=[])
    // gives back the bin. coeff need for vetav2

    if #primr eq 0 then
        primr:=GenPrimRoots(p,k+1); // to convert to bin. to Witt vec
    end if;

    R:=ResidueClassRing(p^(k+1));
    ZZ:=Integers();

    par:=(p mod 2); // partity of p

    res:=[ ];
    for i in [1..k] do
        num:=p^i; 
        tmp:=[ R!1 ];  // 1/p*Biniomial(p,1)
        for j in [2..(num-1)] do
            a:=num-(j-1);
            a:=R!(a div p^(Valuation(j-1,p)));
            b:=R!(j div p^(Valuation(j,p)));
            Append(~tmp,tmp[j-1]*a*b^(-1)); // add next binomial
            delete a, b;
        end for;
        Append(~res,[ IntToWitt(-ZZ!(tmp[j]),p,Valuation(j,p)+1 : primr:=primr)[Valuation(j,p)+1] : j in [1..#tmp] ]); // conv to Witt vec and take right coord

        delete tmp;

        // now add the symmetric part at the end
        bg:=((num-par) div 2)+1;
        for j in [bg..(num-1)] do
            res[i][j]:=res[i][num-j];
        end for;
    end for;

    return res;
end function;




/*
function BinTab2(p,k : primr:=[])

    if #primr eq 0 then
        primr:=GenPrimRoots(p,k+1);
    end if;

    R:=ResidueClassRing(p^(k+1));
    ZZ:=Integers();

    res:=[* *];
    prev:=[R!1,R!1];

    for i in [2..p^k] do
        tmp:=[* R!1 *];
        for j in [1..(#prev-1)] do
            Append(~tmp,prev[j]+prev[j+1]);
        end for;
        if (i mod 2) eq 1 then
            Append(~tmp,tmp[#tmp]);
        end if;
        prev:=tmp;
        if i in [ p^j : j in [1..k] ] then
            Append(~res,Remove(tmp,1));
        end if;
        delete tmp;
    end for;
    //print res;

    par:=(p mod 2);
    for i in [1..k] do
        num:=p^i;
        for j in [1..(#res[i]-par)] do
            val:=Valuation(j,p);
            tmp:=-ZZ!(res[i][j]) div p^(i-val);
            //print i,j,val, tmp;
            tmp:=IntToWitt(tmp,p,val+1)[val+1];
            res[i][j]:=tmp;
            delete tmp;
        end for;
        bg:=((num-par) div 2)+1;
        for j in [bg..(num-1)] do
            res[i][j]:=res[i][num-j];
        end for;
    end for;

    return [ [ X : X in Y] : Y in res];
end function;



function BinTab(p,k : primr:=[])

    if #primr eq 0 then
        primr:=GenPrimRoots(p,k+1);
    end if;

    res:=[ ];

    if p ne 2 then
        for j in [1..k] do
            tmp:=[ IntToWitt(-Binomial(p^j,i) div p^(j-Valuation(i,p)),p,Valuation(i,p)+1 : primr:=primr)[Valuation(i,p)+1] : i in [1..((p^j-1) div 2)] ];
            for i in [1..((p^j-1) div 2)] do
                Append(~tmp,tmp[((p^j-1) div 2)-i+1]);
            end for;
            Append(~res,tmp);
        end for;
    else
        for j in [1..k] do
            tmp:=[ IntToWitt(-Binomial(p^j,i) div p^(j-Valuation(i,p)),p,Valuation(i,p)+1 : primr:=primr)[Valuation(i,p)+1] : i in [1..((p^j) div 2)] ];
            for i in [2..((p^j) div 2)] do
                Append(~tmp,tmp[((p^j) div 2)-i+1]);
            end for;
            Append(~res,tmp);
        end for;
    end if;

    return res;
end function;

*/

/* Comparing the two BinTab

> p:=11; k:=3;
> time v1:=BinTab(p,k);
Time: 0.200
> time v2:=BinTab2(p,k);
Time: 0.560
> p:=11; k:=4;
> ResetMaximumMemoryUsage(); time v1:=BinTab(p,k); print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((102\
4.0)^2);
Time: 116.980
Memory Usage (in MB):  9.03125000000000000000000000000
> ResetMaximumMemoryUsage(); time v2:=BinTab2(p,k); print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((10\
24.0)^2);
Time: 51.950
Memory Usage (in MB):  9.03125000000000000000000000000
ResetMaximumMemoryUsage(); time v3:=BinTab3(p,k); print "Memory Usage (in MB): ", GetMaximumM\
emoryUsage()/((1024.0)^2);
Time: 0.510
Memory Usage (in MB):  9.03125000000000000000000000000

*/



vetav2 := function(p,k,v : bintab:=[] )
    // the input v is a vector!
    // the output is a vector of resulting
    //    polynomials eta_i(v) for i=1, ... , k

    // firts, remove zeros!
    v:=vRemoveZeros(v);
    
    lgt:=#v;
    // print "length = ", lgt;
    if lgt eq 1 then
        P:=Parent(v[1]);
        return [ P!0 : i in [1..k] ];
    end if;

    if #bintab eq 0 then
        // bintab:=BinTab(p,k);
        bintab:=BinTab(p,k);
    end if;

    //print bintab;

    if lgt eq 2 then

        x:=v[1];
        y:=v[2];


        res:=[];

        for j in [1..k] do
           res[j]:=[ bintab[j][i]*x^i*y^(p^j-i) : i in [1..(p^j-1)]];
        end for;

        for j in [1..(k-1)] do
            v:=$$(p,k-j,res[j] : bintab:=bintab );
            for i in [1..(k-j)] do
                if v[i] ne 0 then 
                    Append(~res[j+i],v[i]);
                end if;
            end for;
        end for;
        
        return [ &+term : term in res];


    else
        // //////////////////////////////////////////////////////////////////////
        // now if #v>2
        // split v in two
        v1:=[ v[i] : i in [1..(lgt div 2)] ];
        v2:=[ v[i] : i in [((lgt div 2)+1)..lgt] ];
        x1:=$$(p,k,v1 : bintab:=bintab );
        x2:=$$(p,k,v2 : bintab:=bintab );
        x3:=$$(p,k,[ &+v1, &+v2 ] : bintab:=bintab );
        res:=[ ];
        for i in [1..k] do
            // print "i = ", i;
            // print res;
            // print Parent(x1[i]), Parent(x2[i]), Parent(x3[i]);
            // print x1[i], x2[i], x3[i];
            res[i]:=vRemoveZeros([ x1[i],x2[i],x3[i] ]);
        end for;
        delete x1,x2,x3,v1,v2;
        for i in [2..k] do
            // lim1:=#res[i-1]; // should be 2^(i-2)*3-(2^(i-2)-1);
            // // print lim1 eq 2^(i-2)*3-(2^(i-2)-1);
            // for j in [1..(lim1-1)] do
            //     temp:=$$(p,k-i+1,[res[i-1][j],&+( res[i-1][(j+1)..lim1] )] : pols:=pols[1..(k-i+1)]);
            //     for t in [1..(k-i+1)] do
            //         Append(~res[i-1+t],temp[t]);
            //     end for;
            // end for;
            // delete temp;
            tmp:=$$(p,k-i+1,res[i-1] : bintab:=bintab );
            for t in [i..k] do
                if tmp[t-i+1] ne 0 then 
                    Append(~res[t],tmp[t-i+1]);
                end if;
            end for;
            delete tmp;
        end for;

        return [ &+x : x in res ];
    end if;
end function;


/* TEST


> p:=5;F:=GF(p);P<X0,X1,X2,X3,X4,Y0,Y1,Y2,Y3,Y4>:=PolynomialRing(F,10);
> v1:=vetav(p,3,[X0,Y0]);
> v2:=vetav2(p,3,[X0,Y0]);
> v1 eq v2;
true
> v2:=vetav2(p,2,[X0,X1,X2,X3]);
> v1:=vetav(p,2,[X0,X1,X2,X3]);
> v1 eq v2;
true
> time v1:=vetav(p,3,[X0,X1,X2,X3]);
Time: 350.010
> time v2:=vetav2(p,3,[X0,X1,X2,X3]);
Time: 654.070
> v1 eq v2;
true
> ResetMaximumMemoryUsage();
> time v1:=vetav(p,3,[X0,X1,X2,X3]);
Time: 420.150
> print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
Memory Usage (in MB):  123.593750000000000000000000000
> ResetMaximumMemoryUsage();
> time v2:=vetav2(p,3,[X0,X1,X2,X3]);
Time: 704.390
> print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
Memory Usage (in MB):  140.812500000000000000000000000



*/



// in this one we have first a procedure which will store ALL
// precomputed terms in the "associative array" pre
vetav3p := procedure(p,k,v,~pre : bintab:=[] )
    // the input v is a vector!
    // the output is a vector of resulting
    //    polynomials eta_i(v) for i=1, ... , k


    if (not IsDefined(pre,<k,v>)) and (not IsDefined(pre,<k,vRemoveZeros(v)>)) then 

        // firts, remove zeros!
        v:=vRemoveZeros(v);

        lgt:=#v;
        // print "length = ", lgt;
        if lgt eq 1 then
            P:=Parent(v[1]);
            for i in [1..k] do
                pre[<i,v>]:= P!0;
            end for;
        else

            if #bintab eq 0 then
                // bintab:=BinTab(p,k);
                bintab:=BinTab(p,k);
            end if;

            //print bintab;

            if lgt eq 2 then

                x:=v[1];
                y:=v[2];

                res:=[];

                for j in [1..k] do
                    res[j]:=vRemoveZeros([ bintab[j][i]*x^i*y^(p^j-i) : i in [1..(p^j-1)]]);
                end for;


                for j in [1..(k-1)] do
                    if #res[j] ne 1 then
                        // print "res[",j,"] =", res[j];
                        if not IsDefined(pre,<k-j,res[j]>) then
                            $$(p,k-j,res[j],~pre : bintab:=bintab );
                        end if;
                        for i in [1..(k-j)] do
                            // print "i = ", i;
                            // if we have <k-j,res[j]>, we have all previous!
                            // print res[j];
                            if pre[<i,res[j]>] ne 0 then 
                                Append(~res[j+i],pre[<i,res[j]>]);
                            end if;
                        end for;
                    end if;
                end for;

                for t in [1..k] do
                    if not IsDefined(pre,<t,v>) then
                        pre[<t,v>]:=&+res[t];
                    end if;
                end for;

            else
                // //////////////////////////////////////////////////////////////////////
                // now if #v>2
                // split v in two
                v1:=[ v[i] : i in [1..(lgt div 2)] ];
                v2:=[ v[i] : i in [((lgt div 2)+1)..lgt] ];
                v3:=vRemoveZeros([ &+v1, &+v2 ]);
                
                if (#v1 ne 1) and (not IsDefined(pre,<k,v1>)) then
                    $$(p,k,v1,~pre : bintab:=bintab);
                end if;
                
                if not IsDefined(pre,<k,v2>) then
                    $$(p,k,v2,~pre : bintab:=bintab);
                end if;

                if (#v3 ne 1) and (not IsDefined(pre,<k,v3>)) then
                    $$(p,k,v3,~pre : bintab:=bintab);
                end if;


                // print "keys = ", Keys(pre);
                // print "v1 = ", v1;
                // print "v2 = ", v2;
                // print "v3 = ", v3;
                // print "k = ", k;
                res:=[ ];
                for i in [1..k] do
                    // print "i = ", i;
                    // print res;
                    // print Parent(x1[i]), Parent(x2[i]), Parent(x3[i]);
                    // print x1[i], x2[i], x3[i];
                    if #v1 ne 1 then
                        if #v3 ne 1 then
                            res[i]:=vRemoveZeros([ pre[<i,v1>],pre[<i,v2>],pre[<i,v3>] ]);
                        else
                            res[i]:=vRemoveZeros([ pre[<i,v1>],pre[<i,v2>] ]);
                        end if;
                    else
                        if #v3 ne 1 then
                            res[i]:=vRemoveZeros([ pre[<i,v2>],pre[<i,v3>] ]);
                        else
                            res[i]:=[ pre[<i,v2>] ];
                        end if;
                    end if;
                end for;
                delete v1,v2,v3;
                for i in [2..k] do
                    // lim1:=#res[i-1]; // should be 2^(i-2)*3-(2^(i-2)-1);
                    // // print lim1 eq 2^(i-2)*3-(2^(i-2)-1);
                    // for j in [1..(lim1-1)] do
                    //     temp:=$$(p,k-i+1,[res[i-1][j],&+( res[i-1][(j+1)..lim1] )] : pols:=pols[1..(k-i+1)]);
                    //     for t in [1..(k-i+1)] do
                    //         Append(~res[i-1+t],temp[t]);
                    //     end for;
                    // end for;
                    // delete temp;
                    if not IsDefined(pre,<k-i+1,res[i-1]>) then
                        $$(p,k-i+1,res[i-1],~pre : bintab:=bintab );
                    end if;
                    for t in [i..k] do
                        if pre[<t-i+1,res[i-1]>] ne 0 then
                            Append(~res[t],pre[<t-i+1,res[i-1]>]);
                        end if;
                    end for;
                end for;

                for t in [1..k] do
                    if not IsDefined(pre,<t,v>) then
                        pre[<t,v>]:=&+res[t];
                    end if;
                end for;

            end if;

        end if;
    end if;
end procedure;


function vetav3(p,k,v : bintab:=[] )

    // firts, remove zeros!
    v:=vRemoveZeros(v);

    pre:=AssociativeArray();
    vetav3p(p,k,v,~pre : bintab:=bintab);
    return [ pre[<i,v>] : i in [1..k] ];
end function;




// the new one is faster, but not by much...

vetavold := function(p,k,v : pols:=[ ])
    // the input v is a vector!
    // the output is a vector of resulting
    //    polynomials eta_i(v) for i=1, ... , k

    // firts, remove zeros!
    v:=vRemoveZeros(v);

    lgt:=#v;
    // print "length = ", lgt;
    if lgt eq 1 then
        P:=Parent(v[1]);
        return [ P!0 : i in [1..k] ];
    end if;
    if pols eq [ ] then
        pols:=etapols(p,k);
    end if;
    if lgt eq 2 then
        // //////////////////////////////////////////////////////////////////////
        // now, the special case of #v=2
        // //////////////////////////////////////////////////////////////////////
        return [ Evaluate(pol,v) : pol in pols ]; // needs to have pols precomputed
    else
        // //////////////////////////////////////////////////////////////////////
        // now if #v>2
        // split v in two
        v1:=[ v[i] : i in [1..(lgt div 2)] ];
        v2:=[ v[i] : i in [((lgt div 2)+1)..lgt] ];
        x1:=$$(p,k,v1 : pols:=pols);
        x2:=$$(p,k,v2 : pols:=pols);
        x3:=$$(p,k,[ &+v1, &+v2 ] : pols:=pols);
        res:=[ ];
        for i in [1..k] do
            // print "i = ", i;
            // print res;
            // print Parent(x1[i]), Parent(x2[i]), Parent(x3[i]);
            // print x1[i], x2[i], x3[i];
            res[i]:=[ x1[i],x2[i],x3[i] ];
        end for;
        delete x1,x2,x3,v1,v2;
        for i in [2..k] do
            lim1:=#res[i-1]; // should be 2^(i-2)*3-(2^(i-2)-1);
            // print lim1 eq 2^(i-2)*3-(2^(i-2)-1);
            for j in [1..(lim1-1)] do
                temp:=$$(p,k-i+1,[res[i-1][j],&+( res[i-1][(j+1)..lim1] )] : pols:=pols[1..(k-i+1)]);
                for t in [1..(k-i+1)] do
                    Append(~res[i-1+t],temp[t]);
                end for;
            end for;
            delete temp;
        end for;
        return [ &+x : x in res ];
    end if;
end function;

/*

// test old vs. new

p:=5;
n:=3;
k:=6;

P:=PolynomialRing(GF(p),k);
AssignNames(~P,
  [ "X" cat IntegerToString(i) : i in [1..k] ]);


pols:=etapols(p,n);

time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);
time v2:=vetavold(p, n, [P.i : i in [1..k] ] : pols:=pols);

v1 eq v2;



> p:=3;
> n:=3;
> k:=6;
>
> P:=PolynomialRing(GF(p),k);
> AssignNames(~P,
>   [ "X" cat IntegerToString(i) : i in [1..k] ]);
>
>
> pols:=etapols(p,n);
>
> time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 4.760
> time v2:=vetavold(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 4.870
>
> v1 eq v2;
true


> p:=2;
> n:=5;
> k:=6;
>
> P:=PolynomialRing(GF(p),k);
> AssignNames(~P,
>   [ "X" cat IntegerToString(i) : i in [1..k] ]);
>
>
> pols:=etapols(p,n);
>
> time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 65.060
> time v2:=vetavold(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 75.050
>
> v1 eq v2;
true



> p:=2;
> n:=6;
> k:=4;
>
> P:=PolynomialRing(GF(p),k);
> AssignNames(~P,
>   [ "X" cat IntegerToString(i) : i in [1..k] ]);
>
>
> pols:=etapols(p,n);
>
> time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 16.900
> time v2:=vetavold(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 21.080
>
> v1 eq v2;
true


> p:=2;
> n:=4;
> k:=8;
>
> P:=PolynomialRing(GF(p),k);
> AssignNames(~P,
>   [ "X" cat IntegerToString(i) : i in [1..k] ]);
>
>
> pols:=etapols(p,n);
>
> time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 4.290
> time v2:=vetavold(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 4.330
>
> v1 eq v2;
true


> load 'etas.m';
Loading "etas.m"
>
> p:=3;
> n:=4;
> k:=4;
>
> P:=PolynomialRing(GF(p),k);
> AssignNames(~P,
>   [ "X" cat IntegerToString(i) : i in [1..k] ]);
>
>
> pols:=etapols(p,n);
>
> time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 81.420
> time v2:=vetavold(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 107.840
>
> v1 eq v2;
true



> load 'etas.m';
Loading "etas.m"
>
> p:=2;
> n:=4;
> k:=4;
>
> P:=PolynomialRing(GF(p),k);
> AssignNames(~P,
>   [ "X" cat IntegerToString(i) : i in [1..k] ]);
>
>
> pols:=etapols(p,n);
>
> time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 0.010
> time v2:=vetavold(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 0.010
>
> v1 eq v2;
true
>
> p:=2;
> n:=4;
> k:=5;
>
> P:=PolynomialRing(GF(p),k);
> AssignNames(~P,
>   [ "X" cat IntegerToString(i) : i in [1..k] ]);
>
>
> pols:=etapols(p,n);
>
> time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 0.040
> time v2:=vetavold(p, n, [P.i : i in [1..k] ] : pols:=pols);
Time: 0.020
>
> v1 eq v2;
true

> p:=5;
> n:=3;
> k:=6;
>
> P:=PolynomialRing(GF(p),k);
> AssignNames(~P,
>   [ "X" cat IntegerToString(i) : i in [1..k] ]);
>
>
> pols:=etapols(p,n);
>
> time v1:=vetav(p, n, [P.i : i in [1..k] ] : pols:=pols);

Current total memory usage: 27356.8MB, failed memory request: 22006.0MB
System error: Out of memory.
All virtual memory has been exhausted so Magma cannot perform this statement.
[Interrupt twice in half a second; exiting]

Total time: 3244.510 seconds, Total memory usage: 27356.75MB


*/


/* not tested!!!!  (may not be necessary!)

vvetav := function(p,k,v : pols:=[])
    // the input v is a vector!
    // the output is a vector of vectors, which sum to
    //    the polynomials eta_i(v) for i=1, ... , k
    lgt:=#v;
    // print "length = ", lgt;
    if lgt eq 1 then
        P:=Parent(v[1]);
        return [ [] : i in [1..k] ];
    end if;
    if pols eq [] then
        pols:=etapols(p,k);
    end if;
    if lgt eq 2 then
        // //////////////////////////////////////////////////////////////////////
        // now, the special case of #v=2
        // //////////////////////////////////////////////////////////////////////
        return [ [Evaluate(pol,v)] : pol in pols ]; // needs to have pols precomputed
    else
        // //////////////////////////////////////////////////////////////////////
        // now if #v>2
        // split v in two
        v1:=[ v[i] : i in [1..(lgt div 2)] ];
        v2:=[ v[i] : i in [((lgt div 2)+1)..lgt] ];
        x1:=$$(p,k,v1 : pols:=pols);
        x2:=$$(p,k,v2 : pols:=pols);
        x3:=$$(p,k,[ &+v1, &+v2 ] : pols:=pols);
        res:=[ ];
        for i in [1..k] do
            // print "i = ", i;
            // print res;
            // print Parent(x1[i]), Parent(x2[i]), Parent(x3[i]);
            // print x1[i], x2[i], x3[i];
            res[i]:=[ x1[i] cat x2[i]  cat x3[i] ];
        end for;
        delete x1,x2,x3,v1,v2;
        for i in [2..k] do
            lim1:=#res[i-1]; // should be 2^(i-2)*3-(2^(i-2)-1);
            // print lim1 eq 2^(i-2)*3-(2^(i-2)-1);
            for j in [1..(lim1-1)] do
                temp:=$$(p,k-i+1,[res[i-1][j],&+( res[i-1][(j+1)..lim1] )] : pols:=pols);
                for t in [1..(k-i+1)] do
                    Append(~res[i-1+t],temp[t]);
                end for;
            end for;
            delete temp;
        end for;
        return res;
    end if;
end function;

*/





// raise terms to p^r
pterms := function(f,p,r)
    return &+[ x^(p^r) : x in Terms(f) ];
end function;


// eta of a polynonial done in char p^k (for testing purposes!)
vetav_t := function(p,k,f)
    P:=Parent(f); // needs to be a pol. ring (over ZZ)
    PP:=ChangeRing(P,ResidueClassRing(p^(k+1)));
    PPP:=ChangeRing(P,GF(p));
    res:=[ (pterms(f,p,1) - f^p) div p ];
    for i in [2..k] do
        x:= pterms(f,p,i) - f^(p^i);
        x-:= &+[ p^j*(res[j])^(p^(i-j)) : j in [1..(i-1)] ];
        x:=x div p^i;
        Append(~res,x);
    end for;
    return [ PPP!x : x in res];
end function;



/*
// Testing!

p:=2;
P<X,Y,Z,W,T>:=PolynomialRing(Integers(),5);
PP<x,y,z,w,t>:=PolynomialRing(GF(p),5);
f:=X+Y+Z+T;
k:=6;
time v1:=[ PP!t : t in vetav_t(p,k,f) ];
time v2:=vetav(p,k,[ PP!t : t in Terms(f) ]);
v1 eq v2;


// example:

> p:=3;
> f:=X+Y+Z+W+T;
> k:=3;
> time v1:=[ PP!t : t in vetav_t(p,k,f) ];
Time: 2.240
> time v2:=vetav(p,k,[ PP!t : t in Terms(f) ]);
Time: 0.440
> v1 eq v2;
true
> k:=4;
> time v1:=[ PP!t : t in vetav_t(p,k,f) ];
Time: 31358.080 (8.71057777777778 hours)
> time v2:=vetav(p,k,[ PP!t : t in Terms(f) ]);
Time: 6326.320 (1.75731111111111 hours)
> v1 eq v2;
true


> p:=2;
> P<X,Y,Z,W,T>:=PolynomialRing(Integers(),5);
> PP<x,y,z,w,t>:=PolynomialRing(GF(p),5);
> f:=X+Y+Z+T+W;
> k:=4;
> time v1:=[ PP!t : t in vetav_t(p,k,f) ];
Time: 0.110
> time v2:=vetav(p,k,[ PP!t : t in Terms(f) ]);
Time: 0.040
> v1 eq v2;
true
> k:=5;
> time v1:=[ PP!t : t in vetav_t(p,k,f) ];
Time: 17.120
> time v2:=vetav(p,k,[ PP!t : t in Terms(f) ]);
Time: 3.790
> v1 eq v2;
true
> f:=X+Y+Z+T;
> k:=6;
> time v1:=[ PP!t : t in vetav_t(p,k,f) ];
Time: 52.220
> time v2:=vetav(p,k,[ PP!t : t in Terms(f) ]);
Time: 21.900
> v1 eq v2;
true



*/


function newWittSum(v1,v2 : pols:=[])
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    n:=#v1-1;
    // initialize
    res:=[ vRemoveZeros([v1[i], v2[i]]) : i in [1..(n+1)] ];
    // print res;
    // precompute the nec. pol.
    if #pols eq 0 then
        // print "Computing etapols";
        pols:=etapols(p,n);
    end if;

    if #pols gt n then
        pols:=pols[1..n];
    end if;

    for i in [2..(n+1)] do
        // print "iteration ", i, " of ", n+1;
        l:=#res[i-1];
        for j in [2..l] do
            // print "    iteration", j, " of ", l;
            // add the terms coming from res[i-1]
            temp:=vetav(p,n+2-i,[res[i-1][j],&+(res[i-1][1..(j-1)])] : pols:=pols); // note that here it just eval. the pols.
            for t in [i..(n+1)] do
                if temp[t-i+1] ne 0 then 
                    Append(~res[t],temp[t-i+1]);
                end if;
            end for;
        end for;
        // res[i] is done!
        Remove(~pols,n+2-i); // won't need the last anymore..
        // print res;
    end for;
    return [ &+x : x in res ];
end function;


function newWittSum2(v1,v2 : bintab:=[])
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    n:=#v1-1;
    // initialize
    res:=[ vRemoveZeros([v1[i], v2[i]]) : i in [1..(n+1)] ];
    // print res;
    // precompute the nec. pol.
    if #bintab eq 0 then
        // print "Computing bintab";
        bintab:=BinTab(p,n);
    end if;

    for i in [2..(n+1)] do
        // print "iteration ", i, " of ", n+1;
        l:=#res[i-1];
        for j in [2..l] do
            // print "    iteration", j, " of ", l;
            // add the terms coming from res[i-1]
            temp:=vetav2(p,n+2-i,[res[i-1][j],&+(res[i-1][1..(j-1)])] : bintab:=bintab); // note that here it just eval. the pols.
            for t in [i..(n+1)] do
                if temp[t-i+1] ne 0 then 
                    Append(~res[t],temp[t-i+1]);
                end if;
            end for;
        end for;
        // res[i] is done!
        // print res;
    end for;
    return [ &+x : x in res ];
end function;



function newWittSum3(v1,v2 : bintab:=[])
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    n:=#v1-1;
    // initialize
    res:=[ vRemoveZeros([v1[i], v2[i]]) : i in [1..(n+1)] ];
    // print res;
    // precompute the nec. pol.
    if #bintab eq 0 then
        // print "Computing bintab";
        bintab:=BinTab(p,n);
    end if;

    for i in [2..(n+1)] do
        // print "iteration ", i, " of ", n+1;
        l:=#res[i-1];
        for j in [2..l] do
            // print "    iteration", j, " of ", l;
            // add the terms coming from res[i-1]
            temp:=vetav3(p,n+2-i,[res[i-1][j],&+(res[i-1][1..(j-1)])] : bintab:=bintab); // note that here it just eval. the pols.
            for t in [i..(n+1)] do
                if temp[t-i+1] ne 0 then 
                    Append(~res[t],temp[t-i+1]);
                end if;
            end for;
        end for;
        // res[i] is done!
        // print res;
    end for;
    return [ &+x : x in res ];
end function;






/* TEST newWittSum2

> load 'etas.m';
Loading "etas.m"
Loading "witt_gen.m"
> p:=3;PR<X0,X1,X2,X3,X4,Y0,Y1,Y2,Y3,Y4>:=PolynomialRing(GF(p),10);
> ResetMaximumMemoryUsage(); time w1:=newWittSum([X0,X1,X2],[Y0,Y1,Y2]); print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
Computing etapols
iteration  2  of  3
    iteration 2  of  2
iteration  3  of  3
    iteration 2  of  3
    iteration 3  of  3
Time: 0.000
Memory Usage (in MB):  12.7500000000000000000000000000
> ResetMaximumMemoryUsage(); time w2:=newWittSum2([X0,X1,X2],[Y0,Y1,Y2]); print "Memory \
Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
iteration  2  of  3
    iteration 2  of  2
iteration  3  of  3
    iteration 2  of  3
    iteration 3  of  3
Time: 0.000
Memory Usage (in MB):  12.7500000000000000000000000000
> w1 eq w2;
true
> ResetMaximumMemoryUsage(); time w1:=newWittSum([X0,X1,X2,X3],[Y0,Y1,Y2,Y3]); print "Me\
mory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
Computing etapols
iteration  2  of  4
    iteration 2  of  2
iteration  3  of  4
    iteration 2  of  3
    iteration 3  of  3
iteration  4  of  4
    iteration 2  of  5
    iteration 3  of  5
    iteration 4  of  5
    iteration 5  of  5
Time: 0.000
Memory Usage (in MB):  12.7500000000000000000000000000

> ResetMaximumMemoryUsage(); time w2:=newWittSum2([X0,X1,X2,X3],[Y0,Y1,Y2,Y3]); print "M\
emory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
iteration  2  of  4
    iteration 2  of  2
iteration  3  of  4
    iteration 2  of  3
    iteration 3  of  3
iteration  4  of  4
    iteration 2  of  5
    iteration 3  of  5
    iteration 4  of  5
    iteration 5  of  5
Time: 0.000
Memory Usage (in MB):  12.7500000000000000000000000000
> w1 eq w2;
true
> ResetMaximumMemoryUsage(); time w1:=newWittSum([X0,X1,X2,X3,X4],[Y0,Y1,Y2,Y3,Y4]); pri\
nt "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
Computing etapols
iteration  2  of  5
    iteration 2  of  2
iteration  3  of  5
    iteration 2  of  3
    iteration 3  of  3
iteration  4  of  5
    iteration 2  of  5
    iteration 3  of  5
    iteration 4  of  5
    iteration 5  of  5
iteration  5  of  5
    iteration 2  of  9
    iteration 3  of  9
    iteration 4  of  9
    iteration 5  of  9
    iteration 6  of  9
    iteration 7  of  9
    iteration 8  of  9
    iteration 9  of  9
Time: 0.700
Memory Usage (in MB):  18.3437500000000000000000000000
> ResetMaximumMemoryUsage(); time w2:=newWittSum2([X0,X1,X2,X3,X4],[Y0,Y1,Y2,Y3,Y4]); pr\
int "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
iteration  2  of  5
    iteration 2  of  2
iteration  3  of  5
    iteration 2  of  3
    iteration 3  of  3
iteration  4  of  5
    iteration 2  of  5
    iteration 3  of  5
    iteration 4  of  5
    iteration 5  of  5
iteration  5  of  5
    iteration 2  of  9
    iteration 3  of  9
    iteration 4  of  9
    iteration 5  of  9
    iteration 6  of  9
    iteration 7  of  9
    iteration 8  of  9
    iteration 9  of  9
Time: 0.860
Memory Usage (in MB):  21.0937500000000000000000000000
> w1 eq w2;
true

*/





/*

// testing

load 'witt_gen.m';
p:=17;
n:=3;

P:=PolynomialRing(GF(p),(2*n+2));
AssignNames(~P,
  [ "x" cat IntegerToString(i) : i in [0..n] ] cat
  [ "y" cat IntegerToString(i) : i in [0..n] ] );

// time v1:=newWittSumPol(n,p);
time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);

[ P!v1[i] : i in [1..(n+1)] ] eq v2;






> p:=3;
> n:=2;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.000
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.010
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=5;
> n:=2;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.010
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.000
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=7;
> n:=2;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.000
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.000
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=11;
> n:=2;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.030
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.010
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=2;
> n:=3;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.000
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.010
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=2;
> n:=4;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.000
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.010
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=2;
> n:=5;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.050
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.040
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=2;
> n:=6;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 33.920
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 19.470
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=3;
> n:=3;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.010
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.060
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=3;
> n:=4;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 1.760
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.540
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
> p:=5;
> n:=3;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 0.590
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 0.210
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=7;
> n:=3;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Time: 92.570
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Time: 33.160
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true




///////////////////////////////////////////////////////////////////////

Connection broken with boole after ~46 hour, still incomplete:

> p:=3;
> n:=5;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);

oole(scrn)[~/tmp]$ magma
Magma V2.16-6     Fri Mar 26 2010 15:43:53 on boole    [Seed = 1522278204]
Type ? for help.  Type <Ctrl>-D to quit.
>
> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=3;
> n:=5;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);

>> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
            ^
User error: Identifier 'newWittSum' has not been declared or assigned
> load 'etas.m';
Loading "etas.m"
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Computing etapols
iteration  2  of  6
    iteration 2  of  2
iteration  3  of  6
    iteration 2  of  3
    iteration 3  of  3
iteration  4  of  6
    iteration 2  of  5
    iteration 3  of  5
    iteration 4  of  5
    iteration 5  of  5
iteration  5  of  6
    iteration 2  of  9
    iteration 3  of  9
    iteration 4  of  9
    iteration 5  of  9
    iteration 6  of  9
    iteration 7  of  9
    iteration 8  of  9
    iteration 9  of  9
iteration  6  of  6
    iteration 2  of  17
    iteration 3  of  17
    iteration 4  of  17
    iteration 5  of  17
    iteration 6  of  17
    iteration 7  of  17
    iteration 8  of  17
    iteration 9  of  17
    iteration 10  of  17
    iteration 11  of  17
    iteration 12  of  17
    iteration 13  of  17
    iteration 14  of  17
    iteration 15  of  17
    iteration 16  of  17
    iteration 17  of  17

Current total memory usage: 13003.8MB, failed memory request: 1339.0MB
System error: Out of memory.
All virtual memory has been exhausted so Magma cannot perform this statement.

[running with others...]









> p:=5;
> n:=4;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=newWittSumPol(n,p);
Shared connection to boole closed.

> p:=11;
> n:=3;
> time v1:=WittSumPol(n,p);
Shared connection to boole closed.

(NOT newWittSumPol in the one above!!!!)




Magma V2.16-11    Sat Aug  7 2010 08:07:06 on boole    [Seed = 1180124747]
Type ? for help.  Type <Ctrl>-D to quit.
> load 'etas.m';
Loading "etas.m"
Loading "witt_gen.m"
> p:=3;F:=GF(p);P<X0,X1,X2,X3,X4,X5,Y0,Y1,Y2,Y3,Y4,Y5>:=PolynomialRing(F,12);
> ResetMaximumMemoryUsage(); time w2:=newWittSum2([X0,X1,X2,X3,X4,X5],[Y0,Y1,Y2,Y3,Y4,Y5]); print "Memory Usage\
 (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
Time: 133883.290
Memory Usage (in MB):  14574.5625000000000000000000000

> ResetMaximumMemoryUsage(); time w1:=newWittSum([X0,X1,X2,X3,X4,X5],[Y0,Y1,Y2,Y3,Y4,Y5]); print "Memory Usage \
(in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
Time: 123250.260
Memory Usage (in MB):  16652.2500000000000000000000000
> w1 eq w2;
true





> p:=11;
> n:=3;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Computing etapols
iteration  2  of  4
    iteration 2  of  2
iteration  3  of  4
    iteration 2  of  3
    iteration 3  of  3
iteration  4  of  4
    iteration 2  of  5
    iteration 3  of  5
    iteration 4  of  5
    iteration 5  of  5
Time: 25912.350
> 25912.350/360;
71.9787500000000000000000000000
> 25912.350/3600;
7.19787500000000000000000000000
> time v1:=newWittSumPol(n,p);
iteration  2  of  4
    iteration 2  of  3
    iteration 3  of  3
    iteration 4  of  3
iteration  3  of  4
    iteration 3  of  2
    iteration 4  of  2
iteration  4  of  4
    iteration 4  of  1
Time: 470015.360
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true
> 470015.360/3600;
130.559822222222222222222222222
> 470015.360/(3600*24);
5.43999259259259259259259259260
>
> [ #Terms(x) : x in v2 ];
[ 2, 12, 2368, 25819940 ]



> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=13;
> n:=3;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> // time v1:=newWittSumPol(n,p);
> time v2:=newWittSum([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Computing etapols
iteration  2  of  4
    iteration 2  of  2
iteration  3  of  4
    iteration 2  of  3
    iteration 3  of  3
iteration  4  of  4
    iteration 2  of  5
    iteration 3  of  5
    iteration 4  of  5
    iteration 5  of  5

Current total memory usage: 8025.2MB, failed memory request: 4341.7MB
System error: Out of memory.
All virtual memory has been exhausted so Magma cannot perform this statement.









> load 'etas.m';
Loading "etas.m"
Loading "witt_gen.m"
> p:=3;F:=GF(p);P<X0,X1,X2,X3,X4,X5,Y0,Y1,Y2,Y3,Y4,Y5>:=PolynomialRing(F,12);
> ResetMaximumMemoryUsage(); time w2:=newWittSum2([X0,X1,X2,X3,X4,X5],[Y0,Y1,Y2,Y3,Y4,Y5\
\
time> ]); print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);
iteration  2  of  6
    iteration 2  of  2
iteration  3  of  6
    iteration 2  of  3
    iteration 3  of  3
iteration  4  of  6
    iteration 2  of  5
    iteration 3  of  5
    iteration 4  of  5
    iteration 5  of  5
iteration  5  of  6
    iteration 2  of  9
    iteration 3  of  9
    iteration 4  of  9
    iteration 5  of  9
    iteration 6  of  9
    iteration 7  of  9
    iteration 8  of  9
    iteration 9  of  9
iteration  6  of  6
    iteration 2  of  17
    iteration 3  of  17
    iteration 4  of  17
    iteration 5  of  17
    iteration 6  of  17
    iteration 7  of  17
    iteration 8  of  17
    iteration 9  of  17
    iteration 10  of  17
    iteration 11  of  17
    iteration 12  of  17
    iteration 13  of  17
    iteration 14  of  17
    iteration 15  of  17
    iteration 16  of  17
    iteration 17  of  17
Time: 133848.410
Memory Usage (in MB):  14533.9062500000000000000000000




*/


/*

// Now let's test adding over finite fileds.

p:=11;
k:=10;
n:=5;

F<a>:=GF(p^k);

ntest:=20;

time pols:=etapols(p,n);

t:=Cputime();
for i in [1..ntest] do
    v1:=[ Random(F) : i in [0..n] ];
    v2:=[ Random(F) : i in [0..n] ];
    time x:=newWittSum(v1,v2 : pols:=pols);
end for;

print "average time: "; Cputime(t)/ntest;

//////////////////

 load 'etas.m';
Loading "etas.m"
Loading "witt_gen.m"
>
> p:=11;
> k:=10;
> n:=5;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time pols:=etapols(p,n);
Time: 5353.500
>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum(v1,v2 : pols:=pols);
for> end for;
>
> print "average time: "; Cputime(t)/ntest;
average time:
1.480
Total memory usage: 61.31MB






// now with newWittSum2
p:=11;
k:=10;
n:=5;

F<a>:=GF(p^k);

ntest:=20;

time bintab:=BinTab2(p,n);

t:=Cputime();
for i in [1..ntest] do
    v1:=[ Random(F) : i in [0..n] ];
    v2:=[ Random(F) : i in [0..n] ];
    time x:=newWittSum2(v1,v2 : bintab:=bintab);
end for;

print "average time: "; Cputime(t)/ntest;



> load 'etas.m';
Loading "etas.m"
Loading "witt_gen.m"
> p:=11;
> k:=10;
> n:=5;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time bintab:=BinTab2(p,n);
Time: 5840.940
>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum2(v1,v2 : bintab:=bintab);
for> end for;
Time: 28.230
Time: 28.540
Time: 26.020
Time: 26.370
Time: 25.870
Time: 25.780
Time: 25.760
Time: 25.740
Time: 25.970
Time: 27.020
Time: 28.020
Time: 27.880
Time: 25.920
Time: 26.320
Time: 25.890
Time: 26.120
Time: 26.180
Time: 26.140
Time: 28.190
Time: 26.490
>
> print "average time: "; Cputime(t)/ntest;
average time:


// Try BinTab3

// now with newWittSum2
p:=11;
k:=10;
n:=5;

F<a>:=GF(p^k);

ntest:=20;

time bintab:=BinTab3(p,n);

t:=Cputime();
for i in [1..ntest] do
    v1:=[ Random(F) : i in [0..n] ];
    v2:=[ Random(F) : i in [0..n] ];
    time x:=newWittSum2(v1,v2 : bintab:=bintab);
end for;

print "average time: "; Cputime(t)/ntest;

> p:=11;
> k:=10;
> n:=5;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time bintab:=BinTab3(p,n);
Time: 5.750
>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum2(v1,v2 : bintab:=bintab);
for> end for;
Time: 26.300
Time: 26.580
Time: 25.780
Time: 25.720
Time: 25.650
Time: 25.920
Time: 26.010
Time: 26.050
Time: 25.960
Time: 26.300
Time: 26.810
Time: 26.290
Time: 26.030
Time: 25.830
Time: 25.900
Time: 25.860
Time: 25.890
Time: 25.830
Time: 25.870
Time: 25.930
>
> print "average time: "; Cputime(t)/ntest;
average time:
26.026




// now with newWittSum3
p:=11;
k:=10;
n:=5;

F<a>:=GF(p^k);

ntest:=20;

time bintab:=BinTab(p,n);

t:=Cputime();
for i in [1..ntest] do
    v1:=[ Random(F) : i in [0..n] ];
    v2:=[ Random(F) : i in [0..n] ];
    time x:=newWittSum3(v1,v2 : bintab:=bintab);
end for;

print "average time: "; Cputime(t)/ntest;











p:=11; k:=30; n:=5; F<a>:=GF(p^k)
time bintab:=BinTab(p,n);
v1:=[ Random(F) : i in [0..n] ];
v2:=[ Random(F) : i in [0..n] ];
time x:=newWittSum2(v1,v2 : bintab:=bintab);
time y:=newWittSum3(v1,v2 : bintab:=bintab);
x eq y;

Time: 63.260 
Time: 68.800
true




*/


/*

Markov:


> p:=5;
> k:=10;
> n:=6;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time pols:=etapols(p,n);
Time: 59.630
>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum(v1,v2 : pols:=pols);
for> end for;
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.110
Time: 0.100
Time: 0.100
Time: 0.100
Time: 0.100
>
> print "average time: "; Cputime(t)/ntest;
average time:
0.100





> p:=5;
> k:=10;
> n:=7;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time pols:=etapols(p,n);
Time: 3964.380
>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum(v1,v2 : pols:=pols);
for> end for;
Time: 0.600
Time: 0.590
Time: 0.590
Time: 0.590
Time: 0.590
Time: 0.590
Time: 0.600
Time: 0.590
Time: 0.590
Time: 0.590
Time: 0.590
Time: 0.600
Time: 0.590
>
> print "average time: "; Cputime(t)/ntest;
average time:
0.592



> p:=5;
> k:=10;
> n:=8;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time pols:=etapols(p,n);
^ZTime: 128884.020

[1]+  Stopped                 magma
markov(scrn)[~/tmp]$ fg
magma
>
> pols[1];
4*x^4*y + 3*x^3*y^2 + 3*x^2*y^3 + 4*x*y^4
> 128884.020/3600
> ;
35.8011166666666666666666666666

>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum(v1,v2 : pols:=pols);
for> end for;
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.350
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
Time: 3.410
>
> print "average time: "; Cputime(t)/ntest;
average time:
3.408





> p:=7;
> k:=10;
> n:=5;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time pols:=etapols(p,n);
Time: 77.870
>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum(v1,v2 : pols:=pols);
for> end for;
Time: 0.090
Time: 0.100
Time: 0.100
Time: 0.090
Time: 0.100
Time: 0.090
Time: 0.100
Time: 0.100
Time: 0.090
Time: 0.100
Time: 0.100
Time: 0.090
Time: 0.100
Time: 0.090
Time: 0.100
Time: 0.100
Time: 0.090
Time: 0.100
Time: 0.100
Time: 0.090
>
> print "average time: "; Cputime(t)/ntest;
average time:
0.097


p=7,k=10,n=6
> time pols:=etapols(p,n);
Time: 8403.560
(....)
> print "average time: "; Cputime(t)/ntest;
average time:
0.860
> for x in pols do
for> print #Terms(x);
for> end for;
6
36
288
2084
14490
100912




> p:=11;
> k:=10;
> n:=4;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time pols:=etapols(p,n);
Time: 49.380
>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum(v1,v2 : pols:=pols);
for> end for;
Time: 0.080
Time: 0.070
Time: 0.080
Time: 0.070
Time: 0.080
Time: 0.080
Time: 0.080
Time: 0.070
Time: 0.080
Time: 0.080
Time: 0.070
Time: 0.080
Time: 0.080
Time: 0.080
Time: 0.070
Time: 0.080
Time: 0.080
Time: 0.070
Time: 0.080
Time: 0.080
>
> print "average time: "; Cputime(t)/ntest;
average time:
0.077

> p:=11;
> k:=10;
> n:=5;
>
> F<a>:=GF(p^k);
>
> ntest:=20;
>
> time pols:=etapols(p,n);
Time: 12993.010
>
> t:=Cputime();
> for i in [1..ntest] do
for>     v1:=[ Random(F) : i in [0..n] ];
for>     v2:=[ Random(F) : i in [0..n] ];
for>     time x:=newWittSum(v1,v2 : pols:=pols);
for> end for;

(...)
> print "average time: "; Cputime(t)/ntest;
average time:
1.064

*/



// Run over extension of Zp

// convert a Witt vector over a finite field to p-adic power series
function WittVToSeries(v : Zq:=0)

    FF:=Parent(v[1]);
    p:=Characteristic(FF);
    n:=#v-1;

    // if ring is not given, create it
    if Zq cmpeq pAdicRing(2) then
        k:=Degree(FF);
        Zp:=pAdicRing(p : Precision:=n+1);
        Zq<aa>:= ext<Zp | k>;
    end if;

    return &+[ Zq!TeichmuellerLift(Root(v[i],p^(i-1)),quo<Zq | p^(n+2-i)>)*p^(i-1) : i in [1..(n+1)] ];

end function;


// convert p-adic power series to a Witt vector
function SeriesToWittV(s : F:=0)

    Zq:=Parent(s);
    n:=Zq`DefaultPrecision-1;
    k:=Degree(Zq);
    p:=Prime(Zq);

    // if finite field not given, create it
    if F cmpeq 0 then
        F<a>:=GF(p^k);
    end if;

    v:=[];

    a:=s;
    for i in [0..n] do
        t:=F!a;
        Append(~v,t);
        tt:= Zq!TeichmuellerLift(t,quo<Zq | p^(n+1-i)>);
        a:=(a - tt) div p;
    end for;
    return [ v[i+1]^(p^i) : i in [0..n] ];

end function;

/*

// test!

p:=7;
k:=20;
n:=7;

F<a>:=GF(p^k);
Zp:=pAdicRing(p);
Zq<aa>:=ext< Zp | k >;
Zq`DefaultPrecision:=n+1;

ntest:=20;

for j in [1..ntest] do
    v:=[ Random(F) : i in [0..n] ];
    s:=WittVToSeries(v : Zq:=Zq);
    SeriesToWittV(s : F:=F) eq v;
end for;

for j in [1..ntest] do
    s:=Random(Zq);
    v:=SeriesToWittV(s : F:=F);
    WittVToSeries(v : Zq:=Zq) eq s;
end for;

*/


/*
p:=5;
k:=50;
n:=5;

F<a>:=GF(p^k);
Zp:=pAdicRing(p);
Zq<aa>:= ext<Zp | k>;
Zq`DefaultPrecision:=n+1;

ntest:=20;

tt1:=0;
tt2:=0;

t:=Cputime();
pols:=etapols(p,n);
tt3:=Cputime(t);
print "Time 3 = ", tt3;


for i in [1..ntest] do
    v1:=[ Random(F) : i in [0..n] ];
    v2:=[ Random(F) : i in [0..n] ];

    t:=Cputime();
    // convert to Zq
    x1:=WittVToSeries(v1 : Zq:=Zq);
    x2:=WittVToSeries(v2 : Zq:=Zq);

    x1:=SeriesToWittV(x1 + x2 : F:=F);

    tt:=Cputime(t);

    print "Time 1 = ", tt;
    tt1+:=tt;

    t:=Cputime();
    x2:=newWittSum(v1,v2 : pols:=pols);
    tt:=Cputime(t);

    print "Time 2 = ", tt;
    tt2+:=tt;

    x1 eq x2;

end for;

print "average time 1: "; tt1/ntest;
print "average time 2: "; tt2/ntest;
print "time 3: ", tt3;

*/







function newWittProd(v1,v2 : pols:=[])
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    n:=#v1-1;
    // initialize
    res:=[ vRemoveZeros([ v1[j+1]^(p^(i-j))*v2[i-j+1]^(p^j) : j in [0..i] ]) : i in [0..n] ];

    // print "initial res = ", res;

    if n eq 1 then
        // done
        return [ &+x : x in res];
    end if;
    // print res;
    // precompute the nec. pol.
    if #pols eq 0 then
        print "Computing etapols";
        pols:=etapols(p,n-1); // NOTE: n-1 only!
    end if;

    if #pols gt (n-1) then
        pols:=pols[1..(n-1)];
    end if;

    for i in [3..(n+1)] do
        // print "iteration ", i, " of ", n+1;
        l:=#res[i-1];
        for j in [2..l] do
            // print "    iteration", j, " of ", l;
            // add the terms coming from res[i-1]
            temp:=vetav(p,n+2-i,[res[i-1][j],&+(res[i-1][1..(j-1)])] : pols:=pols); // note that here it just eval. the pols.
            // print "temp = ", temp;
            for t in [i..(n+1)] do
                if temp[t-i+1] ne 0 then
                    Append(~res[t],temp[t-i+1]);
                end if;
            end for;
        end for;
        // res[i] is done!
        Remove(~pols,n+2-i); // won't need the last anymore..
        // print res;
    end for;
    return [ &+x : x in res ];
end function;



function newWittProd2(v1,v2 : bintab:=[])
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    n:=#v1-1;
    // initialize
    res:=[ vRemoveZeros([ v1[j+1]^(p^(i-j))*v2[i-j+1]^(p^j) : j in [0..i] ]) : i in [0..n] ];

    // print "initial res = ", res;

    if n eq 1 then
        // done
        return [ &+x : x in res];
    end if;
    // print res;
    // precompute the nec. binomials
    if #bintab eq 0 then
        // print "Computing bintab";
        bintab:=BinTab(p,n-1);  // note only need n-1
    end if;

    for i in [3..(n+1)] do
        // print "iteration ", i, " of ", n+1;
        l:=#res[i-1];
        for j in [2..l] do
            // print "    iteration", j, " of ", l;
            // add the terms coming from res[i-1]
            temp:=vetav2(p,n+2-i,[res[i-1][j],&+(res[i-1][1..(j-1)])] : bintab:=bintab); // note that here it just eval. the pols.
            // print "temp = ", temp;
            for t in [i..(n+1)] do
                if temp[t-i+1] ne 0 then
                    Append(~res[t],temp[t-i+1]);
                end if;
            end for;
        end for;
        // res[i] is done!

        // print res;
    end for;
    return [ &+x : x in res ];
end function;


function newWittProd3(v1,v2 : bintab:=[])
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    n:=#v1-1;
    // initialize
    res:=[ vRemoveZeros([ v1[j+1]^(p^(i-j))*v2[i-j+1]^(p^j) : j in [0..i] ]) : i in [0..n] ];

    // print "initial res = ", res;

    if n eq 1 then
        // done
        return [ &+x : x in res];
    end if;
    // print res;
    // precompute the nec. binomials
    if #bintab eq 0 then
        // print "Computing bintab";
        bintab:=BinTab(p,n-1);  // note only need n-1
    end if;

    for i in [3..(n+1)] do
        // print "iteration ", i, " of ", n+1;
        l:=#res[i-1];
        for j in [2..l] do
            // print "    iteration", j, " of ", l;
            // add the terms coming from res[i-1]
            temp:=vetav3(p,n+2-i,[res[i-1][j],&+(res[i-1][1..(j-1)])] : bintab:=bintab); // note that here it just eval. the pols.
            // print "temp = ", temp;
            for t in [i..(n+1)] do
                if temp[t-i+1] ne 0 then
                    Append(~res[t],temp[t-i+1]);
                end if;
            end for;
        end for;
        // res[i] is done!

        // print res;
    end for;
    return [ &+x : x in res ];
end function;






/*

// testing

load 'witt_gen.m';
p:=7;
n:=4;

P:=PolynomialRing(GF(p),(2*n+2));
AssignNames(~P,
  [ "x" cat IntegerToString(i) : i in [0..n] ] cat
  [ "y" cat IntegerToString(i) : i in [0..n] ] );

time v1:=WittProdPol(n,p);
time v2:=newWittProd([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);

[ P!v1[i] : i in [1..(n+1)] ] eq v2;




> load 'witt_gen.m';
Loading "witt_gen.m"
> p:=7;
> n:=4;
>
> P:=PolynomialRing(GF(p),(2*n+2));
> AssignNames(~P,
>   [ "x" cat IntegerToString(i) : i in [0..n] ] cat
>   [ "y" cat IntegerToString(i) : i in [0..n] ] );
>
> time v1:=WittProdPol(n,p);
Time: 10850.030
> time v2:=newWittProd([ P.i : i in [1..(n+1)] ], [P.i : i in [(n+2)..(2*n+2)] ]);
Computing etapols
Time: 604.780
>
> [ P!v1[i] : i in [1..(n+1)] ] eq v2;
true



*/



/*
// greenberg transform

p:=3;
n:=4;
// k:=4;

// F<a>:=GF(p^k);

F:=GF(p);

P:=PolynomialRing(F,(2*n+2));
AssignNames(~P,
  [ "x" cat IntegerToString(i) : i in [0..n] ] cat
  [ "y" cat IntegerToString(i) : i in [0..n] ] );

x:=[ P.i : i in [1..(n+1)] ];
y:=[ P.i : i in [(n+2)..(2*n+2)] ];

a:=[F!2,F!1,F!0,F!1,F!1];

pols:=etapols(p,n);
res1:=newWittProd(x,x : pols:=pols);
res1:=newWittProd(x,res1 : pols:=pols);
res1:=newWittProd(a,res1 : pols:=pols);

res1;

*/

// not efficient -- just to test!
function GTaux1(v1,v2)
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    res:=[];
    n:=#v1-1;
    for i in [0..n] do
        res[i+1]:=P!0;
        for j in [0..i] do
            res[i+1]+:=v1[j+1]^(p^(i-j))*v2[i-j+1]^(p^j);
        end for;
    end for;
    return res;
end function;


// not efficient -- just to test!
function GTaux2(v1,v2,i,j)
    P:=Parent(v1[1]);
    p:=Characteristic(P);
    res:=[ P!1 ] cat [ P!0: ii in [2..#v1] ];
    for ii in [1..i] do
        res:=GTaux1(res,v1);
    end for;
    for ii in [1..j] do
        res:=GTaux1(res,v2);
    end for;
    return res;
end function;





/*
p:=3;
n:=4;
// k:=4;

// F<a>:=GF(p^k);

F:=GF(p);

P:=PolynomialRing(F,(2*n+2));
AssignNames(~P,
  [ "x" cat IntegerToString(i) : i in [0..n] ] cat
  [ "y" cat IntegerToString(i) : i in [0..n] ] );

x:=[ P.i : i in [1..(n+1)] ];
y:=[ P.i : i in [(n+2)..(2*n+2)] ];

GTaux(x,y);
*/
