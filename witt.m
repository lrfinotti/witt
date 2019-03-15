// /////////////////////////////////////////////////
// WITT VECTORS
// Witt Vector Operations using the etas
// /////////////////////////////////////////////////

load 'etas.m';


// ///////////////////////////////////////
// SUMS
// ///////////////////////////////////////

// using vetav
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


// using vetav2
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



// using vetav3
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



// ///////////////////////////////////////
// PRODUCTS
// ///////////////////////////////////////


// using vetav
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
        // print "Computing etapols";
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



// using vetav2
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


// using vetav3
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




// ///////////////////////////////////////
// NEGATIVES
// ///////////////////////////////////////

function newWittNeg(v : pols:=[])
    p := Characteristic(Parent(v[1]));
    if p ne 2 then
        return [ -x : x in v ];
    end if;
    n := #v;
    P := Parent(v[1]);
    vnone := [ P!1 : i in [1..n] ];
    return newWittProd(vnone,v : pols:=pols);
end function;


function newWittNeg2(v : bintab:=[])
    p := Characteristic(Parent(v[1]));
    if p ne 2 then
        return [ -x : x in v ];
    end if;
    n := #v;
    P := Parent(v[1]);
    vnone := [ P!1 : i in [1..n] ]; 
    return newWittProd2(vnone,v : bintab:=bintab);
end function;


function newWittNeg3(v : bintab:=[])
    p := Characteristic(Parent(v[1]));
    if p ne 2 then
        return [ -x : x in v ];
    end if;
    n := #v;
    P := Parent(v[1]);
    vnone := [ P!1 : i in [1..n] ]; 
    return newWittProd3(vnone,v : bintab:=bintab);
end function;




// ///////////////////////////////////////
// INVERSES
// ///////////////////////////////////////

function newWittInv(v : pols:=[])
    P := Parent(v[1]);
    p := Characteristic(P);
    n := #v-1;

    if #pols eq 0 then
        //print "Computing etapols";
        pols:=etapols(p,n-1); // NOTE: n-1 only!
    end if;

    res := [ v[1]^(-1) ];
    PR := PolynomialRing(P);
    x := PR.1; // variable

    for i in [1..n] do
        w1 := [ PR!v[j] : j in [1..(i+1)] ];
        w2 := [ PR!res[j] : j in [1..i] ] cat [x];
        coord := newWittProd(w1,w2 : pols:=pols)[i+1];
        Append(~res,-Evaluate(coord,0)/v[1]^(p^i));
    end for;
    return res;
end function;
    


function newWittInv2(v :  bintab:=[])
    P := Parent(v[1]);
    p := Characteristic(P);
    n := #v-1;

    if #bintab eq 0 then
        // print "Computing bintab";
        bintab:=BinTab(p,n-1);  // note only need n-1
    end if;

    res := [ v[1]^(-1) ];
    PR := PolynomialRing(P);
    x := PR.1; // variable

    for i in [1..n] do
        w1 := [ PR!v[j] : j in [1..(i+1)] ];
        w2 := [ PR!res[j] : j in [1..i] ] cat [x];
        coord := newWittProd2(w1,w2 : bintab:=bintab)[i+1];
        Append(~res,-Evaluate(coord,0)/v[1]^(p^i));
    end for;
    return res;
end function;
    


function newWittInv3(v :  bintab:=[])
    P := Parent(v[1]);
    p := Characteristic(P);
    n := #v-1;

    if #bintab eq 0 then
        // print "Computing bintab";
        bintab:=BinTab(p,n-1);  // note only need n-1
    end if;

    res := [ v[1]^(-1) ];
    PR := PolynomialRing(P);
    x := PR.1; // variable

    for i in [1..n] do
        w1 := [ PR!v[j] : j in [1..(i+1)] ];
        w2 := [ PR!res[j] : j in [1..i] ] cat [x];
        coord := newWittProd3(w1,w2 : bintab:=bintab)[i+1];
        Append(~res,-Evaluate(coord,0)/v[1]^(p^i));
    end for;
    return res;
end function;
    


// ///////////////////////////////////////
// POWERS
// ///////////////////////////////////////



// powers in gt.m!


// /////////////////////////////////////////////////
// WRAPPERS
// Use a single function to perform the operations
// /////////////////////////////////////////////////

function WittSum(v,w : choice:=1, pols:=[], bintab:=[])
    // choice (of algorithms): 1, 2 or 3
    // for vetav, vetav2 or vetav3 respectively
    // add sanity check?
    P := Parent(v[1]);

    // if over a finite field, perform in Zq
    if IsFinite(P) and IsField(P) then
        p:=Characteristic(P);
        n:=#v-1;
        k:=Degree(P);
        Zp:=pAdicRing(p : Precision:=n+1);
        Zq<aa>:= ext<Zp | k>;
        v1 := WittVToSeries(v : Zq:=Zq);
        w1 := WittVToSeries(w : Zq:=Zq);
        PX := Parent(v1);
        tmp := SeriesToWittV(v1 + PX!w1);
        return [ P!x : x in tmp ];
    end if;
    
    if choice eq 2 then
        return newWittSum2(v,w :  bintab:=bintab);
    elif choice eq 3 then
        return newWittSum3(v,w :  bintab:=bintab);
    else
        return newWittSum(v,w :  pols:=pols);
    end if;
end function;



function WittProd(v,w : choice:=1, pols:=[], bintab:=[])
    // choice (of algorithms): 1, 2 or 3
    // for vetav, vetav2 or vetav3 respectively
    // add sanity check?
    P := Parent(v[1]);

    // if over a finite field, perform in Zq
    if IsFinite(P) and IsField(P) then
        p:=Characteristic(P);
        n:=#v-1;
        k:=Degree(P);
        Zp:=pAdicRing(p : Precision:=n+1);
        Zq<aa>:= ext<Zp | k>;
        v1 := WittVToSeries(v : Zq:=Zq);
        w1 := WittVToSeries(w : Zq:=Zq);
        PX := Parent(v1);
        tmp := SeriesToWittV(v1 * PX!w1);
        return [ P!x : x in tmp ];
    end if;

    if choice eq 2 then
        return newWittProd2(v,w :  bintab:=bintab);
    elif choice eq 3 then
        return newWittProd3(v,w :  bintab:=bintab);
    else
        return newWittProd(v,w :  pols:=pols);
    end if;
end function;


function WittNeg(v : choice:=1, pols:=[], bintab:=[])
    // choice (of algorithms): 1, 2 or 3
    // for vetav, vetav2 or vetav3 respectively
    // add sanity check?

    // char p diff from 2
    P := Parent(v[1]);
    p:=Characteristic(P);
    if p ne 2 then
        return [ -x : x in v ];
    end if;

    // char p = 2
    // if over a finite field, perform in Zq
    if IsFinite(P) and IsField(P) then
        v1 := WittVToSeries(v);
        tmp := SeriesToWittV(-v1);
        return [ P!x : x in tmp ];
    end if;

    if choice eq 2 then
        return newWittNeg2(v :  bintab:=bintab);
    elif choice eq 3 then
        return newWittNeg3(v :  bintab:=bintab);
    else
        return newWittNeg(v :  pols:=pols);
    end if;
end function;


function WittNeg(v : choice:=1, pols:=[], bintab:=[])
    // choice (of algorithms): 1, 2 or 3
    // for vetav, vetav2 or vetav3 respectively
    // add sanity check?

    P := Parent(v[1]);
    p:=Characteristic(P);

    // if over a finite field, perform in Zq
    if IsFinite(P) and IsField(P) then
        v1 := WittVToSeries(v);
        tmp := SeriesToWittV(1/v1);
        return [ P!x : x in tmp ];
    end if;

    if choice eq 2 then
        return newWittInv2(v :  bintab:=bintab);
    elif choice eq 3 then
        return newWittInv3(v :  bintab:=bintab);
    else
        return newWittInv(v :  pols:=pols);
    end if;
end function;

