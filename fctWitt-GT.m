// computes canonical lifting and elliptic Teichmuller lift
// using new methods (etas, formular for GT, etc.)
// WORKS FOR n>2!!!!

load "gt.m"; // which loads etas.m

function my_int(f)

    deg:=Degree(f);
    P:=Parent(f);
    p:=Characteristic(P);
    res:=0;
    for d in [0..deg] do
        if ((d+1) mod p) ne 0 then
            res+:=Coefficient(f,d)/(d+1)*P.1^(d+1);
        end if;
    end for;
    return res;

end function;


function lift(a0,b0,n : pols:=[])
    // lenght (n+1) -- up to a_n, b_n
    
    F:=Parent(a0);
    p:=Characteristic(F);

    if #pols eq 0 then
        pols:=etapols(p,n);
    end if;

    P1:=PolynomialRing(F,2*n+2*(n+1));
    // a_i's, b_i's, x_i's, y_i's

    Pres<x0>:=PolynomialRing(F);
    // ring for the resulting polynomials in x0

    va:=[a0] cat [ P1.j : j in [1..n]];
    vb:=[b0] cat [ P1.j : j in [n+1..2*n]];
    vx:=[ P1.j : j in [2*n+1..3*n+1]];
    vy:=[ P1.j : j in [3*n+2..4*n+2]];
    vone:=[1] cat [ P1!0 : j in [1..n]];
    
    GTx:=GT( [[* vone,3,0 *], [* va,1,0 *], [* vb,0,0 *]] : pols:=pols);
    GTx:=[ Evaluate(term,vx cat vy) : term in GTx];
    
    // PGT:=Parent(GTx[1]);

    GTy:=GT( [[* vone,0,2 *]] : pols:=pols);
    GTy:=[ Evaluate(term,vx cat vy) : term in GTy];

    
    // don't need the first entry
    Remove(~GTx,1);
    Remove(~GTy,1);

    // results
    resa:=[a0];
    resb:=[b0];
    resF:=[x0];
    resH:=[Pres!1];


    f:=x0^3+a0*x0+b0;
    ff:=f^((p-1) div 2);
    
    // Hasse Invariant
    HI:=F!(Coefficient(ff,p-1));

    
    // to convert from P1 to Pres
    vP1Pres:=[ 0 : j in [1..2*n] ] cat [x0] cat
             [ 0 : j in [2*n+2..4*n+2]];

    
    // main loop
    for i in [1..n] do

        //polynomial ring for this coordinate
        M:=(3*p^(i-1)-1) div 2; // c_i's
        N:=((i+3)*p^i - i*p^(i-1) - 3) div 2; // d_i's
        Pi:=PolynomialRing(F,2 + 2 + (M+1) + (N+1));
        // x_0, y_0, a_n, b_n, c_i's, d_i's

        // a temporary power of f to help compute F_i's
        if i eq 1 then
           tmppf:=ff;
        else
           tmppf:=tmppf^p*ff;
        end if;
           
        Fi:= HI^(-(p^i-1) div (p-1))*tmppf - x0^(p^i-1);
        if i gt 1 then
            Fi-:=&+[ resF[j+1]^(p^(i-j)-1)*Derivative(resF[j+1]) :
                     j in [1..(i-1)] ];
        end if;
        Fi:=my_int(Fi);
        Fi:=Evaluate(Fi,Pi.1);

        // will make c_(p^(n-1))=0
        Fi+:=&+[ Pi.(5+j)*(Pi.1)^(p*j) : j in [0..M] | j ne p^(i-1) ];

        // add condition if not minimal!!!!!!!

        Hi:=&+[ Pi.(6+M+j)*(Pi.1)^j : j in [0..N]];

        print Fi, Hi;
        
        
        // vector to convert from P1 to Pi
        vP1Pi:=[ 0 : j in [1..(i-1)] ] cat [Pi.3] cat
               [ 0 : j in [i+1..n+i-1] ] cat [Pi.4] cat
               [ 0 : j in [n+i+1..2*n] ] cat
               [ Pi.1 ] cat [ 0 : j in [1..(i-1)] ] cat
               [ Fi ] cat [ 0 : j in [i+1..n]] cat
               [ Pi.2 ] cat [ 0 : j in [1..(i-1)] ] cat
               [ Pi.2*Hi ] cat [ 0 : j in [i+1..n]];
        
        // assumes that we replaced previsous values of
        // a_j, b_j, F_j, and H_j already and removed
        // all previous entries
        LHS:=Evaluate(GTy[1],vP1Pi);
        RHS:=Evaluate(GTx[1],vP1Pi);

        // now, replace y0^2 by f
        deg:=Degree(LHS,2);   
        tmppf2:=1;
        tmpLHS:=Coefficient(LHS,2,0);
        for d in [j : j in [2..deg] | (j mod 2) eq 0 ] do
           tmppf2*:=Evaluate(f,Pi.1);
           tmpLHS+:=Coefficient(LHS,2,d)*tmppf2;
        end for;

        eqts:=Coefficients(tmpLHS-RHS,1);

        delete tmpLHS, LHS, RHS, tmppf2;

        // matrix of coefficients
        neqts:=#eqts;
        Mat:=Matrix(F,2+M+(N+1),neqts,
             [[F!(Coefficient(eqts[j],k+2,1)) : j in [1..neqts]]
                 : k in [ ii : ii in [1..(2+(M+1)+(N+1))] | ii ne 3+p^(i-1) ]]);

        vec:=Vector(F,[ -Evaluate(eqts[j],[ 0 : k in [1..(2+2+(M+1)+(N+1))]]) : j in [1..neqts]]);

        vsol:=Solution(Mat,vec);

        // to convert solutoions to Pres/F
        evalvec:=[x0,0] cat [ vsol[j] : j in [1..2+p^(i-1)] ] cat [0] cat
                 [vsol[j] : j in [2+p^(i-1)+1..2+M+(N+1)]]; 

        print evalvec;
        
        Append(~resa,vsol[1]);
        Append(~resb,vsol[2]);
        Append(~resF,Evaluate(Fi,evalvec));
        Append(~resH,Evaluate(Hi,evalvec));

        
        // don't need the first entry
        Remove(~GTx,1);
        Remove(~GTy,1);
        
        // replace by the values found all the next coordinates
        evalvec:=[ 0 : j in [1..(i-1)] ] cat [P1!resa[i+1]] cat
                 [ P1.j : j in [i+1..n] ] cat
                 [ 0 : j in [n+1..(n+i-1)] ] cat  [P1!resb[i+1]] cat
                 [ P1.j : j in [n+i+1..2*n] ] cat [P1.(2*n+1)] cat
                 [ 0 : j in [2*n+2..2*n+i]] cat
                 [Evaluate(resF[i+1],P1.(2*n+1))] cat
                 [ P1.j : j in [2*n+i+2..3*n+2]] cat
                 [ 0 : j in [3*n+3..3*n+1+i]] cat
                 [P1.(3*n+2)*Evaluate(resH[i+1],P1.(2*n+1))] cat
                 [ P1.j : j in [3*n+i+3..4*n+2]];

        tmpGTx:=[ Evaluate(term,evalvec) : term in GTx];
        tmpGTy:=[ Evaluate(term,evalvec) : term in GTy];

        GTx:=tmpGTx;
        GTy:=tmpGTy;

        delete tmpGTx, tmpGTy;
    
    end for;
        
    return resa, resb, resF, resH;

end function;
