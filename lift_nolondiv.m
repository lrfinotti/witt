// computes canonical lifting and elliptic Teichmuller lift
// using new methods (etas, formular for GT, etc.)
//
// This was created to allow to compute more than 2 coordinates.
// Unfortunately, it is *MUCH* slower than fctWitt5.m for 2 coordinates!!!
// I need to check the times of parts to pinpoint the slowdown.
//
// This is very similar to minlift.m, which seems to be efficient...
//
// THIS VERSION TRIES TO COMPUTE IT WITHOUT DOING THE LONG DIVSION!
// IT IS INEFFICIENT!!!!!!


load "gt.m"; // which loads etas.m

function my_int(f)
    // formal integral of polynomials
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


function lift(a0,b0,n : pols:=[], minimal:=false)
    // lenght (n+1) -- up to a_n, b_n
    // returns vectors va, vb, vF, vH
    // pols are etapols from etas.m

    F:=Parent(a0);
    p:=Characteristic(F);

    if #pols eq 0 then
        pols:=etapols(p,n);
    end if;

    // ring for the resulting polynomials in x0
    Pres<x0>:=PolynomialRing(F);


    // results
    resa:=[a0];
    resb:=[b0];
    resF:=[x0];
    resH:=[Pres!1];


    f:=x0^3+a0*x0+b0;
    ff:=f^((p-1) div 2);
    
    // Hasse Invariant
    HI:=F!(Coefficient(ff,p-1));


    // x/y to compute can. lift.
    // load "witt_gen.m"; loaded by gt.m
    if n ge 3 then
        tmpF:=RationalFunctionField(GF(p),2*n+2);
        vx:=[tmpF.j : j in [1..(n+1)]];
        vy:=[tmpF.j : j in [(n+2)..(2*n+2)]];
        vyi:=WittInv(vy);
        delete vy;
        vxoy:=newWittProd(vx,vyi : pols:=pols);
        delete vx;
    end if;

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

        ///////////////////////////////////////////////////////////////
        // ADD CONDITION IF NOT MINIMAL!!!!!!!

        if (not minimal) and (i eq 2) then
            tmp:= F!(3/4)*resF[2]^2;
            Fi+:= &+[(Coefficient(tmp,Integers()!(i/p+p)))^p*Pi.1^i :
                     i in [Integers()!((3*p^2+p)/2) .. 2*p^2-p by p]];
        end if;

        // Condition for i ge 3
        // if (not minimal) and (i ge 3) then
        //
        // end if;

        ///////////////////////////////////////////////////////////////


        Hi:=&+[ Pi.(6+M+j)*(Pi.1)^j : j in [0..N]];

        //print Fi, Hi;
        
        va:=[ Pi!x : x in resa ] cat [Pi.3];
        vb:=[ Pi!x : x in resb ] cat [Pi.4];
        vF:=[ Evaluate(x,Pi.1) : x in resF ] cat [Fi];
        vG:=[ Pi.2*Evaluate(x,Pi.1) : x in resH ] cat [Pi.2*Hi];
        vone:=[Pi!1] cat [ Pi!0 : j in [1..i]];

        vvars:=vF cat vG;

        GTx:=GT( [[* vone,3,0 *], [* va,1,0 *], [* vb,0,0 *]] : pols:=pols, vvars:=vvars);

        GTy:=GT( [[* vone,0,2 *]] : pols:=pols, vvars:=vvars);
        
        RHS:=GTx[i+1];
        delete GTx;
        LHS:=GTy[i+1];
        delete GTy;
        

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

        print "Size = ", NumberOfRows(Mat), NumberOfColumns(Mat);
        ptime:=Cputime();
        vsol:=Solution(Mat,vec);
        print "Partial time = ", Cputime(ptime);

        // to convert solutoions to Pres/F
        evalvec:=[x0,0] cat [ vsol[j] : j in [1..2+p^(i-1)] ] cat [0] cat
                 [vsol[j] : j in [2+p^(i-1)+1..2+M+(N+1)]]; 

        //print evalvec;
        
        Append(~resa,vsol[1]);
        Append(~resb,vsol[2]);
        Append(~resF,Evaluate(Fi,evalvec));
        Append(~resH,Evaluate(Hi,evalvec));

    
    end for;
        
    return resa, resb, resF, resH;

end function;
