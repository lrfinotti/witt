// Computes the canonical lifting of an Elliptic Curve in Edwards Form
// by way of computing the coefficients

// we need the functions for the Greenberg transform
load "gt.m"; // which loads etas.m, and witt_gen.m

function my_int(f, var)
    // formal integral of polynomials
    // var: number for the variable of integration
    deg:=Degree(f);
    P:=Parent(f);
    p:=Characteristic(P);
    res:=0;
    for d in [0..deg] do
        if ((d+1) mod p) ne 0 then
            res+:=Coefficient(f,var,d)/(d+1)*(P.var)^(d+1);
        end if;
    end for;
    return res;

end function;


function my_x_int(f)
    return my_int(f, 1);
end function;

function my_y_int(f)
    return my_int(f, 2);
end function;


// NOT WORKING: the polynomial ring is not only on x0 and y0, so Monomial Coefficient
// does not give the other variables!!!
function reduce_pol(f, a0, b0)
    // reduce polynomial f module b0^2*x0^2 + y0^2 - 1 - a0*x0^2y0^2
    // the result should be of the form f1(x0) + f2(y0) + x0*f3(y0) + y0*f4(x0)
    //
    // NOTE: this seems to be the same as reducing modulo the equation in Magma,
    // which should be faster!!!
    x0 := Parent(f).1;
    y0 := Parent(f).2;

    repl_pol := b0^2/a0*x0^2 + y0^2/a0 - 1/a0;

    res := 0;
    g := f;
    // let's add the "good" terms to the result
    while g ne 0 do
        good_terms := 0;
        for mon in Monomials(g) do
            if (Degree(mon,1) eq 0) or (Degree(mon,2) eq 0) or (Degree(mon,1) eq 1)
                or (Degree(mon,2) eq 1) then
                good_terms +:= MonomialCoefficient(g, mon)*mon;
            end if;
        end for;

        g -:= good_terms; // "bad terms"
        res +:= good_terms;

        // let's remove some bad terms from g
        newg := 0;
        for mon in Monomials(g) do
            m := Degree(mon,1);
            n := Degree(mon,2);
            minexp := Minimum(m div 2, n div 2);
            newg +:= MonomialCoefficient(g, mon)*x0^(m - 2*minexp)*y0^(n - 2*minexp)*repl_pol^minexp;
        end for;

        g:=newg;
    end while;

    return res;
end function;



function lift_ed(a0,b0,n : pols:=[], verb:=true, cstr:=[], check:=true, minimal:=false)

    F:=Universe([a0,b0]);
    p:=Characteristic(F);

    if #pols eq 0 then
            pols:=etapols(p,n);
    end if;

    Pres<x0,y0>:=PolynomialRing(F,2);

    /* CurveEq := b0^2*x0^2 + y0^2 - 1 - a0*x0^2*y0^2; */

    /* // introduce variables for the quotient */
    /* C<xx0, yy0>:=quo<Pres | CurveEq>; */

    /* // let's get the conversion maps */
    /* PtoC := hom< Pres -> C | [xx0, yy0]>; */
    /* CtoP := hom< C -> Pres | [x0, y0]>; */

    resa:=[a0];
    resb:=[b0];
    resx:=[x0];
    resy:=[y0];


    // FIXME: do we need two?
    xf:=(b0^2*x0^2-1)*(a0*x0^2-1); // Important Equation used in derivative of x_n
    yg:=(b0^2-a0*y0^2)*(1-y0^2);   // Important Equation used in derivative of y_n
                                 // and computing Hasse Invariant

    ff:=xf^((p-1) div 2);
    gg:=yg^((p-1) div 2); //Polynomial for Hasse Invaraint

    HI:=F!(Coefficient(gg, 2, p-1)); //Hasse Invaraint

    if HI eq 0 then
            return "This Curve is not Ordinary";
    end if;


    //main loop
    for i in [1..n] do
        if verb then
            print "Working on i equals: ", i;
        end if;

        // necessary variables
        // we only need the coefficients below the degree of the derivative
        // we do not need "+1" since constant terms are zero
        M:=2*p^(i-1)-1;

        // print "The value of M is: ", M;

        Pi:=PolynomialRing(F,2 + 2 + M + M + 1);
        // x_0, y_0, a_n, b_n, epsilon_i's, delta_i's (recall epsilon_0=delta_0=0)
        // plus one additional variable to help with inverses in the computation
        // of higher degree terms

        // a temporary power of f to help compute x_i's and y_i's:
        if i eq 1 then
                tmppf:=ff;
                tmppg:=gg;
        else
                tmppf:=tmppf^p*ff;
                tmppg:=tmppg^p*gg;
        end if;

        xi:= HI^(-(p^i-1) div (p-1))*tmppf - x0^(p^i-1);
        if i gt 1 then
                xi-:=&+[ resx[j+1]^(p^(i-j)-1)*Derivative(resx[j+1],1) :
                        j in [1..(i-1)] ];
        end if;

        xi:=my_x_int(xi);
        xi:=Evaluate(xi,[Pi.1,Pi.2]);


        yi:= HI^(-(p^i-1) div (p-1))*tmppg - y0^(p^i-1);
        if i gt 1 then
                yi-:=&+[ resy[j+1]^(p^(i-j)-1)*Derivative(resy[j+1],2) :
                        j in [1..(i-1)] ];
        end if;

        yi:=my_y_int(yi);
        yi:=Evaluate(yi,[Pi.1,Pi.2]);

        // NOT YET: will make epsilon_(p^(n-1))=0
        // we have epsilon_0 = 0
        // xi+:=&+[ Pi.(5+j)*(Pi.1)^(p*j) : j in [1..M] | j ne p^(i-1) ];
        xi+:=&+[ Pi.(4+j)*(Pi.1)^(p*j) : j in [1..M] ];
        yi+:=&+[ Pi.(4+M+j)*(Pi.2)^(p*j) : j in [1..M] ];


        // find higher degree terms
        if (i gt 1) and (not minimal) then

            if verb then
                print "Computing higher degree terms.";
            end if;

            // to invert
            Fi := FieldOfFractions(Pi);

            vx:=[ Fi!Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [Fi.(2*M+5)];
            vy:=[ Fi!Evaluate(y,[Pi.1,Pi.2]) : y in resy ] cat [Fi.(2*M+5)];


            // xi
            vinvi:=WittInv(vx : pols:=pols)[i+1];
            num:=Pi!Numerator(vinvi);
            den:=Pi!Denominator(vinvi);

            rem_terms:=Coefficient(num,Pi.(2*M+5),0); // terms without xi

            degn:=Degree(Coefficient(num,Pi.(2*M+5),1),1);
            degd := Degree(den,Pi.1);

            for k in [((degd-degn) div p)..((i+1)*p^(i-1)-i*p^(i-2))] do
                xi +:= Coefficient(rem_terms,1,k*p+degn)*(Pi.1)^(k*p);
            end for;


            // yi
            vinvi:=WittInv(vy : pols:=pols)[i+1];
            num:=Pi!Numerator(vinvi);
            den:=Pi!Denominator(vinvi);

            rem_terms:=Coefficient(num,Pi.(2*M+5),0); // terms without xi

            degn:=Degree(Coefficient(num,Pi.(2*M+5),1),2);
            degd := Degree(den,Pi.2);

            for k in [((degd-degn) div p)..((i+1)*p^(i-1)-i*p^(i-2))] do
                yi +:= Coefficient(rem_terms,2,k*p+degn)*(Pi.2)^(k*p);
            end for;


        end if;


        va:=[ Pi!x : x in resa ] cat [Pi.3];
        vb:=[ Pi!x : x in resb ] cat [Pi.4];
        vb2:=WittPower(vb, 2 : pols:=pols);
        vx:=[ Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [xi];
        vy:=[ Evaluate(y,[Pi.1,Pi.2]) : y in resy ] cat [yi];
        vone:=[Pi!1] cat [ Pi!0 : j in [1..i]];

        vvars:=vx cat vy;


        // Compute the Greenberg Transform
        if verb then
            print "Computing the GT.";
        end if;

        LeftGT:=GT( [[* vb2,2,0 *], [* vone,0,2 *]] : pols:=pols, vvars:=vvars);
        RightGT:=GT( [[* va,2,2 *], [* vone,0,0 *]] : pols:=pols, vvars:=vvars);

        difference:=RightGT[i+1] - LeftGT[i+1];

        // We need to add the extra conditions from tau(0,1)=(0,1)
        maxd := Degree(difference,1);

        difference +:= Evaluate(yi, Pi.2, 1)*Pi.1^(maxd + 1);

        // Now tau((1/b0,0)) = (1/b, 0)
        extra_term := WittInv(vb : pols:=pols)[i+1];
        difference +:= (Evaluate(xi, Pi.1, 1/b0) - extra_term)*Pi.1^(maxd + 2);

        // mod out by the curve
        CurveEq := b0^2*(Pi.1)^2 + (Pi.2)^2 - 1 - a0*(Pi.1)^2*(Pi.2)^2;
        Qi:=quo<Pi | CurveEq>;

        // put in reduced form
        // difference:=reduce_pol(difference, a0, b0);

        difference := Pi!(Qi!difference);

       // collect the coefficients, which then have to be zero
        eqts:=[];
        for coefx in Coefficients(difference, Pi.1) do
            for coefxy in Coefficients(coefx, Pi.2) do
                if coefxy ne 0 then
                    Append(~eqts, coefxy);
                end if;
            end for;
        end for;

        neqts:=#eqts;

        // print eqts;

        Mat:=Matrix(F,2 + 2*M,neqts,
                [[Coefficient(eqts[j],k+2,1) : j in [1..neqts]]
                        : k in [1..(2+2*M)]]);

        vec:=Vector(F,
                 [ -Evaluate(eqts[j], [ 0 : k in [1..(4+2*M+1)] ]) : j in [1..neqts] ]);

        if verb then
            print "Solving the system";
        end if;

        vsol:=Solution(Mat,vec);

        //print vsol;

        // to convert solutions to Pres/F
        evalvec:=[x0,y0] cat [ vsol[j] : j in [1..2+2*M]] cat [0];

       // print #evalvec;

        Append(~resa,vsol[1]);
        Append(~resb,vsol[2]);
        Append(~resx,Evaluate(xi,evalvec));
        Append(~resy,Evaluate(yi,evalvec));

    end for;

    return resa, resb, resx, resy;

end function;






function lift_ed_1(a0,n : pols:=[], verb:=true, cstr:=[], check:=true, minimal:=false)

    F:=Parent(a0);
    b0:=F!1;
    p:=Characteristic(F);

    if #pols eq 0 then
            pols:=etapols(p,n);
    end if;

    Pres<x0,y0>:=PolynomialRing(F,2);

    /* CurveEq := b0^2*x0^2 + y0^2 - 1 - a0*x0^2*y0^2; */

    /* // introduce variables for the quotient */
    /* C<xx0, yy0>:=quo<Pres | CurveEq>; */

    /* // let's get the conversion maps */
    /* PtoC := hom< Pres -> C | [xx0, yy0]>; */
    /* CtoP := hom< C -> Pres | [x0, y0]>; */

    resa:=[a0];
    resx:=[x0];
    resy:=[y0];


    // FIXME: do we need two?
    xf:=(x0^2-1)*(a0*x0^2-1); // Important Equation used in derivative of x_n

    ff:=xf^((p-1) div 2);

    HI:=F!(Coefficient(ff, 1, p-1)); //Hasse Invaraint

    if HI eq 0 then
            return "This Curve is not Ordinary";
    end if;


    //main loop
    for i in [1..n] do
        if verb then
            print "Working on i equals: ", i;
        end if;

        // necessary variables
        // we only need the coefficients below the degree of the derivative
        // we do not need "+1" since constant terms are zero
        M:=2*p^(i-1)-1;

        // print "The value of M is: ", M;

        Pi:=PolynomialRing(F,2 + 1 + M + 1);
        // x_0, y_0, a_n, epsilon_i' (recall epsilon_0=delta_0=0)
        // plus one additional variable to help with inverses in the computation
        // of higher degree terms

        // a temporary power of f to help compute x_i's and y_i's:
        if i eq 1 then
                tmppf:=ff;
        else
                tmppf:=tmppf^p*ff;
        end if;

        xi:= HI^(-(p^i-1) div (p-1))*tmppf - x0^(p^i-1);
        if i gt 1 then
                xi-:=&+[ resx[j+1]^(p^(i-j)-1)*Derivative(resx[j+1],1) :
                        j in [1..(i-1)] ];
        end if;

        xi:=my_x_int(xi);
        yi:=Evaluate(xi,[Pi.2,Pi.1]);
        xi:=Evaluate(xi,[Pi.1,Pi.2]);

        // NOT YET: will make epsilon_(p^(n-1))=0
        // we have epsilon_0 = 0
        // xi+:=&+[ Pi.(5+j)*(Pi.1)^(p*j) : j in [1..M] | j ne p^(i-1) ];
        xi+:=&+[ Pi.(3+j)*(Pi.1)^(p*j) : j in [1..M] ];
        yi+:=&+[ Pi.(3+j)*(Pi.2)^(p*j) : j in [1..M] ];



        // find higher degree terms
        if (i gt 1) and (not minimal) then

            if verb then
                print "Computing higher degree terms.";
            end if;

            // to invert
            Fi := FieldOfFractions(Pi);

            vx:=[ Fi!Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [Fi.(M+4) ];

            // xi
            vinvi:=WittInv(vx : pols:=pols)[i+1];
            num:=Pi!Numerator(vinvi);
            den:=Pi!Denominator(vinvi);

            rem_terms:=Coefficient(num,Pi.(M+4),0); // terms without xi

            degn:=Degree(Coefficient(num,Pi.(M+4),1),1);
            degd := Degree(den,Pi.1);

            for k in [((degd-degn) div p)..((i+1)*p^(i-1)-i*p^(i-2))] do
                xi +:= Coefficient(rem_terms,1,k*p+degn)*(Pi.1)^(k*p);
                yi +:= Coefficient(rem_terms,1,k*p+degn)*(Pi.2)^(k*p);
            end for;

        end if;

        va:=[ Pi!x : x in resa ] cat [Pi.3];
        vx:=[ Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [xi];
        vy:=[ Evaluate(y,[Pi.1,Pi.2]) : y in resy ] cat [yi];
        vone:=[Pi!1] cat [ Pi!0 : j in [1..i]];

        vx2:=WittPower(vx,2 : pols:=pols);
        vy2:=[ Evaluate(x,[Pi.2,Pi.1] cat [Pi.j : j in [3..(M+4)]]) : x in vx2 ];

        // vvars:=vx cat vy;
        vvars:=vx2 cat vy2;

        // Compute the Greenberg Transform
        if verb then
            print "Computing the GT.";
        end if;

        // LeftGT:=GT( [[* vone,2,0 *], [* vone,0,2 *]] : pols:=pols, vvars:=vvars);
        // RightGT:=GT( [[* va,2,2 *], [* vone,0,0 *]] : pols:=pols, vvars:=vvars);
        LeftGT:=GT( [[* vone,1,0 *], [* vone,0,1 *]] : pols:=pols, vvars:=vvars);
        RightGT:=GT( [[* va,1,1 *], [* vone,0,0 *]] : pols:=pols, vvars:=vvars);

        difference:=RightGT[i+1] - LeftGT[i+1];

        // We need to add the extra conditions from tau(0,1)=(0,1)
        maxd := Degree(difference,1);

        // this adds that xi(1) = yi(1) = 0
        difference +:= Evaluate(yi, Pi.2, 1)*Pi.1^(maxd + 1);

        // mod out by the curve
        CurveEq := (Pi.1)^2 + (Pi.2)^2 - 1 - a0*(Pi.1)^2*(Pi.2)^2;
        Qi:=quo<Pi | CurveEq>;

        // put in reduced form
        // difference:=reduce_pol(difference, a0, b0);

        difference := Pi!(Qi!difference);

       // collect the coefficients, which then have to be zero
        eqts:=[];
        for coefx in Coefficients(difference, Pi.1) do
            for coefxy in Coefficients(coefx, Pi.2) do
                if coefxy ne 0 then
                    Append(~eqts, coefxy);
                end if;
            end for;
        end for;

        neqts:=#eqts;

        // print eqts;

        Mat:=Matrix(F,1 + M,neqts,
                [[Coefficient(eqts[j],k+2,1) : j in [1..neqts]]
                        : k in [1..(1+M)]]);

        vec:=Vector(F,
                 [ -Evaluate(eqts[j], [ 0 : k in [1..(3+M+1)] ]) : j in [1..neqts] ]);

        if verb then
            print "Solving the system";
        end if;

        vsol:=Solution(Mat,vec);

        // to convert solutions to Pres/F
        evalvec:=[x0,y0] cat [ vsol[j] : j in [1..1+M]] cat [0];

       // print #evalvec;

        Append(~resa,vsol[1]);
        Append(~resx,Evaluate(xi,evalvec));
        Append(~resy,Evaluate(yi,evalvec));

    end for;

    return resa, resx, resy;

end function;



function lift_ed_1a(a0,n : pols:=[], verb:=true, cstr:=[], check:=true, minimal:=false)

    F:=Parent(a0);
    b0:=F!1;
    p:=Characteristic(F);

    if #pols eq 0 then
            pols:=etapols(p,n);
    end if;

    Pres<x0,y0>:=PolynomialRing(F,2);

    /* CurveEq := b0^2*x0^2 + y0^2 - 1 - a0*x0^2*y0^2; */

    /* // introduce variables for the quotient */
    /* C<xx0, yy0>:=quo<Pres | CurveEq>; */

    /* // let's get the conversion maps */
    /* PtoC := hom< Pres -> C | [xx0, yy0]>; */
    /* CtoP := hom< C -> Pres | [x0, y0]>; */

    resa:=[a0];
    resx:=[x0];
    resy:=[y0];


    // FIXME: do we need two?
    xf:=(x0^2-1)*(a0*x0^2-1); // Important Equation used in derivative of x_n

    ff:=xf^((p-1) div 2);

    HI:=F!(Coefficient(ff, 1, p-1)); //Hasse Invaraint

    if HI eq 0 then
            return "This Curve is not Ordinary";
    end if;


    //main loop
    for i in [1..n] do
        if verb then
            print "Working on i equals: ", i;
        end if;

        // necessary variables
        // we only need the coefficients below the degree of the derivative
        // we do not need "+1" since constant terms are zero
        M:=2*p^(i-1)-1;

        // print "The value of M is: ", M;

        Pi:=PolynomialRing(F,2 + 1 + M + 1);
        // x_0, y_0, a_n, epsilon_i' (recall epsilon_0=delta_0=0)
        // plus one additional variable to help with inverses in the computation
        // of higher degree terms

        // a temporary power of f to help compute x_i's and y_i's:
        if i eq 1 then
                tmppf:=ff;
        else
                tmppf:=tmppf^p*ff;
        end if;

        xi:= HI^(-(p^i-1) div (p-1))*tmppf - x0^(p^i-1);
        if i gt 1 then
                xi-:=&+[ resx[j+1]^(p^(i-j)-1)*Derivative(resx[j+1],1) :
                        j in [1..(i-1)] ];
        end if;

        xi:=my_x_int(xi);
        yi:=Evaluate(xi,[Pi.2,Pi.1]);
        xi:=Evaluate(xi,[Pi.1,Pi.2]);

        // NOT YET: will make epsilon_(p^(n-1))=0
        // we have epsilon_0 = 0
        // xi+:=&+[ Pi.(5+j)*(Pi.1)^(p*j) : j in [1..M] | j ne p^(i-1) ];
        xi+:=&+[ Pi.(3+j)*(Pi.1)^(p*j) : j in [1..M] ];
        yi+:=&+[ Pi.(3+j)*(Pi.2)^(p*j) : j in [1..M] ];



        // find higher degree terms
        if (i gt 1) and (not minimal) then

            if verb then
                print "Computing higher degree terms.";
            end if;

            // to invert
            Fi := FieldOfFractions(Pi);

            vx:=[ Fi!Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [Fi.(M+4) ];

            // xi
            vinvi:=WittInv(vx : pols:=pols)[i+1];
            num:=Pi!Numerator(vinvi);
            den:=Pi!Denominator(vinvi);

            rem_terms:=Coefficient(num,Pi.(M+4),0); // terms without xi

            degn:=Degree(Coefficient(num,Pi.(M+4),1),1);
            degd := Degree(den,Pi.1);

            for k in [((degd-degn) div p)..((i+1)*p^(i-1)-i*p^(i-2))] do
                xi +:= Coefficient(rem_terms,1,k*p+degn)*(Pi.1)^(k*p);
                yi +:= Coefficient(rem_terms,1,k*p+degn)*(Pi.2)^(k*p);
            end for;

        end if;


        // mod out by the curve
        CurveEq := (Pi.1)^2 + (Pi.2)^2 - 1 - a0*(Pi.1)^2*(Pi.2)^2;
        Qi:=quo<Pi | CurveEq>;

        va:=[ Qi!x : x in resa ] cat [Qi.3];
        vx:=[ Evaluate(x,[Qi.1,Qi.2]) : x in resx ] cat [Qi!xi];
        vy:=[ Evaluate(y,[Qi.1,Qi.2]) : y in resy ] cat [Qi!yi];
        vone:=[Qi!1] cat [ Qi!0 : j in [1..i]];

        vx2:=WittPower(vx,2 : pols:=pols);
        vy2:=[ Evaluate(x,[Qi.2,Qi.1] cat [Qi.j : j in [3..(M+4)]]) : x in vx2 ];

        // vvars:=vx cat vy;
        vvars:=vx2 cat vy2;

        // Compute the Greenberg Transform
        if verb then
            print "Computing the GT.";
        end if;

        // LeftGT:=GT( [[* vone,2,0 *], [* vone,0,2 *]] : pols:=pols, vvars:=vvars);
        // RightGT:=GT( [[* va,2,2 *], [* vone,0,0 *]] : pols:=pols, vvars:=vvars);
        LeftGT:=GT( [[* vone,1,0 *], [* vone,0,1 *]] : pols:=pols, vvars:=vvars);
        RightGT:=GT( [[* va,1,1 *], [* vone,0,0 *]] : pols:=pols, vvars:=vvars);

        difference:=RightGT[i+1] - LeftGT[i+1];

        // We need to add the extra conditions from tau(0,1)=(0,1)
        maxd := Degree(difference,1);

        // this adds that xi(1) = yi(1) = 0
        difference +:= Qi!(Evaluate(yi, Pi.2, 1)*Pi.1^(maxd + 1));


        // put in reduced form
        // difference:=reduce_pol(difference, a0, b0);

        difference := Pi!difference;

       // collect the coefficients, which then have to be zero
        eqts:=[];
        for coefx in Coefficients(difference, Pi.1) do
            for coefxy in Coefficients(coefx, Pi.2) do
                if coefxy ne 0 then
                    Append(~eqts, coefxy);
                end if;
            end for;
        end for;

        neqts:=#eqts;

        // print eqts;

        Mat:=Matrix(F,1 + M,neqts,
                [[Coefficient(eqts[j],k+2,1) : j in [1..neqts]]
                        : k in [1..(1+M)]]);

        vec:=Vector(F,
                 [ -Evaluate(eqts[j], [ 0 : k in [1..(3+M+1)] ]) : j in [1..neqts] ]);

        if verb then
            print "Solving the system";
        end if;

        vsol:=Solution(Mat,vec);

        // to convert solutions to Pres/F
        evalvec:=[x0,y0] cat [ vsol[j] : j in [1..1+M]] cat [0];

       // print #evalvec;

        Append(~resa,vsol[1]);
        Append(~resx,Evaluate(xi,evalvec));
        Append(~resy,Evaluate(yi,evalvec));

    end for;

    return resa, resx, resy;

end function;





/* function lift_ed_2(a0,n : pols:=[], verb:=true, cstr:=[], check:=true, minimal:=false) */

/*     F:=Parent(a0); */
/*     b0:=F!1; */
/*     p:=Characteristic(F); */

/*     if #pols eq 0 then */
/*             pols:=etapols(p,n); */
/*     end if; */

/*     Pres<x0,y0>:=PolynomialRing(F,2); */

/*     /\* CurveEq := b0^2*x0^2 + y0^2 - 1 - a0*x0^2*y0^2; *\/ */

/*     /\* // introduce variables for the quotient *\/ */
/*     /\* C<xx0, yy0>:=quo<Pres | CurveEq>; *\/ */

/*     /\* // let's get the conversion maps *\/ */
/*     /\* PtoC := hom< Pres -> C | [xx0, yy0]>; *\/ */
/*     /\* CtoP := hom< C -> Pres | [x0, y0]>; *\/ */

/*     resa:=[a0]; */
/*     resx:=[x0]; */
/*     resy:=[y0]; */


/*     // FIXME: do we need two? */
/*     xf:=(x0^2-1)*(a0*x0^2-1); // Important Equation used in derivative of x_n */

/*     ff:=xf^((p-1) div 2); */

/*     HI:=F!(Coefficient(ff, 1, p-1)); //Hasse Invaraint */

/*     if HI eq 0 then */
/*             return "This Curve is not Ordinary"; */
/*     end if; */


/*     //main loop */
/*     for i in [1..n] do */
/*         if verb then */
/*             print "Working on i equals: ", i; */
/*         end if; */

/*         // necessary variables */
/*         // we only need the coefficients below the degree of the derivative */
/*         // we do not need "+1" since constant terms are zero */
/*         M:=2*p^(i-1)-1; */

/*         // print "The value of M is: ", M; */

/*         Pi:=PolynomialRing(F,1 + 1 + M + 1); */
/*         // x_0, a_n, epsilon_i' (recall epsilon_0=delta_0=0) */
/*         // plus one additional variable to help with inverses in the computation */
/*         // of higher degree terms */

/*         // a temporary power of f to help compute x_i's and y_i's: */
/*         if i eq 1 then */
/*                 tmppf:=ff; */
/*         else */
/*                 tmppf:=tmppf^p*ff; */
/*         end if; */

/*         xi:= HI^(-(p^i-1) div (p-1))*tmppf - x0^(p^i-1); */
/*         if i gt 1 then */
/*                 xi-:=&+[ resx[j+1]^(p^(i-j)-1)*Derivative(resx[j+1],1) : */
/*                         j in [1..(i-1)] ]; */
/*         end if; */

/*         xi:=my_x_int(xi); */
/*         xi:=Evaluate(xi,[Pi.1,0]); */

/*         // NOT YET: will make epsilon_(p^(n-1))=0 */
/*         // we have epsilon_0 = 0 */
/*         // xi+:=&+[ Pi.(5+j)*(Pi.1)^(p*j) : j in [1..M] | j ne p^(i-1) ]; */
/*         xi+:=&+[ Pi.(2+j)*(Pi.1)^(p*j) : j in [1..M] ]; */




/*         // find higher degree terms */
/*         if i gt 1 then */

/*             if verb then */
/*                 print "Computing higher degree terms."; */
/*             end if; */

/*             // to invert */
/*             Fi := FieldOfFractions(Pi); */

/*             vx:=[ Fi!Evaluate(x,[Pi.1,0]) : x in resx ] cat [Fi.(M+3)]; */

/*             // xi */
/*             vinvi:=WittInv(vx : pols:=pols)[i+1]; */
/*             num:=Pi!Numerator(vinvi); */
/*             den:=Pi!Denominator(vinvi); */

/*             rem_terms:=Coefficient(num,Pi.(M+3),0); // terms without xi */

/*             degn:=Degree(Coefficient(num,Pi.(M+3),1),1); */
/*             degd := Degree(den,Pi.1); */

/*             for k in [((degd-degn) div p)..((i+1)*p^(i-1)-i*p^(i-2))] do */
/*                 xi +:= Coefficient(rem_terms,1,k*p+degn)*(Pi.1)^(k*p); */
/*             end for; */

/*         end if; */

/*         va:=[ Pi!x : x in resa ] cat [Pi.2]; */
/*         vx:=[ Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [xi]; */
/*         vx2:=WittPower(vx, 2 : pols:=pols); */
/*         vzero:=[ Pi!0 : j in [1..(i+1)]]; */
/*         vone:=[Pi!1] cat [ Pi!0 : j in [1..i]]; */
/*         vtwo:=[ Pi!x : x in IntToWitt(2,p,i) ]; */

/*         vvars:=vx2 cat vzero; */


/*         // Compute the Greenberg Transform */
/*         if verb then */
/*             print "Computing the GT."; */
/*         end if; */

/*         // LeftGT:=GT( [[* vtwo,2,0 *]] : pols:=pols, vvars:=vvars ); */
/*         LeftGT:=WittProd(vtwo,vx2 : pols:=pols); */
/*         RightGT:=GT( [[* va,2,0 *], [* vone,0,0 *]] : pols:=pols, vvars:=vvars ); */

/*         difference:=RightGT[i+1] - LeftGT[i+1]; */

/*         // We need to add the extra conditions from tau(0,1)=(0,1) */
/*         maxd := Degree(difference,1); */

/*         // this adds that xi(1) = yi(1) = 0 */
/*         difference +:= Evaluate(xi, Pi.1, 1)*Pi.1^(maxd + 1); */

/*         // No reduction? */

/*        // collect the coefficients, which then have to be zero */
/*         eqts:=[]; */
/*         for coefx in Coefficients(difference, Pi.1) do */
/*             if coefx ne 0 then */
/*                 Append(~eqts, coefx); */
/*             end if; */
/*         end for; */

/*         neqts:=#eqts; */

/*         // print eqts; */

/*         Mat:=Matrix(F,1 + M,neqts, */
/*                 [[Coefficient(eqts[j],k+1,1) : j in [1..neqts]] */
/*                         : k in [1..(1+M)]]); */

/*         vec:=Vector(F, */
/*                  [ -Evaluate(eqts[j], [ 0 : k in [1..(2+M+1)] ]) : j in [1..neqts] ]); */

/*         if verb then */
/*             print "Solving the system"; */
/*         end if; */

/*         vsol:=Solution(Mat,vec); */

/*         // to convert solutions to Pres/F */
/*         evalvec:=[x0] cat [ vsol[j] : j in [1..1+M]] cat [0]; */

/*        // print #evalvec; */

/*         xxi:=Evaluate(xi,evalvec); */

/*         Append(~resa,vsol[1]); */
/*         Append(~resx,xxi); */
/*         Append(~resy,Evaluate(xxi,[y0,0])); */

/*     end for; */

/*     return resa, resx, resy; */

/* end function; */
