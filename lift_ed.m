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
            res+:=Coefficient(f,var,d)/(d+1)*P.1^(d+1);
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


function reduce_pol(f, a0, b0)
    // reduce polynomial f module b0*x0^2 + y0^2 - 1 - a0*x0^2y0^2
    // the result should be of the form f1(x0) + f2(y0) + x0*f3(y0) + y0*f4(x0)
    //
    // NOTE: this seems to be the same as reducing modulo the equation in Magma,
    // which should be faster!!!
    x0 := Parent(f).1;
    y0 := Parent(f).2;

    repl_pol := b0/a0*x0^2 + y0^2/a0 - 1/a0;

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

	F:=Parent(a0);
	p:=Characteristic(F);

    	if #pols eq 0 then
        	pols:=etapols(p,n);
    	end if;

	Pres<x0,y0>:=PolynomialRing(F,2);

	CurveEq := b0*x0^2 + y0^2 - 1 - a0*x0^2*y0^2;

        // introduce variables for the quotient
	C<xx0, yy0>:=quo<Pres | CurveEq>;

        // let's get the conversion maps
        PtoC := hom< Pres -> C | [xx0, yy0]>;
        CtoP := hom< C -> Pres | [x0, y0]>;

	resa:=[a0];
	resb:=[b0];
	resx:=[x0];
	resy:=[y0];


        // FIXME: do we need two?
	xf:=(b0*x0^2-1)*(a0*x0^2-1); // Important Equation used in derivative of x_n
	yg:=(b0-a0*y0^2)*(1-y0^2);   // Important Equation used in derivative of y_n
                                     // and computing Hasse Invariant

	ff:=xf^((p-1) div 2);
	gg:=yg^((p-1) div 2); //Polynomial for Hasse Invaraint

	HI:=F!(Coefficient(gg, 2, p-1)); //Hasse Invaraint

	if HI eq 0 then
		return "This Curve is not Ordinary";
	else
		HI;
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
            if i gt 1 then

                if verb then
                    print "Coputing higher degree terms.";
                end if;

                vx:=[ Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [Pi.(2*M+5)];
                vy:=[ Evaluate(y,[Pi.1,Pi.2]) : y in resy ] cat [Pi.(2*M+5)];

                // xi
                vinvi:=WittInv(vx, pols:=pols)[i+1];
                num:=Numerator(vinvi);
                rem_terms:=Coefficient(num,Pi.(2*M+5),0); // terms without xi

                for k in [2*p^(i-1)..((i+1)*p^(i-1)-i*p^(i-2))] do
                    xi += -Coefficient(rem_terms,1,k*p+(i-1)*p^i)*(P.1)^(k*p);
                end for;

                // yi
                vinvi:=WittInv(vy, pols:=pols)[i+1];
                num:=Numerator(vinvi);
                rem_terms:=Coefficient(num,Pi.(2*M+5),0); // terms without xi

                for k in [2*p^(i-1)..((i+1)*p^(i-1)-i*p^(i-2))] do
                    yi += -Coefficient(rem_terms,2,k*p+(i-1)*p^i)*(P.2)^(k*p);
                end for;

            end if;


            va:=[ Pi!x : x in resa ] cat [Pi.3];
            vb:=[ Pi!x : x in resb ] cat [Pi.4];
            vx:=[ Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [xi];
            vy:=[ Evaluate(y,[Pi.1,Pi.2]) : y in resy ] cat [yi];
            vone:=[Pi!1] cat [ Pi!0 : j in [1..i]];

            vvars:=vx cat vy;


            // Compute the Greenberg Transform
            if verb then
                print "Computing the GT.";
            end if;

            LeftGT:=GT( [[* vb,2,0 *], [* vone,0,2 *]] : pols:=pols, vvars:=vvars);
            RightGT:=GT( [[* va,2,2 *], [* vone,0,0 *]] : pols:=pols, vvars:=vvars);

            difference:=RightGT[i+1] - LeftGT[i+1];

            // put in reduced form
            difference:=reduce_pol(difference, a0, b0);
            eqts:=[ MonomialCoefficient(mon) : mon in Monomials(difference) ];

            neqts:=#eqts;

            // STOPPED HERE!


            Mat:=Matrix(Parent(vrem[1]),1+M+M,neqts,
                    [[Coefficient(vrem[j],k+2,1) : j in [1..neqts]]
                            : k in [ ii : ii in [1..(2+M+M)] | ii ne 2+p^(i-1) ]]);

            vec:=Vector(Parent(vrem[1]),
                     [ -Evaluate(vrem[j], [ 0 : k in [1..(2+2+M+M)]]) : j in [1..neqts]]);

            //Mat;
            //vec;
            //vrem;

            //return ":)";

            vsol:=Solution(Mat,vec);

            //print vsol;

            // to convert solutions to Pres/F
            evalvec:=[Pi.1,Pi.2] cat [ vsol[j] : j in [1..2+p^(i-1)] ] cat [0] cat
                     [vsol[j] : j in [2+p^(i-1)+1..(1+M+M)]];

           // print #evalvec;



            Append(~resa,vsol[1]);
            Append(~resb,vsol[2]);


            //These are the polynomials of interest. Due to Magma Multivariate stuff, need to reindex everything to make $.1=x0 and $.2=y0.
            //If desiring to see unintelligible Magma output, un-comment the two lines below. Call it a right of passage haha.
            //Evaluate(xi,evalvec);
            //Evaluate(yi,evalvec);

            //Appropriate reindexing done below

            Append(~resx,Evaluate(Evaluate(xi,evalvec),[x0,y0] cat [ 0 : i in [1..2 + M + M] ]));

            Append(~resy,Evaluate(Evaluate(yi,evalvec),[x0,y0] cat [ 0 : i in [1..2 + M + M] ] ));

	end for;

	return resa, resb, resx, resy;

end function;
