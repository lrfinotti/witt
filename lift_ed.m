// Computes the canonical lifting of an Elliptic Curve in Edwards Form by way of computing the coefficients

load "Finminlift.m";
load "polyDiv.m";

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


function my_x_int(f)
    // formal integral of polynomials
    deg:=Degree(f);
    P:=Parent(f);
    p:=Characteristic(P);
    res:=0;
    for d in [0..deg] do
        if ((d+1) mod p) ne 0 then
            res+:=Coefficient(f,1,d)/(d+1)*P.1^(d+1);
        end if;
    end for;
    return res;

end function;

function my_y_int(f)
    // formal integral of polynomials
    deg:=Degree(f);
    P:=Parent(f);
    p:=Characteristic(P);
    res:=0;
    for d in [0..deg] do
        if ((d+1) mod p) ne 0 then
            res+:=Coefficient(f,2,d)/(d+1)*P.2^(d+1);
        end if;
    end for;
    return res;

end function;






function LiftEd(a0,b0,n : pols:=[], verb:=true, cstr:=[], check:=true, minimal:=false)

	F:=Parent(a0);
	p:=Characteristic(F);

    	if #pols eq 0 then
        	pols:=etapols(p,n);
    	end if;

	Pres<x0,y0>:=PolynomialRing(F,2);

	CurveEq:=b0*x0^2+y0^2-1-a0*x0^2*y0^2;

	C:=quo<Pres | CurveEq>;

	resa:=[a0];
	resb:=[b0];
	resx:=[x0];
	resy:=[y0];

	xf:=(b0*x0^2-1)*(a0*x0^2-1); //Important Equation used in derivative of x_n

	yg:=(b0-a0*y0^2)*(1-y0^2); //Important Equation used in derivative of y_n and computing Hasse Invariant

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

	print "Working on i equals: ", i;

	M:=2*p^(i-1)-1;

	print "The value of M is: ", M;

	Pi:=PolynomialRing(F,2 + 2 + M + M); // x_0, y_0, a_n, b_n, epsilon_i's, delta_i's (recall epsilon_0=delta_0=0)

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

	// will make epsilon_(p^(n-1))=0
        xi+:=&+[ Pi.(5+j)*(Pi.1)^(p*j) : j in [0..M] | j ne p^(i-1) ];




        // //////////////////////////////////////////////////////////
        // Non-Minimal Conditions.

        if (not minimal) and (i eq 2) then
            tmp:= resx[2]^2;
		for i in [2*p^2 .. 3*p^2-2*p by p] do
			xi+:=Pi.1^i*Pi!(Coefficient(tmp,1, Integers()!(p + i/p)));
		end for;
        end if;

        if (not minimal) and (i eq 2) then
            tmp:= resy[2]^2;
		for i in [2*p^2 .. 3*p^2-2*p by p] do
			yi+:=Pi.2^i*Pi!(Coefficient(tmp,2, Integers()!(p + i/p)));
		end for;
        end if;

        // Condition for i ge 3
        // if (not minimal) and (i ge 3) then
        //
        // end if;



        // //////////////////////////////////////////////////////////

	va:=[ Pi!x : x in resa ] cat [Pi.3];
        vb:=[ Pi!x : x in resb ] cat [Pi.4];
        vx:=[ Evaluate(x,[Pi.1,Pi.2]) : x in resx ] cat [xi];
        vy:=[ Evaluate(y,[Pi.1,Pi.2]) : y in resy ] cat [yi];
        vone:=[Pi!1] cat [ Pi!0 : j in [1..i]];

	vvars:=vx cat vy;


	///Get all the Greenberg Transform Stuff

	LeftGT:=GT( [[* vb,2,0 *], [* vone,0,2 *]] : pols:=pols, vvars:=vvars);

	RightGT:=GT( [[* va,2,2 *], [* vone,0,0 *]] : pols:=pols, vvars:=vvars);

	RHS:=RightGT[i+1];

	LHS:=LeftGT[i+1];

	tmpLHS:=Coefficient(LHS,2,0);


	LHS:= LHS - tmpLHS; //LHS now only has terms in that are multiples of y0

	RHS:= RHS - tmpLHS; //Since RHS must be equal to LHS, it must be divisible by y0.

	SystemEq:=RHS;

	vrem:=Coefficients(SystemEq,2); //Gather all the terms in y0

        // matrix of coefficients

	neqts:=#vrem;

	//print "neqts equals: ", neqts;


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
