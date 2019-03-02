load "fctWitt-GTnew.m";

n:=3;
p:=5;

F<a0,b0>:=RationalFunctionField(GF(p),2);

D0:=4*a0^3 + 27*b0^2;

a1 := (a0^3*b0^2 + b0^4)/a0;

b1 := 4*a0^6*b0 + a0^3*b0^3 + b0^5;


a2 := (2*a0^36 + a0^33*b0^2 + a0^30*b0^4 + 3*a0^27*b0^6 + 2*a0^24*b0^8
    + a0^18*b0^12 + 4*a0^12*b0^16 + 3*a0^9*b0^18 + 4*a0^6*b0^20 +
      4*a0^3*b0^22 + 4*b0^24)/a0^11;

b2 := a0^36*b0 + 4*a0^33*b0^3 + 3*a0^27*b0^7 + 4*a0^21*b0^11 +
4*a0^15*b0^15 + a0^12*b0^17 + 3*a0^6*b0^21 + b0^25;


a3 := (2*a0^225 + 2*a0^219*b0^4 + a0^216*b0^6 + 3*a0^213*b0^8 +
        a0^210*b0^10 + a0^207*b0^12 + a0^204*b0^14 + 3*a0^198*b0^18 +
        4*a0^195*b0^20 + 4*a0^192*b0^22 + 4*a0^189*b0^24 +
        3*a0^186*b0^26 + 4*a0^183*b0^28 + a0^180*b0^30 +
        3*a0^177*b0^32 + 3*a0^174*b0^34 + 3*a0^171*b0^36 +
        4*a0^168*b0^38 + 2*a0^165*b0^40 + 2*a0^162*b0^42 +
        2*a0^159*b0^44 + 3*a0^156*b0^46 + 3*a0^153*b0^48 +
        a0^150*b0^50 + 2*a0^147*b0^52 + 3*a0^141*b0^56 + a0^138*b0^58
        + 2*a0^132*b0^62 + 4*a0^129*b0^64 + 3*a0^126*b0^66 +
        4*a0^123*b0^68 + 2*a0^120*b0^70 + 3*a0^117*b0^72 +
        a0^114*b0^74 + 3*a0^111*b0^76 + 4*a0^108*b0^78 + a0^102*b0^82
        + 4*a0^99*b0^84 + a0^96*b0^86 + 4*a0^93*b0^88 + a0^90*b0^90 +
        2*a0^87*b0^92 + 2*a0^84*b0^94 + 4*a0^81*b0^96 + 4*a0^78*b0^98
        + a0^75*b0^100 + 2*a0^72*b0^102 + a0^63*b0^108 + a0^60*b0^110
        + a0^57*b0^112 + a0^54*b0^114 + 2*a0^51*b0^116 +
        2*a0^48*b0^118 + 3*a0^45*b0^120 + a0^42*b0^122 + a0^39*b0^124
      + 2*a0^30*b0^130 + a0^15*b0^140 + 3*b0^150)/a0^100;

b3 := (4*a0^261*b0 + a0^258*b0^3 + 4*a0^255*b0^5 + 2*a0^252*b0^7 +
        4*a0^249*b0^9 + 2*a0^246*b0^11 + 3*a0^243*b0^13 +
        4*a0^240*b0^15 + 2*a0^237*b0^17 + 3*a0^228*b0^23 +
        2*a0^225*b0^25 + a0^222*b0^27 + 2*a0^216*b0^31 +
        4*a0^213*b0^33 + 3*a0^207*b0^37 + a0^204*b0^39 +
        2*a0^201*b0^41 + 4*a0^198*b0^43 + a0^192*b0^47 +
        3*a0^189*b0^49 + 3*a0^186*b0^51 + 4*a0^183*b0^53 +
        a0^180*b0^55 + 4*a0^177*b0^57 + 2*a0^174*b0^59 + a0^168*b0^63
        + 2*a0^165*b0^65 + 3*a0^162*b0^67 + 2*a0^159*b0^69 +
        3*a0^156*b0^71 + a0^153*b0^73 + 2*a0^147*b0^77 +
        3*a0^144*b0^79 + 3*a0^141*b0^81 + 2*a0^135*b0^85 +
        3*a0^132*b0^87 + a0^129*b0^89 + 3*a0^126*b0^91 +
        2*a0^123*b0^93 + a0^117*b0^97 + 3*a0^111*b0^101 +
        3*a0^105*b0^105 + 2*a0^99*b0^109 + 4*a0^96*b0^111 +
        4*a0^93*b0^113 + 3*a0^90*b0^115 + 3*a0^87*b0^117 +
        2*a0^75*b0^125 + a0^66*b0^131 + a0^63*b0^133 + 3*a0^57*b0^137
        + a0^54*b0^139 + 4*a0^51*b0^141 + 4*a0^48*b0^143 +
        a0^45*b0^145 + 4*a0^42*b0^147 + 2*a0^39*b0^149 +
        2*a0^36*b0^151 + a0^33*b0^153 + 2*a0^30*b0^155 +
        2*a0^27*b0^157 + a0^24*b0^159 + 2*a0^21*b0^161 +
      3*b0^175)/a0^75;

// VERY SLOW!!!  Better copy and paste...
// va, vb:=lift(a0,b0,n);


va:=[F.1,a1,a2,a3];
vb:=[F.2,b1,b2,b3];
v4:=IntToWitt(4,p,n+1);
v27:=IntToWitt(27,p,n+1);

D:=GT([[* v4,3,0 *], [* v27,0,2 *]] : vvars:=(va cat vb));

for i in [1..n] do
    Di:=F!Numerator(D[i+1]);
    //print Di;
    print Denominator(Di/F!D0);
end for;
    
