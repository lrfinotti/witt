load 'lift_ed.m';

p:=11;

n:=2;
F<a0,b0>:=RationalFunctionField(GF(p), 2);

pols := etapols(p,n);


va, vb, vx, vy := lift_ed(a0, b0, n);

v14 := IntToWitt(14, p, n);
v16 := IntToWitt(16, p, n);
va2 := WittPower(va, 2: pols:=pols);
vb2 := WittPower(vb, 2: pols:=pols);
vb4 := WittPower(vb, 4: pols:=pols);

num:=WittProd(va, vb2: pols:=pols);
num:=WittProd(v14, num: pols:=pols);
num:=WittSum(va2, num: pols:=pols);
num:=WittSum(vb4, num: pols:=pols);
num:=WittPower(num, 3: pols:=pols);
num:=WittProd(v16, num: pols:=pols);

den:=WittDiff(va, vb2: pols:=pols);
den:=WittPower(den, 4: pols:=pols);
den:=WittProd(va, den: pols:=pols);
den:=WittProd(vb2, den: pols:=pols);

j:=WittDiv(num,den: pols:=pols);

j0 := j[1];


if p eq 5 then
    j1 := 3*j0^3 + j0^4;
    j2 := 3*j0^5 + 2*j0^10 + 2*j0^13 + 4*j0^14 + 4*j0^15 + 4*j0^16 +
    j0^17 + 4*j0^18 + j0^19 + j0^20 + 3*j0^23 + j0^24;
elif p eq 7 then
    j1 := 3*j0^5+5*j0^6;
    j2 := (3*j0^21 + 6*j0^28 + 3*j0^33 + 5*j0^34 + 4*j0^35 + 2*j0^36 +
           3*j0^37 + 6*j0^38 + 3*j0^39 + 5*j0^40 + 5*j0^41 + 5*j0^42 + 2*j0^43
           + 3*j0^44 + 6*j0^45 + 3*j0^46 + 5*j0^47 + 5*j0^48 + 3*j0^49 +
           3*j0^54 + 5*j0^55)/(1 + j0^7);
elif p eq 11 then
    j1 := 8*j0^7 + 5*j0^8 + 8*j0^9 + 4*j0^10;
    j2 :=(7*j0^33 + 10*j0^44 + 3*j0^55 + 4*j0^66 + 4*j0^73 + 8*j0^74 +
          4*j0^75 + 2*j0^76 + 7*j0^77 + 3*j0^78 + 4*j0^79 + 10*j0^80 +
          6*j0^81 + 10*j0^82 + 8*j0^83 + 3*j0^84 + 2*j0^86 + 4*j0^88 +
          7*j0^89 + j0^90 + 3*j0^91 + j0^92 + 3*j0^93 + 3*j0^95 + 7*j0^96 +
          3*j0^97 + 10*j0^98 + 4*j0^99 + 5*j0^100 + 2*j0^101 + 5*j0^102 +
          6*j0^103 + 9*j0^104 + 10*j0^105 + 7*j0^106 + 4*j0^108 + 9*j0^109 +
          7*j0^111 + 4*j0^112 + 4*j0^113 + 9*j0^114 + 4*j0^116 + 9*j0^117 +
          4*j0^118 + 2*j0^119 + 3*j0^120 + 10*j0^121 + 7*j0^128 + 3*j0^129 +
          7*j0^130 + 9*j0^131)/(6 + 5*j0^11);
elif p eq 13 then
    j1 :=(9*j0^9 + 10*j0^10 + j0^11 + 3*j0^12 + 6*j0^13)/(10 + 11*j0);


    j2 :=(10*j0^195 + 5*j0^194 + 6*j0^193 + 8*j0^192 + 2*j0^191 + 6*j0^183 +
          j0^182 + 3*j0^181 + 7*j0^180 + 3*j0^179 + 11*j0^178 + 9*j0^177 +
          6*j0^176 + 3*j0^175 + 5*j0^174 + 6*j0^173 + 3*j0^172 + 9*j0^171 +
          7*j0^170 + 4*j0^169 + 2*j0^167 + 7*j0^166 + 7*j0^165 + 4*j0^164 +
          8*j0^163 + 4*j0^161 + 11*j0^160 + 4*j0^159 + 5*j0^158 + 7*j0^157 +
          12*j0^156 + 11*j0^155 + 11*j0^154 + 2*j0^153 + 2*j0^152 + 6*j0^151
          + 8*j0^150 + 9*j0^149 + 12*j0^148 + 4*j0^147 + 11*j0^146 +
          9*j0^145 + 11*j0^144 + 7*j0^143 + 3*j0^142 + 10*j0^140 + 5*j0^139
          + 5*j0^138 + 4*j0^137 + 5*j0^136 + 8*j0^135 + 7*j0^134 + 9*j0^133
          + 10*j0^132 + 4*j0^131 + 9*j0^130 + j0^129 + 7*j0^128 + 3*j0^127 +
          3*j0^125 + 9*j0^123 + 2*j0^122 + 11*j0^120 + 5*j0^119 + 5*j0^116 +
          6*j0^115 + 8*j0^114 + 2*j0^113 + 4*j0^105 + 6*j0^104 + 12*j0^92 +
          5*j0^91 + 2*j0^79 + 3*j0^78 + 10*j0^66 + 2*j0^65)/(j0^27 + 8*j0^26
                                                             + 3*j0^14 + 11*j0^13 + 12*j0 + 5);
end if;


jj := [j0, j1, j2];

if j1 eq j[2] then
    print "j1 is correct";
else
    print "j1 FAILS!!!!!";
end if;
if j2 eq j[3] then
    print "j2 is correct";
else
    print "j2 FAILS!!!!!";
end if;

P := Parent(vx[1]);
F := FieldOfFractions(P);

vvx := [ F!x : x in vx ];
vvy := [ F!x : x in vy ];

xninv := WittInv(vvx: pols:=pols)[n+1];
yninv := WittInv(vvy: pols:=pols)[n+1];

if Degree(Numerator(xninv),P.1) lt Degree(Denominator(xninv),P.1) then
    print "1/x works";
else
    print "1/x FAILS!!!!";
end if;

if Degree(Numerator(yninv),P.2) lt Degree(Denominator(yninv),P.2) then
    print "1/y works";
else
    print "1/y FAILS!!!!";
end if;

x1 := vx[2]; x2 := vx[3];
y1 := vy[2]; y2 := vy[3];

if [0, Evaluate(x1, [0,1]), Evaluate(x2, [0,1])] eq [0, 0 , 0] then
    print "x(0, 1) = 0 works";
else
    print "x(0, 1) = 0  FAILS!!!!";
end if;

if [1, Evaluate(y1, [0,1]), Evaluate(y2, [0,1])] eq [1, 0 , 0] then
    print "y(0, 1) = 1 works";
else
    print "y(0, 1) = 1  FAILS!!!!";
end if;


if [1/b0, Evaluate(x1, [1/b0,0]), Evaluate(x2, [1/b0,0])] eq WittInv(vb : pols:=pols) then
    print "x(1/b0, 0) = 1/b works";
else
    print "x(1/b0, 0) = 1/b  FAILS!!!!";
end if;

if [0, Evaluate(y1, [1/b0,0]), Evaluate(y2, [1/b0,0])] eq [0, 0 , 0] then
    print "y(1/b0, 0) = 0 works";
else
    print "y(1/b0, 0) = 0  FAILS!!!!";
end if;
