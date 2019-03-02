p:=31;
F<a0,b0>:=RationalFunctionField(GF(p),2);

load "fctWitt5.m";

time a1, b1, F1, H1, a2, b2, F2, H2 := lift(a0,b0 : ncoord:=3, minimal:=true, tm:=true);
