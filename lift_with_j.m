 load 'witt_gen.m';

 //choose p
 p:=13;


 F:=GF(p);
 FF<a0,b0>:=RationalFunctionField(F,2);


j0:=1728*(4*a0^3)/(4*a0^3+27*b0^2);

j1fct := function(p);
    if p eq 5 then
        return 3*j0^3 + j0^4;
    end if;
    if p eq 7 then
        return 3*j0^5+5*j0^6;
    end if;    
    if p eq 11 then
        return 8*j0^7 + 5*j0^8 + 8*j0^9 + 4*j0^10;
    end if;    
    if p eq 13 then
        return (9*j0^9 + 10*j0^10 + j0^11 + 3*j0^12 + 6*j0^13)/(10 + 11*j0);
    end if;    
    if p eq 17 then
        return (15*j0^11 + 16*j0^12 + 12*j0^13 + 15*j0^14 + 5*j0^15 + 10*j0^16 +
    10*j0^17)/(14 + 11*j0);
    end if;    
    if p eq 19 then
        return (j0^13 + 6*j0^14 + 8*j0^15 + 18*j0^16 + 17*j0^17 + 4*j0^18 +
    7*j0^19)/(10 + 4*j0);
    end if;    
end function;

j1:=j1fct(p);

j:=[j0,j1];

// compute the coefficients
num:=WittProd(IntToWitt(27,p,2),j);
den:=WittSum(IntToWitt(4*1728,p,2),WittNeg(WittProd(IntToWitt(4,p,2),j)));

jj:=WittQuot(num,den); // the coefficients!

lambda2:=[b0/a0,0]; // factor to switch back to a0 and b0
A1:=WittProd(WittPwr(lambda2,2),jj);
B1:=WittProd(WittPwr(lambda2,3),jj);

print(A1);
print(B1);
