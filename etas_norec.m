function rec_v(n)
    // retunrs a vector that splits [1,2,...,n]
    // into "halves" to replace the recursion
    v:=[[ i : i in [1..n] ]]; // v contains the pieces
    maxl:=n; // maximal length
    while maxl gt 2 do 
        // continue until no piece has length > 2
        i:=1; 
        while i le #v do
            // go over the pieces in v
            vv:=v[i];
            lvv:=#vv;
            if lvv gt 2 then
                // split
                v1:=[ vv[j] : j in [1..(lvv div 2)] ];
                v2:=[ vv[j] : j in [((lvv div 2)+1)..lvv] ];
                v[i]:=v1; // replace vv with v1
                Insert(~v,i+1,v2); // add v2 after v1
                i+:=1; // add one to the index as we added another
                       // piece
            end if;
            i+:=1; 
        end while;
        maxl:=Ceiling(maxl/2); // max length after splitting
    end while;
    return v;
end function;



