
R<x,y> := PolynomialRing(GF(31), 2);
Xmin := (x-2)^2;
Ymin := (y-3)^3;

//Get an ideal contained in (x-2, y-3)
S<z,w> := quo<R | Xmin, Ymin>;
J := ideal<S | 26*z + 28, 17*w^2 + 16*w + 26>;
var := true;
while var do
    temp1<z1> := PolynomialRing(GF(31));
    s1 := quo<temp1 | (z1-2)^2>;
    f1 := Random(s1);
    temp2<z2> := PolynomialRing(GF(31));
    s2 := quo<temp2 | (z2-3)^3>;
    f2 := Random(s2);

    J := ideal<S | Evaluate(f1, z), Evaluate(f2, w)>;
    if not ((z-2)in J) or not ((w-3) in J) then
        var := true;
    end if;
end while;
