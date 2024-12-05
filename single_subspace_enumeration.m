// Sanity test: Enumerate the subspaces when operators are single nilpotent Jordan blocks.
Jordan_blocks := function(e_dim, R)
    blocks := <>;
    for eigdim in e_dim do
        eig := eigdim[1];
        dim := eigdim[2];
        eigM := Matrix(R, dim, dim, [<i,i,eig>: i in [1..dim]]) + Matrix(R, dim, dim, [<i,i+1,1>: i in [1..dim-1]]); //[<i,i,eig>: i in [1..dim], <i,i+1,1>: i in [1..dim]]
        Append(~blocks, eigM);
    end for;
    JNF := DiagonalJoin(blocks);
    return JNF;
end function;

VariablesInMonomial := function(m, R)
    exponents := Exponents(m);
    vars := [];
    for i in [1..#exponents] do
        if exponents[i] ne 0 then
            Append(~vars, R.i);
        end if;
    end for;
    return vars;
end function;

ComputeMinors := function(matrix, int, rvars, R)
    minors := { m: m in Minors(matrix, int)}; //matrix has rank <= r iff r+1 times r+1 minors all vanish
    I := Ideal(minors);
    G := GroebnerBasis(I);
    elimvars := []; //Key: 3 = a, 6 = d, 8 = f.
    for gen in G do
        if #Monomials(gen) eq 1 then
            vars := VariablesInMonomial(Monomials(gen)[1], R);
            if #vars eq 1 then
                Append(~elimvars, vars[1]);
            end if;
        end if;
    end for;
    return EliminationIdeal(I, {i: i in rvars | not (i in elimvars)}), elimvars; //Eliminates the x^k variables and returns these coordinates.
end function;

R := AlgebraicClosure(Rationals());

ROWS:=2; //Dimensions
COLS:=3;

// coordinates := ["t" * IntegerToString(i) * IntegerToString(j): j in [0..COLS-1], i in [0..ROWS-1]]; //General
coordinates := ["a", "b", "c", "d", "e", "f"]; //Specific to this 2x3 matrix small example
scalarvars := ["l" * IntegerToString(i) * IntegerToString(j): j in [0..COLS-1], i in [0..ROWS-1]]; // | (i ne 0) or (j ne 0)];
othervars := ["x","y"];
allvars := othervars cat coordinates cat scalarvars;

S := PolynomialRing(R, #allvars);
AssignNames(~S, allvars);

XEigs:=[<0,ROWS>]; //Eigenvalues and Block sizes of each operator in JNF
YEigs:=[<0,COLS>];

X := Jordan_blocks(XEigs, R);
Y := Jordan_blocks(YEigs, R);

Xmin := Evaluate(MinimalPolynomial(X), S.1);
Ymin := Evaluate(MinimalPolynomial(Y), S.2);

X := Jordan_blocks(XEigs, S);
Y := Jordan_blocks(YEigs, S);

Xdeg := Degree(Xmin);
Ydeg := Degree(Ymin);


//Goal is to repeat this calculation now producing the submodule generators.
// This will become a GetGens function:


M := Matrix(S, ROWS, COLS, [S.i: i in [1+#othervars..(ROWS*COLS + #othervars)]]); //want to make this print the submodule generators.
print(M);

mat := KMatrixSpace(R, ROWS, COLS);

for m in Basis(mat) do
    orbitlist := [];
    Append(~orbitlist, [X^i * m * Y^j: i in [0..Xdeg-1], j in [0..Ydeg-1]]);
end for;

orbit := [X^i * M * Y^j: i in [0..Xdeg-1], j in [0..Ydeg-1]];

flattened_orbit := [Eltseq(mat) : mat in orbit];
flattened_matrix := Matrix(S,flattened_orbit);

svars := [S.i: i in [1..Rank(S)]];

// aff := Spec(S);
// asols := Scheme(aff, G);

// To be able to get the rest, ask for generators of these submodules and then have subspaces paramaterized by an arbitrary linear combination of basis vectors

dimeqs := [];

ksubmodules := function(k)
    minors := { m: m in Minors(flattened_matrix, k+1)}; //matrix has rank <= r iff r+1 times r+1 minors all vanish
    eqtns, elimvars := ComputeMinors(flattened_matrix, k, svars, S);
    return eqtns, elimvars;
end function;

allsubmodules := function()
    for k in [1..#orbit] do //start with the 2x2 minors. Ignore the entire space.
        print "dimension", k, ":";
        eqtns, elimvars:= ksubmodules(k);
        print "equations:", eqtns, "eliminated variables", elimvars;
    end for;
end function;
