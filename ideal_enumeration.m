Attach("example-3D.m");
// Want to enumerate paramaterizations of ideals of k[x,y,z] containing a certain ideal and contained in a certain ideal


//Functions

//e_dim is a list of eigenvalues, size of this Jordan block
//R is a ring
Jordan_blocks := function(e_dim, R)
    blocks := [];
    for eigdim in e_dim do
        eig := eigdim[1];
        dim := eigdim[2];
        eigM := Matrix(R, dim, dim, [<i,i,eig>: i in [1..dim]]) + Matrix(R, dim, dim, [<i,i+1,1>: i in [1..dim-1]]); //[<i,i,eig>: i in [1..dim], <i,i+1,1>: i in [1..dim]]
        Append(~blocks, eigM);
    end for;
JNF := DiagonalJoin(blocks);
return JNF;
end function;


// Want to enumerate paramaterizations of ideals of k[x,y,z] containing a certain ideal and contained in a certain ideal
// Begin by enumerating principle ideals and then can combine to get ideals with more generators. Just need to handle containment
// Then reverse engineer to find the corresponding submodules.
R := RationalField();
S := PolynomialRing(R, 3);

ROWS:=2; //Dimensions
COLS:=3;
DEPTH:=2;

XEigs:=[<0,ROWS>]; //Eigenvalues and Block sizes of each operator in JNF
YEigs:=[<0,COLS>];
ZEigs:=[<0,DEPTH>];

X := Jordan_blocks(XEigs, R);
Y := Jordan_blocks(YEigs, R);
Z := Jordan_blocks(ZEigs, R);

Xmin := Evaluate(MinimalPolynomial(X), S.1);
Ymin := Evaluate(MinimalPolynomial(Y), S.2);
Zmin := Evaluate(MinimalPolynomial(Z), S.3);

Xdeg := Degree(Xmin);
Ydeg := Degree(Ymin);
Zdeg := Degree(Zmin);


// axis := Ideal([Xmin, Ymin, Zmin]); //Quotient is not needed yet
// S<x,y,z> := quo<S | axis>;

//Now adding in linear dependencies
scalarvars := ["l" * IntegerToString(i) * IntegerToString(j) * IntegerToString(k): i in [0..Xdeg-1], j in [0..Ydeg-1], k in [0..Zdeg-1] | (i ne 0) or (j ne 0) or (k ne 0)];
othervars := ["x","y","z"];
allvars := othervars cat scalarvars;

S := PolynomialRing(R, #allvars);
AssignNames(~S, allvars);


// Enumerate all the monomials
Shoms := [];
for i in [0..Xdeg-1] do
    for j in [0..Ydeg-1] do
        for k in [0..Zdeg-1] do
            if (i ne 0) or (j ne 0) or (k ne 0) then
                Append(~Shoms, S.1^i * S.2^j * S.3^k);
            end if;
        end for;
    end for;
end for;

values := [0,1];
indexlist := [elts : elts in CartesianPower(values, #scalarvars)];


//Attach the products to every element in Shoms
//Create list of monomial ideals of form l_ijk^indexlist[i,j,k]
//Turn list of monomials into poset by declaring equality of ideals
monomiallist := [];
for index in indexlist do
    poly := &+[S.(3+i) * (index[i]) * Shoms[i]: i in [1..#index]];
    Append(~monomiallist, ideal<S | poly>);
end for;

// Rel := [<i, j> : i in monomiallist, j in monomiallist | i ne j and (i subset j)]; //Magma doesn't seem to have a good package for posets, this is slow, and this is not needed yet.
// P := Poset(monomiallist, Rel);


// Want to enumerate submodules of k^2 \otimes k^3 \otimes k^2 of which are acted on by k[x,y,z]/(min polys) via the Jordan operators
// Begin by enumerating the images under the operators and then take linear combinations. 
// Afterwards test this with the 2x3 example I already have.

coordinates := ["t" * IntegerToString(i) * IntegerToString(j) * IntegerToString(k): i in [0..Xdeg-1], j in [0..Ydeg-1], k in [0..Zdeg-1] | (i ne 0) or (j ne 0) or (k ne 0)];
allvars := coordinates cat scalarvars;

P := PolynomialRing(R, #allvars);
AssignNames(~P, allvars);

//Create a tensor in these symbols
// slice1 := Matrix(ROWS, COLS, [0]cat[P.i : i in [1..((#coordinates - 1)/2)]]); //adding 0 to the top left
// slice2 := Matrix(ROWS, COLS, [P.i : i in [((#coordinates + 1)/2)..#coordinates]]); //adding the rest of the entries
MS := RMatrixSpace(P, ROWS,COLS);

values := [0,1];
indexlist := [elts : elts in CartesianPower(values, #coordinates)];

tensorlist := [];
for index in indexlist do
    slice1 := Matrix(ROWS, COLS, [0]cat[index[i] * P.i * P.(#coordinates + i) : i in [1..((#coordinates - 1)/2)]]); //adding 0 to the top left
    slice2 := Matrix(ROWS, COLS, [index[i] * P.i * P.(#coordinates + i) : i in [((#coordinates + 1)/2)..#coordinates]]); //adding the rest of the entries
    M1 := MS!slice1;
    M2 := MS!slice2;
    T := [M1, M2];
    Append(~tensorlist, T);
end for;




 // Can alternatively inductively create these tensors. From a dimension i subspace with no related entries, to get a dimension i+1 subspace we can either:
 //     add a single cube with new color
 //     add a pair of cubes to different axes with a single relation
 //     add a triple of cubes to all axes with all pairwise related

 // If a dimension i subspace has related entries, then cubes can only be added outside with the induced relations.
