// Goal: Want to enumerate paramaterizations of the orbits

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

R := AlgebraicClosure(Rationals());

ROWS:=2; //Dimensions
COLS:=3;
//DEPTH:=2;


// coordinates := ["t" * IntegerToString(i) * IntegerToString(j): j in [0..COLS-1], i in [0..ROWS-1]]; //General
coordinates := ["a", "b", "c", "d", "e", "f"];
scalarvars := ["l" * IntegerToString(i) * IntegerToString(j): j in [0..COLS-1], i in [0..ROWS-1]]; // | (i ne 0) or (j ne 0)];
othervars := ["x","y"];
allvars := othervars cat coordinates cat scalarvars;

S := PolynomialRing(R, #allvars);
AssignNames(~S, allvars);

XEigs:=[<0,ROWS>]; //Eigenvalues and Block sizes of each operator in JNF
YEigs:=[<0,2>, <0,1>];
// ZEigs:=[<0,DEPTH>];

X := Jordan_blocks(XEigs, R);
Y := Jordan_blocks(YEigs, R);
// Z := Jordan_blocks(ZEigs, R);

Xmin := Evaluate(MinimalPolynomial(X), S.1);
Ymin := Evaluate(MinimalPolynomial(Y), S.2);
// Zmin := Evaluate(MinimalPolynomial(Z), S.3);

X := Jordan_blocks(XEigs, S); //Magma does not like minimal polynomials of matrices over polynomial rings
Y := Jordan_blocks(YEigs, S);

Xdeg := Degree(Xmin);
Ydeg := Degree(Ymin);
// Zdeg := Degree(Zmin);

//Compute X action on M
M := Matrix(S, ROWS, COLS, [S.i: i in [1+#othervars..(ROWS*COLS + #othervars)]]);

orbit := [X^i * M * Y^j: i in [0..Xdeg-1], j in [0..Ydeg-1]];

flattened_orbit := [Eltseq(mat) : mat in orbit];
flattened_matrix := Matrix(S,flattened_orbit);

minors := { m: m in Minors(flattened_matrix, 3)}; //matrix has rank <= r iff r+1 times r+1 minors all vanish
I := Ideal(minors);
G := GroebnerBasis(I);
print(G);
print(M);

aff := Spec(S);
asols := Scheme(aff, G);

elimvars := {}; //Key: 3 = a, 6 = d, 8 = f.
EliminationIdeal(I, {S.i: i in {1..#allvars} | not (i in elimvars)}); //Eliminates only the elimvars vars.

//M2 := M - Matrix(S, ROWS, COLS, [<2,1,S.6>, <2,3,S.8>]); //subtracting d and f = 0
//orbit2 := [X^i * M2 * Y^j: i in [0..Xdeg-1], j in [0..Ydeg-1]];                 
//flattened_orbit2 := [Eltseq(mat) : mat in orbit2];
//flattened_matrix2 := Matrix(S,flattened_orbit2);
//Rank(flattened_matrix2);


// proj := Proj(S);
// psols := Scheme(proj, G);


// // axis := Ideal([Xmin, Ymin, Zmin]); //Quotient is not needed yet
// // S<x,y,z> := quo<S | axis>;

// //Now adding in linear dependencies
// scalarvars := ["l" * IntegerToString(i) * IntegerToString(j): i in [0..Xdeg-1], j in [0..Ydeg-1]| (i ne 0) or (j ne 0)];
// othervars := ["x","y"];
// allvars := othervars cat scalarvars;

// S := PolynomialRing(R, #allvars);
// AssignNames(~S, allvars);

// coordinates := ["t" * IntegerToString(i) * IntegerToString(j): i in [0..Xdeg-1], j in [0..Ydeg-1]| (i ne 0) or (j ne 0)];
// allvars := coordinates cat scalarvars;

// P := PolynomialRing(R, #allvars);
// AssignNames(~P, allvars);

// // Enumerate all the monomials
// Shoms := [];
// for i in [0..Xdeg-1] do
//     for j in [0..Ydeg-1] do
//         if (i ne 0) or (j ne 0) then
//             Append(~Shoms, S.1^i * S.2^j);
//         end if;
//     end for;
// end for;

// InitialTensor := P.(ROWS*COLS -1);

// // Then want to define a function which inductively builds the tensors out by adding a single block
// // Then want to have the function include whether entries are linearly independent or not.








// // Grassmannian;
// // INPUT: two natural numbers k,n, and a field FF;
// // OUTPUT: a pair (Gr,pl), where Gr is the Grassmannian of k-planes in P^n over FF, and pl is its Plücker embedding.

// Grassmannian := function(k,n,FF)
// A<[x]> := AffineSpace(FF,(k+1)*(n-k));
// RA := CoordinateRing(A);
// P<[z]> := ProjectiveSpace(FF,Binomial(n+1,k+1)-1);
// Id := DiagonalMatrix(RA,k+1,[1: i in [1 .. k+1]]);
// M2 := Matrix(RA,k+1,n-k,&cat[[x[i+(k+1)*(j-1)]: i in [1 .. k+1]]: j in [1 .. n-k]]);
// M := HorizontalJoin(Id,M2);
// pl :=map<A->P|Minors(M,k+1)>;
// Gr :=pl(A);
// return Gr,pl;
// end function;


// // Example (Grassmannian of lines in P^3)

// Gr,pl := Grassmannian(1,3,Rationals()); 