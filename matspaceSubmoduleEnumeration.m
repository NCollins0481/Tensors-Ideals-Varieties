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


QQ := AlgebraicClosure(Rationals());
R<x,y> := PolynomialRing(QQ, 2);

ROWS:=2; //Dimensions
COLS:=3;

XEigs:=[<2,ROWS>]; //Eigenvalues and Block sizes of each operator in JNF
YEigs:=[<3,COLS>];

X := Jordan_blocks(XEigs, QQ);
Y := Jordan_blocks(YEigs, QQ);

Xmin := Evaluate(MinimalPolynomial(X), x);
Ymin := Evaluate(MinimalPolynomial(Y), y);

Xdeg := Degree(Xmin);
Ydeg := Degree(Ymin);

Mat := KMatrixSpace(QQ, ROWS, COLS);


S<z,w> := quo<R | Xmin, Ymin>;

// orbitlist := [];
// for m in Basis(Mat) do
//     orbit := [X^i * m * Y^j: i in [0..Xdeg-1], j in [0..Ydeg-1]];
//     orbitlist := orbitlist cat orbit;
// end for;

// Orbitspace := sub<Mat | orbitlist>; //Will be the whole space by definition.

// Computing orbit for a specific matrix and operators
M := Matrix(QQ, ROWS, COLS, [1,0,0, 1,0,0]);

orbit := [X^i * M * Y^j: i in [0..Xdeg-1], j in [0..Ydeg-1]];

flattened_orbit := [Eltseq(mat) : mat in orbit];
flattened_matrix := Matrix(QQ,flattened_orbit);
Rank(flattened_matrix);


// M1 := Matrix(QQ, ROWS, COLS, [0,0,0, 0,0,0]);

// orbit1 := [X^i * M1 * Y^j: i in [0..Xdeg-1], j in [0..Ydeg-1]];

// flattened_orbit1 := [Eltseq(mat) : mat in orbit1];
// flattened_matrix1 := Matrix(QQ,flattened_orbit1);

// M2 := Matrix(QQ, ROWS, COLS, [0,1,0, 27,-29,6]);

// orbit2 := [X^i * M2 * Y^j: i in [0..Xdeg-1], j in [0..Ydeg-1]];

// flattened_orbit2 := [Eltseq(mat) : mat in orbit2];
// flattened_matrix2 := Matrix(QQ,flattened_orbit2);
