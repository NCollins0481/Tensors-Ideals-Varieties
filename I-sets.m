// TIZ set example computation for 3-tensors

R := Integers();
P<x,y, z> := PolynomialRing(R, 3);

// S: list of matricies.
// z: |S| x |S| matrix.
// Produces a contraction on the slices of the matricies. (i.e S'[i] = \sum_j S[j]z_{ij})
outer_action := function (S, z)
     T := [ ];
     for i in [1..#S] do 
          // Magma notation: &+ sums every element in a list.
          F := &+ [ z[i][j] * S[j] : j in [1..#S] ];
          Append (~T, F);
     end for;
return T;
end function;

nullspace_mem_to_polynomial := function (nullspace_row, n_rows, n_cols, n_depth)
  // A nullspace row having non-zero value corresponds to a triple (i,j,k) being part of a recurrence relation.
  // Translate the position into a monomial term.

  // The data of the rows of the dense matrix is organized first ranging over columns, then rows, then depths, hence the necessary striding.
  i_stride := (n_depth+1) * (n_cols + 1);
  j_stride := n_cols + 1;
  trait := 0;
  for i in [0..n_rows] do
    for j in [0..n_cols] do
      for k in [0..n_depth] do
        monomial := nullspace_row[k + j*j_stride + i*i_stride + 1] * x^i * y^j * z^k;
        trait := trait + monomial;
      end for;
    end for;
  end for;
  return trait;
end function;

// Creates a matrix by computing the result of outer_action(X^i * T * Y^j, Z^k) as (i,j,k) ranges 
// Each value `outer_action(X^i * T * Y^j, Z^k)` is flattened as a single row, and we first range through rows, then cols, then slices - this striding is then used above in `nullspace_mem_to_polynomial`.
compute_recurrence_ideal := function (S, X, Y, Z)
  n_rows := #Rows(X);
  n_cols := #Rows(Y);
  n_slices := #Rows(Z);
  flattened_row_collection := [];
  for i in [0..n_rows] do
    for j in [0..n_cols] do
      before_z := [X^i * M * Y^j :  M in S];
      for k in [0..n_slices] do
        // print <i,j,k>;
        after_z := outer_action(before_z, Z^k);
        flattened_row := Eltseq(VerticalJoin(after_z));
        Append(~flattened_row_collection, flattened_row);
      end for;
    end for;
  end for;

  flattened_as_matrix := Matrix(R, (n_rows+1)*(n_cols+1)*(n_slices+1), n_rows*n_cols*n_slices, flattened_row_collection);
  // print flattened_as_matrix;
  nullspace := Nullspace(flattened_as_matrix);
  nullspace_polys := [nullspace_mem_to_polynomial(nullspace.i, n_rows, n_cols, n_slices) : i in [1..Dimension(nullspace)]];
  J := ideal<P | nullspace_polys>;
  _ := GroebnerBasis(J);
  return J, flattened_as_matrix, nullspace;
end function;

//Computes the set I(S,Ops) given a finite set of systems of forms S, operators embedded into a block diagonal matrix, and tuple of dimensions.
Compute_I_set := function(S, Ops, dims)
  a := dims[1];
  b := dims[2];
  c := dims[3];
  offset := [1,a+1,a+b+1];
  blocks := [* [ ExtractBlock(X,offset[i],offset[i],dims[i],dims[i]) : X in Ops ] : i in [1..3] *]; //grabs each block from diagonal. Stores in [X,Y,Z] where e.g X is the top left block of each

  I := ideal<P | 1>;
  for t in S do
    for i in [1..#blocks[1]] do
      I := I meet compute_recurrence_ideal(t, blocks[1][i], blocks[2][i], blocks[3][i]);
    end for;
  end for;

  return I;
end function;

ROWS:=2;
COLS:=2;
DEPTH:=2;

// Example from Figure 2.2 of https://arxiv.org/pdf/1911.02518
slice1 := [1,0, 0,0];
slice2 := [0,0, 0,1];
MS := RMatrixSpace(R, ROWS,COLS);
M1 := MS!slice1;
M2 := MS!slice2;
GHZ := [M1, M2];

W_slice1 := MS![0,1, 1,0];
W_slice2 := MS![1,0, 0,0];
W := [W_slice1, W_slice2];

X_MS := RMatrixSpace(R, ROWS, ROWS);
X := Matrix(R, ROWS,ROWS, [0,1, 1,0]);

Y_MS := RMatrixSpace(R, COLS, COLS);
Y := X;

Z_MS := RMatrixSpace(R, COLS, COLS);
Z := X;

GHZ_trait := compute_recurrence_ideal(GHZ, X, Y, Z);
GHZ_trait eq ideal<P | x^2-1, y^2-1,z^2-1,x*y-z,x-y*z,y-x*z>;

W_trait := compute_recurrence_ideal(W, X, Y, Z);
W_trait eq ideal<P | x^2-1, y^2-1,z^2-1>;

Ops := [DiagonalJoin(X,DiagonalJoin(Y,Z))];
I := Compute_I_set([GHZ, W], Ops, [ROWS,COLS,DEPTH]);