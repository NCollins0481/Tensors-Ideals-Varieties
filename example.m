// TIZ set example computation.

Q := RationalField();
Z := Integers();
// Some other rings I've tried around with
// Z := NumberField(elt<PolynomialAlgebra(Rationals()) | [-8, 0, 1]>);
// Z := GF(3);
// Z := GF(997);
P<x,y> := PolynomialRing(Z, 2);

nullspace_mem_to_polynomial := function (nullspace_row, mat_n_rows, mat_n_cols)
  // Keep track of number of rows/cols to translate back to monomials.
  // Ex ( 0  0  0  0  1  0 -1  0  0  0  0  0)
  // Translates to "1" at position (1,1) and (0,2), which is xy -y^2
  trait := 0;
  for i in [0..mat_n_cols] do
    for j in [0..mat_n_rows] do
      trait := trait + (nullspace_row[i*(mat_n_rows+1)+j+1] * x^j * y^i);
      // print [i,j, i*(mat_n_rows+1)+j+1, monomial];
    end for;
  end for;
  return trait;
end function;

compute_recurrence_ideal := function (M, row_op, col_op)
  n_rows := Dimension(Domain(M));
  n_cols := Dimension(Codomain(M));
  Xi_M_yj_seq_of_vectors := [ Eltseq(row_op^i * M * col_op^j) : i in [0..n_rows], j in [0..n_cols]];
  Xi_M_yj_matrix := Matrix(Z, (n_rows+1)*(n_cols+1), n_rows*n_cols, Xi_M_yj_seq_of_vectors);
  nullspace := Nullspace(Xi_M_yj_matrix);
  nullspace_polys := [nullspace_mem_to_polynomial(nullspace.i, n_rows, n_cols) : i in [1..Dimension(nullspace)]];
  J := ideal<P | nullspace_polys>;
  _ := GroebnerBasis(J);
  return J;
end function;

compute_recurrence_ideal_without_axis_polynomials := function(M, row_op, col_op)
  J := compute_recurrence_ideal(M, row_op, col_op);
  mp_row := MultivariatePolynomial(P, MinimalPolynomial(row_op), 1);
  mp_col := MultivariatePolynomial(P, MinimalPolynomial(col_op), 2);

  axes_polynomials := ideal<P | mp_row, mp_col>;
  J_mod_axes := J / axes_polynomials;
  return J, axes_polynomials, J_mod_axes;
end function;

ROWS:=2;
COLS:=3;
// Example from Section 7.3 of https://arxiv.org/pdf/1911.02518
v := [1, 2, 3, 2, 3, 0];
MS := RMatrixSpace(Z, ROWS,COLS);
M := MS!v;

X_MS := RMatrixSpace(Z, ROWS, ROWS);
X := Matrix(Z, ROWS,ROWS, [0,1, 0,0]);

Y_MS := RMatrixSpace(Z, COLS, COLS);
Y := Matrix(Z, COLS,COLS, [0,0,0, 1,0,0, 0,1,0]);

J, axes_ideal, J_mod_axes := compute_recurrence_ideal_without_axis_polynomials(M, X, Y);
J_mod_axes_empty := IsZero(J_mod_axes);
