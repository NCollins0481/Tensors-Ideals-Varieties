R := GF(997);

outer_action := function (S, z)
     T := [ ];
     for i in [1..#S] do 
          // Magma notation: &+ sums every element in a list.
          F := &+ [ z[i][j] * S[j] : j in [1..#S] ];
          Append (~T, F);
     end for;
return T;
end function;

seq_of_matrix_to_seq := function(S) 
  acc := [];
  for slice in S do
    for v in Eltseq(slice) do
      Append(~acc, v);
    end for;
  end for;
  return acc;
end function;

seq_of_matrix_to_col_vec := function(S)
  return Matrix(R, ROWS * COLS * DEPTH, 1, seq_of_matrix_to_seq(S));
end function;

CUBE_SIDE := 3; DEG_CUTOFF := 3;
ROWS := CUBE_SIDE; COLS := CUBE_SIDE; DEPTH := CUBE_SIDE;

// Compute for a cubic tensor (CUBE_SIDE^3 entries), fixed operators X,Y, and tensors s,t, whether there exists a polynomial p, such that p(X,Y) \cdot s = t where p is at most DEG_CUTOFF.
compute_coeffs := function(X,Y,Z,s,t)
  solution_vector := seq_of_matrix_to_col_vec(t);
  col_vecs := [];
  for k in [0..DEG_CUTOFF-1] do
    cur_s := outer_action(s, Z^k);
    for j in [0..DEG_CUTOFF-1] do
      for i in [0..DEG_CUTOFF-1] do
        // j is power of Y, i is power of X
        Append(~col_vecs, seq_of_matrix_to_col_vec([X^i * s_slice * Y^j : s_slice in cur_s]));
      end for;
    end for;
  end for;
  // Stack the flattened vectors of (X^i t Y^j)^(Z^k) as columns, and then create an augmented matrix for which taking its RREF allows us to read off solutions or to verify no solution exists.
  coeff_matrix := HorizontalJoin(col_vecs);
  augmented := HorizontalJoin(coeff_matrix, solution_vector);
  print "has solution if rank equal... ", Rank(augmented) = Rank(coeff_matrix);
  rrefed := EchelonForm(augmented);
  return rrefed;
end function;


MS := MatrixAlgebra(R, CUBE_SIDE);
// s and t are conjugate, s and r are not.
R_op := Random(MS); C_op := Random(MS); D_op := Random(MS);
s := [Random(MS) : i in [1..CUBE_SIDE]];
t := outer_action([R_op * s_slice * C_op: s_slice in s], D_op);
r := [Random(MS) : i in [1..CUBE_SIDE]];

// In general, for random X,Y,Z triple, there is a polynomial solution - but it's not very meaningful. 
X := Random(MS);Y := Random(MS);Z := Random(MS);
rrefed := compute_coeffs(X,Y,Z,s,t);

// Sanity check solution - this does work
proposed_solution := ColumnSubmatrix(rrefed, DEG_CUTOFF^3+1, 1); // This gets the last column.
t_prime := [ZeroMatrix(R, CUBE_SIDE) : i in [1..CUBE_SIDE]];
for k in [0..DEG_CUTOFF-1] do
  outer_actioned := outer_action(s, Z^k);
  for j in [0..DEG_CUTOFF-1] do
    for i in [0..DEG_CUTOFF-1] do
      idx := k * DEG_CUTOFF^2 + j*DEG_CUTOFF + i + 1;
      output := [proposed_solution[idx][1] * X^i * slice * Y^j : slice in outer_actioned];

      for l in [1..#output] do
        t_prime[l] := t_prime[l] + output[l];
      end for;

    end for;
  end for;
end for;