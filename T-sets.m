// Trying to implement with a nullspace computation.
R := GF(997);
P<x,y, z> := PolynomialRing(R, 3);

// Given a single operator, a polynomial ideal, and dims [a,b,c], compute the tensors t (flat arrays) such that Ops is a subset of Der(t).
Compute_T_relations := function( K, Polys, Ops, dims )
  /* Assumption:
   *   dims = [a, b, c]
   *   T : K^a x K^b >-> K^c
   *   #Ngens(Ops) = d
   * Solves abc linear equations with abcd variables. 
   */

  // Ops := [D : D in Generators(Delta_Lie)];
  // K := BaseRing(Delta_Lie);
  d := 1; // d := #Ops;
  a := dims[1];
  b := dims[2];
  c := dims[3];

  vprint Densor, 1 : "Setting linear system: " cat IntegerToString(a*b*c) cat " by " cat IntegerToString(a*b*c*d);

  offset := [1,a+1,a+b+1];
  blocks := [* [ ExtractBlock(X,offset[i],offset[i],dims[i],dims[i]) : X in Ops ] : i in [1..3] *]; //grabs each block from diagonal. Stores in [X,Y,Z] where e.g X is the top left block of each
  Z := ZeroMatrix(K, a*b*c, a*b*c*d);
  Y := ZeroMatrix(K, a*b*c*d, a*b*c);
  X := ZeroMatrix(K, a*b*c*d, a*b*c);
  
  // Z Block. Writing the Z block as id \otimes id \otimes Z
  jpos := 1;
  for i in [1..d] do
    ipos := 1;
    for j in [1..a*b] do
      InsertBlock(~Z,blocks[3][i],ipos,jpos);
      ipos +:= c;
      jpos +:= c;
    end for;
  end for;

  // Y Block. Writing Y as id \otimes Y \otimes id
  ipos := 1;
  for i in [1..d] do
    jpos := 1;
    for j in [1..a] do
      Yblock := [];
      for k in [1..b] do
        vec := &cat[ [x] cat [0 : i in [1..c-1]] : x in Eltseq(blocks[2][i][k]) ];
        for l in [1..c] do
          Append(~Yblock,vec);
          Remove(~vec,#vec);  
          vec := [0] cat vec;
        end for;      
      end for;
      InsertBlock(~Y,Matrix(K,Yblock),ipos,jpos);
      ipos +:= b*c;
      jpos +:= b*c;
    end for;
  end for;

  // X Block. Writing X as X \otimes id \otimes id
  ipos := 1;
  for i in [1..d] do
    for j in [1..a] do
      vec := &cat[[x] cat [0 : i in [1..b*c-1]] : x in Eltseq(blocks[1][i][j])];
      for k in [1..b*c] do
        InsertBlock(~X,Matrix(K,1,a*b*c,vec),ipos,1);
        Remove(~vec,#vec);
        vec := [0] cat vec;
        ipos +:= 1;
      end for;
    end for;
  end for;

  // Compute the each matrix of operators for a polynomial
  polylist := [];
  for poly in Generators(Polys) do
    matrixlist := [];
    for mon in Monomials(poly) do
      powers := Exponents(mon);
      coeff := MonomialCoefficient(poly, mon);
      Append(~matrixlist, coeff * Transpose(X)^powers[1] * Transpose(Y)^powers[2] * Z^powers[3]);
    end for;
    matrixpoly := 0;
    for M in matrixlist do
      matrixpoly := matrixpoly + M;
    end for;
    Append(~polylist, matrixpoly);
  end for;

  delete X, Y, Z;

  vprint Densor, 1 : "Solving linear system: " cat IntegerToString(a*b*c) cat " by " cat IntegerToString(a*b*c*d);

  Null := Nullspace(polylist[1]);
  for i in [2..#polylist] do
    Null := Null meet Nullspace(polylist[i]);
  end for;

  return Null;
end function;

Compute_T_set := function( T, Polys, Ops )
  // remove redundant dims if coords are fused.
  ess_dims := [ Dimension(X) : X in T`Frame ];
  RC := {};
  for P in T`Cat`Repeats do
    m := Maximum(P);
    Include(~RC, m);
    for p in P diff {m} do
      ess_dims[Valence(T) - p] := 0;
    end for;
  end for;
  ess_dims := [d : d in ess_dims | d ne 0] cat [0];
  dims := [ Dimension(X) : X in T`Frame ] cat [0];

  matched_fused := Nrows(Ops[1]) eq &+ess_dims and Ncols(Ops[1]) eq &+ess_dims;
  matched_not_fused := Nrows(Ops[1]) eq &+dims and Ncols(Ops[1]) eq &+dims;

  if matched_not_fused then
    RC := {0..2};
  end if;

  // Package the Delta into a Lie algebra 
  Delta_Lie := sub< MatrixLieAlgebra(BaseRing(T), Nrows(Ops[1])) | Ops >;
  DerivedFrom(~Delta_Lie, T!0, {0..2}, RC : Fused := matched_fused); // Really just need the tensor category info

  Ops := [D : D in Generators(Delta_Lie)];
  Null := Compute_T_relations(BaseRing(Delta_Lie), Polys, [Ops[1]], dims);
  for index in [2..#Ops] do
    Null := Null meet Compute_T_relations(BaseRing(Delta_Lie), Polys, [Ops[index]], dims);
  end for;
  N := BasisMatrix(Null);
  S := TensorSpace(T`Frame, T`Cat);
  S`Mod := sub< T`Mod | [ T`Mod!N[i] : i in [1..Nrows(N)] ] >;

  return S;
end function;

ROWS:=2;
COLS:=2;
DEPTH:=2;
dims := [ROWS,COLS,DEPTH];
T := RTensorSpace(R, dims);
MS := RMatrixSpace(R, ROWS,COLS);


// Example from Figure 2.2 of https://arxiv.org/pdf/1911.02518
// slice1 := [1,0, 0,0];
// slice2 := [0,0, 0,1];
// MS := RMatrixSpace(R, ROWS,COLS);
// M1 := MS!slice1;
// M2 := MS!slice2;
// GHZ := [M1, M2];

// W_slice1 := MS![0,1, 1,0];
// W_slice2 := MS![1,0, 0,0];
// W := [W_slice1, W_slice2];

// GHZ_ideal := ideal<P | x^2-1, y^2-1,z^2-1,x*y-z,x-y*z,y-x*z>;
// W_ideal := ideal<P | x^2-1, y^2-1,z^2-1>;
// Derivation_ideal := ideal<P | x+y-z>;
// Group_ideal := ideal<P | y-x*z, x-z>;

// X_MS := RMatrixSpace(R, ROWS, ROWS);
// X := Matrix(R, ROWS,ROWS, [0,1, 1,0]);

// Y_MS := RMatrixSpace(R, COLS, COLS);
// Y := X;

// Z_MS := RMatrixSpace(R, DEPTH, DEPTH);
// Z := X;


// ops := [DiagonalJoin(X,DiagonalJoin(Y,Z))];

// // ops := Basis(DerivationAlgebra(T![0,1,1,0,1,0,0,0]));

// // Compute_T_set(T, Derivation_ideal, ops);
// for b in Basis(Compute_T_set(T, Group_ideal, ops)) do
//   SystemOfForms(b);
// end for;

//Matrix example
// X := Matrix(R, ROWS,ROWS, [0,1, 1,0]);

// polys := ideal<P | x-y>;
// ops := [DiagonalJoin(X,DiagonalJoin(-X,IdentityMatrix(R,2)))];
// for b in Basis(Compute_T_set(T, polys, ops)) do
//   SystemOfForms(b);
// end for;

//More w example with solvable part of the derivation
W_slice1 := MS![0,1, 1,0];
W_slice2 := MS![1,0, 0,0];
W := [W_slice1, W_slice2];
flattened_W := Eltseq(W[1]) cat Eltseq(W[2]);
Wb := T!flattened_W;

Derivation_ideal := ideal<P | x+y-z>;

X_MS := RMatrixSpace(R, ROWS, ROWS);
X1 := Matrix(R, ROWS,ROWS, [1,0, 0,0]);
X2 := Matrix(R, ROWS,ROWS, [0,1, 0,0]);
X3 := Matrix(R, ROWS,ROWS, [0,0, 0,0]);

Y_MS := RMatrixSpace(R, COLS, COLS);
Y1 := Matrix(R, ROWS,ROWS, [0,0, 0,-1]);
Y2 := Matrix(R, ROWS,ROWS, [0,0, 0,0]);
Y3 := Matrix(R, ROWS,ROWS, [0,1, 0,0]);

Z_MS := RMatrixSpace(R, DEPTH, DEPTH);
Z1 := Matrix(R, ROWS,ROWS, [0,0, 0,1]);
Z2 := Matrix(R, ROWS,ROWS, [0,0, 1,0]);
Z3 := Matrix(R, ROWS,ROWS, [0,0, 1,0]);

X4 := Matrix(R, ROWS,ROWS, [1,0, 0,1]);
Y4 := Matrix(R, ROWS,ROWS, [0,0, 0,0]);
Z4 := Matrix(R, ROWS,ROWS, [1,0, 0,1]);

X5 := Matrix(R, ROWS,ROWS, [0,0, 0,0]);
Y5 := Matrix(R, ROWS,ROWS, [1,0, 0,1]);
Z5 := Matrix(R, ROWS,ROWS, [1,0, 0,1]);

//DiagonalJoin(X4,DiagonalJoin(Y4,Z4)),DiagonalJoin(X5,DiagonalJoin(Y5,Z5))
ops := [DiagonalJoin(X1,DiagonalJoin(Y1,Z1)),DiagonalJoin(X2,DiagonalJoin(Y2,Z2)),DiagonalJoin(X3,DiagonalJoin(Y3,Z3))];

// ops := Basis(DerivationAlgebra(T![0,1,1,0,1,0,0,0]));

// Compute_T_set(T, Derivation_ideal, ops);
for b in Basis(Compute_T_set(T, Derivation_ideal, ops)) do
  SystemOfForms(b);
end for;

X := Matrix(R, ROWS,ROWS, [0,1, 0,0]);
Y := Matrix(R, ROWS,ROWS, [0,1, 0,0]);
Z := Matrix(R, ROWS,ROWS, [0,0, 1,0]);

ops := [IdentityMatrix(R,6), DiagonalJoin(X,DiagonalJoin(Y,Z))];
Centroid_ideal := ideal<P | x-y, y-z>;
for b in Basis(Compute_T_set(T, Centroid_ideal, ops)) do
  SystemOfForms(b);
end for;