// Computing equations for a Z-set

//returns the sum of every tensor in tensorlist.
tensorsum := function(tensorlist, valence)
    if #tensorlist eq 1 then
        return tensorlist[1];
    elif valence eq 2 then
        newtensor := tensorlist[1] + tensorlist[2];
        for i in [3..#tensorlist] do
            newtensor := newtensor + tensorlist[i];
        end for;
        return newtensor;
    else
        newtensor := [];
        for i in [1..#tensorlist[1]] do
            entrytensor := &+ [tensorlist[j][i] : j in [1..#tensorlist]];
            Append(~newtensor, entrytensor);
        end for;
        return newtensor;
    end if;
end function;


// Valence 2 case

//poly: polynomial in 2 variables
//return the image of this polynomial action on T via X and Y and the permuted version.
compute_2_valent_polynomial_action := function(poly,T,X,Y)
    tensormons := [];
    for mon in Monomials(poly) do
        powers := Exponents(mon);
        coeff := MonomialCoefficient(poly, mon);
        montensor := coeff * (Transpose(X)^(powers[1]) * T * Y^powers[2]);
        Append(~tensormons, montensor);
    end for;
    ptensor := tensorsum(tensormons,2);

    return ptensor;
end function;

// compute_2_valent_operators := function(poly,T,X,Y)
//     LHS := [];
//     for T in Basis do
//         polytensor := compute_2_valent_polynomial_action(poly,T,X,Y);
//         Append(~LHS,polytensor);
//     end for;

//     eqtnlist := [];
//     for index in [1..#LHS] do
//         for d in [1..dim] do
//             LHS2 := Eltseq(LHS[index][d]);
//         end for;
//         for entry in [1..#LHS2] do
//             Append(~eqtnlist, LHS2[entry]);
//         end for;
//     end for;

//     S := Scheme(affinespace,eqtnlist);
//     return S;

// end function;


// Valence 3 case


outer_action := function (S, z)
     T := [ ];
     for i in [1..#S] do 
          // Magma notation: &+ sums every element in a list.
          F := &+ [ z[i][j] * S[j] : j in [1..#S] ];
          Append (~T, F);
     end for;
return T;
end function;

//S : System of forms
//returns the transverse action of a monomial and operator on T
compute_action := function(T, X, Y, Z, coeff)
    before_z := [coeff * Transpose(X) * M * Y : M in T]; //can transpose X if needed.
    // print <i,j,k>;
    after_z := outer_action(before_z, Z);
    flattened_row := Eltseq(VerticalJoin(after_z));
return after_z;
end function;

compute_3_valent_polynomial_action := function(poly,T,X,Y,Z)
    tensormons := [];
    for mon in Monomials(poly) do
        powers := Exponents(mon);
        coeff := MonomialCoefficient(poly, mon);
        montensor := compute_action(T, X^powers[1], Y^powers[2], Z^powers[3], coeff);

        Append(~tensormons, montensor);
    end for;
    ptensor := tensorsum(tensormons,3);

    return ptensor;
end function;

tensors_over_coordinate_ring := function(tensors, coordinatering, isform)
    newscalars := [];
    for bas in tensors do
        newbas := [];
        if isform then
            sys := bas;
        else
            sys := SystemOfForms(bas);
        end if;
        for mat in sys do
            Append(~newbas,ChangeRing(mat, coordinatering));
        end for;
        Append(~newscalars, newbas);
    end for;

    return newscalars;
end function;


//Basis: Basis for the tensor space.
//returns the affine scheme of operators carved by preserving symmetry under a collection of polynomials and a matrix algebra.
compute_3_valent_operators := function(B_forms, affinespace, valence, poly, X, Y, Z, frame)
    LHS := [];
    for T in B_forms do //now this is a symmetric basis
        polytensor := compute_3_valent_polynomial_action(poly,T,X,Y,Z);
        Append(~LHS,polytensor);

        // print("tensor");
        // print(polytensor);
    end for;

    eqtnlist := [];
    for index in [1..#LHS] do
        for d in [1..#LHS[index]] do
            eqtnlist := eqtnlist cat Eltseq(LHS[index][d]);
        end for;
    end for;

    S := Scheme(affinespace,eqtnlist);
    return S;
end function;


// User code
compute_Z_set := function(k, frame, valence, TS, tensors, polys, isSystemOfForms);
    if valence eq 3 then
        aff := AffineSpace(k, frame[1]^2 +frame[2]^2 +frame[3]^2 );
        coord := CoordinateRing(aff);
    else
        aff := AffineSpace(k, (frame[1]+frame[2])*valence);
        coord := CoordinateRing(aff);
    end if;

    X := Matrix(coord, frame[1], frame[1], [coord.i : i in [1..frame[1]^2]]);
    Y := Matrix(coord, frame[2], frame[2], [coord.(i + frame[1]^2) : i in [1..frame[2]^2]]);
    if valence eq 3 then
        Z := Matrix(coord, frame[3], frame[3], [coord.(i + (frame[1]^2 + frame[2]^2)) : i in [1..frame[3]^2]]);
        P<x,y,z> := PolynomialRing(k,3);
        tensors_as_forms := tensors_over_coordinate_ring(tensors, coord, isSystemOfForms);

        S := aff;
        for poly in polys do
            S := Intersection(S,compute_3_valent_operators(tensors_as_forms, aff, valence, poly, X, Y, Z, frame));
        end for;
    end if;
    return S;
end function;


generate_bowtie_t := function(K,links)
  VS := VectorSpace(K, 2 * links + 1);
  V1 := VectorSpace(K, 1);
  nonzeros := {};
  for i in [1..links] do
    nonzeros := nonzeros join Permutations({1,2*i,2*i+1});
  end for;
  t_bowtie_action := function(x)
    arg1 := x[1];
    arg2 := x[2];
    arg3 := x[3];
    acc := 0;
    for nonzero_inds in nonzeros do
      i := nonzero_inds[1]; j := nonzero_inds[2]; k := nonzero_inds[3];
      acc := acc + (arg1[i] * arg2[j] * arg3[k]);
    end for;
    return [acc];
  end function;

  t_bowtie := Tensor([VS,VS,VS,V1], t_bowtie_action);
  return Compress(Shuffle(t_bowtie,[1,0,2,3]));
end function;


ROWS:=2;
COLS:=2;
DEPTH:=2;

// k := RationalField();
k := GF(997);
frame := [ROWS,COLS,DEPTH];
valence := #frame; //only works for valence 2 or 3.

P<x,y,z> := PolynomialRing(k,valence);
MS := RMatrixSpace(k, ROWS,COLS);
T := KTensorSpace(k,frame);

//3 - valent example(s)

//bowtie tensor
// bowtie := generate_bowtie_t(k,2);
// MS := KMatrixSpace(k,5,5);
// t1 := [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0 ]; //(3,1,2) & (2,1,3) & (5,1,4) & (4,1,5)
// t2 := [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]; //(3,2,1) & (1,2,3)
// t3 := [ 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]; //(2,3,1) & (1,3,2)
// t4 := [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ]; //(5,4,1) & (1,4,5)
// t5 := [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]; //(4,5,1) & (1,5,4)
// full_t := t1 cat t2 cat t3 cat t4 cat t5;
// BT := T!full_t;
// SystemOfForms(BT);
// D := DerivationAlgebra(BT);
// Dimension(D);
// DirectSumDecomposition(D);

// polys := [y-x*z, x-z];

// S := compute_Z_set(k,frame,valence, T, tensors, polys, true);

//w tensor
// W_slice1 := MS![0,1, 1,0];
// W_slice2 := MS![1,0, 0,0];
W_slice1 := [0,1, 1,0];
W_slice2 := [1,0, 0,0];
full_t := W_slice1 cat W_slice2;
W := T!full_t;
// W := [W_slice1, W_slice2];
Derivation_ideal := ideal<P | x+y-z>;
D := DerivationAlgebra(W);
Dimension(D);
DirectSumDecomposition(D);


// S := compute_Z_set(k,frame,valence,T,[W], Generators(Derivation_ideal), true);