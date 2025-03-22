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
//returns the transverse action of X,Y,Z on t
compute_action := function(S, X, Y, Z, coeff)
    before_z := [coeff * Transpose(X) * M * Y : M in S]; //can transpose X if needed.
    // print <i,j,k>;
    after_z := outer_action(before_z, Z);
    flattened_row := Eltseq(VerticalJoin(after_z));
return after_z;
end function;

//S1,S2: Systems of forms
//return the sum of S1 and S2 added together by sheets
sumTwoTensors := function(S1,S2)
    sumoftensors := [];
    for i in [1..#S1] do
        Append(~sumoftensors, S1[i] + S2[i]);
    end for;

    return sumoftensors;
end function;

//S1,S2,S3: Systems of forms
//return the sum of S1 and S2 added together by sheets
sumThreeTensors := function(S1,S2,S3)
    sumoftensors := [];
    for i in [1..#S1] do
        Append(~sumoftensors, S1[i] + S2[i] + S3[i]);
    end for;

    return sumoftensors;
end function;

compute_2_valent_e1_action := function(k,dim,T,X,Y)
    e1tensor := (Transpose(X) * T) + (T * Y);
    e1with12actiontensor := (Transpose(Y) * Y) + (T * X);

    return e1tensor,e1with12actiontensor;
end function;

compute_3_valent_e1_action := function(k,dim,T,X,Y,Z)
    XLtensor := compute_action(T,X,IdentityMatrix(k,dim),IdentityMatrix(k,dim),1);
    XRtensor := compute_action(SystemOfForms(T),IdentityMatrix(k,dim), Y, IdentityMatrix(k,dim),1);
    XDtensor := compute_action(SystemOfForms(T),IdentityMatrix(k,dim),IdentityMatrix(k,dim),Z,1);
    e1tensor := sumThreeTensors(XLtensor,XRtensor,XDtensor);

    XL12tensor := compute_action(SystemOfForms(T),Y,IdentityMatrix(k,dim),IdentityMatrix(k,dim),1);
    XR12tensor := compute_action(SystemOfForms(T),IdentityMatrix(k,dim), X, IdentityMatrix(k,dim),1);
    XD12tensor := compute_action(SystemOfForms(T),IdentityMatrix(k,dim),IdentityMatrix(k,dim),Z,1);
    e1with12actiontensor := sumThreeTensors(XL12tensor,XR12tensor,XD12tensor);

    XL23tensor := compute_action(SystemOfForms(T),X,IdentityMatrix(k,dim),IdentityMatrix(k,dim),1);
    XR23tensor := compute_action(SystemOfForms(T),IdentityMatrix(k,dim), Z, IdentityMatrix(k,dim),1);
    XD23tensor := compute_action(SystemOfForms(T),IdentityMatrix(k,dim),IdentityMatrix(k,dim),Y,1);
    e1with23actiontensor := sumThreeTensors(XLtensor,XRtensor,XDtensor);

    return e1tensor,e1with12actiontensor,e1with23actiontensor;
end function;

k := RationalField();
valence := 2; //can be 2 or 3
dim := 2;
if valence eq 2 then
    frame := [2,2];
else 
    frame := [2,2,2];
end if;

aff := AffineSpace(k, dim^2 *valence);
coord := CoordinateRing(aff);
if valence eq 3 then
    TS := KTensorSpace(k,frame);
else
    TS := RMatrixSpace(k,dim,dim);
end if;

X := Matrix(coord, dim, dim, [coord.i : i in [1..dim^2]]);
Y := Matrix(coord, dim, dim, [coord.(i + dim^2) : i in [1..dim^2]]);
if valence eq 3 then
    Z := Matrix(coord, dim, dim, [coord.(i + 2*dim^2) : i in [1..dim^2]]);
else
    Z := IdentityMatrix(k,dim);
end if;

LHS := [];
RHS := [];
for T in Basis(TS) do //e1 action
    if valence eq 2 then
        e1tensor, e1with12actiontensor := compute_2_valent_e1_action(k,dim,T,X,Y);
        Append(~LHS,e1tensor);
        Append(~RHS, e1with12actiontensor);
    else 
        e1tensor, e1with12actiontensor, e1with23actiontensor := compute_3_valent_e1_action(k,dim,T,X,Y,Z);
        Append(~LHS,e1tensor);
        Append(~RHS, e1with12actiontensor);
        Append(~LHS,e1tensor);
        Append(~RHS, e1with23actiontensor);
    end if;
end for;

// for M in Basis(A) do //e2 action
//     Append(~LHS,Transpose(X) * M * Y);
//     Append(~RHS, Transpose(Y) * M * X);
// end for;

// for M in Basis(A) do //e3 action
//     Append(~LHS,Transpose(X) * M * Y);
//     Append(~RHS, Transpose(Y) * M * X);
// end for;

eqtnlist := [];
for index in [1..#LHS] do
    for d in [1..dim] do
        LHS2 := Eltseq(LHS[index][d]);
        RHS2 := Eltseq(RHS[index][d]);
    end for;
    for entry in [1..#LHS2] do
        Append(~eqtnlist, LHS2[entry] - RHS2[entry]);
    end for;
end for;

S := Scheme(aff,eqtnlist);
S;
IrreducibleComponents(S);