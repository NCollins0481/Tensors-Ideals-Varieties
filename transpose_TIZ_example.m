// Want to test whether matrix of your choice is in the span of powers of matrices X^i Y^j.

R := RationalField(); //Ring is Q to start
R := AlgebraicClosure(R); //If you wish to work over the algebraic closure

dimX := 2; //X is 2x2
dimY := 2; //Y is 3x3

tdim := dimX*dimY;
target := Matrix(R, tdim, tdim, [1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1]);




xvars := ["x" * IntegerToString(j) * IntegerToString(i): i in [1..dimX], j in [1..dimX]];
yvars := ["y" * IntegerToString(j) * IntegerToString(i): i in [1..dimY], j in [1..dimY]];
vars := xvars cat yvars;
scalarvars := ["l" * IntegerToString(j) * IntegerToString(i): i in [0..(dimX-1)], j in [0..(dimY-1)]];
allvars := vars cat scalarvars;

S := PolynomialRing(R, #allvars);
AssignNames(~S, allvars);




X := Matrix(S, dimX, dimX, [S.i: i in [1..dimX^2]]);
Y := Matrix(S, dimY, dimY, [S.i: i in [(dimX^2)+1..(dimX^2)+(dimY^2)]]);
// X; Y;

Matrixpowers := [];
counter := 1;
for i in [0..(dimX-1)] do
    for j in [0..(dimX-1)] do
        Append(~Matrixpowers, S.((dimX^2)+(dimY^2)+counter) * TensorProduct(X^i, Y^j));
        counter := counter + 1;
    end for;
end for;
// Matrixpowers;

polysystem := [];

for i in [1..tdim] do
    for j in [1..tdim] do
        Append(~polysystem, &+[M[i,j] : M in Matrixpowers] - target[i,j]);
    end for;
end for;
// polysystem;
// #polysystem;




sys := ideal<S | Set(polysystem)>;
// sys;
GroebnerBasis(sys); // 1 in the Groebner basis <=> no solutions to your system