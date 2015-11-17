function n = MV3norm(A)
%expects 3 x n set of vectors, returns array of n vector norms
n = sqrt(A(1,:).^2 + A(2,:).^2 + A(3,:).^2);
