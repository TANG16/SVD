function T = tridiag(A)

m = size(A,1);

for k = 1:m-2
    v = A(k+1:m,k);
    v(1) = sign(v(1)) * norm(v) + v(1);
    v = v / norm(v);
    A(k+1:m,k:m) = A(k+1:m,k:m) - 2 * v * (v' * A(k+1:m,k:m));
    A(1:m,k+1:m) = A(1:m,k+1:m) - 2 * (A(1:m,k+1:m) * v) * v';
end

% Forcing to be symmetric
v = diag(A,-1); 
T = spdiags([[v; 0] diag(A) [0; v]],[-1 0 1],m,m);
end