function [ B, U, V ] = bidighh( A );

M = size( A, 1 ); % rows
N = size( A, 2 ); % cols
if M < N
    error( 'A must have M >= N' );
end
B = A; % B will be M x N
betas = zeros( N, 1 );
gammas = zeros( N-2, 1 );

for k = 1:N
    [ v, beta, mu ] = house( B(k:M,k) );
    betas( k ) = beta;
    B(k:M,k:N) = B(k:M,k:N) - beta *v * v' * B(k:M,k:N);
    B((k+1):M,k) = v(2:end);
    if k < N - 1
        [ v, beta, mu ] = house( B(k,(k+1):N)' );
        gammas( k ) = beta;
        B(k:M,(k+1):N) = B(k:M,(k+1):N) - beta * B(k:M,(k+1):N) * v *v';
        B(k,(k+2):N) = v(2:end)';
    end

end
if nargout > 1
    % compute U and V
    v = zeros( M, 1 );
    U = eye( M );
    for k = N:-1:1
        v(k) = 1;
        v(k+1:M) = B((k+1):M,k);
        U(k:M,k:M) = U(k:M,k:M) - betas( k ) * v(k:M) * v(k:M)' * U(k:M,k:M);
        B(k+1:M,k) = 0;
    end
    V = eye( N );
    for k = N-1:-1:2
        v(k) = 1;
        v(k+1:N) = B(k-1,(k+1):N)';
        V(k:N,k:N) = V(k:N,k:N) - gammas( k-1 ) * v(k:N) * v(k:N)' * V(k:N,k:N);
        B(k-1,k+1:N)= 0;
    end
end
B = [ diag(B) [ diag(B,1) ; 0 ] ];

end