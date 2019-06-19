function [ A, Q ] = symtridhh( A )

if any(any(A ~= A'))
    error('Input must be symmetric!')
end

N = size( A, 1 );
for k = 1:(N-2)
    [ v, beta, mu ] = house( A( (k+1):N, k ) );
    p = A( (k+1):N, (k+1):N )*(beta*v); 
    w = p - (beta/2*p'*v)*v;
    A( k, (k+2):N ) = 0 ; A(k,k+1) = mu;
    if nargout == 2
        A( k+1:N, k ) = v ; A( k+1, k ) = beta;
    else
        A( (k+2):N, k ) = 0 ; A(k+1,k) = mu;
    end
    A( (k+1):N, (k+1):N ) = A( (k+1):N, (k+1):N ) - v*w' - w*v';
end

if nargout == 2
    Q = eye( N );
    for k = (N-1):-1:1
        v = A( k+1:N, k );
        beta = v( 1 ) ;
        v( 1 ) = 1;
        Q(k+1:N, k+1:N) = Q( k+1:N, k+1:N ) - beta * (v * v') * Q( k+1:N, k+1:N );

        A(k+1:N, k) = A(k, k+1:N )';
    end
end

end