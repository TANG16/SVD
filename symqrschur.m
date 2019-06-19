function [ D, Q ] = symqrschur( A )

tol = 10^-15;

[ T, Q ] = symtridhh( A );
if nargin == 1
    tol = 1e-14;
end
%TS = stridl2s ( T ) ;
TS=T;
N = size( TS, 1 );
q = 0;
while q < N
    q = 0; p = N-1;
    for k = N-1:-1:1
        if abs( TS( k, 2 ) ) <= tol * ( abs( TS( k,1 ) ) + abs( TS( k+1, 1 ) ) )
            TS( k, 2 ) = 0;
            if q == N-k-1
                q = N-k;
                p = k-1;
            end
            else
            if p == k
                p = k-1;
            end
        end
    end

    if q == N-1
        q = N;
    else
        [ TS( p+1:N-q, : ), Q(:, p+1:N-q ) ] = qrexstep( TS( p+1:N-q, : ), Q(:, p+1:N-q ) );
    end
end
D = TS( :,1 );

end