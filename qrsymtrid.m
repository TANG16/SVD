function [ Q, R ] = qrsymtrid( T )

N = size( T, 1 );
if N ~= size( T, 2 ) || any(any(T ~= T'))
    error( 'input T must be symmetric' );
end;
R = T;
Q = eye( N );
for k=1:(N-1)
    [ c , s ] = givens( R( k,k ) , R( k+1, k ) ) ;
    G = [ c s; -s c ];
    K = min( N, k+2 );
    R( k:k+1, k:K ) = G' * R( k:k+1, k:K ); 
    R( k+1, k ) = 0; 
    Q(1:k+1, k:k+1) = Q(1:k+1, k:k+1) * G;
end

end