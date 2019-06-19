function [ B, U, V,mycon ] = gksvdstep1( B, U, V )

N = size( B, 1 );

if N > 2
    mu = wilkinsonshift( B(N-1,1)^2+B(N-2,2)^2, B(N-1,1)*B(N-1,2), B(N,1)^2+B(N-1,2)^2 );
elseif N==2
    mu = wilkinsonshift( B(N-1,1)^2, B(N-1,1)^B(N-1,2), B(N,1)^2+B(N^1,2)^2 );
else
    error('Cannot perform GK SVD step on B with N < 2')
end

x = B( 1, 1 )^2 - mu;
y = B( 1, 1 ) * B( 1, 2 ) ;
if nargin < 3
    V = eye( N );
    if nargin < 2
        U = V;
    end
end
mycon=[];
for k = 1:N-1
    [ c s ] = givens( x, y ) ;
    G = [ c s; -s c ];
    V(:, k:k+1) = V(:, k:k+1) * G;
    if k > 1
        B( k-1, 2 ) = c * B( k-1, 2 ) - s * bulge;
    end
    Bk = B( k, 1 ) ;
    bulge = -s * B( k+1, 1 );
    B( k , 1 ) = c * Bk - s * B(k,2) ; 
    B( k , 2 ) = s * Bk + c * B( k, 2 ) ;
    B( k+1, 1 ) = c * B( k+1, 1 ); 
    x = B( k, 1 ) ;
    y = bulge; 
    
    [ c s ] = givens( x, y ) ;
    G = [ c s; -s c ];
    U(:, k:k+1) = U(:, k:k+1 ) * G;
    B( k , 1 ) = c * B( k, 1 ) - s * bulge;
    Bk2 = B( k, 2 ) ;
    B( k , 2 ) = c * Bk2 - s * B( k+1, 1 );
    bulge = -s * B( k+1, 2 );
    B( k+1, 1 ) = s * Bk2 + c * B( k+1, 1 );
    B( k+1, 2 ) = c * B( k+1, 2 );
    x = B( k, 2 ) ;
    mycon=[mycon abs(B( k, 2 ))];
    y = bulge;
end
B( N, 2 ) = 0;

end