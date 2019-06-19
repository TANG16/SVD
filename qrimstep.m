function [ TS, Q ] = qrimstep( TS, Q )

N = size( TS, 1 );

mu = wilkinsonshift( TS( N-1, 1 ), TS( N-1, 2 ), TS( N, 1 ) );

x = TS(1, 1) - mu;
y = TS(1, 2);

if nargin < 2
    Q = eye( N );
end;
for k = 1:(N-1)
    [ c s ] = givens( x, y ) ; 

    G = [ c s; -s c ];
    K = min( N, k+2 );
    Q(:, k:k+1 ) = Q(:, k:k+1 ) * G;
    if k > 1
        TS( k-1, 2 ) = c * TS( k-1, 2 ) - s * bulge;
    end
    Tk1 = TS( k, 1 );
    Tk2 = TS( k, 2 );
    bulge = -s * TS( k+1, 2 );
    c2 = c^2; s2 = s^2; sc2 = 2*c*s;
    TS( k, 1 ) = c2 * Tk1 + s2 * TS( k+1, 1 ) - sc2*Tk2;
    TS( k, 2 ) = c*s*( Tk1 - TS( k+1, 1 ) ) + ( c2 - s2 ) * Tk2;
    TS( k+1, 1 ) = s2 * Tk1 + c2 * TS( k+1, 1 ) + sc2*Tk2;
    TS( k+1, 2 ) = c * TS( k+1, 2 );
    x = TS( k, 2 );
y = bulge;
end
TS( N, 2 ) = 0;

end