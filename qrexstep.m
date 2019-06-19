function T = qrexstep( T )

ident = eye( size( T ) );
mu = wilkinsonshift ( T( end-1, end-1 ), T( end, end-1 ), T( end, end ) );
[Q,R] = qrsymtrid( T - mu * ident ); 
T = R*Q + mu * ident

end