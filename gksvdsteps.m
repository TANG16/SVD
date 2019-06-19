function [ D, U, V ] = gksvdsteps( A )
M = size( A, 1 );
N = size( A, 2 );
transposed=0;
if N>M
    A=A';
    transposed=1;
end

[ B, U, V ] = bidighh( A );

tol = 10^-15;

M = size( U, 1 ); 
N = size( V, 1 );

q = 0;
while q < N
    q = 0; p = N-1;
    for k = N-1:-1:1
        if abs( B(k,2) ) <= tol*( abs( B(k,1) ) + abs( B(k+1,1) ) )
            B(k,2) = 0;
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
        q=N; 
    else 
        k = p+1;
        smalldiag = tol * norm( B, 'inf' ) ;
        while abs( B( k, 1 ) ) > smalldiag && k < N-q
            k = k+1;
        end

        if abs( B( k, 1 ) ) <= smalldiag
        B( k, 1 ) = 0;

            if k < N-q
                bulge = B(k,2);
                B(k,2)= 0;
                for j=k+1:N-q
                    [ c s ] = givens ( B(j,1) , bulge ) ;

                    B(j,1) = -s * bulge + c * B(j,1);
                    bulge = s * B(j,2) ;
                    B(j,2) = c * B(j,2) ;
                    G = [ c s; -s c ];
                    U(:,[ k j ]) = U(:,[k j ]) * G'; 
                end
            else 

            bulge = B(N-q-1,2);
            B(N-q-1,2)= 0;
                for j=N-q-1:-1:p+1
                    [ c s ] = givens( B(j,1) , bulge ) ;
                    B(j,1) = c * B(j,1) - s * bulge;
                    if j>p+1
                        bulge = s * B(j-1,2);
                        B(j-1,2)= c * B(j-1,2);
                    end

                    G = [ c s; -s c ];
                    V(:,[ j k]) = V(:,[ j k]) * G;
                end
            end 
        else 
            [ B(p+1:N-q,:), U(:,p+1:N-q), V(:,p+1:N-q) ] = gksvdstep( B(p+1:N-q,:), U(:,p+1:N-q), V(:,p+1:N-q) );
        end
    end
end 
D = B( :,1 );
for k = 1:N
    if D(k) < 0
        D(k) = -D(k);
        V(:,k) = -V(:,k);
    end
end
if transposed==1
    VV=V;
    V=U';
    U=VV';
end
V=real(V);
U=real(U);
D=real(D);
end