clear all
A1= rand(30,30)
[ D, U, V,c ] = gksvdsteps1( A1 );
plot(real(c))
ylabel("lower off-diagonal entries of submatrices")
xlabel("Number of iterations")
title("Golub-Kahan Iterations - Convergence")

[U,Sjacobi,V, c] = jacobi_SVD(A1);
figure;
plot(c)
ylabel("lower off-diagonal entries of submatrices")
xlabel("Number of iterations")
title("Jacobi Iterations - Convergence")

[Sqr,c] = qrdriverprog(diag(30:-1:1)+ones(30,30),0);
figure;
semilogy(c(1:end-1))
ylabel("lower off-diagonal entries of submatrices")
xlabel("Number of iterations")
title("QR without shift iterations- Convergence")

[Sqr,c] = qrdriverprog(A1*A1',1);
figure;
semilogy(c(1:end-1))
ylabel("lower off-diagonal entries of submatrices")
xlabel("Number of iterations")
title("QR with shift iterations- Convergence")