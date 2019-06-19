clear all

m=150;
Ejacobi=zeros(m,1);
Egk=zeros(m,1);
Eqr=zeros(m,1);
Eqrshift=zeros(m,1);
for j=1:10
    A1=rand(m,m);
    [U,S,V] = svd(A1);
    [U,Sjacobi,V, convergence] = jacobi_SVD(A1);
    Ejacobi = (j*Ejacobi+abs(diag(S)-sort(diag(Sjacobi),'descend')))./(j+1);
    [Sgk,U,V] = gksvdsteps(A1);
    Egk = (j*Egk+abs(diag(S)-sort((Sgk),'descend')))./(j+1);
    [Sqr, tmm] = qrdriverprog(A1*A1',0);
    Eqr = (j*Eqr+abs(diag(S)-sort((Sqr),'descend')))./(j+1);
    [Sqrshift, tmm] = qrdriverprog(A1*A1',1);
    Eqrshift = (j*Eqrshift+abs(diag(S)-sort((Sqrshift),'descend')))./(j+1);
end

sval=1:m;
%%
semilogy(sval,Ejacobi)
 hold on
 semilogy(sval,Egk)
 hold on
 semilogy(sval,Eqr)
 hold on
 semilogy(sval,Eqrshift)
legend("Jacobi","Golub-Kahan","QR without Shift", "QR with Shift", 'Location','southwest')
xlabel("Number of singular value (largest to smallest)")
ylabel("Absolute error in each singular values")