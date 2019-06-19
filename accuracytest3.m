clear all

Ematlab=zeros(10,1);
Ejacobi=zeros(10,1);
Egk=zeros(10,1);
Eqr=zeros(10,1);
Eqrshift=zeros(10,1);
for j=0:10
i=0;
size=[];
for m=10:10:100
    i=i+1
    A1=rand(m,m);
    [U,S,V] = svd(A1);
    Ematlab(i) = (j*Ematlab(i)+norm(U'*U-eye(m))/m)/(j+1);
    [U,Sjacobi,V, convergence] = jacobi_SVD(A1);
    Ejacobi(i) = (j*Ejacobi(i)+norm(U'*U-eye(m))/m)/(j+1);
    [Sgk,U,V] = gksvdsteps(A1);
    Egk(i) = (j*Egk(i)+norm(U'*U-eye(m))/m)/(j+1);
    [U,S,V,tmm] = qrdriverprog1(A1*A1',1);
    Eqr(i) = (j*Eqr(i)+norm(U'*U-eye(m))/m)/(j+1);
    [U,S,V,tmm] = qrdriverprog1(A1*A1',1);
    Eqrshift(i) = (j*Eqrshift(i)+norm(U'*U-eye(m))/m)/(j+1);
    size=[size m];
    
end
end
%%
 semilogy(size,Ematlab)
 hold on
 semilogy(size,Ejacobi)
 hold on
 semilogy(size,Egk)
 hold on
 semilogy(size,Eqr)
 hold on
 semilogy(size,Eqrshift)
legend("Maltab","Jacobi","Golub-Kahan","QR without Shift", "QR with Shift", 'Location','northeast')
xlabel("Size of matrix")
ylabel("norm( U^TU ? I_{m x m} )/M")