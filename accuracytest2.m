
Ejacobi=zeros(5,1);
Egk=zeros(5,1);
Eqr=zeros(5,1);
Eqrshift=zeros(5,1);
for j=0:10
j
i=0;
size=[];
for m=10:10:50
    i=i+1
    A1=magic(m);
    [U,S,V] = svd(A1);
    [U,Sjacobi,V, convergence] = jacobi_SVD(A1);
    Ejacobi(i) = (j*Ejacobi(i)+norm(diag(S)-sort(diag(Sjacobi),'descend'))/m)/(j+1);
    %[Sgk,U,V] = gksvdsteps(A1);
    %Egk(i) = (j*Egk(i)+norm(diag(S)-sort(Sgk,'descend'))/m)/(j+1);
    [Sqr, tmm] = qrdriverprog(A1*A1',0);
    Eqr(i) = (j*Eqr(i)+norm(diag(S)-sort(Sqr,'descend'))/m)/(j+1);
    [Sqrshift, tmm] = qrdriverprog(A1*A1',1);
    Eqrshift(i) = (j*Eqrshift(i)+norm(diag(S)-sort(Sqrshift,'descend'))/m)/(j+1);
    size=[size m];
    
end
end
%%
 semilogy(size,Ejacobi)
 hold on
 %semilogy(size,Ejacobi/10)
 %hold on
 semilogy(size,Eqr)
 hold on
 semilogy(size,Eqrshift)
legend("Jacobi","QR without Shift", "QR with Shift", 'Location','southeast')
xlabel("Size of matrix")
ylabel("Mean squared error")