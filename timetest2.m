clear all
m=100;
n=50;

Tmatlab=zeros(10,1);
Tjacobi=zeros(10,1);
Tgk=zeros(10,1);
Tqr=zeros(10,1);
Tqrshift=zeros(10,1);

for j=1:10
    
    i=0;
    size=[];
    for s=0.1:0.1:1
        i=i+1;
        A1=full(sprand(m,m,s));
        

        tic;
        [U,S,V] = svd(A1);
        time=toc;
        Tmatlab(i)=(j*Tmatlab(i)+time)./(j+1);

        tic;
        [U,S,V, convergence] = jacobi_SVD(A1);
        time=toc;
        Tjacobi(i)=(j*Tjacobi(i)+time)./(j+1);

        tic;
        [D,U,V] = gksvdsteps(A1);
        time=toc;
        Tgk(i)=(j*Tgk(i)+time)./(j+1);


        tic;
        [Eval, tmm] = qrdriverprog(A1*A1',0);
        time=toc;
        Tqr(i)=(j*Tqr(i)+time)./(j+1);

        tic;
        [Eval, tmm] = qrdriverprog(A1*A1',1);
        time=toc;
        Tqrshift(i)=(j*Tqrshift(i)+time)./(j+1);
        
        size=[size m];
    end
end


%%
s=0.1:0.1:1
semilogy(s,Tmatlab)
 hold on
 semilogy(s,Tjacobi)
 hold on
 semilogy(s,Tgk)
 hold on
 semilogy(s,Tqr)
 hold on
 semilogy(s,Tqrshift)
 
 legend("Matlab","Jacobi","Golub-Kahan","QR without Shift", "QR with Shift", 'Location','northwest')
xlabel("Size of matrix")
ylabel("Time to convergence (seconds)")
%axis([0 inf 0 1])