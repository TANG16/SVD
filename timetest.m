clear all
m=100;
n=50;

Tmatlab=zeros(15,1);
Tjacobi=zeros(15,1);
Tgk=zeros(15,1);
Tqr=zeros(15,1);
Tqrshift=zeros(15,1);

for j=1:10
    
    i=0;
    size=[];
    for m=10:10:150
        i=i+1;
        A1=rand(m,m);

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

plot(size,Tmatlab)
 hold on
 plot(size,Tjacobi)
 hold on
 plot(size,Tgk)
 hold on
 plot(size,Tqr)
 hold on
 plot(size,Tqrshift)
 
 legend("Matlab","Jacobi","Golub-Kahan","QR without Shift", "QR with Shift", 'Location','northwest')
xlabel("Size of matrix")
ylabel("Time to convergence (seconds)")
axis([0 inf 0 1])