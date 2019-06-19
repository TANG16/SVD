function [Tnew,tmm1,U] = qralg1(A)
    m=length(A);
    tmm1(1)=abs(A(m,m-1));
    U=eye(m);
    %wilkinson shift
    while tmm1(end)>=10^-15
        u=0;
        [Q,R] = qr(A-u.*eye(m));
        A=R*Q+u*eye(m);
        U=U*Q;
        tmm1=[tmm1,abs(A(m,m-1))];
    end
    
    Tnew=A;
end