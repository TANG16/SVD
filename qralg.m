function [Tnew,tmm1] = qralg(A)
    m=length(A);
    tmm1(1)=abs(A(m,m-1));
    
     while tmm1(end)>=10^-15
             [Q,R] = qr(A);
             A=R*Q;
             tmm1=[tmm1,abs(A(m,m-1))];
     end
    Tnew=A;
end