function [Tnew,tmm1] = qralgwilkinson(A)
    m=length(A);
    tmm1(1)=abs(A(m,m-1));
    
    %wilkinson shift
    while tmm1(end)>=10^-12
        sigma=(A(m-1,m-1)-A(m,m))/2;
        if sign(sigma)==0
            sgn=-1
        else
            sgn=sign(sigma);
        end
        u=A(m,m)-sgn*(A(m-1,m)^2)/(abs(sigma)+sqrt(sigma^2+A(m-1,m)^2));
        [Q,R] = qr(A-u.*eye(m));
        A=R*Q+u*eye(m);
        tmm1=[tmm1,abs(A(m,m-1))];
    end
    
    Tnew=A;
end