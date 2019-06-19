function [Eval, tmm] = qrdriverprog(A,wilkinson)

T=tridiag(A);
m=length(T);
Tnew=T;
tmv=[];
Eval=[];
for i=m:-1:2 
   if wilkinson==1
       [Tnew,tmm1]=qralgwilkinson(full(Tnew));
   else
       [Tnew,tmm1]=qralg(full(Tnew));
   end
   Eval=[Eval Tnew(i,i)];
   tmv=[tmv tmm1];
   Tnew=Tnew(1:i-1,1:i-1);
end

Eval = [Eval Tnew(1,1)];
tmm=tmv;

Eval=Eval';
Eval=Eval.^(1/2);

end