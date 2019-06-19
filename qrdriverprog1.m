function [U, Eval, V, tmm] = qrdriverprog1(A,wilkinson)
[h,N] = size(A);
[ T, QQ ] = symtridhh( A );
m=length(T);
Tnew=T;
tmv=[];
Eval=[];
U=QQ;
for i=m:-1:2 
   if wilkinson==1
       [Tnew,tmm1,Q]=qralgwilkinson1(full(Tnew));
   else
       [Tnew,tmm1,Q]=qralg1(full(Tnew));
   end
   Eval=[Eval Tnew(i,i)];
   U(:,1:i)=U(:,1:i)*Q;
   tmv=[tmv tmm1];
   Tnew=Tnew(1:i-1,1:i-1);
end

Eval = [Eval Tnew(1,1)];
tmm=tmv;

Eval=Eval';
Eval=sort(Eval.^(1/2),'descend');
V=zeros(N,N);
% for i=1:m
% V(1:i)=A\(Eval(i)*U(:,i));
% end
%V=A\(diag(Eval)*U);
end