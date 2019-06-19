function [U,S,V, convergence] = jacobi_SVD(A)
  M = size( A, 1 );
  N = size( A, 2 );
  transposed=0;
  if M>N
     A=A';
     transposed=1;
  end

  TOL=1.e-15;
  n=size(A,1);
  U=A';
  V=eye(n);
  convergence=[];
  converge=TOL+1;
  while converge>TOL
  converge=0;
  for i=1:n-1
    for j=i+1:n
      % compute [alpha gamma;gamma beta]=(i,j) submatrix of U*U'
      alpha=sum(sum(U(i, :).*U(i, :))); 
      beta=sum(sum(U(j, :).*U(j, :)));
      gamma=sum(U(i, :).* U(j, :));
      converge=max(converge,abs(gamma)/sqrt(alpha*beta));
      convergence=[convergence abs(gamma)/sqrt(alpha*beta)];
      % compute Jacobi rotation that diagonalizes 
      %    [alpha gamma;gamma beta]
      zeta=(beta-alpha)/(2*gamma);
      t=sign(zeta)/(abs(zeta)+sqrt(1+zeta^2));
      c= 1.0 / (sqrt(1 + t * t));
      s= c * t;

      % update columns i and j of U
      t=U(i, :);
      U(i, :)=c*t-s*U(j, :);
      U(j, :)=s*t+c*U(j, :);

      % update matrix V of right singular vectors
      t=V(i, :);
      V(i, :)=c*t-s*V(j, :);
      V(j, :)=s*t+c*V(j, :);
    end
  end
end

% the singular values are the norms of the columns of U
% the left singular vectors are the normalized columns of U
for j=1:n
  singvals(j)=norm(U(j, :));
  U(j, :)=U(j, :)/singvals(j);
end
S=diag(singvals);
U = U';
V = V'; %return V, not V'

if transposed==1
    VV=V;
    V=U';
    U=VV';
end

end