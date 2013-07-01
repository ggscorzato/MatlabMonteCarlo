function [a] = sun_gen(N,norm)
% if a = sun_gen(N,norm)
% a(i).m is the i-th (i=1,N^2-1) su(N) generator in the adjoint rep. (generalized Gell-Mann)
% a(N^2).m is proportional to the identity matrix.
% they are all normalized as Tr ( a(i).m * a(j).m ) = norm* delta_ij 
% gell-mann: a=sun_gen(3,2); lambda_3= a(7).m, lambda_8= a(8).m, 

for n = 1:N^2
 a(:,:,n) = zeros(N);
end

n=1;
for h=1:N
  for k=1:h-1
    a(h,k,n) = 1;
    a(k,h,n) = 1;
    n=n+1;
    a(h,k,n) = i;
    a(k,h,n) =-i;
    n=n+1;
  end
end
% cartan sub-algebra
for k=1:N-1
 a([1:k],[1:k],n) = eye(k); 
 a(k+1,k+1,n) = -k; 
 n=n+1;
end

% u(1) gen
a(:,:,N^2) = eye(N);

% normalization
for n = 1:N^2
 no=trace(a(:,:,n) * a(:,:,n));
 a(:,:,n) = a(:,:,n) * sqrt(norm/no);
end
