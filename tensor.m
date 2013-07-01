function c= tensor(a,b)
% TENSOR Tensor product
% if X and Y are 2 dimensional arrays, 
% Z=TENSOR(X,Y) is the tensor product of X and Y
% i.e. Z(i,j,k,l) = X(i,j) * Y(k,l)

c=permute(reshape(kron(a,b),size(b,1),size(a,1),size(b,2),size(a,2)),[2 4 1 3]);


