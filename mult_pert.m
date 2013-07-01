function C= mult_pert(A,B)

[d d p]=size(A);
C=zeros(d,d,p);

for n=0:p-1
  for k=0:p-1
    C(:,:,n +1)=C(:,:,n +1) + A(:,:,k +1) * B(:,:,n-k +1);
  end
end