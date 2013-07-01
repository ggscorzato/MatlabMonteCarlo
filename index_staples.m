function res=index_staples(L)

D=length(L); send=1; recv=2; up=1; dn=2;
for mu=1:D
  temp(1:2,mu,up,send)=[mu,dn];
  temp(1:2,mu,dn,recv)=[mu,dn];
  temp(1:2,mu,dn,send)=[mu,up];
  temp(1:2,mu,up,recv)=[mu,up];
end

res==reshape(temp,2*D,2,2);