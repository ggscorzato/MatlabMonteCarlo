function [gi, gv, soi, sov, G, SO] = cliff_gen(t,s)
% gamma matrices and so(t,s) generators with space-time metric (t(-),s(+))
% they are used on a vector v by:
% (Gamma_mu * v)(j) = gv(mu,j) v(gi(mu,j)).
% (SO_mu_nu * v)(j) = sov(mu,nu,j) v(soi(mu,nu,j)).

% space time metric
d=s+t;
dm=floor(d/2);
eta=eye(d);
eta(1:t,1:t)=-eta(1:t,1:t);

% pauli matrices
s0=[1,0;0,1];s1=[0,1;1,0];s2=[0,-i;i,0];s3=[1,0;0,-1];

% d dimensional euclidean gamma matrices
tu = s0;
for j=1:dm-1
  tu = kron(tu,s0);
end
one.m = tu;  % unfortunately matlab indices strat with 1.
             % Mind that in mathematica code I have Gamma_0=one.

for l=1:dm
  to = s1;
  te = s2;
  for j=1:l-1
    to = kron(s3,to);
    te = kron(s3,te);
  end
  for j=1:dm-l
    to = kron(to,s0);
    te = kron(te,s0);
  end
  G(2*l-1).m = to;
  G(2*l).m = te;
end
tc=s3;
for j=1:dm-1
  tc = kron(tc,s3);
end
G(2*dm+1).m=tc;

% !!!! CHECK that it is possible !!! 
% economical gamma matrices: gv and gi are such that 
% (G * v)(j) = gv(j) v(gi(j)).

for l=1:2*dm+1
  for j=1:2^dm
    gi(l,j) = find(G(l).m(j,:)); 
    gv(l,j) = G(l).m(gi(l,j),j);
  end
end

% gamma matrices in (t,s) metric

for l=1:t
  G(l).m  = i * G(l).m;
  gv(l,:) = i * gv(l,:);
end

% so(t,s) generators

for mu=1:d
  for nu=mu+1:d
    SO(mu,nu).m = -1/4 * (G(mu).m * G(nu).m - G(nu).m * G(mu).m);
    SO(nu,mu).m = -SO(mu,nu).m;
    for j=1:2^dm
      soi(mu,nu,j) = find(SO(mu,nu).m(j,:)); 
      sov(mu,nu,j) = SO(mu,nu).m(soi(mu,nu,j),j);
      soi(nu,mu,j) =  soi(mu,nu,j);
      sov(nu,mu,j) = -sov(mu,nu,j);
    end
  end
end

% A, B, C ... maybe