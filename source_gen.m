function eta = source_gen(L,Nc,DomainType,color_flag,spin_flag,spacetime_flags,Rep)
% eta = source_gen()
%
% generate a set of N spin-color fields eta(c,s,x,i) (where c=1:Nc, s=2:S, x=...L, i=1:N).
% with value in DomainType where
% L=lattice sizes vector;
% Nc = number of colors;
% S = number of spin components;
% DomainType can be: 1 (one), n (random-ZN), 0 (random-gauss), 0.5 (Z2+iZ2);
% color_flag can be: 
% c0 (delta for c=c0), Nc+1 (repeated delta forall c), Nc+2 (indep. foreach c), Nc+3 (constant forall c);
% spin_flag can be: 
% s0 (delta for s=s0), S+1 (repeated delta forall s), S+2 (indep. foreach s), S+3 (constant forall s);
% spacetime_flag can be: 
% x0 (delta for x(i)=x0(i)), x0(i)=L(i)+1 (repeated for all x0(i)), x0(i)=L(i)+2 (indep. foreach x0(i)), 
% x0(i)=L(i)+3 (constant for all x0(i)).
% 